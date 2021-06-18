use crate::dashboard::{StatusUpdate, StatusUpdateSender};
use crate::squared_topologies::CutkoskyCut;
use crate::{FloatLike, JetSliceSettings, ObservableMode, SelectorMode, Settings};
use itertools::Itertools;
use libc::{c_double, c_int, c_void};
use lorentz_vector::LorentzVector;
use num::Complex;
use smallvec::SmallVec;
use std::fs::File;
use std::io::{BufWriter, Write};

mod fjcore {
    use libc::{c_double, c_int, c_void};

    #[link(name = "fjcore", kind = "static")]
    extern "C" {
        pub fn fastjet_workspace() -> *mut c_void;
        //pub fn fastjet_free(workspace: *mut c_void);
        pub fn fastjetppgenkt_(
            workspace: *mut c_void,
            p: *const c_double,
            npart: *const c_int,
            R: *const c_double,
            ptjet_min: *const c_double,
            palg: *const c_double,
            jets: *mut c_double,
            njets: *mut c_int,
            whichjets: *mut c_int,
        );
    }
}

#[derive(Debug, Default, Clone)]
pub struct EventInfo {
    pub accepted_event_counter: usize,
    pub rejected_event_counter: usize,
    pub no_phase_space_counter: usize,
    pub zero_eval_counter: usize
}

#[derive(Debug, Default, Clone)]
pub struct AverageAndErrorAccumulator {
    sum: f64,
    sum_sq: f64,
    weight_sum: f64,
    avg_sum: f64,
    avg: f64,
    err: f64,
    guess: f64,
    chi_sq: f64,
    chi_sum: f64,
    chi_sq_sum: f64,
    num_samples: usize,
    cur_iter: usize,
}

impl AverageAndErrorAccumulator {
    pub fn add_sample(&mut self, sample: f64) {
        self.sum += sample;
        self.sum_sq += sample * sample;
        self.num_samples += 1;
    }

    pub fn merge_samples(&mut self, other: &mut AverageAndErrorAccumulator) {
        self.sum += other.sum;
        self.sum_sq += other.sum_sq;
        self.num_samples += other.num_samples;

        // reset the other
        other.sum = 0.;
        other.sum_sq = 0.;
        other.num_samples = 0;
    }

    pub fn update_iter(&mut self) {
        // TODO: we could be throwing away events that are very rare
        if self.num_samples < 2 {
            self.cur_iter += 1;
            return;
        }

        let n = self.num_samples as f64;
        self.sum /= n;
        self.sum_sq /= n * n;
        let mut w = (self.sum_sq * n).sqrt();

        w = ((w + self.sum) * (w - self.sum)) / (n - 1.);
        if w == 0. {
            w = std::f64::EPSILON;
        }
        w = 1. / w;

        self.weight_sum += w;
        self.avg_sum += w * self.sum;
        let sigsq = 1. / self.weight_sum;
        self.avg = sigsq * self.avg_sum;
        self.err = sigsq.sqrt();
        if self.cur_iter == 0 {
            self.guess = self.sum;
        }
        w *= self.sum - self.guess;
        self.chi_sum += w;
        self.chi_sq_sum += w * self.sum;
        self.chi_sq = self.chi_sq_sum - self.avg * self.chi_sum;

        // reset
        self.sum = 0.;
        self.sum_sq = 0.;
        self.num_samples = 0;
        self.cur_iter += 1;
    }
}

#[derive(Default)]
pub struct EventManager {
    pub event_selector: Vec<Selectors>,
    pub observables: Vec<Observables>,
    pub event_buffer: Vec<Event>,
    pub track_events: bool,
    pub accepted_event_counter: usize,
    pub rejected_event_counter: usize,
    pub no_phase_space_counter: usize,
    pub zero_eval_counter: usize,
    pub time_integrand_evaluation: bool,
    pub integrand_evaluation_timing: u128,
    pub integrand_evaluation_timing_start: Option<std::time::Instant>,
    pub event_group_counter: usize,
    pub status_update_sender: Option<StatusUpdateSender>,
}

impl EventManager {
    pub fn new(
        track_events: bool,
        settings: Settings,
        status_update_sender: StatusUpdateSender,
    ) -> EventManager {
        let mut observables = vec![];
        for o in &settings.observables.active_observables {
            match o {
                ObservableMode::Jet1PT => {
                    if track_events {
                        observables.push(Observables::Jet1PT(Jet1PTObservable::new(
                            settings.observables.Jet1PT.x_min,
                            settings.observables.Jet1PT.x_max,
                            settings.observables.Jet1PT.n_bins,
                            settings.observables.Jet1PT.dR,
                            settings.observables.Jet1PT.min_jpt,
                            settings.observables.Jet1PT.write_to_file,
                            settings.observables.Jet1PT.filename.clone(),
                            settings.observables.Jet1PT.use_fastjet,
                        )));
                    }
                }
                ObservableMode::CrossSection => {
                    observables.push(Observables::CrossSection(CrossSectionObservable::new(
                        status_update_sender.clone(),
                    )));
                }
            }
        }

        let mut selectors = vec![];
        for s in &settings.selectors.active_selectors {
            match s {
                SelectorMode::Jet => {
                    selectors.push(Selectors::Jet(JetSelector::new(&settings.selectors.jet)));
                }
            }
        }

        EventManager {
            event_selector: selectors,
            observables,
            event_buffer: Vec::with_capacity(100),
            track_events,
            accepted_event_counter: 0,
            rejected_event_counter: 0,
            no_phase_space_counter: 0,
            zero_eval_counter : 0,
            time_integrand_evaluation: false,
            integrand_evaluation_timing: 0,
            integrand_evaluation_timing_start: None,
            event_group_counter: 0,
            status_update_sender: Some(status_update_sender),
        }
    }

    pub fn add_event<T: FloatLike>(
        &mut self,
        orig_incoming_momenta: &[LorentzVector<f64>],
        cut_momenta: &[LorentzVector<T>],
        cut_info: &[CutkoskyCut],
        rot_matrix: &[[f64; 3]; 3],
    ) -> bool {
        if !self.track_events && self.event_selector.is_empty() {
            return true;
        }

        let mut incoming_momenta = SmallVec::new();

        // rotate all momenta with the inverse of the rotation matrix
        for e in orig_incoming_momenta {
            incoming_momenta.push(LorentzVector::from_args(
                e.t,
                rot_matrix[0][0] * e.x + rot_matrix[1][0] * e.y + rot_matrix[2][0] * e.z,
                rot_matrix[0][1] * e.x + rot_matrix[1][1] * e.y + rot_matrix[2][1] * e.z,
                rot_matrix[0][2] * e.x + rot_matrix[1][2] * e.y + rot_matrix[2][2] * e.z,
            ));
        }

        let mut outgoing_momenta = SmallVec::new();
        let mut final_state_particle_ids = SmallVec::new();
        for (cut_mom, cut) in cut_momenta[..cut_info.len()].iter().zip(cut_info.iter()) {
            // make sure all momenta are outgoing
            let e = cut_mom.cast::<f64>() * cut.sign as f64;

            outgoing_momenta.push(LorentzVector::from_args(
                e.t,
                rot_matrix[0][0] * e.x + rot_matrix[1][0] * e.y + rot_matrix[2][0] * e.z,
                rot_matrix[0][1] * e.x + rot_matrix[1][1] * e.y + rot_matrix[2][1] * e.z,
                rot_matrix[0][2] * e.x + rot_matrix[1][2] * e.y + rot_matrix[2][2] * e.z,
            ));
            final_state_particle_ids.push(cut.particle_id);
        }

        let mut e = Event {
            kinematic_configuration: (incoming_momenta, outgoing_momenta),
            final_state_particle_ids,
            integrand: Complex::new(1., 0.),
            weights: SmallVec::from_slice(&[0.]),
        };

        // run the event through the selectors
        for f in &mut self.event_selector {
            if !f.process_event(&mut e) {
                if self.track_events {
                    self.rejected_event_counter += 1;
                }
                return false;
            }
        }

        if self.track_events {
            self.event_buffer.push(e);
            self.accepted_event_counter += 1;
        }
        true
    }

    pub fn merge_samples(&mut self, other: &mut EventManager) {
        for (o1, o2) in self
            .observables
            .iter_mut()
            .zip(other.observables.iter_mut())
        {
            o1.merge_samples(o2);
        }

        self.accepted_event_counter += other.accepted_event_counter;
        self.rejected_event_counter += other.rejected_event_counter;
        self.no_phase_space_counter += other.no_phase_space_counter;
        self.zero_eval_counter += other.zero_eval_counter;
        self.event_group_counter += other.event_group_counter;
        self.integrand_evaluation_timing = other.integrand_evaluation_timing;
        other.accepted_event_counter = 0;
        other.rejected_event_counter = 0;
        other.event_group_counter = 0;
        other.no_phase_space_counter = 0;
        other.zero_eval_counter = 0;
        other.integrand_evaluation_timing = 0;
    }

    pub fn update_result(&mut self) {
        for o in &mut self.observables {
            o.update_result(self.event_group_counter);
        }

        self.status_update_sender
            .as_mut()
            .unwrap()
            .send(StatusUpdate::EventInfo(EventInfo {
                accepted_event_counter: self.accepted_event_counter,
                rejected_event_counter: self.rejected_event_counter,
                no_phase_space_counter: self.no_phase_space_counter,
                zero_eval_counter: self.zero_eval_counter
            }))
            .unwrap();
    }

    pub fn update_live_result(&mut self) {
        for o in &mut self.observables {
            // for now, only live update the cross section
            match o {
                Observables::CrossSection(c) => c.update_live_result(),
                _ => {}
            }
        }

        self.status_update_sender
            .as_mut()
            .unwrap()
            .send(StatusUpdate::EventInfo(EventInfo {
                accepted_event_counter: self.accepted_event_counter,
                rejected_event_counter: self.rejected_event_counter,
                no_phase_space_counter: self.no_phase_space_counter,
                zero_eval_counter : self.zero_eval_counter
            }))
            .unwrap();
    }

    pub fn clear(&mut self, count_as_rejected: bool) {
        self.accepted_event_counter -= self.event_buffer.len();
        if count_as_rejected {
            self.rejected_event_counter += self.event_buffer.len();
        }
        self.event_buffer.clear();
    }

    pub fn process_events(&mut self, integrand_result: Complex<f64>, integrand_weight: f64) {
        if self.track_events {
            // give the events to an observable function
            for o in &mut self.observables {
                o.process_event_group(&self.event_buffer, integrand_weight);
            }
            self.event_buffer.clear();
        } else {
            for o in &mut self.observables {
                if let Observables::CrossSection(c) = o {
                    c.add_sample(integrand_result, integrand_weight);
                    break;
                }
            }
        }

        self.event_group_counter += 1;
    }
}

#[derive(Default, Debug, Clone)]
pub struct Event {
    pub kinematic_configuration: (
        SmallVec<[LorentzVector<f64>; 2]>,
        SmallVec<[LorentzVector<f64>; 4]>,
    ),
    pub final_state_particle_ids: SmallVec<[isize; 5]>,
    pub integrand: Complex<f64>,
    pub weights: SmallVec<[f64; 1]>,
}

pub enum Selectors {
    All(NoEventSelector),
    Jet(JetSelector),
}

impl Selectors {
    #[inline]
    fn process_event(&mut self, event: &mut Event) -> bool {
        match self {
            Selectors::All(f) => f.process_event(event),
            Selectors::Jet(f) => f.process_event(event),
        }
    }
}

pub trait EventSelector {
    /// Process a group of events and return a new integrand value that will be returned to the integrator.
    fn process_event(&mut self, event: &mut Event) -> bool;
}

#[derive(Default)]
pub struct NoEventSelector {}

impl EventSelector for NoEventSelector {
    #[inline]
    fn process_event(&mut self, _event: &mut Event) -> bool {
        true
    }
}

#[derive(Debug, Clone)]
pub struct JetClustering {
    //ordered_jets: Vec<LorentzVector<f64>>,
    ordered_pt: Vec<f64>,
    //jet_structure: Vec<usize>, // map from original momenta to jets
    d_r: f64,
    min_jpt: f64,
    use_fastjet: bool,
    fastjet_jets_in: Vec<f64>,
    fastjet_jets_out: Vec<f64>,
    fastjet_jets_map: Vec<c_int>,
    fastjet_workspace: *mut c_void,
}

unsafe impl std::marker::Send for JetClustering {}

impl JetClustering {
    pub fn new(use_fastjet: bool, d_r: f64, min_jpt: f64) -> JetClustering {
        let fastjet_workspace = unsafe { fjcore::fastjet_workspace() }; // TODO: free
        JetClustering {
            //ordered_jets: vec![],
            //jet_structure: vec![],
            ordered_pt: vec![],
            fastjet_jets_in: vec![],
            fastjet_jets_out: vec![],
            fastjet_jets_map: vec![],
            d_r,
            min_jpt,
            fastjet_workspace,
            use_fastjet,
        }
    }

    pub fn cluster_fastjet(&mut self, event: &Event) {
        self.fastjet_jets_in.clear();

        let mut len: c_int = 0;
        for (e, id) in event
            .kinematic_configuration
            .1
            .iter()
            .zip_eq(&event.final_state_particle_ids)
        {
            // filter for jet particles: u, d, c, s, d, g, QCD_ghost
            //TODO make it a hyperparam of the observable!
            if id.abs() < 6 || *id == 21 || id.abs() == 82 {
                self.fastjet_jets_in.extend(&[e.t, e.x, e.y, e.z]);
                len += 1;
            }
        }

        self.fastjet_jets_out.clear();
        self.fastjet_jets_out.resize(self.fastjet_jets_in.len(), 0.);

        self.fastjet_jets_map.clear();
        self.fastjet_jets_map.resize(self.fastjet_jets_in.len(), 0);

        let mut actual_len: c_int = 0;
        let palg = -1.0;
        let clustering_ptjet_min = 0.;

        if len > 0 {
            unsafe {
                fjcore::fastjetppgenkt_(
                    self.fastjet_workspace,
                    &self.fastjet_jets_in[0] as *const c_double,
                    &len as *const c_int,
                    &self.d_r as *const c_double,
                    &clustering_ptjet_min as *const c_double,
                    &palg as *const c_double,
                    &mut self.fastjet_jets_out[0] as *mut c_double,
                    &mut actual_len as *mut c_int,
                    &mut self.fastjet_jets_map[0] as *mut c_int,
                );
            }
        }

        self.ordered_pt.clear();
        for i in 0..actual_len as usize {
            let jet = LorentzVector::from_slice(&self.fastjet_jets_out[i * 4..(i + 1) * 4]);
            self.ordered_pt.push(jet.pt());
        }

        // Filter order_pt jets by removing those below the necessary min_jpt
        let min_pt = self.min_jpt;
        self.ordered_pt.retain(|&pt| pt >= min_pt);
    }

    pub fn cluster(&mut self, event: &Event) {
        if self.use_fastjet {
            self.cluster_fastjet(event)
        } else {
            self.cluster_custom(event)
        }
    }

    pub fn cluster_custom(&mut self, event: &Event) {
        //self.ordered_jets.clear();
        self.ordered_pt.clear();
        //self.jet_structure.clear();

        // group jets based on their dR
        let ext = &event.kinematic_configuration.1;

        let index_offset = 2;
        if ext.len() == 5 {
            let pt1 = ext[index_offset + 0].pt();
            let pt2 = ext[index_offset + 1].pt();
            let pt3 = ext[index_offset + 2].pt();
            let dr12 = ext[index_offset + 0].delta_r(&ext[index_offset + 1]);
            let dr13 = ext[index_offset + 0].delta_r(&ext[index_offset + 2]);
            let dr23 = ext[index_offset + 1].delta_r(&ext[index_offset + 2]);

            // we have at least 2 jets
            if dr12 < self.d_r {
                let m12 = ext[index_offset + 0] + ext[index_offset + 1];
                let pt12 = m12.pt();
                self.ordered_pt.push(pt12);
                self.ordered_pt.push(pt3);
            } else if dr13 < self.d_r {
                let m13 = ext[index_offset + 0] + ext[index_offset + 2];
                let pt13 = m13.pt();
                self.ordered_pt.push(pt13);
                self.ordered_pt.push(pt2);
            } else if dr23 < self.d_r {
                let m23 = ext[index_offset + 1] + ext[index_offset + 2];
                let pt23 = m23.pt();
                self.ordered_pt.push(pt23);
                self.ordered_pt.push(pt1);
            } else {
                self.ordered_pt.push(pt1);
                self.ordered_pt.push(pt2);
                self.ordered_pt.push(pt3);
            }

            // sort from large pt to small
            self.ordered_pt
                .sort_unstable_by(|a, b| b.partial_cmp(a).unwrap());
        } else {
            // treat every external momentum as its own jet for now
            for m in ext {
                self.ordered_pt.push(m.pt());
            }
            self.ordered_pt
                .sort_unstable_by(|a, b| b.partial_cmp(a).unwrap());
        }

        // Filter order_pt jets by removing those below the necessary min_jpt
        let min_pt = self.min_jpt;
        self.ordered_pt.retain(|&pt| pt >= min_pt);
    }
}

pub struct JetSelector {
    jet_selector_settings: JetSliceSettings,
    clustering: JetClustering,
}

impl JetSelector {
    pub fn new(settings: &JetSliceSettings) -> JetSelector {
        JetSelector {
            jet_selector_settings: settings.clone(),
            clustering: JetClustering::new(settings.use_fastjet, settings.dR, settings.min_jpt),
        }
    }
}

impl EventSelector for JetSelector {
    #[inline]
    fn process_event(&mut self, event: &mut Event) -> bool {
        self.clustering.cluster(event);
        //println!("Event: {:#?}",event);
        //println!("clustering: {:#?}",self.clustering);
        self.clustering.ordered_pt.len() >= self.jet_selector_settings.min_jets
            && self.clustering.ordered_pt.len() <= self.jet_selector_settings.max_jets
            && (self.clustering.ordered_pt.len() == 0
                || (self.clustering.ordered_pt[0] >= self.jet_selector_settings.min_j1pt
                    && (self.jet_selector_settings.max_j1pt < 0.0
                        || (self.clustering.ordered_pt[0] <= self.jet_selector_settings.max_j1pt))))
    }
}

#[derive(Debug, Clone)]
pub enum Observables {
    CrossSection(CrossSectionObservable),
    Jet1PT(Jet1PTObservable),
}

impl Observables {
    #[inline]
    pub fn process_event_group(&mut self, event: &[Event], integrator_weight: f64) {
        use self::Observables::*;
        match self {
            CrossSection(o) => o.process_event_group(event, integrator_weight),
            Jet1PT(o) => o.process_event_group(event, integrator_weight),
        }
    }

    #[inline]
    pub fn merge_samples(&mut self, other: &mut Observables) {
        use self::Observables::*;
        match (self, other) {
            (CrossSection(o1), CrossSection(o2)) => o1.merge_samples(o2),
            (Jet1PT(o1), Jet1PT(o2)) => o1.merge_samples(o2),
            (o1, o2) => panic!(
                "Cannot merge observables of different types: {:?} vs {:?}",
                o1, o2
            ),
        }
    }

    /// Produce the result (histogram, etc.) of the observable from all processed event groups.
    #[inline]
    pub fn update_result(&mut self, total_events: usize) {
        use self::Observables::*;
        match self {
            CrossSection(o) => o.update_result(total_events),
            Jet1PT(o) => o.update_result(total_events),
        }
    }
}

pub trait Observable {
    /// Process a group of events and return a new integrand value that will be returned to the integrator.
    fn process_event_group(&mut self, event: &[Event], integrator_weight: f64);

    fn merge_samples(&mut self, other: &mut Self)
    where
        Self: Sized;

    /// Produce the result (histogram, etc.) of the observable from all processed event groups.
    fn update_result(&mut self, total_events: usize);
}

#[derive(Debug, Clone)]
pub struct CrossSectionObservable {
    re: AverageAndErrorAccumulator,
    im: AverageAndErrorAccumulator,
    status_update_sender: StatusUpdateSender,
}

impl CrossSectionObservable {
    pub fn new(status_update_sender: StatusUpdateSender) -> CrossSectionObservable {
        CrossSectionObservable {
            re: AverageAndErrorAccumulator::default(),
            im: AverageAndErrorAccumulator::default(),
            status_update_sender,
        }
    }

    pub fn add_sample(&mut self, integrand: Complex<f64>, integrator_weight: f64) {
        self.re.add_sample(integrand.re * integrator_weight);
        self.im.add_sample(integrand.im * integrator_weight);
    }

    /// Give a live update on a copy of the statistics
    pub fn update_live_result(&self) {
        let mut re = self.re.clone();
        re.update_iter();
        let mut im = self.im.clone();
        im.update_iter();

        self.status_update_sender
            .send(StatusUpdate::NewPoint(
                re.cur_iter,
                re.avg,
                re.err,
                re.chi_sq / re.cur_iter as f64,
                im.avg,
                im.err,
                im.chi_sq / im.cur_iter as f64,
                true,
            ))
            .unwrap();
    }
}

impl Observable for CrossSectionObservable {
    fn process_event_group(&mut self, events: &[Event], integrator_weight: f64) {
        let mut integrand = Complex::<f64>::default();
        for e in events {
            integrand += e.integrand;
        }

        self.re.add_sample(integrand.re * integrator_weight);
        self.im.add_sample(integrand.im * integrator_weight);
    }

    fn merge_samples(&mut self, other: &mut CrossSectionObservable) {
        self.re.merge_samples(&mut other.re);
        self.im.merge_samples(&mut other.im);
    }

    fn update_result(&mut self, _total_events: usize) {
        self.re.update_iter();
        self.im.update_iter();

        self.status_update_sender
            .send(StatusUpdate::NewPoint(
                self.re.cur_iter,
                self.re.avg,
                self.re.err,
                self.re.chi_sq / self.re.cur_iter as f64,
                self.im.avg,
                self.im.err,
                self.im.chi_sq / self.im.cur_iter as f64,
                false,
            ))
            .unwrap();
    }
}

#[derive(Debug, Clone)]
pub struct Jet1PTObservable {
    x_min: f64,
    x_max: f64,
    jet_clustering: JetClustering,
    index_event_accumulator: Vec<(isize, f64)>,
    bins: Vec<AverageAndErrorAccumulator>,
    write_to_file: bool,
    filename: String,
    num_events: usize,
}

impl Jet1PTObservable {
    pub fn new(
        x_min: f64,
        x_max: f64,
        num_bins: usize,
        d_r: f64,
        min_jpt: f64,
        write_to_file: bool,
        filename: String,
        use_fastjet: bool,
    ) -> Jet1PTObservable {
        Jet1PTObservable {
            x_min,
            x_max,
            index_event_accumulator: Vec::with_capacity(20),
            bins: vec![AverageAndErrorAccumulator::default(); num_bins],
            jet_clustering: JetClustering::new(use_fastjet, d_r, min_jpt),
            write_to_file,
            filename,
            num_events: 0,
        }
    }
}

impl Observable for Jet1PTObservable {
    fn process_event_group(&mut self, events: &[Event], integrator_weight: f64) {
        // add events in a correlated manner such that cancellations are realized in the error
        self.index_event_accumulator.clear();
        for e in events {
            self.jet_clustering.cluster_fastjet(e);
            let max_pt = self.jet_clustering.ordered_pt[0];

            let index = ((max_pt - self.x_min) / (self.x_max - self.x_min) * self.bins.len() as f64)
                as isize;
            if index >= 0 && index < self.bins.len() as isize {
                let mut new = true;
                for i in self.index_event_accumulator.iter_mut() {
                    if i.0 == index {
                        *i = (index, i.1 + e.integrand.re * integrator_weight);
                        new = false;
                        break;
                    }
                }
                if new {
                    self.index_event_accumulator
                        .push((index, e.integrand.re * integrator_weight));
                }
            }
        }

        for i in &self.index_event_accumulator {
            self.bins[i.0 as usize].add_sample(i.1);
        }
    }

    fn merge_samples(&mut self, other: &mut Jet1PTObservable) {
        // TODO: check if bins are compatible?
        for (b, bo) in self.bins.iter_mut().zip_eq(&mut other.bins) {
            b.merge_samples(bo);
        }
    }

    fn update_result(&mut self, total_events: usize) {
        let diff = total_events - self.num_events;
        self.num_events = total_events;

        // rescale the entries in each bin by the number of samples
        // per bin over the complete number of events (including rejected) ones
        // this is not the same as recaling the average after recombining
        // with previous iterations
        for b in &mut self.bins {
            let scaling = b.num_samples as f64 / diff as f64;
            b.sum *= scaling;
            b.sum_sq *= scaling * scaling;
            b.update_iter();
        }

        if self.write_to_file {
            let mut f =
                BufWriter::new(File::create(&self.filename).expect("Could not create log file"));

            writeln!(f, "##& xmin & xmax & central value & dy &\n").unwrap();
            writeln!(
                f,
                "<histogram> {} \"j1 pT |X_AXIS@LIN |Y_AXIS@LOG |TYPE@LOALPHALOOP\"",
                self.bins.len()
            )
            .unwrap();

            for (i, b) in self.bins.iter().enumerate() {
                let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
                let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64
                    + self.x_min;

                writeln!(
                    f,
                    "  {:.8e}   {:.8e}   {:.8e}   {:.8e}",
                    c1, c2, b.avg, b.err,
                )
                .unwrap();
            }

            writeln!(f, "<\\histogram>").unwrap();
        } else {
            for (i, b) in self.bins.iter().enumerate() {
                let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
                let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64;
                println!("{}={}: {} +/ {}", c1, c2, b.avg, b.err,);
            }
        }
    }
}
