use dashboard::{StatusUpdate, StatusUpdateSender};
use itertools::Itertools;
use num::Complex;
use squared_topologies::CutkoskyCut;
use std::fs::File;
use std::io::{BufWriter, Write};
use vector::LorentzVector;
use {FloatLike, JetSliceSettings, ObservableMode, SelectorMode, Settings};

#[derive(Debug, Default, Clone)]
pub struct EventInfo {
    pub accepted_event_counter: usize,
    pub rejected_event_counter: usize,
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
            return;
        }

        let n = self.num_samples as f64;
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
                            settings.observables.Jet1PT.write_to_file,
                            settings.observables.Jet1PT.filename.clone(),
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

        // TODO: recycle vectors from the buffer
        let mut incoming_momenta = Vec::with_capacity(orig_incoming_momenta.len());

        // rotate all momenta with the inverse of the rotation matrix
        for e in orig_incoming_momenta {
            incoming_momenta.push(LorentzVector::from_args(
                e.t,
                rot_matrix[0][0] * e.x + rot_matrix[1][0] * e.y + rot_matrix[2][0] * e.z,
                rot_matrix[0][1] * e.x + rot_matrix[1][1] * e.y + rot_matrix[2][1] * e.z,
                rot_matrix[0][2] * e.x + rot_matrix[1][2] * e.y + rot_matrix[2][2] * e.z,
            ));
        }

        let mut outgoing_momenta = Vec::with_capacity(cut_info.len());
        for (cut_mom, cut) in cut_momenta[..cut_info.len()].iter().zip(cut_info.iter()) {
            if cut.level == 0 {
                // make sure all momenta are outgoing
                let e = cut_mom.cast::<f64>() * cut.sign as f64;

                outgoing_momenta.push(LorentzVector::from_args(
                    e.t,
                    rot_matrix[0][0] * e.x + rot_matrix[1][0] * e.y + rot_matrix[2][0] * e.z,
                    rot_matrix[0][1] * e.x + rot_matrix[1][1] * e.y + rot_matrix[2][1] * e.z,
                    rot_matrix[0][2] * e.x + rot_matrix[1][2] * e.y + rot_matrix[2][2] * e.z,
                ));
            }
        }

        let mut e = Event {
            kinematic_configuration: (incoming_momenta, outgoing_momenta),
            integrand: Complex::new(1., 0.),
            weights: vec![0.],
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
        other.accepted_event_counter = 0;
        other.rejected_event_counter = 0;
    }

    pub fn update_result(&mut self) {
        for o in &mut self.observables {
            o.update_result();
        }

        self.status_update_sender
            .as_mut()
            .unwrap()
            .send(StatusUpdate::EventInfo(EventInfo {
                accepted_event_counter: self.accepted_event_counter,
                rejected_event_counter: self.rejected_event_counter,
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
    }
}

#[derive(Default, Debug, Clone)]
pub struct Event {
    pub kinematic_configuration: (Vec<LorentzVector<f64>>, Vec<LorentzVector<f64>>),
    pub integrand: Complex<f64>,
    pub weights: Vec<f64>,
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
}

impl JetClustering {
    pub fn new(d_r: f64) -> JetClustering {
        JetClustering {
            //ordered_jets: vec![],
            //jet_structure: vec![],
            ordered_pt: vec![],
            d_r,
        }
    }

    pub fn cluster(&mut self, event: &Event) {
        //self.ordered_jets.clear();
        self.ordered_pt.clear();
        //self.jet_structure.clear();

        // group jets based on their dR
        let ext = &event.kinematic_configuration.1;

        if ext.len() == 3 {
            let pt1 = (ext[0].x * ext[0].x + ext[0].y * ext[0].y).sqrt();
            let pt2 = (ext[1].x * ext[1].x + ext[1].y * ext[1].y).sqrt();
            let pt3 = (ext[2].x * ext[2].x + ext[2].y * ext[2].y).sqrt();
            let dr12 = ext[0].deltaR(&ext[1]);
            let dr13 = ext[0].deltaR(&ext[2]);
            let dr23 = ext[1].deltaR(&ext[2]);

            // we have at least 2 jets
            if dr12 < self.d_r {
                let m12 = ext[0] + ext[1];
                let pt12 = (m12.x * m12.x + m12.y * m12.y).sqrt();
                self.ordered_pt.push(pt12);
                self.ordered_pt.push(pt3);
            } else if dr13 < self.d_r {
                let m13 = ext[0] + ext[2];
                let pt13 = (m13.x * m13.x + m13.y * m13.y).sqrt();
                self.ordered_pt.push(pt13);
                self.ordered_pt.push(pt2);
            } else if dr23 < self.d_r {
                let m23 = ext[1] + ext[2];
                let pt23 = (m23.x * m23.x + m23.y * m23.y).sqrt();
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
                let pt = (m.x * m.x + m.y * m.y).sqrt();
                self.ordered_pt.push(pt);
            }
            self.ordered_pt
                .sort_unstable_by(|a, b| b.partial_cmp(a).unwrap());
        }
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
            clustering: JetClustering::new(settings.dR),
        }
    }
}

impl EventSelector for JetSelector {
    #[inline]
    fn process_event(&mut self, event: &mut Event) -> bool {
        self.clustering.cluster(event);
        self.clustering.ordered_pt.len() >= self.jet_selector_settings.min_jets
            && self.clustering.ordered_pt.len() <= self.jet_selector_settings.max_jets
            && self.clustering.ordered_pt[0] >= self.jet_selector_settings.min_j1pt
            && self.clustering.ordered_pt[0] <= self.jet_selector_settings.max_j1pt
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
    pub fn update_result(&mut self) {
        use self::Observables::*;
        match self {
            CrossSection(o) => o.update_result(),
            Jet1PT(o) => o.update_result(),
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
    fn update_result(&mut self);
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

    fn update_result(&mut self) {
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
            ))
            .unwrap();
    }
}

#[derive(Debug, Clone)]
pub struct Jet1PTObservable {
    x_min: f64,
    x_max: f64,
    jet_clustering: JetClustering,
    bins: Vec<AverageAndErrorAccumulator>,
    write_to_file: bool,
    filename: String,
}

impl Jet1PTObservable {
    pub fn new(
        x_min: f64,
        x_max: f64,
        num_bins: usize,
        d_r: f64,
        write_to_file: bool,
        filename: String,
    ) -> Jet1PTObservable {
        Jet1PTObservable {
            x_min,
            x_max,
            bins: vec![AverageAndErrorAccumulator::default(); num_bins],
            jet_clustering: JetClustering::new(d_r),
            write_to_file,
            filename,
        }
    }
}

impl Observable for Jet1PTObservable {
    fn process_event_group(&mut self, events: &[Event], integrator_weight: f64) {
        for e in events {
            self.jet_clustering.cluster(e);
            let max_pt = self.jet_clustering.ordered_pt[0];

            // insert into histogram
            let index = ((max_pt - self.x_min) / (self.x_max - self.x_min) * self.bins.len() as f64)
                as isize;
            if index >= 0 && index < self.bins.len() as isize {
                self.bins[index as usize].add_sample(e.integrand.re * integrator_weight);
            }
        }
    }

    fn merge_samples(&mut self, other: &mut Jet1PTObservable) {
        // TODO: check if bins are compatible?
        for (b, bo) in self.bins.iter_mut().zip_eq(&mut other.bins) {
            b.merge_samples(bo);
        }
    }

    fn update_result(&mut self) {
        for b in &mut self.bins {
            b.update_iter();
        }

        if self.write_to_file {
            let mut f =
                BufWriter::new(File::create(&self.filename).expect("Could not create log file"));

            writeln!(f, "##& xmin & xmax & central value & dy &\n").unwrap();
            writeln!(
                f,
                "<histogram> {} \"j1 pT |X_AXIS@LIN |Y_AXIS@LOG |TYPE@NLOALPHALOOP\"",
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
                    c1, c2, b.avg, b.err
                )
                .unwrap();
            }

            writeln!(f, "<\\histogram>").unwrap();
        } else {
            for (i, b) in self.bins.iter().enumerate() {
                let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
                let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64;
                println!("{}={}: {} +/ {}", c1, c2, b.avg, b.err);
            }
        }
    }
}
