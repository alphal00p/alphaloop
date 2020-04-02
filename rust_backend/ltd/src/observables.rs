use num::Complex;
use std::fs::File;
use std::io::{BufWriter, Write};
use vector::LorentzVector;

#[derive(Default, Clone)]
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

#[derive(Default, Debug, Clone)]
pub struct Event {
    pub kinematic_configuration: (Vec<LorentzVector<f64>>, Vec<LorentzVector<f64>>),
    pub integrand: Complex<f64>,
    pub weights: Vec<f64>,
}

pub trait EventFilter {
    /// Process a group of events and return a new integrand value that will be returned to the integrator.
    fn process_event_group(
        &mut self,
        event: &mut Vec<Event>,
        integrator_weight: f64,
    ) -> Complex<f64>;
}

#[derive(Default)]
pub struct NoEventFilter {}

impl EventFilter for NoEventFilter {
    #[inline]
    fn process_event_group(
        &mut self,
        events: &mut Vec<Event>,
        _integrator_weight: f64,
    ) -> Complex<f64> {
        let mut integrand = Complex::default();
        for e in events {
            integrand += e.integrand;
        }
        integrand
    }
}

pub trait Observable {
    /// Process a group of events and return a new integrand value that will be returned to the integrator.
    fn process_event_group(&mut self, event: &[Event], integrator_weight: f64);

    /// Produce the result (histogram, etc.) of the observable from all processed event groups.
    fn update_result(&mut self);
}

#[derive(Default)]
pub struct CrossSectionObservable {
    re: AverageAndErrorAccumulator,
    im: AverageAndErrorAccumulator,
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

    fn update_result(&mut self) {
        self.re.update_iter();
        self.im.update_iter();

        println!("Iteration {}", self.re.cur_iter);
        println!(
            " re: {} +- {} chisq {}",
            self.re.avg, self.re.err, self.re.chi_sq
        );
        println!(
            " im: {} +- {} chisq {}",
            self.im.avg, self.im.err, self.im.chi_sq
        );
    }
}

#[derive(Default)]
pub struct Jet1PTObservable {
    x_min: f64,
    x_max: f64,
    bins: Vec<AverageAndErrorAccumulator>,
    d_r: f64,
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
            d_r,
            write_to_file,
            filename,
        }
    }
}

impl Observable for Jet1PTObservable {
    fn process_event_group(&mut self, events: &[Event], integrator_weight: f64) {
        for e in events {
            let mut max_pt = 0.;
            let mut three_jets = false;

            // group jets based on their dR
            let ext = &e.kinematic_configuration.1;
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
                    max_pt = pt12.max(pt3);
                } else if dr13 < self.d_r {
                    let m13 = ext[0] + ext[2];
                    let pt13 = (m13.x * m13.x + m13.y * m13.y).sqrt();
                    max_pt = pt13.max(pt2);
                } else if dr23 < self.d_r {
                    let m23 = ext[1] + ext[2];
                    let pt23 = (m23.x * m23.x + m23.y * m23.y).sqrt();
                    max_pt = pt23.max(pt1);
                } else {
                    max_pt = pt1.max(pt2).max(pt3);

                    //max_pt = pt1.min(pt2).min(pt3);
                    //three_jets = true;
                }
            } else {
                // treat every external momentum as its own jet
                for m in ext {
                    let pt = (m.x * m.x + m.y * m.y).sqrt();
                    if pt > max_pt {
                        max_pt = pt;
                    }
                }
            }

            if !three_jets {
                //continue;
            }

            // insert into histogram
            let index = ((max_pt - self.x_min) / (self.x_max - self.x_min) * self.bins.len() as f64)
                as isize;
            if index >= 0 && index < self.bins.len() as isize {
                self.bins[index as usize].add_sample(e.integrand.re * integrator_weight);
            }
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
