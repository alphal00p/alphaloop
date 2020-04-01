use num::Complex;
use std::fs::File;
use std::io::{BufWriter, Write};
use vector::LorentzVector;

#[derive(Default, Debug, Clone)]
pub struct Event {
    pub kinematic_configuration: (Vec<LorentzVector<f64>>, Vec<LorentzVector<f64>>),
    pub weights: Vec<f64>,
}

pub trait EventFilter {
    /// Process a group of events and return a new integrand value that will be returned to the integrator.
    fn process_event_group(
        &mut self,
        event: &mut Vec<Event>,
        integrator_weight: f64,
        integrand: Complex<f64>,
    ) -> Complex<f64>;
}

#[derive(Default)]
pub struct NoEventFilter {}

impl EventFilter for NoEventFilter {
    #[inline]
    fn process_event_group(
        &mut self,
        _events: &mut Vec<Event>,
        _integrator_weight: f64,
        integrand: Complex<f64>,
    ) -> Complex<f64> {
        integrand
    }
}

pub trait Observable {
    /// Process a group of events and return a new integrand value that will be returned to the integrator.
    fn process_event_group(
        &mut self,
        event: &[Event],
        integrator_weight: f64,
        integrand: Complex<f64>,
    );

    /// Produce the result (histogram, etc.) of the observable from all processed event groups.
    fn update_result(&mut self);
}

#[derive(Default)]
pub struct Jet1PTObservable {
    x_min: f64,
    x_max: f64,
    bins: Vec<f64>,
    bins_sq: Vec<f64>,
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
            bins: vec![0.; num_bins],
            bins_sq: vec![0.; num_bins],
            d_r,
            write_to_file,
            filename,
        }
    }
}

impl Observable for Jet1PTObservable {
    fn process_event_group(
        &mut self,
        events: &[Event],
        integrator_weight: f64,
        integrand: Complex<f64>,
    ) {
        for e in events {
            let mut max_pt = 0.;

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

            // insert into histogram
            let index = ((max_pt - self.x_min) / (self.x_max - self.x_min) * self.bins.len() as f64)
                as isize;
            if index >= 0 && index < self.bins.len() as isize {
                let r = integrand.re * integrator_weight;
                self.bins[index as usize] += r;
                self.bins_sq[index as usize] += r * r;
            }
        }
    }

    fn update_result(&mut self) {
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

            for (i, (b, bsq)) in self.bins.iter().zip(self.bins_sq.iter()).enumerate() {
                let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
                let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64
                    + self.x_min;

                let n = 100000.; // FIXME: this varies
                let w = (bsq * n).sqrt();
                let err = ((w + b) * (w - b) / (n - 1.)).sqrt();
                writeln!(f, "  {:.8e}   {:.8e}   {:.8e}   {:.8e}", c1, c2, b, err).unwrap();
            }

            writeln!(f, "<\\histogram>").unwrap();
        } else {
            println!("{:?}", self.bins);
        }

        // TODO: merge results with the results from previous iterations
        for (b, bsq) in self.bins.iter_mut().zip(self.bins_sq.iter_mut()) {
            *b = 0.;
            *bsq = 0.;
        }
    }
}
