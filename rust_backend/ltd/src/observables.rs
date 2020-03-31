use num::Complex;
use vector::LorentzVector;

#[derive(Default, Debug, Clone)]
pub struct Event {
    pub kinematic_configuration: (Vec<LorentzVector<f64>>, Vec<LorentzVector<f64>>),
    pub weights: Vec<f64>,
}

pub trait Observable {
    /// Process a group of events and return a new integrand value that will be returned to the integrator.
    fn process_event_group(
        &mut self,
        event: &mut Vec<Event>,
        integrator_weight: f64,
        integrand: Complex<f64>,
    ) -> Complex<f64>;

    /// Produce the result (histogram, etc.) of the observable from all processed event groups.
    fn generate_result(&mut self);
}

#[derive(Default)]
pub struct NoObservable {}

impl Observable for NoObservable {
    #[inline]
    fn process_event_group(
        &mut self,
        events: &mut Vec<Event>,
        integrator_weight: f64,
        integrand: Complex<f64>,
    ) -> Complex<f64> {
        events.clear();
        integrand
    }

    #[inline]
    fn generate_result(&mut self) {}
}

#[derive(Default)]
pub struct Jet1PTObservable {}

impl Observable for Jet1PTObservable {
    #[inline]
    fn process_event_group(
        &mut self,
        events: &mut Vec<Event>,
        integrator_weight: f64,
        integrand: Complex<f64>,
    ) -> Complex<f64> {
        events.clear();
        // TODO
        integrand
    }

    #[inline]
    fn generate_result(&mut self) {}
}
