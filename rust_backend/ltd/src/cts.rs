use topologies::Topology;
use vector::LorentzVector;
type Complex = num::Complex<f64>;

impl Topology {
    pub fn counterterm(&self) -> Complex {
        match self.n_loops {
            1 => {
                match (self.on_shell_flag, self.loop_lines[0].propagators.len()) {
                    // 4-point 1-loop with on-shell externals
                    (0b1111, 4) => {
                        // TODO
                        Complex::new(0., 0.)
                    }
                    _ => Complex::new(0., 0.),
                }
            }
            _ => Complex::new(0., 0.),
        }
    }
}
