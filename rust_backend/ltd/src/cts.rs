use topologies::Topology;
use vector::LorentzVector;
use Complex;

impl Topology {
    pub fn counterterm(&self) -> Complex {
        match self.n_loops {
            1 => {
                match (self.on_shell_flag, self.loop_lines[0].propagators.len()) {
                    // 4-point 1-loop with on-shell externals
                    (0b1111, 4) => {
                        // TODO
                        Complex::default()
                    }
                    _ => Complex::default(),
                }
            }
            _ => Complex::default(),
        }
    }
}
