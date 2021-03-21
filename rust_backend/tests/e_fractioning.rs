use lorentz_vector::LorentzVector;
use ltd::partial_fractioning::{PartialFractioning, PartialFractioningMultiLoops};
use ltd::topologies::{LoopLine, Propagators};

// Testing the one-loop version of E-fractioning
mod partial_fractioning {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_partial_fractioning() {
        let start = std::time::Instant::now();
        let pf_expr = PartialFractioning::new(10, 9);

        println!("{:?}", start.elapsed());
        println!(
            "{} -> size: {}",
            10,
            pf_expr.partial_fractioning_element.len()
        );
        assert_eq!(pf_expr.partial_fractioning_element.len(), 92378);
        //for expr in pf_expr.partial_fractioning_element.iter() {
        //    println!("{:?} :: num({:?})", expr.ellipsoids_product, expr.numerator);
        //}
    }
}

// Testing the multi loop version of E-fractioning
mod partial_fractioning_multi_loops {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_pf_multi_loops_decagon() {
        //Dummy prop
        let propagator = Propagators {
            name: "dummy".to_owned(),
            uv: false,
            id: 0, // global id
            m_squared: 0.0,
            q: LorentzVector::default(),
            parametric_shift: (Vec::new(), Vec::new()),
            signature: Vec::new(),
            power: 1,
        };
        // Replicate one loop looplines
        let mut looplines_1loop = Vec::new();
        looplines_1loop.push(LoopLine {
            start_node: 0,
            end_node: 0,
            signature: vec![1],
            propagators: vec![propagator.clone(); 10],
        });
        // TEST 1 LOOP
        let start = std::time::Instant::now();
        let pf_expr = PartialFractioningMultiLoops::new(&looplines_1loop, 9);
        print!("1LOOP -> #{}", pf_expr.partial_fractioning_element.len());
        println!(" :: {:?}", start.elapsed());
        assert_eq!(pf_expr.partial_fractioning_element.len(), 92378);
    }

    #[test]
    fn test_pf_multi_loops_2loop() {
        //Dummy prop
        let propagator = Propagators {
            name: "dummy".to_owned(),
            uv: false,
            id: 0, // global id
            m_squared: 0.0,
            q: LorentzVector::default(),
            parametric_shift: (Vec::new(), Vec::new()),
            signature: Vec::new(),
            power: 1,
        };
        // Replicate two loops looplines
        let signatures = [vec![1, 0], vec![1, -1], vec![0, 1]];
        let multiplicity = [6, 1, 4];
        let mut looplines_2loop = Vec::new();
        for (s, m) in signatures.iter().zip(multiplicity.iter()) {
            looplines_2loop.push(LoopLine {
                start_node: 0,
                end_node: 0,
                signature: s.clone(),
                propagators: vec![propagator.clone(); m.clone()],
            });
        }
        // TEST 2 LOOPS
        let start = std::time::Instant::now();
        let pf_expr = PartialFractioningMultiLoops::new(&looplines_2loop, 0);
        print!("2LOOPS -> #{}", pf_expr.partial_fractioning_element.len());
        println!(" :: {:?}", start.elapsed());
        assert_eq!(pf_expr.partial_fractioning_element.len(), 48620);
    }

    #[test]
    fn test_pf_multi_loops_4loop() {
        //Dummy prop
        let propagator = Propagators {
            name: "dummy".to_owned(),
            uv: false,
            id: 0, // global id
            m_squared: 0.0,
            q: LorentzVector::default(),
            parametric_shift: (Vec::new(), Vec::new()),
            signature: Vec::new(),
            power: 1,
        };

        // Replicate 2x2 FISHNET looplines
        let signatures = [
            vec![1, 0, 0, 0],
            vec![0, 1, 0, 0],
            vec![1, -1, 0, 0],
            vec![-1, 0, -1, 0],
            vec![0, -1, 0, -1],
            vec![0, 0, 1, 0],
            vec![0, 0, -1, 1],
            vec![0, 0, 0, -1],
        ];
        let multiplicity = [1, 1, 1, 1, 1, 1, 1, 1];
        let mut looplines_2x2_fishnet = Vec::new();
        for (s, m) in signatures.iter().zip(multiplicity.iter()) {
            looplines_2x2_fishnet.push(LoopLine {
                start_node: 0,
                end_node: 0,
                signature: s.clone(),
                propagators: vec![propagator.clone(); m.clone()],
            });
        }
        // TEST 2x2 FISHNET
        let start = std::time::Instant::now();
        let pf_expr = PartialFractioningMultiLoops::new(&looplines_2x2_fishnet, 0);
        print!("FISHNET - >#{}", pf_expr.partial_fractioning_element.len());
        println!(" :: {:?}", start.elapsed());
        assert_eq!(pf_expr.partial_fractioning_element.len(), 98);
    }
}
