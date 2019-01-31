use vector::LorentzVector;
use ltd::LoopLine;

pub fn build_loop_line(
    loop_momenta: &[(usize, bool)],
    ext: Vec<LorentzVector<f64>>,
    mass: Vec<f64>,
    compute_qs: bool,
) -> LoopLine {
    let qs = if compute_qs {
        // convention of Dario
        let mut qs = vec![ext[0]];
        if ext.len() > 1 {
            for e in &ext[1..ext.len() - 1] {
                let n = qs.last().unwrap() + e;
                qs.push(n);
            }
        }
        qs.push(LorentzVector::default());
        qs
    } else {
        // externals are already summed
        ext
    };

    let q_and_mass: Vec<_> = qs.into_iter().zip(mass.into_iter()).collect();

    LoopLine::new(loop_momenta, q_and_mass)
}

pub fn create_topology(topology: &str) -> (usize, f64, Vec<LoopLine>) {
    match topology {
        "P1" => {
            // does not need deformation
            let p1 = LorentzVector::from_args(5.23923, -4.18858, 0.74966, -3.05669);
            let p2 = LorentzVector::from_args(6.99881, -2.93659, 5.03338, 3.87619);
            let m = 7.73358;
            let e_cm = (p1 + p2).square().abs().sqrt();
            let ll = build_loop_line(&[(0, true)], vec![p1, p2, -p1 - p2], vec![m, m, m], true);
            (1, e_cm, vec![ll])
        }
        "P5" => {
            // does not need deformation
            let p1 = LorentzVector::from_args(31.54872, -322.40325, 300.53015, -385.58013);
            let p2 = LorentzVector::from_args(103.90430, 202.00974, -451.27794, -435.12848);
            let p3 = LorentzVector::from_args(294.76653, 252.88958, 447.09194, 311.71630);
            let m = 4.68481;
            let e_cm = (p1 + p2).square().abs().sqrt();
            let ll = build_loop_line(
                &[(0, true)],
                vec![p1, p2, p3, -p1 - p2 - p3],
                vec![m; 4],
                true,
            );
            (1, e_cm, vec![ll])
        }
        "P3" => {
            //P3 in https://arxiv.org/pdf/1510.00187.pdf
            let p1 = LorentzVector::from_args(10.51284, 6.89159, -7.40660, -2.85795);
            let p2 = LorentzVector::from_args(6.45709, 2.46635, 5.84093, 1.22257);
            let m = 0.52559;
            let e_cm = (p1 + p2).square().abs().sqrt();
            let ll = build_loop_line(&[(0, true)], vec![p1, p2, -p1 - p2], vec![m; 3], true);
            (1, e_cm, vec![ll])
        }
        "P4" => {
            //P4 in https://arxiv.org/pdf/1510.00187.pdf
            let p1 = LorentzVector::from_args(95.77004, 31.32025, -34.08106, -9.38565);
            let p2 = LorentzVector::from_args(94.54738, -53.84229, 67.11107, 45.56763);
            let m1 = 83.02643;
            let m2 = 76.12873;
            let m3 = 55.00359;
            let e_cm = (p1 + p2).square().abs().sqrt();
            let ll = build_loop_line(&[(0, true)], vec![p1, p2, -p1 - p2], vec![m1, m2, m3], true);
            (1, e_cm, vec![ll])
        }
        "P7" => {
            let p1 = LorentzVector::from_args(62.80274, -49.71968, -5.53340, -79.44048);
            let p2 = LorentzVector::from_args(48.59375, -1.65847, 34.91140, 71.89564);
            let p3 = LorentzVector::from_args(76.75934, -19.14334, -17.10279, 30.22959);
            let m = 9.82998;

            let e_cm = (p1 + p2).square().abs().sqrt();
            let ll = build_loop_line(
                &[(0, true)],
                vec![p1, p2, p3, -p1 - p2 - p3],
                vec![m, m, m, m],
                true,
            );
            (1, e_cm, vec![ll])
        }
        "double-triangle" => {
            let p = LorentzVector::from_args(1.0, 0., 0., 0.);
            let m = 0.;
            let e_cm = p.square().abs().sqrt();

            let ll1 = build_loop_line(
                &[(0, true)],
                vec![LorentzVector::default(), -p],
                vec![m, m],
                false,
            );
            let ll2 = build_loop_line(
                &[(1, true)],
                vec![LorentzVector::default(), p],
                vec![m, m],
                false,
            );
            let ll3 = build_loop_line(
                &[(0, true), (1, true)],
                vec![LorentzVector::default()],
                vec![m],
                false,
            );
            (2, e_cm, vec![ll1, ll2, ll3])
        }
        "triangle-box" => {
            let p1 = LorentzVector::from_args(1.0, 0.0, 0.0, 0.0);
            let p2 = LorentzVector::from_args(-0.59375, -0.0, -0.0, -0.32021721143623744);
            let _p3 = LorentzVector::from_args(-0.40625, -0.0, -0.0, 0.32021721143623744);
            let m = 0.;
            let e_cm = 1.0;

            let ll1 = build_loop_line(
                &[(0, true)],
                vec![LorentzVector::default(), -p1],
                vec![m, m],
                false,
            );
            let ll2 = build_loop_line(
                &[(1, true)],
                vec![LorentzVector::default(), -p2, p1],
                vec![m, m, m],
                false,
            );
            let ll3 = build_loop_line(
                &[(0, true), (1, true)],
                vec![-p1],
                vec![m],
                false,
            );
            (2, e_cm, vec![ll1, ll2, ll3])
        }
        x => unimplemented!("Unknown topology {}", x),
    }
}
