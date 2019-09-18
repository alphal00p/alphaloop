/// A propagator or external momentum
#[derive(Debug, Clone)]
struct Edge {
    index: usize, // index in the edge list, for easy access
    momentum_map: Vec<(usize, bool)>,
    vertices: (usize, usize), // directed edge map (a,b) from vertex a to vertex b
    mass: f64,
}

impl Edge {
    fn from_edge(index: usize, v1: usize, v2: usize) -> Edge {
        Edge {
            index,
            momentum_map: vec![],
            vertices: (v1, v2),
            mass: 0.,
        }
    }

    fn from_edge_list(edges: &[(usize, usize)]) -> Vec<Edge> {
        edges
            .iter()
            .enumerate()
            .map(|(i, (v1, v2))| Edge {
                index: i,
                momentum_map: vec![],
                vertices: (*v1, *v2),
                mass: 0.,
            })
            .collect()
    }
}

#[derive(Debug, Clone)]
struct Topology {
    num_vertices: usize,
    externals: Vec<usize>, // edge index of external
    external_values: Vec<f64>,
    loop_momenta: Vec<usize>,
    edges: Vec<Edge>, // Edges parameterized with loop momenta
    momentum_map: Option<Vec<Vec<(usize, bool)>>>,
}

impl Topology {
    pub fn new(
        num_vertices: usize,
        edges: Vec<Edge>,
        loop_momenta: Vec<usize>,
        external_values: Vec<f64>,
    ) -> Topology {
        let mut vertex_count = vec![0; num_vertices];

        for x in &edges {
            vertex_count[x.vertices.0] += 1;
            vertex_count[x.vertices.1] += 1;
        }

        let externals = edges
            .iter()
            .enumerate()
            .filter_map(|(i, x)| {
                if vertex_count[x.vertices.0] == 1 || vertex_count[x.vertices.1] == 1 {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();

        Topology {
            num_vertices,
            edges,
            loop_momenta,
            externals,
            external_values,
            momentum_map: None,
        }
    }

    /// Find all paths from edge with index `start` to edge with index `dest`.
    /// Returns a path that includes a direction boolean.
    pub fn find_paths(&self, start: usize, dest: usize) -> Vec<Vec<(usize, bool)>> {
        assert!(start < self.edges.len() && dest < self.edges.len());

        let is_loop = start == dest;
        let start_check = if is_loop { 1 } else { 0 };
        let mut paths = vec![vec![(start, true)]];
        if !is_loop {
            paths.push(vec![(start, false)]);
        }

        let mut res = vec![];
        loop {
            let mut newpaths = vec![];
            for p in paths {
                if p.len() > 1 && p.last().unwrap().0 == dest {
                    res.push(p);
                    continue;
                }
                let last_vertex = p
                    .last()
                    .map(|&(v, dir)| {
                        if dir {
                            self.edges[v].vertices.1
                        } else {
                            self.edges[v].vertices.0
                        }
                    })
                    .unwrap();
                for x in self.edges.iter() {
                    if p.iter().skip(start_check).all(|pp| x.index != pp.0) {
                        if is_loop && x.index == start {
                            // if we have a loop, we need to enter from the right direction
                            if x.vertices.0 == last_vertex {
                                let mut n = p.clone();
                                n.push((x.index, true));
                                newpaths.push(n);
                            }
                            continue;
                        }

                        // vertex check
                        let mut has_vertex = [false, false];
                        if p.len() > start_check {
                            for &(prop, _) in &p[start_check..p.len() - 1] {
                                has_vertex[0] |= self.edges[prop].vertices.0 == x.vertices.0
                                    || self.edges[prop].vertices.1 == x.vertices.0;
                                has_vertex[1] |= self.edges[prop].vertices.0 == x.vertices.1
                                    || self.edges[prop].vertices.1 == x.vertices.1;
                            }
                        }

                        if x.vertices.0 == last_vertex && !has_vertex[0] {
                            let mut n = p.clone();
                            n.push((x.index, true));
                            newpaths.push(n);
                        }
                        if x.vertices.1 == last_vertex && !has_vertex[1] {
                            let mut n = p.clone();
                            n.push((x.index, false));
                            newpaths.push(n);
                        }
                    }
                }
            }
            paths = newpaths;
            if paths.is_empty() {
                break;
            }
        }
        res
    }

    ///Generate momentum flow where `loop_momenta` are a
    ///list of all edge indices that go like 1/k^2
    fn generate_momentum_flow(&mut self) {
        // first we fix the loop momentum flows
        let mut flows = vec![];
        for &l in &self.loop_momenta {
            let mut paths = self.find_paths(l, l);
            // make sure other loop momenta propagators are not in the path
            // TODO: prevent clone
            paths = paths
                .iter()
                .filter(|x| {
                    x[1..x.len() - 1]
                        .iter()
                        .all(|y| !self.loop_momenta.contains(&y.0))
                })
                .cloned()
                .collect();

            // TODO: take shortest?
            let mut e = paths.swap_remove(0);
            e.pop();
            flows.push(e)
        }

        // now route the external loop_momenta to the sink
        let mut ext_flows = vec![];
        for &e in &self.externals[..self.externals.len() - 1] {
            let mut paths = self.find_paths(e, *self.externals.last().unwrap());
            paths = paths
                .into_iter()
                .filter(|x| {
                    x[1..x.len() - 1]
                        .iter()
                        .all(|y| !self.loop_momenta.contains(&y.0))
                })
                .map(|mut x| {
                    // if the external momentum is outgoing, we need to flip all the signs
                    if !x[0].1 {
                        for y in &mut x {
                            y.1 = !y.1;
                        }
                    }
                    x
                })
                .collect();
            assert!(paths.len() > 0);
            ext_flows.push(paths.swap_remove(0))
        }

        // propagator momenta
        let mut mom = vec![];
        let mut s = vec!["1".to_owned()];
        for (i, x) in self.edges.iter().enumerate() {
            if self.externals.contains(&i) {
                mom.push(vec![]);
                continue;
            }
            if self.loop_momenta.contains(&i) {
                mom.push(vec![(i, true)]);
                s.push(format!("/k{}^2", i));
            } else {
                let mut newmom = vec![];
                let mut s1 = "/(".to_owned();
                for (j, y) in flows.iter().enumerate() {
                    for yy in y {
                        if yy.0 == i {
                            s1 += &format!(
                                "{}k{}",
                                if yy.1 { "+" } else { "-" },
                                self.loop_momenta[j]
                            );
                            newmom.push((self.loop_momenta[j], yy.1));
                            break;
                        }
                    }
                }
                for (j, y) in ext_flows.iter().enumerate() {
                    for yy in y {
                        if yy.0 == i {
                            s1 +=
                                &format!("{}p{}", if yy.1 { "+" } else { "-" }, self.externals[j]);
                            newmom.push((self.externals[j], yy.1));
                            break;
                        }
                    }
                }
                mom.push(newmom);
                s.push(s1 + ")^2");
            }
        }

        println!("Mom map: {:?}", mom);
        println!("Propagator structure: {}", s.join(""));
        self.momentum_map = Some(mom);
    }

    /// Evaluate propagators
    pub fn evaluate_propagators(&self, loop_momenta_values: &[f64]) -> f64 {
        let mut res = 1.;
        for (i, x) in self.momentum_map.as_ref().unwrap().iter().enumerate() {
            if x.is_empty() {
                continue;
            }
            res *= self.evaluate(i, loop_momenta_values).powi(2);
        }
        1. / res
    }

    /// Evaluate a propagator momentum
    fn evaluate(&self, index: usize, loop_momenta_values: &[f64]) -> f64 {
        let mom_map = self.momentum_map.as_ref().unwrap();

        let mut r = 0.;
        for (y, sign) in &mom_map[index] {
            let val = match self.loop_momenta.iter().position(|a| a == y) {
                Some(x) => loop_momenta_values[x],
                None => self.external_values[self.externals.iter().position(|a| a == y).unwrap()],
            };
            r += if *sign { 1. } else { -1. } * val;
        }
        r
    }

    /// Construct the cycles used for Weinzierls' deformation
    pub fn construct_cycles_per_loop_momentum(
        &self,
        index: usize,
        loop_momenta_values: &[f64],
    ) -> Vec<Vec<f64>> {
        let mut flows = vec![];
        let paths = self.find_paths(self.loop_momenta[index], self.loop_momenta[index]);

        let mom = self.evaluate(self.loop_momenta[index], loop_momenta_values);

        for mut path in paths {
            path.pop();
            let mut path_eval = vec![];
            for (link, sign) in path {
                // pick up a minus sign if the direction is opposite
                let s = if sign { 1. } else { -1. };
                path_eval.push(s * self.evaluate(link, loop_momenta_values) - mom);
            }
            flows.push(path_eval);
        }

        flows
    }
}

#[test]
fn test_oneloop() {
    let edges = Edge::from_edge_list(&[(0, 1), (1, 2), (1, 3), (2, 3), (3, 4), (2, 5)]);
    let mut t = Topology::new(6, edges, vec![1], vec![1.0, 2.0, 3.0]);
    println!("{:?}", t.find_paths(1, 4));
    println!("{:?}", t.generate_momentum_flow());
    println!("{}", t.evaluate_propagators(&[5.0]));
    println!("{:?}", t.construct_cycles_per_loop_momentum(0, &[5.0]));
    assert!(false);
}

#[test]
fn test_doubletriangle() {
    let edges = Edge::from_edge_list(&[(1, 0), (4, 1), (2, 3), (1, 2), (2, 4), (3, 4), (3, 5)]);
    let mut t = Topology::new(6, edges, vec![1, 2], vec![1.0, -1.0]);
    println!("{:?}", t.find_paths(0, 6));
    println!("{:?}", t.generate_momentum_flow());
    println!("{}", t.evaluate_propagators(&[5.0, 6.0]));
    println!("{:?}", t.construct_cycles_per_loop_momentum(0, &[5.0, 6.0]));
    assert!(false);
}
