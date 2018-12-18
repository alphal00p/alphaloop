class Topology(object):
    def __init__(self, num_vertices, edge_map):
        self.num_vertices = num_vertices
        self.edge_map = edge_map

        vertices = [y for x in edge_map for y in x]

        self.ext = [i for i, x in enumerate(edge_map) if vertices.count(x[0]) == 1 or vertices.count(x[1]) == 1]
        self.loop_momenta = None
        self.momentum_map = None

    def find_path(self, start, dest):
        # find all paths from source to dest
        loop = start == dest
        start_check = 1 if loop else 0
        paths = [[(start, True)]] # store direction
        if not loop:
            paths.append([(start, False)])
        res = []
        while True:
            newpaths = []
            for p in paths:
                if len(p) > 1 and p[-1][0] == dest:
                    res.append(p)
                    continue
                last_vertex = self.edge_map[p[-1][0]][1] if p[-1][1] else self.edge_map[p[-1][0]][0]
                for i, x in enumerate(self.edge_map):
                    if all(i != pp[0] for pp in p[start_check:]):
                        if loop and i == start:
                            # if we need a loop, we need to enter from the right direction
                            if x[0] == last_vertex:
                                newpaths.append(p + [(i, True)])
                            continue

                        if x[0] == last_vertex and all(x[1] not in self.edge_map[pp[0]] for pp in p[start_check:-1]):
                            newpaths.append(p + [(i, True)])
                        if x[1] == last_vertex and all(x[0] not in self.edge_map[pp[0]] for pp in p[start_check:-1]):
                            newpaths.append(p + [(i, False)])
            paths = newpaths
            if len(paths) == 0:
                break
        return res

    def generate_momentum_flow(self, loop_momenta):
        """Generate momentum flow where `loop_momenta` are a
        list of all edge indices that go like 1/k^2
        """
        self.loop_momenta = loop_momenta

        # first we fix the loop momentum flows
        flows = []
        for l in loop_momenta:
            paths = self.find_path(l, l)
            # make sure other loop momenta propagators are not in the path
            paths = [x for x in paths if all(y[0] not in self.loop_momenta for y in x[1:-1])]
            # TODO: take shortest?
            flows.append(paths[0][:-1])

        # now route the external loop_momenta to the sink
        ext_flows = []
        for e in self.ext[:-1]:
            paths = self.find_path(e, self.ext[-1])
            paths = [x for x in paths if all(y[0] not in self.loop_momenta for y in x[1:-1])]
            ext_flows.append(paths[0])

        # propagator momenta
        mom = []
        s = []
        for i, x in enumerate(self.edge_map):
            if i in self.ext:
                mom.append(())
                continue
            if i in loop_momenta:
                mom.append([(i, True)])
                s.append("1/k{}^2".format(i))
            else:
                newmom = []
                s1 = "1/("
                for j, y in enumerate(flows):
                    for yy in y:
                        if yy[0] == i:
                            s1 += ("+" if yy[1] else "-") + "k{}".format(loop_momenta[j])
                            newmom.append((loop_momenta[j], yy[1]))
                            break
                for j, y in enumerate(ext_flows):
                    for yy in y:
                        if yy[0] == i:
                            s1 += ("+" if yy[1] else "-") + "p{}".format(self.ext[j])
                            newmom.append((self.ext[j], yy[1]))
                            break
                mom.append(tuple(newmom))
                s.append(s1 + ")^2")

        self.momentum_map = mom
        print(mom)
        print("Propagator structure: {}".format('*'.join(s)))


    def evaluate(self, loop_mom, ext_mom):
        res = 1.
        for x in self.momentum_map:
            if len(x) == 0:
                continue
            r = 0.
            for y, sign in x:
                val = loop_mom[self.loop_momenta.index(y)] if y in self.loop_momenta else ext_mom[self.ext.index(y)]
                r += (1. if sign else -1.) * val
            res *= r
        return 1. / res


t = Topology(6, [(0,1), (1,2), (1,3), (2,3), (3,4), (2,5)])
t.generate_momentum_flow([1])
t.evaluate([1], [2,3])
