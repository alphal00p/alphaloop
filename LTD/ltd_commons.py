import vectors
zero_lv = vectors.LorentzVector([0.,0.,0.,0.])

class LoopTopology(object):
    """ A simple container for describing a loop topology."""

    def __init__(self, ltd_cut_structure, loop_lines, n_loops=1, name=None, **opts):
        """
            loop_lines          : A tuple of loop lines instances corresponding to each edge of the directed
                                  graph of this topology.
            ltd_cut_structure   : A tuple of tuples instructing what cuts result from the iterative
                                  application of LTD. Example: ((0,1,-1),(1,-1,0),...)
                                  Each entry of the tuple is either 0 (NO_CUT), -1 (NEGATIVE_CUT sol.)
                                  or +1 (POSTIVIVE_CUT sol.).
            n_loops             : Number of loops in this topology
            name                : Name of this topology
        """
        self.loop_lines         = loop_lines
        self.ltd_cut_structure  = ltd_cut_structure
        self.n_loops            = n_loops 
        self.name               = name
        
    def print_topology(self):
        """ Print the topology using graphviz."""

        try:
            from graphviz import Digraph
        except ImportError:
            print "The print function of the LoopTopology requires the package graphviz."

        dot = Digraph(comment=self.name if self.name else 'Topology',
                        graph_attr={'fontsize':'8','sep':'10', 'splines':'true'},
                        edge_attr={'fontsize':'8'},
                        node_attr={'fontsize':'8'})

        # Define node and 
        node_properties = {'width':'0.1','height':'0.1'}
        edge_properties = {'constraint':'true', 'arrowsize':'.5'}
        
        node_added = []
        for i_ll, loop_line in enumerate(self.loop_lines):
            if loop_line.start_node not in node_added:
                node_added.append('%d'%loop_line.start_node)
                dot.node('%d'%loop_line.start_node,'', style='filled', **node_properties)
            last_node = '%d'%loop_line.start_node
            for i,propagator in enumerate(loop_line.propagators[:-1]):
                new_node = 'LL%dP%d'%(i_ll,i+1)
                node_added.append(new_node)
                dot.node(new_node,'', **node_properties)
                dot.edge(last_node, new_node, label=('' if propagator.q!=zero_lv else str(propagator.signature)),**edge_properties)
                last_node = new_node
            if loop_line.end_node not in node_added:
                node_added.append('%d'%loop_line.end_node)
                dot.node('%d'%loop_line.end_node,'', style='filled', **node_properties)
            dot.edge(last_node, '%d'%loop_line.end_node, label=
                     ('' if loop_line.propagators[-1].q!=zero_lv else str(loop_line.propagators[-1].signature)),**edge_properties)

        filename = 'topology%s.gv'%('_%s'%self.name if self.name else '')
        dot.render(filename, view=True)

class LoopLine(object):
    """ A simple container for describing a loop line."""

    NO_CUT          = 0
    POSITIVE_CUT    = 1
    NEGATIVE_CUT    = -1

    def __init__(self, signature, propagators, start_node, end_node, **opts):
        """
            signature           : The signature of the loop line specifies its dependence on the loop
                                  momenta. It is a tuple of length n_loops and with entries which are
                                  either 0 (no dependence), +1 (positive dependence) or -1 (negative dependence).
                                  For instance (1,0,-1) means a loop momenta dependency reading 1*k + 0*l -1*m.
            propagators         : List of Propagator instances specifying all propagators making up this loop line.
            start_node          : Integer labeling the starting node of this directed edge.
            end_node            : Integer labeling the end node of this directed edge (where the arrow points).
        """
        self.signature      = signature
        self.propagators    = propagators
        # Automatically forward the signature attribute to the loop propagators if not specified by the user.
        for propagator in self.propagators:
            propagator.set_signature(self.signature, force=False)
        self.start_node     = start_node
        self.end_node       = end_node

class Propagator(object):
    """ A simple container for describing a loop propagator."""

    def __init__(self, q, m_squared, signature = None, **opts):
        self.q          = q
        self.m_squared  = m_squared
        # Note that this signature member is not strictly speaking necessary as it should always be the
        # same as the signature attribute of the LoopLine containing it. We forward it here for convenience however.
        self.signature  = signature

    def evaluate_inverse(self, loop_momenta):
        """ Evaluates the inverse propagator with the provided list loop momenta, given as a list of LorentzVector.""" 
        return (sum(wgt*loop_momentum for wgt, loop_momentum in zip(self.signature, loop_momenta)) + q).square() - self.m_squared

    def set_signature(self, signature, force=True):
        """ Overwrite (if force=True) the signature attribute of this propagator."""
        if (self.signature is None) or force:
            self.signature = signature

class TopologyCollection(dict):
    """ Collection of hard-coded topologies."""

    def add_topology(self, topology_to_add):

        assert(isinstance(topology_to_add,LoopTopology))
        if topology_to_add.name is None:
            raise BaseException("Specify the name of the hard-coded topology to add to the collection.")
        self[topology_to_add.name] = topology_to_add


hard_coded_topology_collection = TopologyCollection()

p       = vectors.LorentzVector([1.,0.,0.,0.])
# Add the double-triangle to the hard-coded topology collection:
hard_coded_topology_collection.add_topology(
    LoopTopology(
        name    = 'DoubleTriange',
        n_loops = 2,
        ltd_cut_structure = (
            (LoopLine.NEGATIVE_CUT  , LoopLine.NO_CUT           , LoopLine.POSITIVE_CUT ),
            (LoopLine.POSITIVE_CUT  , LoopLine.POSITIVE_CUT     , LoopLine.NO_CUT       ),
            (LoopLine.NO_CUT        , LoopLine.POSITIVE_CUT     , LoopLine.POSITIVE_CUT ),
        ),
        loop_lines = (
            LoopLine(
                start_node  = 1, 
                end_node    = 2,
                signature   = (1,0),
                propagators = (
                    Propagator(q=-p, m_squared=0.),
                    Propagator(q=zero_lv, m_squared=0.),
                )
            ),
            LoopLine(
                start_node  = 1, 
                end_node    = 2,
                signature   = (0,1),
                propagators = (
                    Propagator(q=p, m_squared=0.),
                    Propagator(q=zero_lv, m_squared=0.),
                )
            ),
            LoopLine(
                start_node  = 2, 
                end_node    = 1,
                signature   = (1,1),
                propagators = (
                    Propagator(q=zero_lv, m_squared=0.),
                )
            ),
        ) 
    )
)

# Example printout
#hard_coded_topology_collection['DoubleTriange'].print_topology()
