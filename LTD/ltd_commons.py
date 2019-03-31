#!/usr/bin/env python

import os
import vectors
import math
zero_lv = vectors.LorentzVector([0.,0.,0.,0.])

#############################################################################################################
# HyperParameters
#############################################################################################################
class HyperParameters(dict):

    def __init__(self, *args, **opts):
        super(HyperParameters, self).__init__(*args, **opts)

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        
        # For now the hyperparameters dict is already supposed to be directly exportable to yaml
        return dict(self)

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        # Directly get HyperParameters from the flat_dict
        return HyperParameters(flat_dict)

    def export_to(self, output_path, format='yaml'):
        """ Exports these hyperparameters to a given format."""
        
        export_format = format.lower()
        allowed_export_format = ['yaml']
        if export_format not in ['yaml']:
            raise BaseException("Hyperparameters can only be exported in the following formats: %s"%(', '.join(allowed_export_format)))

        if export_format=='yaml':
            try:
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to export hyperparameters to yaml.")

        if output_path is not None:
            open(output_path,'w').write(yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False))
        else:
            return yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False)


    @staticmethod
    def import_from(input_path, format='yaml'):
        """ Imports this topology from a given format."""
        
        import_format = format.lower()
        allowed_import_format = ['yaml']
        if import_format not in ['yaml']:
            raise BaseException("Hyperparameters can only be imported from the following formats: %s"%(', '.join(allowed_import_format)))

        if import_format=='yaml':
            try: 
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to import hyperparameters from yaml.")
 
        if '\n' in input_path:
            return HyperParameters.from_flat_format(yaml.load(input_path, Loader=Loader))
        else:
            return HyperParameters.from_flat_format(yaml.load(open(input_path,'r'), Loader=Loader))


#############################################################################################################
# Define topology structures
#############################################################################################################

class LoopTopology(object):
    """ A simple container for describing a loop topology."""

    def __init__(self, ltd_cut_structure, loop_lines, external_kinematics, n_loops=1, name=None, analytical_result=None, **opts):
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
        self.external_kinematics = external_kinematics
        self.loop_lines          = loop_lines
        self.ltd_cut_structure   = ltd_cut_structure
        self.n_loops             = n_loops 
        self.name                = name
        if callable(analytical_result):
            self.analytical_result   = analytical_result(self.external_kinematics)
        else:
            self.analytical_result   = analytical_result
       
    def get_com_energy(self):
        """ Returns c.o.m energy of the external kinematics."""
        incoming_momenta_sum = vectors.LorentzVector()
        for v in self.external_kinematics:
            if v[0] > 0.:
                incoming_momenta_sum += v
        return math.sqrt(abs(incoming_momenta_sum.square()))

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

        filename = 'topology%s'%('_%s'%self.name if self.name else '')
        dot.render(filename, view=True)

    def export_to(self, output_path, format='yaml'):
        """ Exports this topology to a given format."""
        
        export_format = format.lower()
        allowed_export_format = ['yaml']
        if export_format not in ['yaml']:
            raise BaseException("Topology can only be exported in the following formats: %s"%(', '.join(allowed_export_format)))

        if export_format=='yaml':
            try:
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to export topologies to yaml.")

        if output_path is not None:
            open(output_path,'w').write(yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False))
        else:
            return yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False)

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        
        res={}
        
        res['name'] = self.name
        res['ltd_cut_structure'] = [list(cs) for cs in self.ltd_cut_structure]
        res['n_loops'] = self.n_loops
        res['loop_lines'] = [ll.to_flat_format() for ll in self.loop_lines]
        res['external_kinematics'] = [ [float(v) for v in vec] for vec in self.external_kinematics]
        res['analytical_result_real'] = float(self.analytical_result.real) if self.analytical_result else 0.
        res['analytical_result_imag'] = float(self.analytical_result.imag) if self.analytical_result else 0. 

        return res

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        return LoopTopology(
            name                =   flat_dict['name'],
            ltd_cut_structure   =   tuple([tuple(cs) for cs in flat_dict['ltd_cut_structure']]),
            n_loops             =   flat_dict['n_loops'],
            loop_lines          =   tuple([LoopLine.from_flat_format(ll) for ll in flat_dict['loop_lines']]),
            external_kinematics =   vectors.LorentzVectorList([vectors.LorentzVector(v) for v in flat_dict['external_kinematics']]),
            analytical_result   =   (None if (flat_dict['analytical_result_real']==0. and flat_dict['analytical_result_imag']==0.) 
                                     else complex(flat_dict['analytical_result_real'],flat_dict['analytical_result_imag']))
        ) 

    @staticmethod
    def import_from(input_path, format='yaml'):
        """ Imports this topology from a given format."""
        
        import_format = format.lower()
        allowed_import_format = ['yaml']
        if import_format not in ['yaml']:
            raise BaseException("Topology can only be imported from the following formats: %s"%(', '.join(allowed_import_format)))

        if import_format=='yaml':
            try: 
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to import topologies from yaml.")
 
        if '\n' in input_path:
            return LoopTopology.from_flat_format(yaml.load(input_path, Loader=Loader))
        else:
            return LoopTopology.from_flat_format(yaml.load(open(input_path,'r'), Loader=Loader))

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

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        return LoopLine(
            signature   =   tuple(flat_dict['signature']),
            propagators =   tuple([Propagator.from_flat_format(p) for p in flat_dict['propagators']]),
            start_node  =   flat_dict['start_node'],
            end_node    =   flat_dict['end_node']             
        ) 


    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        
        res={}
        
        res['signature'] = list(self.signature)
        res['propagators'] = [p.to_flat_format() for p in self.propagators]
        res['start_node'] = self.start_node
        res['end_node'] = self.end_node

        return res

class Propagator(object):
    """ A simple container for describing a loop propagator."""

    def __init__(self, q, m_squared, signature = None, **opts):
        self.q          = q
        self.m_squared  = m_squared
        # Note that this signature member is not strictly speaking necessary as it should always be the
        # same as the signature attribute of the LoopLine containing it. We forward it here for convenience however.
        self.signature  = signature

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        
        res={}
        
        res['q'] = [float(v) for v in self.q]
        res['m_squared'] = self.m_squared

        return res

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        return Propagator(
            q           =   vectors.LorentzVector(flat_dict['q']),
            m_squared   =   flat_dict['m_squared'] 
        ) 

    def evaluate_inverse(self, loop_momenta):
        """ Evaluates the inverse propagator with the provided list loop momenta, given as a list of LorentzVector.""" 
        return (sum(wgt*loop_momentum for wgt, loop_momentum in zip(self.signature, loop_momenta)) + self.q).square() - self.m_squared

    def set_signature(self, signature, force=True):
        """ Overwrite (if force=True) the signature attribute of this propagator."""
        if (self.signature is None) or force:
            self.signature = signature

class TopologyCollection(dict):
    """ Collection of hard-coded topologies."""

    def add_topology(self, topology_to_add, entry_name=None):

        assert(isinstance(topology_to_add,LoopTopology))
        if entry_name is None:
            if topology_to_add.name is None:
                raise BaseException("Specify the name of the hard-coded topology to add to the collection.")
            self[topology_to_add.name] = topology_to_add
        else:
            self[entry_name] = topology_to_add

    def export_to(self, output_path, format='yaml'):
        """ Exports this topology to a given format."""
        
        export_format = format.lower()
        allowed_export_format = ['yaml']
        if export_format not in ['yaml']:
            raise BaseException("Topology can only be exported in the following formats: %s"%(', '.join(allowed_export_format)))

        if export_format=='yaml':
            try:
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to export topologies to yaml.")

        flat_record = [self[k].to_flat_format() for k in sorted(self.keys())]
        if output_path is not None:
            open(output_path,'w').write(yaml.dump(flat_record, Dumper=Dumper))
        else:
            return yaml.dump(flat_record, Dumper=Dumper)

    @staticmethod
    def import_from(input_path, format='yaml'):
        """ Imports this topology from a given format."""
        
        import_format = format.lower()
        allowed_import_format = ['yaml']
        if import_format not in ['yaml']:
            raise BaseException("Topology can only be imported from the following formats: %s"%(', '.join(allowed_import_format)))

        if import_format=='yaml':
            try: 
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to import topologies from yaml.")
        
        if '\n' in input_path:
            flat_record = yaml.load(input_path, Loader=Loader)
        else:
            flat_record = yaml.load(open(input_path,'r'), Loader=Loader)

        result = TopologyCollection()
        for topology in flat_record:
            result.add_topology(LoopTopology.from_flat_format(topology))
        
        return result

#############################################################################################################
# Definition of hard-coded topol
#############################################################################################################

def create_hard_coded_topoloogy(topology_type, external_momenta, analytical_result=None, name=None, parameter_values = {}):
    """ Creates a hard-coded topology of the given name with specified kinematics. 
    
    ================================================
    /!\ ALL EXTERNAL MOMENTA ARE DEFINED AS INCOMING 
    ================================================
   
    Convention for loop normalisation in analytical result is 
        (1./ 2.*math.pi)**(4*n_loops)

    """

    if name is None:
        name = topology_type

    if topology_type =='DoubleTriangle':
        p1 = external_momenta[0]
        p2 = external_momenta[1]
        return LoopTopology(
            name    = name,
            n_loops = 2,
            external_kinematics = external_momenta,
            # Analytical result is simple and should be -(6*Zeta[3])/((16*pi^2)^2 s)
            # so we can use the analytical result here.
            analytical_result = (lambda ps: -0.00028922566024/ps[0].square()),
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
                        Propagator(q=-p1, m_squared=0.0),
                        Propagator(q=zero_lv, m_squared=0.0),
                    )
                ),
                LoopLine(
                    start_node  = 1, 
                    end_node    = 2,
                    signature   = (0,1),
                    propagators = (
                        Propagator(q=p1, m_squared=0.0),
                        Propagator(q=zero_lv, m_squared=0.0),
                    )
                ),
                LoopLine(
                    start_node  = 2, 
                    end_node    = 1,
                    signature   = (1,1),
                    propagators = (
                        Propagator(q=zero_lv, m_squared=0.0),
                    )
                ),
            ) 
        )

    if topology_type =='AltDoubleTriangle':
        p1 = external_momenta[0]
        p2 = external_momenta[1]
        return LoopTopology(
            name    = name,
            n_loops = 2,
            external_kinematics = external_momenta,
            # Analytical result is simple and should be -(6*Zeta[3])/((16*pi^2)^2 s)
            # so we can use the analytical result here.
            analytical_result = (lambda ps: -0.00028922566024/ps[0].square()),
            ltd_cut_structure = (
                (LoopLine.POSITIVE_CUT  , LoopLine.NO_CUT           , LoopLine.POSITIVE_CUT ),
                (LoopLine.NEGATIVE_CUT  , LoopLine.POSITIVE_CUT     , LoopLine.NO_CUT       ),
                (LoopLine.NO_CUT        , LoopLine.POSITIVE_CUT     , LoopLine.POSITIVE_CUT ),
            ),
            loop_lines = (
                LoopLine(
                    start_node  = 1, 
                    end_node    = 2,
                    signature   = (1,0),
                    propagators = (
                        Propagator(q=-p1, m_squared=0.),
                        Propagator(q=zero_lv, m_squared=0.),
                    )
                ),
                LoopLine(
                    start_node  = 2, 
                    end_node    = 1,
                    signature   = (0,1),
                    propagators = (
                        Propagator(q=zero_lv, m_squared=0.),
                        Propagator(q=p2, m_squared=0.),
                    )
                ),
                LoopLine(
                    start_node  = 1, 
                    end_node    = 2,
                    signature   = (-1,1),
                    propagators = (
                        Propagator(q=zero_lv, m_squared=0.),
                    )
                ),
            ) 
        )

    elif topology_type =='TriangleBox':
        p1 = external_momenta[0]
        p2 = external_momenta[1]
        p3 = external_momenta[2]
        return LoopTopology(
            name    = name,
            n_loops = 2,
            external_kinematics = external_momenta,
            analytical_result = analytical_result,
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
                        Propagator(q=zero_lv, m_squared=0.),
                        Propagator(q=p3, m_squared=0.),
                    )
                ),
                LoopLine(
                    start_node  = 1, 
                    end_node    = 2,
                    signature   = (0,1),
                    propagators = (
                        Propagator(q=zero_lv, m_squared=0.),
                        Propagator(q=p2, m_squared=0.),                        
                        Propagator(q=p1+p2, m_squared=0.),
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

    elif topology_type =='DoubleBox':
        p1 = external_momenta[0]
        p2 = external_momenta[1]
        p3 = external_momenta[2]
        p4 = external_momenta[3]
        return LoopTopology(
            name    = name,
            n_loops = 2,
            external_kinematics = external_momenta, 
            analytical_result = analytical_result,
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
                        Propagator(q=-p1-p2, m_squared=0.),
                        Propagator(q=-p2, m_squared=0.),                        
                        Propagator(q=zero_lv, m_squared=0.),
                    )
                ),
                LoopLine(
                    start_node  = 1, 
                    end_node    = 2,
                    signature   = (0,1),
                    propagators = (
                        Propagator(q=-p4-p3, m_squared=0.),
                        Propagator(q=-p3, m_squared=0.),                        
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
	
    elif topology_type =='TriangleBoxBox':
        p1,p2,p3 = external_momenta
        return LoopTopology(
            name    = name,
            n_loops = 3,
            external_kinematics = external_momenta,
            analytical_result = analytical_result,
            ltd_cut_structure = (
            	(LoopLine.POSITIVE_CUT	, LoopLine.NO_CUT	    , LoopLine.POSITIVE_CUT	,   LoopLine.NO_CUT		, LoopLine.POSITIVE_CUT		, LoopLine.NO_CUT	),
            	(LoopLine.POSITIVE_CUT	, LoopLine.NEGATIVE_CUT	    , LoopLine.POSITIVE_CUT	,   LoopLine.NO_CUT		, LoopLine.NO_CUT		, LoopLine.NO_CUT	),
            	(LoopLine.POSITIVE_CUT	, LoopLine.NEGATIVE_CUT	    , LoopLine.NO_CUT		,   LoopLine.NO_CUT		, LoopLine.NO_CUT		, LoopLine.NEGATIVE_CUT	),
            	(LoopLine.POSITIVE_CUT	, LoopLine.NO_CUT	    , LoopLine.NO_CUT		,   LoopLine.NO_CUT		, LoopLine.POSITIVE_CUT		, LoopLine.NEGATIVE_CUT	),
            	(LoopLine.POSITIVE_CUT	, LoopLine.NO_CUT	    , LoopLine.NO_CUT		,   LoopLine.NEGATIVE_CUT	, LoopLine.NO_CUT		, LoopLine.NEGATIVE_CUT	),
            	(LoopLine.POSITIVE_CUT	, LoopLine.NO_CUT	    , LoopLine.NEGATIVE_CUT	,   LoopLine.NO_CUT		, LoopLine.NO_CUT		, LoopLine.NEGATIVE_CUT	),
            	(LoopLine.POSITIVE_CUT	, LoopLine.NO_CUT	    , LoopLine.POSITIVE_CUT	,   LoopLine.NEGATIVE_CUT	, LoopLine.NO_CUT		, LoopLine.NO_CUT	),
            	(LoopLine.NO_CUT	, LoopLine.POSITIVE_CUT	    , LoopLine.NO_CUT		,   LoopLine.NO_CUT		, LoopLine.POSITIVE_CUT		, LoopLine.POSITIVE_CUT	),
            	(LoopLine.NO_CUT	, LoopLine.POSITIVE_CUT	    , LoopLine.NEGATIVE_CUT	,   LoopLine.NO_CUT		, LoopLine.POSITIVE_CUT		, LoopLine.NO_CUT	),
            	(LoopLine.NO_CUT	, LoopLine.NO_CUT	    , LoopLine.NO_CUT		,   LoopLine.POSITIVE_CUT	, LoopLine.POSITIVE_CUT		, LoopLine.POSITIVE_CUT	),
            	(LoopLine.NO_CUT	, LoopLine.NO_CUT	    , LoopLine.NEGATIVE_CUT	,   LoopLine.POSITIVE_CUT	, LoopLine.POSITIVE_CUT		, LoopLine.NO_CUT	),
            	(LoopLine.NO_CUT	, LoopLine.NO_CUT	    , LoopLine.POSITIVE_CUT	,   LoopLine.NO_CUT		, LoopLine.POSITIVE_CUT		, LoopLine.POSITIVE_CUT	),
            ),
            loop_lines = (
                LoopLine(
                    start_node  = 1, 
                    end_node    = 2,
                    signature   = (1,0,0),
                    propagators = (
                        Propagator(q=zero_lv, m_squared=0.),
                        Propagator(q=p3, m_squared=0.),
                    )
                ),
                LoopLine(
                    start_node  = 2, 
                    end_node    = 3,
                    signature   = (1,-1,0),
                    propagators = (
                        Propagator(q=p3, m_squared=0.),
                    )
                ),
                LoopLine(
                    start_node  = 3, 
                    end_node    = 4,
                    signature   = (1,-1,-1),
                    propagators = (
                        Propagator(q=p3, m_squared=0.),
                        Propagator(q=p3+p1, m_squared=0.),
                        Propagator(q=zero_lv, m_squared=0.),
                    )
                ),
                LoopLine(
                    start_node  = 4, 
                    end_node    = 1,
                    signature   = (1,-1,0),
                    propagators = (
                        Propagator(q=zero_lv, m_squared=0.),
                    )
                ),
                LoopLine(
                    start_node  = 2, 
                    end_node    = 1,
                    signature   = (0,1,0),
                    propagators = (
                        Propagator(q=zero_lv, m_squared=0.),
                    )
                ),
                LoopLine(
                    start_node  = 3, 
                    end_node    = 4,
                    signature   = (0,0,1),
                    propagators = (
                        Propagator(q=zero_lv, m_squared=0.),
                    )
                ),
            ) 
        )
	
    elif topology_type =='Triangle':
        p1 = external_momenta[0]
        p2 = external_momenta[1]
        p3 = external_momenta[2]
        if parameter_values == {}:
            parameters = {'m1': 0., 'm2': 0., 'm3': 0.}
        else:
            parameters = parameter_values
        return LoopTopology(
            name    = name,
            n_loops = 1,
            external_kinematics = external_momenta,
            analytical_result = analytical_result,
            ltd_cut_structure = (
                (LoopLine.POSITIVE_CUT,),
            ),
            loop_lines = (
                LoopLine(
                    start_node  = 1, 
                    end_node    = 1,
                    signature   = (1,),
                    propagators = (
                        Propagator(q=p1, m_squared=parameters['m1']**2),
                        Propagator(q=p1+p2, m_squared=parameters['m2']**2),
                        Propagator(q=zero_lv, m_squared=parameters['m3']**2),
                    )
                ),
            ) 
        )

    elif topology_type =='Box':
        p1 = external_momenta[0]
        p2 = external_momenta[1]
        p3 = external_momenta[2]
        p4 = external_momenta[3]
        if parameter_values == {}:
            parameters = {'m1': 0., 'm2': 0., 'm3': 0., 'm4': 0.}
        else:
            parameters = parameter_values
        return LoopTopology(
            name    = name,
            n_loops = 1,
            external_kinematics = external_momenta,
            analytical_result = analytical_result,
            ltd_cut_structure = (
                (LoopLine.POSITIVE_CUT,),
            ),
            loop_lines = (
                LoopLine(
                    start_node  = 1, 
                    end_node    = 1,
                    signature   = (1,),
                    propagators = (
                        Propagator(q=p1, m_squared=parameters['m1']**2),
                        Propagator(q=p1+p2, m_squared=parameters['m2']**2),
                        Propagator(q=p1+p2+p3, m_squared=parameters['m3']**2),
                        Propagator(q=zero_lv, m_squared=parameters['m4']**2),
                    )
                ),
            ) 
        )

    elif topology_type =='Pentagon':
        p1 = external_momenta[0]
        p2 = external_momenta[1]
        p3 = external_momenta[2]
        p4 = external_momenta[3]
        p5 = external_momenta[4]
        if parameter_values == {}:
            parameters = {'m1': 0., 'm2': 0., 'm3': 0., 'm4': 0., 'm5': 0.}
        else:
            parameters = parameter_values
        return LoopTopology(
            name    = name,
            n_loops = 1,
            external_kinematics = external_momenta,
            analytical_result = analytical_result,
            ltd_cut_structure = (
                (LoopLine.POSITIVE_CUT,),
            ),
            loop_lines = (
                LoopLine(
                    start_node  = 1, 
                    end_node    = 1,
                    signature   = (1,),
                    propagators = (
                        Propagator(q=p1, m_squared=parameters['m1']**2),
                        Propagator(q=p1+p2, m_squared=parameters['m2']**2),
                        Propagator(q=p1+p2+p3, m_squared=parameters['m3']**2),
                        Propagator(q=p1+p2+p3+p4, m_squared=parameters['m4']**2),                        
                        Propagator(q=zero_lv, m_squared=parameters['m5']**2),
                    )
                ),
            ) 
        )

    elif topology_type =='Decagon':
        p1 = external_momenta[0]
        p2 = external_momenta[1]
        p3 = external_momenta[2]
        p4 = external_momenta[3]
        p5 = external_momenta[4]
        p6 = external_momenta[5]
        p7 = external_momenta[6]
        p8 = external_momenta[7]
        p9 = external_momenta[8]
        p10 = external_momenta[9]
        if parameter_values == {}:
            parameters = {'m1': 0., 'm2': 0., 'm3': 0., 'm4': 0., 'm5': 0.,
                          'm6': 0., 'm7': 0., 'm8': 0., 'm9': 0., 'm10': 0.
                         }
        else:
            parameters = parameter_values
        return LoopTopology(
            name    = name,
            n_loops = 1,
            external_kinematics = external_momenta,
            analytical_result = analytical_result,
            ltd_cut_structure = (
                (LoopLine.POSITIVE_CUT,),
            ),
            loop_lines = (
                LoopLine(
                    start_node  = 1, 
                    end_node    = 1,
                    signature   = (1,),
                    propagators = (
                        Propagator(q=p1, m_squared=parameters['m1']**2),
                        Propagator(q=p1+p2, m_squared=parameters['m2']**2),
                        Propagator(q=p1+p2+p3, m_squared=parameters['m3']**2),
                        Propagator(q=p1+p2+p3+p4, m_squared=parameters['m4']**2),
                        Propagator(q=p1+p2+p3+p4+p5, m_squared=parameters['m5']**2),                        
                        Propagator(q=p1+p2+p3+p4+p5+p6, m_squared=parameters['m6']**2),                        
                        Propagator(q=p1+p2+p3+p4+p5+p6+p7, m_squared=parameters['m7']**2),                        
                        Propagator(q=p1+p2+p3+p4+p5+p6+p7+p8, m_squared=parameters['m8']**2),                        
                        Propagator(q=p1+p2+p3+p4+p5+p6+p7+p8+p9, m_squared=parameters['m9']**2), 
                        Propagator(q=zero_lv, m_squared=parameters['m10']**2),
                    )
                ),
            ) 
        )

    elif topology_type =='Tringigon':
        p1 = external_momenta[0]
        p2 = external_momenta[1]
        p3 = external_momenta[2]
        p4 = external_momenta[3]
        p5 = external_momenta[4]
        p6 = external_momenta[5]
        p7 = external_momenta[6]
        p8 = external_momenta[7]
        p9 = external_momenta[8]
        p10 = external_momenta[9]
        p11 = external_momenta[10]
        p12 = external_momenta[11]
        p13 = external_momenta[12]
        p14 = external_momenta[13]
        p15 = external_momenta[14]
        p16 = external_momenta[15]
        p17 = external_momenta[16]
        p18 = external_momenta[17]
        p19 = external_momenta[18]
        p20 = external_momenta[19]
        p21 = external_momenta[20]
        p22 = external_momenta[21]
        p23 = external_momenta[22]
        p24 = external_momenta[23]
        p25 = external_momenta[24]
        p26 = external_momenta[25]
        p27 = external_momenta[26]
        p28 = external_momenta[27]
        p29 = external_momenta[28]
        p30 = external_momenta[29]
        psum10 = p1+p2+p3+p4+p5+p6+p7+p8+p9
        psum20 = psum10+p10+p11+p13+p13+p14+p15+p16+p17+p18+p19
        if parameter_values == {}:
            parameters = {'m1': 0., 'm2': 0., 'm3': 0., 'm4': 0., 'm5': 0.,
                          'm6': 0., 'm7': 0., 'm8': 0., 'm9': 0., 'm10': 0.,
                          'm11': 0., 'm12': 0., 'm13': 0., 'm14': 0., 'm15': 0.,
                          'm16': 0., 'm17': 0., 'm18': 0., 'm19': 0., 'm20': 0.,
                          'm21': 0., 'm22': 0., 'm23': 0., 'm24': 0., 'm25': 0.,
                          'm26': 0., 'm27': 0., 'm28': 0., 'm29': 0., 'm30': 0.,
                         }
        else:
            parameters = parameter_values
        return LoopTopology(
            name    = name,
            n_loops = 1,
            external_kinematics = external_momenta,
            analytical_result = analytical_result,
            ltd_cut_structure = (
                (LoopLine.POSITIVE_CUT,),
            ),
            loop_lines = (
                LoopLine(
                    start_node  = 1, 
                    end_node    = 1,
                    signature   = (1,),
                    propagators = (
                        Propagator(q=p1, m_squared=parameters['m1']**2),
                        Propagator(q=p1+p2, m_squared=parameters['m2']**2),
                        Propagator(q=p1+p2+p3, m_squared=parameters['m3']**2),
                        Propagator(q=p1+p2+p3+p4, m_squared=parameters['m4']**2),
                        Propagator(q=p1+p2+p3+p4+p5, m_squared=parameters['m5']**2),                        
                        Propagator(q=p1+p2+p3+p4+p5+p6, m_squared=parameters['m6']**2),                        
                        Propagator(q=p1+p2+p3+p4+p5+p6+p7, m_squared=parameters['m7']**2),                        
                        Propagator(q=p1+p2+p3+p4+p5+p6+p7+p8, m_squared=parameters['m8']**2),                        
                        Propagator(q=p1+p2+p3+p4+p5+p6+p7+p8+p9, m_squared=parameters['m9']**2),
                        Propagator(q=psum10+p10, m_squared=parameters['m10']**2),
                        Propagator(q=psum10+p10+p11, m_squared=parameters['m11']**2),
                        Propagator(q=psum10+p10+p11+p12, m_squared=parameters['m12']**2),
                        Propagator(q=psum10+p10+p11+p12+p13, m_squared=parameters['m13']**2),
                        Propagator(q=psum10+p10+p11+p12+p13+p14, m_squared=parameters['m14']**2),                        
                        Propagator(q=psum10+p10+p11+p12+p13+p14+p15, m_squared=parameters['m15']**2),                        
                        Propagator(q=psum10+p10+p11+p12+p13+p14+p15+p16, m_squared=parameters['m16']**2),                        
                        Propagator(q=psum10+p10+p11+p12+p13+p14+p15+p16+p17, m_squared=parameters['m17']**2),                        
                        Propagator(q=psum10+p10+p11+p12+p13+p14+p15+p16+p17+p18, m_squared=parameters['m18']**2),
                        Propagator(q=psum10+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19, m_squared=parameters['m19']**2),                        
                        Propagator(q=psum20+p20, m_squared=parameters['m20']**2),
                        Propagator(q=psum20+p20+p21, m_squared=parameters['m21']**2),
                        Propagator(q=psum20+p20+p21+p22, m_squared=parameters['m22']**2),
                        Propagator(q=psum20+p20+p21+p22+p23, m_squared=parameters['m23']**2),
                        Propagator(q=psum20+p20+p21+p22+p23+p24, m_squared=parameters['m24']**2),                        
                        Propagator(q=psum20+p20+p21+p22+p23+p24+p25, m_squared=parameters['m25']**2),                        
                        Propagator(q=psum20+p20+p21+p22+p23+p24+p25+p26, m_squared=parameters['m26']**2),                        
                        Propagator(q=psum20+p20+p21+p22+p23+p24+p25+p26+p27, m_squared=parameters['m27']**2),                        
                        Propagator(q=psum20+p20+p21+p22+p23+p24+p25+p26+p27+p28, m_squared=parameters['m28']**2),
                        Propagator(q=psum20+p20+p21+p22+p23+p24+p25+p26+p27+p28+p29, m_squared=parameters['m29']**2),                        
                        Propagator(q=zero_lv, m_squared=parameters['m30']**2),
                    )
                ),
            ) 
        )

    else:
        raise BaseException("Unknown hard-coded topology name '%s'"%name)

#############################################################################################################
# Create the collection of hard-coded topologies.
#############################################################################################################

hard_coded_topology_collection = TopologyCollection()

# Add the double-triangle to the hard-coded topology collection:
hard_coded_topology_collection.add_topology(   
    create_hard_coded_topoloogy(
        'DoubleTriangle',
        vectors.LorentzVectorList([
                vectors.LorentzVector([1.,0.1,0.2,0.1]),
                vectors.LorentzVector([-1.,-0.1,-0.2,-0.1])
        ])
        # Analytical is given by its exact function directly for 'DoubleTriangle'        
    ),
)

hard_coded_topology_collection.add_topology(   
    create_hard_coded_topoloogy(
        'AltDoubleTriangle',
        vectors.LorentzVectorList([
                vectors.LorentzVector([1.,0.1,0.2,0.1]),
                vectors.LorentzVector([-1.,-0.1,-0.2,-0.1])
        ])
        # Analytical is given by its exact function directly for 'DoubleTriangle'        
    ),
)

hard_coded_topology_collection.add_topology(   
    create_hard_coded_topoloogy(
        'DoubleTriangle',
        vectors.LorentzVectorList([
                vectors.LorentzVector([1.,1.3,0.5,2.1]),
                vectors.LorentzVector([-1.,-1.3,-0.5,-2.1])
        ]),
        name = 'DoubleTriangle_euclidean',
        # Analytical is given by its exact function directly for 'DoubleTriangle'        
    ),
    entry_name = 'DoubleTriangle_euclidean'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'TriangleBox',
        vectors.LorentzVectorList([
                vectors.LorentzVector([39.7424,-14.1093,0.102709,20.4908]),
                vectors.LorentzVector([50.2576,14.1093,-0.102709,-20.4908]),
                vectors.LorentzVector([-90.,0.,0.,0.]),
        ]),
        # Analytical result known but exact numerical result simply copied here from C^(2) of table 1
        # https://arxiv.org/pdf/1211.0509.pdf
        analytical_result = -1.832e-11
    ),
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'DoubleBox',        
        vectors.LorentzVectorList([
                vectors.LorentzVector([19.6586,-7.15252,-0.206016,8.96383]),
                vectors.LorentzVector([26.874,7.04203,-0.0501295,-12.9055]),
                vectors.LorentzVector([43.4674,0.110491,0.256146,3.9417]),
                vectors.LorentzVector([-90.,0.,0.,0.]),
        ]),
        # Analytical result known but exact numerical result simply copied here from D^(2) of table 1 of
        # https://arxiv.org/pdf/1211.0509.pdf
        analytical_result =  -5.897e-14
    ),
)

"""
hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'TriangleBoxBox',        
       	vectors.LorentzVectorList([
                vectors.LorentzVector([0.,5.,0.,0.]),
                vectors.LorentzVector([0.,0.,6.,0.]),
                vectors.LorentzVector([0.,-5.,-6.,0.]),
        ]),
        analytical_result =  1.76576112187575e-10j
    ),
)
"""

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'TriangleBoxBox',        
        vectors.LorentzVectorList([
                vectors.LorentzVector([0.,5.,0.,0.]),
                vectors.LorentzVector([0.,0.,6.,0.]),
                vectors.LorentzVector([-0.,-5.,-6.,0.]),
        ]),
        analytical_result =  4.52887419598747e-11j,
    ),
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Triangle', #P3 from Rodrigo
        vectors.LorentzVectorList([
                vectors.LorentzVector([10.51284,6.89159,-7.40660,-2.85795]),
                vectors.LorentzVector([6.45709,2.46635,5.84093,1.22257]),
                -vectors.LorentzVector([10.51284,6.89159,-7.40660,-2.85795]) - vectors.LorentzVector([6.45709,2.46635,5.84093,1.22257]),
        ]),
        analytical_result = 6.68103e-4 + 5.37305e-4j,
        parameter_values = {'m1': 0.52559, 'm2': 0.52559, 'm3': 0.52559}
    ),
)
"""
hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Triangle', #P4 from Rodrigo
        vectors.LorentzVectorList([
                vectors.LorentzVector([95.77004,31.32025,-34.08106,-9.38565]),
                vectors.LorentzVector([94.54738,-53.84229,67.11107,45.56763]),
                -vectors.LorentzVector([95.77004,31.32025,-34.08106,-9.38565]) - vectors.LorentzVector([94.54738,-53.84229,67.11107,45.56763]),
        ]),
        analytical_result = 1.01665e-6 - 5.61370-7j,
        parameter_values = {'m1': 83.02643, 'm2': 76.12873, 'm3': 55.00359}
    ),
)
"""

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Box', #P7 from Rodrigo
        vectors.LorentzVectorList([
                vectors.LorentzVector([62.80274,-49.71968,-5.53340,-79.44048]),
                vectors.LorentzVector([48.59375,-1.65847,34.91140,71.89564]),
                vectors.LorentzVector([76.75934,-19.14334,-17.10279,30.22959]),
                 (- vectors.LorentzVector([62.80274,-49.71968,-5.53340,-79.44048]) -vectors.LorentzVector([48.59375,-1.65847,34.91140,71.89564])
                 - vectors.LorentzVector([76.75934,-19.14334,-17.10279,30.22959])),
        ]),
        analytical_result = 3.03080e-10 -2.38766e-10j,
        parameter_values = {'m1': 9.82998, 'm2': 9.82998, 'm3': 9.82998, 'm4': 9.82998}
    ),
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Pentagon', #P15 from Rodrigo
        vectors.LorentzVectorList([
                vectors.LorentzVector([94.79774, -70.04005, -84.77221, 36.09812]),
                vectors.LorentzVector([-42.15872, -36.33754, -14.72331, -41.24018]),
                vectors.LorentzVector([73.77293, 88.37064, 33.47296, -24.17542]),
                vectors.LorentzVector([81.85638, 77.17370, -62.39774, -6.89737]),
                (   -vectors.LorentzVector([94.79774, -70.04005, -84.77221, 36.09812])
                    -vectors.LorentzVector([-42.15872, -36.33754, -14.72331, -41.24018])
                    -vectors.LorentzVector([73.77293, 88.37064, 33.47296, -24.17542])
                    -vectors.LorentzVector([81.85638, 77.17370, -62.39774, -6.89737])
                )
        ]),
        analytical_result = 6.55440e-14 -4.29464e-15j,
        parameter_values = {'m1': 1.30619, 'm2': 1.30619, 'm3': 1.30619, 'm4': 1.26692, 'm5': 1.26692}
    ),
)

# ================================================================================
# Decagons with 10 external particles with masses increasing from 100 GeV to 2 TeV
# ================================================================================

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Decagon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.9998000000000000e+04,     0.0000000000000000e+00,      0.0000000000000000e+00,     0.9997499887471868e+04]),
                vectors.LorentzVector([ 0.1000200000000000e+05,     0.0000000000000000e+00,      0.0000000000000000e+00,    -0.9997499887471868e+04]),
                vectors.LorentzVector([-0.1237893881022244e+04,     0.1267918824409219e+03,     -0.6964708931616466e+03,     0.8838740714592880e+03]),
                vectors.LorentzVector([-0.2105530188905366e+04,     0.7964815951476172e+03,      0.1814849773817702e+04,    -0.1232669601183012e+03]),
                vectors.LorentzVector([-0.1344535096182690e+04,     0.8763606318834720e+03,      0.4476622405756132e+03,    -0.1713627325722238e+03]),
                vectors.LorentzVector([-0.5681047787072152e+04,    -0.2392205859044582e+04,     -0.4819795346193552e+04,     0.1453006506441443e+04]),
                vectors.LorentzVector([-0.1508708755159517e+04,     0.4159161811996088e+03,     -0.3306462124937450e+02,    -0.6419677320029004e+03]),
                vectors.LorentzVector([-0.3225969255726788e+04,     0.3968037749821656e+03,      0.7654076259435024e+03,    -0.2722788197638934e+04]),
                vectors.LorentzVector([-0.2526503858267339e+04,    -0.8560561209020956e+03,      0.1284795475177741e+04,     0.1053418364501306e+04]),
                vectors.LorentzVector([-0.2369811177663906e+04,     0.6359079142928917e+03,      0.1236615745090013e+04,     0.2690866799303233e+03]),
        ]),
        analytical_result = -1.01888646237782899e-063+6.86568168003053015e-063j,
        parameter_values = {'m1': 200.0, 'm2': 400.0, 'm3': 600.0, 'm4': 800.0, 'm5': 1000.0,
                            'm6': 1200.0, 'm7': 1400.0, 'm8': 1600.0, 'm9': 1800.0, 'm10': 2000.0,
                           },
        name = 'Decagon_P1_physical_massive',
    ),
    entry_name = 'Decagon_P1_physical_massive'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Decagon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.9997499887471868e+04,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.9998000000000000e+04]),
                vectors.LorentzVector([-0.9997499887471868e+04,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.1000200000000000e+05]),
                vectors.LorentzVector([ 0.8838740714592880e+03,     0.1267918824409219e+03,    -0.6964708931616466e+03,    -0.1237893881022244e+04]),
                vectors.LorentzVector([-0.1232669601183012e+03,     0.7964815951476172e+03,     0.1814849773817702e+04,    -0.2105530188905366e+04]),
                vectors.LorentzVector([-0.1713627325722238e+03,     0.8763606318834720e+03,     0.4476622405756132e+03,    -0.1344535096182690e+04]),
                vectors.LorentzVector([ 0.1453006506441443e+04,    -0.2392205859044582e+04,    -0.4819795346193552e+04,    -0.5681047787072152e+04]),
                vectors.LorentzVector([-0.6419677320029004e+03,     0.4159161811996088e+03,    -0.3306462124937450e+02,    -0.1508708755159517e+04]),
                vectors.LorentzVector([-0.2722788197638934e+04,     0.3968037749821656e+03,     0.7654076259435024e+03,    -0.3225969255726788e+04]),
                vectors.LorentzVector([ 0.1053418364501306e+04,    -0.8560561209020956e+03,     0.1284795475177741e+04,    -0.2526503858267339e+04]),
                vectors.LorentzVector([ 0.2690866799303233e+03,     0.6359079142928917e+03,     0.1236615745090013e+04,    -0.2369811177663906e+04]),
        ]),
        analytical_result = -8.60853357693996182e-065+1.50530418005552852e-063j,
        parameter_values = {'m1': 200.0, 'm2': 400.0, 'm3': 600.0, 'm4': 800.0, 'm5': 1000.0,
                            'm6': 1200.0, 'm7': 1400.0, 'm8': 1600.0, 'm9': 1800.0, 'm10': 2000.0,
                           },
        name = 'Decagon_P1_euclidean_massive',
    ),
    entry_name = 'Decagon_P1_euclidean_massive'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Decagon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.9998000000000000e+04,     0.0000000000000000e+00,      0.0000000000000000e+00,     0.9997499887471868e+04]),
                vectors.LorentzVector([ 0.1000200000000000e+05,     0.0000000000000000e+00,      0.0000000000000000e+00,    -0.9997499887471868e+04]),
                vectors.LorentzVector([-0.1237893881022244e+04,     0.1267918824409219e+03,     -0.6964708931616466e+03,     0.8838740714592880e+03]),
                vectors.LorentzVector([-0.2105530188905366e+04,     0.7964815951476172e+03,      0.1814849773817702e+04,    -0.1232669601183012e+03]),
                vectors.LorentzVector([-0.1344535096182690e+04,     0.8763606318834720e+03,      0.4476622405756132e+03,    -0.1713627325722238e+03]),
                vectors.LorentzVector([-0.5681047787072152e+04,    -0.2392205859044582e+04,     -0.4819795346193552e+04,     0.1453006506441443e+04]),
                vectors.LorentzVector([-0.1508708755159517e+04,     0.4159161811996088e+03,     -0.3306462124937450e+02,    -0.6419677320029004e+03]),
                vectors.LorentzVector([-0.3225969255726788e+04,     0.3968037749821656e+03,      0.7654076259435024e+03,    -0.2722788197638934e+04]),
                vectors.LorentzVector([-0.2526503858267339e+04,    -0.8560561209020956e+03,      0.1284795475177741e+04,     0.1053418364501306e+04]),
                vectors.LorentzVector([-0.2369811177663906e+04,     0.6359079142928917e+03,      0.1236615745090013e+04,     0.2690866799303233e+03]),
        ]),
        analytical_result = 9.32266615611493377e-064-8.49710641987076969e-063j,
        parameter_values = {'m1': 0.0, 'm2': 0.0, 'm3': 0.0, 'm4': 0.0, 'm5': 0.0,
                            'm6': 0.0, 'm7': 0.0, 'm8': 0.0, 'm9': 0.0, 'm10': 0.0,
                           },
        name = 'Decagon_P1_physical_massless',
    ),
    entry_name = 'Decagon_P1_physical_massless'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Decagon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.9997499887471868e+04,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.9998000000000000e+04]),
                vectors.LorentzVector([-0.9997499887471868e+04,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.1000200000000000e+05]),
                vectors.LorentzVector([ 0.8838740714592880e+03,     0.1267918824409219e+03,    -0.6964708931616466e+03,    -0.1237893881022244e+04]),
                vectors.LorentzVector([-0.1232669601183012e+03,     0.7964815951476172e+03,     0.1814849773817702e+04,    -0.2105530188905366e+04]),
                vectors.LorentzVector([-0.1713627325722238e+03,     0.8763606318834720e+03,     0.4476622405756132e+03,    -0.1344535096182690e+04]),
                vectors.LorentzVector([ 0.1453006506441443e+04,    -0.2392205859044582e+04,    -0.4819795346193552e+04,    -0.5681047787072152e+04]),
                vectors.LorentzVector([-0.6419677320029004e+03,     0.4159161811996088e+03,    -0.3306462124937450e+02,    -0.1508708755159517e+04]),
                vectors.LorentzVector([-0.2722788197638934e+04,     0.3968037749821656e+03,     0.7654076259435024e+03,    -0.3225969255726788e+04]),
                vectors.LorentzVector([ 0.1053418364501306e+04,    -0.8560561209020956e+03,     0.1284795475177741e+04,    -0.2526503858267339e+04]),
                vectors.LorentzVector([ 0.2690866799303233e+03,     0.6359079142928917e+03,     0.1236615745090013e+04,    -0.2369811177663906e+04]),
        ]),
        analytical_result = -3.17196357423582536e-063-1.69292136375179309e-063j,
        parameter_values = {'m1': 0.0, 'm2': 0.0, 'm3': 0.0, 'm4': 0.0, 'm5': 0.0,
                            'm6': 0.0, 'm7': 0.0, 'm8': 0.0, 'm9': 0.0, 'm10': 0.0,
                           },
        name = 'Decagon_P1_euclidean_massless',
    ),
    entry_name = 'Decagon_P1_euclidean_massless'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Decagon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.9997499887471868e+04,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.9998000000000000e+04]),
                vectors.LorentzVector([-0.9997499887471868e+04,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.1000200000000000e+05]),
                vectors.LorentzVector([ 0.8838740714592880e+03,     0.1267918824409219e+03,    -0.6964708931616466e+03,    -0.1237893881022244e+04]),
                vectors.LorentzVector([-0.1232669601183012e+03,     0.7964815951476172e+03,     0.1814849773817702e+04,    -0.2105530188905366e+04]),
                vectors.LorentzVector([-0.1713627325722238e+03,     0.8763606318834720e+03,     0.4476622405756132e+03,    -0.1344535096182690e+04]),
                vectors.LorentzVector([ 0.1453006506441443e+04,    -0.2392205859044582e+04,    -0.4819795346193552e+04,    -0.5681047787072152e+04]),
                vectors.LorentzVector([-0.6419677320029004e+03,     0.4159161811996088e+03,    -0.3306462124937450e+02,    -0.1508708755159517e+04]),
                vectors.LorentzVector([-0.2722788197638934e+04,     0.3968037749821656e+03,     0.7654076259435024e+03,    -0.3225969255726788e+04]),
                vectors.LorentzVector([ 0.1053418364501306e+04,    -0.8560561209020956e+03,     0.1284795475177741e+04,    -0.2526503858267339e+04]),
                vectors.LorentzVector([ 0.2690866799303233e+03,     0.6359079142928917e+03,     0.1236615745090013e+04,    -0.2369811177663906e+04]),
        ]),
        analytical_result = -3.17196357423582536e-063-1.69292136375179309e-063j,
        parameter_values = {'m1': 0.0, 'm2': 0.0, 'm3': 0.0, 'm4': 0.0, 'm5': 0.0,
                            'm6': 0.0, 'm7': 0.0, 'm8': 0.0, 'm9': 0.0, 'm10': 0.0,
                           },
        name = 'Decagon_P1_euclidean_massless',
    ),
    entry_name = 'Decagon_P1_euclidean_massless'
)

# ===========================================================================================
# Decagons with 10 external particles with masses all equal to 500 GeV and no internal masses
# ===========================================================================================

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Decagon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.4999749993749687e+00,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.1750000000000000e+01]),
                vectors.LorentzVector([-0.4999749993749687e+00,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.1750000000000000e+01]),
                vectors.LorentzVector([ 0.5319269421240436e-01,     0.7630500824806391e-02,    -0.4191452656442592e-01,    -0.2391684683852854e+00]),
                vectors.LorentzVector([-0.7418366402851829e-02,     0.4793330102618212e-01,     0.1092200259939211e+00,    -0.4186363970344227e+00]),
                vectors.LorentzVector([-0.1031283270711509e-01,     0.5274052562103488e-01,     0.2694089739963054e-01,    -0.2111265889119278e+00]),
                vectors.LorentzVector([ 0.8744382631132949e-01,    -0.1439660680884015e+00,    -0.2900615690571847e+00,    -0.1174109758803568e+01]),
                vectors.LorentzVector([-0.3863445525252511e-01,     0.2503037814879582e-01,    -0.1989872023808945e-02,    -0.1622166199303486e+00]),
                vectors.LorentzVector([-0.1638609443119325e+00,     0.2388016861961569e-01,     0.4606322903831265e-01,    -0.6018320233011958e+00]),
                vectors.LorentzVector([ 0.6339609085730495e-01,    -0.5151857367262940e-01,     0.7732066709885181e-01,    -0.3940681898996678e+00]),
                vectors.LorentzVector([ 0.1619398729338573e-01,     0.3826976752059603e-01,     0.7442114811470348e-01,    -0.2988419537335841e+00]),
        ]),
        analytical_result = -1.05384532825251748e-2+7.69076155224268383e-3j,
        parameter_values = {'m1': 0.0, 'm2': 0.0, 'm3': 0.0, 'm4': 0.0, 'm5': 0.0,
                            'm6': 0.0, 'm7': 0.0, 'm8': 0.0, 'm9': 0.0, 'm10': 0.0,
                           },
        name = 'Decagon_P1_one_ellipse',
    ),
    entry_name = 'Decagon_P1_one_ellipse'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Decagon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.4999749993749686e+00,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.3499999999999999e+01]),
                vectors.LorentzVector([-0.4999749993749686e+00,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.3499999999999999e+01]),
                vectors.LorentzVector([ 0.1048610088837169e+00,    -0.6590329264207770e-01,     0.3386480688194951e-01,    -0.8994632372733035e+00]),
                vectors.LorentzVector([ 0.5601434019370606e-01,     0.1832339746854182e-01,    -0.6206486061145505e-02,    -0.4163013450570706e+00]),
                vectors.LorentzVector([ 0.3009395901511873e-01,     0.4678686980455751e-01,    -0.1955231641799209e+00,    -0.1423411148974809e+01]),
                vectors.LorentzVector([-0.2469715427257115e-01,    -0.4358936740107071e-01,     0.4428064486243019e-01,    -0.4693528912897491e+00]),
                vectors.LorentzVector([-0.1200693811817634e-01,     0.4276555521667993e-01,     0.9076045429518785e-01,    -0.7081952044823810e+00]),
                vectors.LorentzVector([ 0.3089208388056481e-02,     0.7524676827788172e-01,    -0.5791613252166179e-01,    -0.6659534441677145e+00]),
                vectors.LorentzVector([-0.2411258803768644e-01,     0.7371721423248271e-02,    -0.9750760233496275e-01,    -0.7058725755366821e+00]),
                vectors.LorentzVector([-0.1332418360521643e+00,    -0.8100165214776085e-01,     0.1882474790581234e+00,    -0.1711450153218290e+01]),
        ]),
        analytical_result = -1.29287202675867206e-5+1.61392009112210223e-5j,
        parameter_values = {'m1': 0.0, 'm2': 0.0, 'm3': 0.0, 'm4': 0.0, 'm5': 0.0,
                            'm6': 0.0, 'm7': 0.0, 'm8': 0.0, 'm9': 0.0, 'm10': 0.0,
                           },
        name = 'Decagon_P2_one_ellipse',
    ),
    entry_name = 'Decagon_P2_one_ellipse'
)


hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Decagon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.5000000000000000e+01,     0.0000000000000000e+00,      0.0000000000000000e+00,     0.4974937185533100e+01]),
                vectors.LorentzVector([ 0.5000000000000000e+01,     0.0000000000000000e+00,      0.0000000000000000e+00,    -0.4974937185533100e+01]),
                vectors.LorentzVector([-0.7858249945336048e+00,     0.6787719733389865e-01,     -0.3728510953725848e+00,     0.4731761498589067e+00]),
                vectors.LorentzVector([-0.1174780067368942e+01,     0.4263911645277731e+00,      0.9715678469101043e+00,    -0.6599015343587409e-01]),
                vectors.LorentzVector([-0.7320893245338694e+00,     0.4691538795768828e+00,      0.2396530255526758e+00,    -0.9173790774737208e-01]),
                vectors.LorentzVector([-0.3025359192261747e+01,    -0.1280651616110646e+01,     -0.2580245623965889e+01,     0.7778574421837994e+00]),
                vectors.LorentzVector([-0.6465300920084710e+00,     0.2226579822158695e+00,     -0.1770092673212064e-01,    -0.3436731878120021e+00]),
                vectors.LorentzVector([-0.1608633085477756e+01,     0.2124262817049831e+00,      0.4097559202281392e+00,    -0.1457626688961433e+01]),
                vectors.LorentzVector([-0.1118539827394456e+01,    -0.4582839936494981e+00,      0.6878067769281205e+00,     0.5639405680069816e+00]),
                vectors.LorentzVector([-0.9082434164211539e+00,     0.3404291044007367e+00,      0.6620140764515544e+00,     0.1440537779069934e+00]),
        ]),
        analytical_result =  1.44683214438469199e-011-5.62046076823293754e-010j,
        parameter_values = {'m1': 0.0, 'm2': 0.0, 'm3': 0.0, 'm4': 0.0, 'm5': 0.0,
                            'm6': 0.0, 'm7': 0.0, 'm8': 0.0, 'm9': 0.0, 'm10': 0.0,
                           },
        name = 'Decagon_P2_physical_massless',
    ),
    entry_name = 'Decagon_P2_physical_massless'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Decagon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.4999749993749687e+00,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.2500000000000000e+01]),
                vectors.LorentzVector([-0.4999749993749687e+00,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.2500000000000000e+01]),
                vectors.LorentzVector([ 0.5319269421240436e-01,     0.7630500824806391e-02,    -0.4191452656442592e-01,    -0.3416692405504076e+00]),
                vectors.LorentzVector([-0.7418366402851829e-02,     0.4793330102618212e-01,     0.1092200259939211e+00,    -0.5980519957634609e+00]),
                vectors.LorentzVector([-0.1031283270711509e-01,     0.5274052562103488e-01,     0.2694089739963054e-01,    -0.3016094127313255e+00]),
                vectors.LorentzVector([ 0.8744382631132949e-01,    -0.1439660680884015e+00,    -0.2900615690571847e+00,    -0.1677299655433668e+01]),
                vectors.LorentzVector([-0.3863445525252511e-01,     0.2503037814879582e-01,    -0.1989872023808945e-02,    -0.2317380284719265e+00]),
                vectors.LorentzVector([-0.1638609443119325e+00,     0.2388016861961569e-01,     0.4606322903831265e-01,    -0.8597600332874225e+00]),
                vectors.LorentzVector([ 0.6339609085730495e-01,    -0.5151857367262940e-01,     0.7732066709885181e-01,    -0.5629545569995255e+00]),
                vectors.LorentzVector([ 0.1619398729338573e-01,     0.3826976752059603e-01,     0.7442114811470348e-01,    -0.4269170767622630e+00]),
        ]),
        analytical_result = 1.58733080719071658e-5j,
        parameter_values = {'m1': 0.0, 'm2': 0.0, 'm3': 0.0, 'm4': 0.0, 'm5': 0.0,
                            'm6': 0.0, 'm7': 0.0, 'm8': 0.0, 'm9': 0.0, 'm10': 0.0,
                           },
        name = 'Decagon_P2_euclidean_massless',
    ),
    entry_name = 'Decagon_P2_euclidean_massless'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Decagon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.4999999999999998e+01,     0.0000000000000000e+00,      0.0000000000000000e+00,     0.4999749993749686e+01]),
                vectors.LorentzVector([ 0.4999999999999998e+01,     0.0000000000000000e+00,      0.0000000000000000e+00,    -0.4999749993749686e+01]),
                vectors.LorentzVector([-0.3411557641735493e+01,     0.8034973126969602e+00,     -0.5044645607458048e+00,     0.3276603900126494e+01]),
                vectors.LorentzVector([-0.5035286543834834e+00,    -0.3506951976120900e-01,     -0.3731809945129911e+00,     0.3324863604770198e+00]),
                vectors.LorentzVector([-0.1867802059883001e+01,     0.6723228534567237e+00,      0.7020515084050917e+00,    -0.1594142463887442e+01]),
                vectors.LorentzVector([-0.6892026754834055e+00,     0.4261712129825938e-01,      0.6853484710203325e+00,    -0.3133021123175139e-01]),
                vectors.LorentzVector([-0.6498460792144278e+00,    -0.4454022198840913e+00,      0.4204895170750307e+00,     0.2111997992960589e+00]),
                vectors.LorentzVector([-0.8326129756066623e+00,    -0.4538046926431276e+00,     -0.1686174789943128e+00,    -0.6755544492202383e+00]),
                vectors.LorentzVector([-0.1123212674826061e+01,    -0.5411628118251340e+00,     -0.7789668793605800e+00,    -0.5995499352403172e+00]),
                vectors.LorentzVector([-0.9222372388674636e+00,    -0.4299804333838157e-01,      0.1734041711323346e-01,    -0.9197130003198243e+00]),
        ]),
        analytical_result =  1.19800873245882517e-5-5.68876257310534512e-5j,
        parameter_values = {'m1': 0.0, 'm2': 0.0, 'm3': 0.0, 'm4': 0.0, 'm5': 0.0,
                            'm6': 0.0, 'm7': 0.0, 'm8': 0.0, 'm9': 0.0, 'm10': 0.0,
                           },
        name = 'Decagon_P3_physical',
    ),
    entry_name = 'Decagon_P3_physical'
)

# ===================================================================================
# Tringigons with 30 external particles with masses all equal and no internal masses
# ===================================================================================

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Tringigon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.4497499305169485e+00,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.2250000000000000e+01]),
                vectors.LorentzVector([-0.4497499305169485e+00,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.2250000000000000e+01]),
                vectors.LorentzVector([ 0.1190914841688229e-01,     0.4811227902220810e-02,    -0.8341860844658933e-02,    -0.1071870295667918e+00]),
                vectors.LorentzVector([-0.2953676144156443e-02,     0.1653208835624463e-01,     0.2848734934413556e-01,    -0.1815601292160575e+00]),
                vectors.LorentzVector([-0.3076524095310454e-02,     0.1513401940618041e-01,     0.7861188448553311e-02,    -0.1145973197460721e+00]),
                vectors.LorentzVector([ 0.1826390036566392e-01,    -0.1970736165619368e-01,    -0.6152450853423360e-01,    -0.3439554953185766e+00]),
                vectors.LorentzVector([-0.9504372484440455e-02,     0.7908062483943984e-02,     0.6339866998242896e-03,    -0.9724614230247199e-01]),
                vectors.LorentzVector([-0.3995306507647221e-01,     0.1291832979641057e-01,     0.1465796670892942e-01,    -0.2346799225420517e+00]),
                vectors.LorentzVector([ 0.1397668611845013e-01,    -0.7382426010658950e-02,     0.2042622601469497e-01,    -0.1493384309094558e+00]),
                vectors.LorentzVector([ 0.2922013023842157e-02,     0.1278514318474546e-01,     0.1951856314298546e-01,    -0.1394607788531113e+00]),
                vectors.LorentzVector([ 0.2083295515366754e-01,    -0.1206148100688736e-01,     0.1196453807046953e-01,    -0.1539189328635323e+00]),
                vectors.LorentzVector([ 0.1119494688534474e-01,     0.4323911541453070e-02,     0.9659621201751422e-03,    -0.9617121835580723e-01]),
                vectors.LorentzVector([ 0.5716287720297196e-02,     0.1106557732182106e-01,    -0.3298377382829205e-01,    -0.1915758519961705e+00]),
                vectors.LorentzVector([-0.5369562427290601e-02,    -0.8116977264730580e-02,     0.1165446985452119e-01,    -0.1067174058163273e+00]),
                vectors.LorentzVector([-0.3038203642399710e-02,     0.9971571409418041e-02,     0.2249169047900925e-01,    -0.1448740020111760e+00]),
                vectors.LorentzVector([ 0.3487267645289652e-03,     0.1615148442435443e-01,    -0.8542986948582871e-02,    -0.1182132274187326e+00]),
                vectors.LorentzVector([-0.5117047914552306e-02,     0.2247503633420992e-02,    -0.1651133441512298e-01,    -0.1149847398690462e+00]),
                vectors.LorentzVector([-0.2846706995927437e-01,    -0.1362233500133967e-01,     0.4794467816646061e-01,    -0.2966325271881200e+00]),
                vectors.LorentzVector([ 0.6091584485756565e-01,     0.2195915724549356e-01,    -0.1197142659379113e-01,    -0.3376854970358458e+00]),
                vectors.LorentzVector([ 0.5980168657704429e-02,     0.2629583881848908e-03,    -0.7665766571385161e-02,    -0.8938617642223451e-01]),
                vectors.LorentzVector([-0.3342042485085184e-01,     0.1650655929454859e-01,     0.1312973979301172e-01,    -0.2113516834560396e+00]),
                vectors.LorentzVector([-0.1264861426574454e-02,     0.1915890345560194e-02,     0.1316019919940445e-01,    -0.1004317371666761e+00]),
                vectors.LorentzVector([ 0.3562605710800451e-02,    -0.7763262588943596e-02,     0.7852910865420467e-02,    -0.9481938284651482e-01]),
                vectors.LorentzVector([-0.1421836695496606e-01,    -0.7456157604864838e-02,    -0.3770019054652306e-02,    -0.1114640400271648e+00]),
                vectors.LorentzVector([-0.1310529030415713e-01,    -0.8542658344382956e-02,    -0.1598285983983683e-01,    -0.1346455953795890e+00]),
                vectors.LorentzVector([-0.1913722145628823e-01,     0.8006367042236684e-03,    -0.5783245758276338e-04,    -0.1216426743367967e+00]),
                vectors.LorentzVector([ 0.2334791001580955e-01,     0.5082287378373648e-02,    -0.3533908119628511e-02,    -0.1421656651323283e+00]),
                vectors.LorentzVector([-0.9648375412153016e-03,     0.5048990265441869e-03,    -0.5142668558714967e-02,    -0.7947214518508028e-01]),
                vectors.LorentzVector([ 0.3717510139731589e-03,    -0.4957556091851899e-01,    -0.4157498399677061e-01,    -0.3320901132800368e+00]),
                vectors.LorentzVector([ 0.2475795734193290e-03,    -0.2665308744662153e-01,    -0.3145539144342609e-02,    -0.1537321357581931e+00]),
        ]),
        analytical_result =  -3.17518049311205017e-2-1.00739776983530312e-2j,
        parameter_values = {'m1': 0., 'm2': 0., 'm3': 0., 'm4': 0., 'm5': 0.,
                          'm6': 0., 'm7': 0., 'm8': 0., 'm9': 0., 'm10': 0.,
                          'm11': 0., 'm12': 0., 'm13': 0., 'm14': 0., 'm15': 0.,
                          'm16': 0., 'm17': 0., 'm18': 0., 'm19': 0., 'm20': 0.,
                          'm21': 0., 'm22': 0., 'm23': 0., 'm24': 0., 'm25': 0.,
                          'm26': 0., 'm27': 0., 'm28': 0., 'm29': 0., 'm30': 0.,
                         },
        name = 'Tringigon_P1_physical_few_ellipses',
    ),
    entry_name = 'Tringigon_P1_physical_few_ellipses'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Tringigon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.3000000000000000e+01,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.2998332870112990e+01]),
                vectors.LorentzVector([ 0.3000000000000000e+01,     0.0000000000000000e+00,     0.0000000000000000e+00,    -0.2998332870112990e+01]),
                vectors.LorentzVector([-0.1429160394223891e+00,     0.3207485268147207e-01,    -0.5561240563105956e-01,     0.7939432277921528e-01]),
                vectors.LorentzVector([-0.2420801722880766e+00,     0.1102139223749642e+00,     0.1899156622942371e+00,    -0.1969117429437628e-01]),
                vectors.LorentzVector([-0.1527964263280961e+00,     0.1008934627078694e+00,     0.5240792299035540e-01,    -0.2051016063540302e-01]),
                vectors.LorentzVector([-0.4586073270914355e+00,    -0.1313824110412912e+00,    -0.4101633902282240e+00,     0.1217593357710928e+00]),
                vectors.LorentzVector([-0.1296615230699626e+00,     0.5272041655962656e-01,     0.4226577998828597e-02,    -0.6336248322960303e-01]),
                vectors.LorentzVector([-0.3129065633894023e+00,     0.8612219864273712e-01,     0.9771977805952947e-01,    -0.2663537671764814e+00]),
                vectors.LorentzVector([-0.1991179078792744e+00,    -0.4921617340439300e-01,     0.1361748400979665e+00,     0.9317790745633422e-01]),
                vectors.LorentzVector([-0.1859477051374817e+00,     0.8523428789830303e-01,     0.1301237542865697e+00,     0.1948008682561438e-01]),
                vectors.LorentzVector([-0.2052252438180430e+00,    -0.8040987337924907e-01,     0.7976358713646353e-01,     0.1388863676911169e+00]),
                vectors.LorentzVector([-0.1282282911410763e+00,     0.2882607694302047e-01,     0.6439747467834282e-02,     0.7463297923563160e-01]),
                vectors.LorentzVector([-0.2554344693282273e+00,     0.7377051547880703e-01,    -0.2198918255219470e+00,     0.3810858480198131e-01]),
                vectors.LorentzVector([-0.1422898744217697e+00,    -0.5411318176487053e-01,     0.7769646569680792e-01,    -0.3579708284860401e-01]),
                vectors.LorentzVector([-0.1931653360149014e+00,     0.6647714272945361e-01,     0.1499446031933950e+00,    -0.2025469094933140e-01]),
                vectors.LorentzVector([-0.1576176365583102e+00,     0.1076765628290296e+00,    -0.5695324632388580e-01,     0.2324845096859768e-02]),
                vectors.LorentzVector([-0.1533129864920616e+00,     0.1498335755613994e-01,    -0.1100755627674865e+00,    -0.3411365276368204e-01]),
                vectors.LorentzVector([-0.3955100362508266e+00,    -0.9081556667559777e-01,     0.3196311877764041e+00,    -0.1897804663951624e+00]),
                vectors.LorentzVector([-0.4502473293811277e+00,     0.1463943816366237e+00,    -0.7980951062527424e-01,     0.4061056323837710e+00]),
                vectors.LorentzVector([-0.1191815685629794e+00,     0.1753055921232606e-02,    -0.5110511047590108e-01,     0.3986779105136286e-01]),
                vectors.LorentzVector([-0.2818022446080528e+00,     0.1100437286303240e+00,     0.8753159862007812e-01,    -0.2228028323390123e+00]),
                vectors.LorentzVector([-0.1339089828889014e+00,     0.1277260230373463e-01,     0.8773466132936300e-01,    -0.8432409510496359e-02]),
                vectors.LorentzVector([-0.1264258437953531e+00,    -0.5175508392629064e-01,     0.5235273910280311e-01,     0.2375070473866967e-01]),
                vectors.LorentzVector([-0.1486187200362198e+00,    -0.4970771736576559e-01,    -0.2513346036434871e-01,    -0.9478911303310708e-01]),
                vectors.LorentzVector([-0.1795274605061186e+00,    -0.5695105562921971e-01,    -0.1065523989322456e+00,    -0.8736860202771418e-01]),
                vectors.LorentzVector([-0.1621902324490622e+00,     0.5337578028157789e-02,    -0.3855497172184225e-03,    -0.1275814763752549e+00]),
                vectors.LorentzVector([-0.1895542201764378e+00,     0.3388191585582432e-01,    -0.2355938746419007e-01,     0.1556527334387304e+00]),
                vectors.LorentzVector([-0.1059628602467737e+00,     0.3365993510294580e-02,    -0.3428445705809978e-01,    -0.6432250274768677e-02]),
                vectors.LorentzVector([-0.4427868177067157e+00,    -0.3305037394567933e+00,    -0.2771665599784707e+00,     0.2478340093154393e-02]),
                vectors.LorentzVector([-0.2049761810109241e+00,    -0.1776872496441435e+00,    -0.2097026096228406e-01,     0.1650530489462193e-02]),
        ]),
        analytical_result = 1.71112768603422270e-22-2.99165437321843644e-9j,
        parameter_values = {'m1': 0., 'm2': 0., 'm3': 0., 'm4': 0., 'm5': 0.,
                          'm6': 0., 'm7': 0., 'm8': 0., 'm9': 0., 'm10': 0.,
                          'm11': 0., 'm12': 0., 'm13': 0., 'm14': 0., 'm15': 0.,
                          'm16': 0., 'm17': 0., 'm18': 0., 'm19': 0., 'm20': 0.,
                          'm21': 0., 'm22': 0., 'm23': 0., 'm24': 0., 'm25': 0.,
                          'm26': 0., 'm27': 0., 'm28': 0., 'm29': 0., 'm30': 0.,
                         },
        name = 'Tringigon_P1_physical_many_ellipses',
    ),
    entry_name = 'Tringigon_P1_physical_many_ellipses'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Tringigon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.7495832175282474e+00,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.3750000000000000e+01]),
                vectors.LorentzVector([-0.7495832175282474e+00,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.3750000000000000e+01]),
                vectors.LorentzVector([-0.6817311405082294e-03,     0.1769549021091204e-01,    -0.2261479667294971e-01,    -0.1903959138027442e+00]),
                vectors.LorentzVector([-0.2297126915576033e-02,    -0.2486861508307056e-01,    -0.1218262344611994e-01,    -0.1868917549175284e+00]),
                vectors.LorentzVector([ 0.5552116874509821e-02,    -0.1013041668135510e-01,    -0.1565251934470179e-01,    -0.1583866081588306e+00]),
                vectors.LorentzVector([-0.2004525777719839e-01,     0.7722419512278923e-02,     0.1262863851925462e+00,    -0.6525852282935560e+00]),
                vectors.LorentzVector([ 0.3285470847103632e-02,     0.3170068757872742e-03,     0.1607483555551681e-01,    -0.1495238400877408e+00]),
                vectors.LorentzVector([-0.1318288584333414e-01,    -0.4241211672984990e-01,     0.2070012714958520e-01,    -0.2750486951738632e+00]),
                vectors.LorentzVector([ 0.8956956028529564e-01,    -0.4181923330839848e-01,     0.5297985740059507e-02,    -0.5105052168466390e+00]),
                vectors.LorentzVector([ 0.8596069193012913e-02,     0.2157747630002380e-01,    -0.1065806392368144e-01,    -0.1787508192743766e+00]),
                vectors.LorentzVector([ 0.1171718027633463e-02,     0.3119604228228752e-01,    -0.7927203338651304e-02,    -0.2038630992839849e+00]),
                vectors.LorentzVector([-0.3190379168731026e-01,     0.6827221444062047e-02,    -0.1981376649359754e-01,    -0.2281473328430751e+00]),
                vectors.LorentzVector([-0.1816471305963102e-01,     0.6790578299609332e-02,    -0.1435441309868346e-01,    -0.1737180135729676e+00]),
                vectors.LorentzVector([ 0.2565826506991252e-02,     0.1308838652596334e-01,     0.1466999843306017e-01,    -0.1595382541542142e+00]),
                vectors.LorentzVector([-0.1314912045648699e-01,     0.5885203547639661e-03,    -0.2240074501028731e-01,    -0.1802802749889105e+00]),
                vectors.LorentzVector([ 0.7336654778063456e-02,    -0.1083798819930219e-01,    -0.1633887920578607e-01,    -0.1630373780181989e+00]),
                vectors.LorentzVector([-0.4130262051206410e-01,     0.7410142225023706e-02,    -0.1633009877332162e-01,    -0.2575115917247754e+00]),
                vectors.LorentzVector([-0.1120155391289462e-02,    -0.1021679395338564e-03,    -0.1123363177845045e-02,    -0.1252524581812510e+00]),
                vectors.LorentzVector([-0.3304245278313268e-01,    -0.1507805246966791e-01,    -0.5512450605421093e-02,    -0.2221788954806924e+00]),
                vectors.LorentzVector([-0.2555567017686268e-02,     0.8913510655265457e-01,     0.5875908576282214e-01,    -0.5483891935162849e+00]),
                vectors.LorentzVector([-0.7837524382871369e-03,     0.1388228545616210e-01,     0.2670930318187453e-02,    -0.1436546193249197e+00]),
                vectors.LorentzVector([ 0.1083156678543094e+00,     0.9466799463750012e-02,    -0.8299608557180177e-01,    -0.6952563271027057e+00]),
                vectors.LorentzVector([ 0.1167901446699877e-02,    -0.3534342864696930e-02,     0.2206032585614081e-02,    -0.1268584006917978e+00]),
                vectors.LorentzVector([-0.3533679600653622e-01,    -0.2442264934702566e-01,     0.4833844494643144e-01,    -0.3466540138587216e+00]),
                vectors.LorentzVector([-0.3495602566032728e-01,     0.7406426766330819e-03,    -0.2120952087880850e-01,    -0.2396516251482828e+00]),
                vectors.LorentzVector([-0.1512068736634098e-02,     0.1454286461805387e-01,    -0.7397269291384289e-02,    -0.1494574233797544e+00]),
                vectors.LorentzVector([ 0.3856831198166811e-02,     0.2175800020561430e-01,    -0.8494837411934431e-01,    -0.4563308006632564e+00]),
                vectors.LorentzVector([ 0.1704351011754387e-01,    -0.2693882752091508e-01,     0.2728646031444713e-01,    -0.2442198078483976e+00]),
                vectors.LorentzVector([ 0.3303859085676054e-02,     0.1147321141773022e-02,    -0.4376771802157857e-02,    -0.1281003475862310e+00]),
                vectors.LorentzVector([-0.1731120789003846e-02,    -0.6374189400153728e-01,     0.4354665875627292e-01,    -0.4058120660763001e+00]),
        ]),
        analytical_result =  -3.41393026657494931e-12-4.16637655564940194e-12j,
        parameter_values = {'m1': 0., 'm2': 0., 'm3': 0., 'm4': 0., 'm5': 0.,
                          'm6': 0., 'm7': 0., 'm8': 0., 'm9': 0., 'm10': 0.,
                          'm11': 0., 'm12': 0., 'm13': 0., 'm14': 0., 'm15': 0.,
                          'm16': 0., 'm17': 0., 'm18': 0., 'm19': 0., 'm20': 0.,
                          'm21': 0., 'm22': 0., 'm23': 0., 'm24': 0., 'm25': 0.,
                          'm26': 0., 'm27': 0., 'm28': 0., 'm29': 0., 'm30': 0.,
                         },
        name = 'Tringigon_P2_physical_few_ellipses',
    ),
    entry_name = 'Tringigon_P2_physical_few_ellipses'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Tringigon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.3000000000000000e+01,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.2998332870112990e+01]),
                vectors.LorentzVector([ 0.3000000000000000e+01,     0.0000000000000000e+00,     0.0000000000000000e+00,    -0.2998332870112990e+01]),
                vectors.LorentzVector([-0.3697071195646340e+00,     0.2976193557472875e+00,    -0.5005611927042255e-01,     0.1886808370368684e+00]),
                vectors.LorentzVector([-0.1204801338912180e+00,    -0.5373888539081542e-01,    -0.3856825181225509e-01,    -0.1183574296548784e-01]),
                vectors.LorentzVector([-0.1366080963118978e+00,    -0.8728423034057854e-01,     0.3229902999905747e-01,    -0.8816424193215502e-04]),
                vectors.LorentzVector([-0.1356547589326014e+00,    -0.3297191718091518e-01,    -0.7495174996335303e-01,     0.4119831884808025e-01]),
                vectors.LorentzVector([-0.2502365923769989e+00,    -0.7058110938682028e-02,     0.6818112189080966e-01,    -0.2189060754116708e+00]),
                vectors.LorentzVector([-0.1786459353420698e+00,     0.1046992424065473e+00,    -0.2851267254157882e-01,    -0.1006949172415109e+00]),
                vectors.LorentzVector([-0.1592327297983732e+00,    -0.6421670596738141e-01,    -0.9536426408069075e-01,     0.4622698400381409e-01]),
                vectors.LorentzVector([-0.2019394989686633e+00,    -0.1195593600818853e+00,    -0.6042765215591156e-01,     0.1132855662273399e+00]),
                vectors.LorentzVector([-0.5264146891335212e+00,    -0.2698403820376176e+00,    -0.4406113756930751e+00,     0.1265736019832709e-01]),
                vectors.LorentzVector([-0.2456916193465228e+00,    -0.3896042734023340e-01,     0.2210121784814447e+00,     0.2718114586922751e-03]),
                vectors.LorentzVector([-0.3105628314950897e+00,     0.1144275545892296e-01,     0.1609748083390276e+00,     0.2457751955221848e+00]),
                vectors.LorentzVector([-0.2360632028929489e+00,     0.8579400904826549e-01,     0.1770718816736143e-01,     0.1950683963606339e+00]),
                vectors.LorentzVector([-0.1044324563436895e+00,    -0.2318951842512641e-01,    -0.1366604968458472e-01,     0.1347676738709212e-01]),
                vectors.LorentzVector([-0.1251929550496708e+00,    -0.5796954816989604e-01,     0.3143675208148838e-01,     0.3639420417613268e-01]),
                vectors.LorentzVector([-0.1315091659009482e+00,    -0.4103165102902245e-01,     0.5431993345322663e-01,    -0.5157915431094884e-01]),
                vectors.LorentzVector([-0.1533868966226102e+00,     0.1739147346848744e-01,    -0.5892606379293033e-01,    -0.9875624391390005e-01]),
                vectors.LorentzVector([-0.2538164388828844e+00,    -0.2066102767507780e+00,    -0.1171397986290710e-02,     0.1083217707338358e+00]),
                vectors.LorentzVector([-0.1555565963532279e+00,     0.5008114260575076e-01,     0.8417009009219953e-01,     0.6786110637305087e-01]),
                vectors.LorentzVector([-0.3712553969073230e+00,     0.2342895141799043e+00,     0.1626403220830896e+00,     0.2156087171495329e+00]),
                vectors.LorentzVector([-0.1386366152975742e+00,     0.3891149076364325e-01,    -0.7506591296912074e-02,    -0.8746232374468231e-01]),
                vectors.LorentzVector([-0.3007702658884222e+00,    -0.6430703231338415e-01,     0.2188884858413501e+00,    -0.1685680551104467e+00]),
                vectors.LorentzVector([-0.1384629962891667e+00,     0.5409712846398516e-01,    -0.3898246337517689e-01,    -0.6874496041549427e-01]),
                vectors.LorentzVector([-0.1614890920976794e+00,    -0.1228455600692281e-01,    -0.3539991142999072e-01,     0.1211390227011692e+00]),
                vectors.LorentzVector([-0.2583749513957764e+00,     0.2113054523365527e+00,     0.5890038969164927e-01,    -0.9294280722996563e-01]),
                vectors.LorentzVector([-0.1009843843836903e+00,    -0.3081508234576515e-02,    -0.9861971236881199e-02,     0.9544198220703085e-02]),
                vectors.LorentzVector([-0.1804323179350839e+00,    -0.1358148630147986e+00,    -0.5741109526132627e-01,    -0.2853262134108276e-01]),
                vectors.LorentzVector([-0.2500517724315000e+00,     0.3899722216643952e-01,    -0.8356044850002008e-01,    -0.2098160075063663e+00]),
                vectors.LorentzVector([-0.3044104901662132e+00,     0.7329018657682732e-01,    -0.1555222203930474e-01,    -0.2775831829639689e+00]),
        ]),
        analytical_result = 2.90177568848127685e-17-1.75774788743910202e-7j,
        parameter_values = {'m1': 0., 'm2': 0., 'm3': 0., 'm4': 0., 'm5': 0.,
                          'm6': 0., 'm7': 0., 'm8': 0., 'm9': 0., 'm10': 0.,
                          'm11': 0., 'm12': 0., 'm13': 0., 'm14': 0., 'm15': 0.,
                          'm16': 0., 'm17': 0., 'm18': 0., 'm19': 0., 'm20': 0.,
                          'm21': 0., 'm22': 0., 'm23': 0., 'm24': 0., 'm25': 0.,
                          'm26': 0., 'm27': 0., 'm28': 0., 'm29': 0., 'm30': 0.,
                         },
        name = 'Tringigon_P2_physical_many_ellipses',
    ),
    entry_name = 'Tringigon_P2_physical_many_ellipses'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Tringigon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.4497499305169484e-01,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.2250000000000000e+01]),
                vectors.LorentzVector([-0.4497499305169484e-01,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.2250000000000000e+01]),
                vectors.LorentzVector([ 0.1190914841688229e-02,     0.4811227902220809e-03,    -0.8341860844658933e-03,    -0.1071870295667918e+00]),
                vectors.LorentzVector([-0.2953676144156442e-03,     0.1653208835624463e-02,     0.2848734934413556e-02,    -0.1815601292160574e+00]),
                vectors.LorentzVector([-0.3076524095310454e-03,     0.1513401940618040e-02,     0.7861188448553309e-03,    -0.1145973197460720e+00]),
                vectors.LorentzVector([ 0.1826390036566392e-02,    -0.1970736165619368e-02,    -0.6152450853423359e-02,    -0.3439554953185766e+00]),
                vectors.LorentzVector([-0.9504372484440454e-03,     0.7908062483943983e-03,     0.6339866998242894e-04,    -0.9724614230247197e-01]),
                vectors.LorentzVector([-0.3995306507647221e-02,     0.1291832979641057e-02,     0.1465796670892942e-02,    -0.2346799225420517e+00]),
                vectors.LorentzVector([ 0.1397668611845013e-02,    -0.7382426010658949e-03,     0.2042622601469497e-02,    -0.1493384309094558e+00]),
                vectors.LorentzVector([ 0.2922013023842157e-03,     0.1278514318474545e-02,     0.1951856314298545e-02,    -0.1394607788531113e+00]),
                vectors.LorentzVector([ 0.2083295515366754e-02,    -0.1206148100688736e-02,     0.1196453807046953e-02,    -0.1539189328635323e+00]),
                vectors.LorentzVector([ 0.1119494688534474e-02,     0.4323911541453069e-03,     0.9659621201751422e-04,    -0.9617121835580723e-01]),
                vectors.LorentzVector([ 0.5716287720297196e-03,     0.1106557732182105e-02,    -0.3298377382829204e-02,    -0.1915758519961705e+00]),
                vectors.LorentzVector([-0.5369562427290601e-03,    -0.8116977264730579e-03,     0.1165446985452119e-02,    -0.1067174058163272e+00]),
                vectors.LorentzVector([-0.3038203642399709e-03,     0.9971571409418041e-03,     0.2249169047900924e-02,    -0.1448740020111760e+00]),
                vectors.LorentzVector([ 0.3487267645289651e-04,     0.1615148442435443e-02,    -0.8542986948582869e-03,    -0.1182132274187326e+00]),
                vectors.LorentzVector([-0.5117047914552305e-03,     0.2247503633420991e-03,    -0.1651133441512298e-02,    -0.1149847398690462e+00]),
                vectors.LorentzVector([-0.2846706995927436e-02,    -0.1362233500133966e-02,     0.4794467816646061e-02,    -0.2966325271881199e+00]),
                vectors.LorentzVector([ 0.6091584485756564e-02,     0.2195915724549356e-02,    -0.1197142659379113e-02,    -0.3376854970358457e+00]),
                vectors.LorentzVector([ 0.5980168657704429e-03,     0.2629583881848909e-04,    -0.7665766571385161e-03,    -0.8938617642223451e-01]),
                vectors.LorentzVector([-0.3342042485085184e-02,     0.1650655929454859e-02,     0.1312973979301171e-02,    -0.2113516834560396e+00]),
                vectors.LorentzVector([-0.1264861426574454e-03,     0.1915890345560194e-03,     0.1316019919940445e-02,    -0.1004317371666761e+00]),
                vectors.LorentzVector([ 0.3562605710800450e-03,    -0.7763262588943595e-03,     0.7852910865420466e-03,    -0.9481938284651482e-01]),
                vectors.LorentzVector([-0.1421836695496606e-02,    -0.7456157604864837e-03,    -0.3770019054652306e-03,    -0.1114640400271648e+00]),
                vectors.LorentzVector([-0.1310529030415712e-02,    -0.8542658344382956e-03,    -0.1598285983983683e-02,    -0.1346455953795890e+00]),
                vectors.LorentzVector([-0.1913722145628823e-02,     0.8006367042236683e-04,    -0.5783245758276338e-05,    -0.1216426743367967e+00]),
                vectors.LorentzVector([ 0.2334791001580955e-02,     0.5082287378373648e-03,    -0.3533908119628510e-03,    -0.1421656651323283e+00]),
                vectors.LorentzVector([-0.9648375412153015e-04,     0.5048990265441869e-04,    -0.5142668558714967e-03,    -0.7947214518508026e-01]),
                vectors.LorentzVector([ 0.3717510139731588e-04,    -0.4957556091851899e-02,    -0.4157498399677060e-02,    -0.3320901132800367e+00]),
                vectors.LorentzVector([ 0.2475795734193290e-04,    -0.2665308744662153e-02,    -0.3145539144342609e-03,    -0.1537321357581931e+00]),
        ]),
        analytical_result = 0.35864077119483945j,
        parameter_values = {'m1': 0., 'm2': 0., 'm3': 0., 'm4': 0., 'm5': 0.,
                          'm6': 0., 'm7': 0., 'm8': 0., 'm9': 0., 'm10': 0.,
                          'm11': 0., 'm12': 0., 'm13': 0., 'm14': 0., 'm15': 0.,
                          'm16': 0., 'm17': 0., 'm18': 0., 'm19': 0., 'm20': 0.,
                          'm21': 0., 'm22': 0., 'm23': 0., 'm24': 0., 'm25': 0.,
                          'm26': 0., 'm27': 0., 'm28': 0., 'm29': 0., 'm30': 0.,
                         },
        name = 'Tringigon_P1_euclidean',
    ),
    entry_name = 'Tringigon_P1_euclidean'
)

hard_coded_topology_collection.add_topology(
    create_hard_coded_topoloogy(
        'Tringigon', # Analytical result from MadLoop
        vectors.LorentzVectorList([
                vectors.LorentzVector([ 0.2998332870112990e-01,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.3000000000000000e+01]),
                vectors.LorentzVector([-0.2998332870112990e-01,     0.0000000000000000e+00,     0.0000000000000000e+00,     0.3000000000000000e+01]),
                vectors.LorentzVector([-0.5696908365050267e-02,     0.6034707142760056e-03,     0.4877422735515625e-04,    -0.5815610092502989e+00]),
                vectors.LorentzVector([-0.7673128823705573e-05,     0.4774325983672727e-03,     0.3886279207491904e-04,    -0.1108833206346892e+00]),
                vectors.LorentzVector([ 0.3571184600905674e-03,    -0.3557010403250827e-03,     0.6924145550597886e-03,    -0.1316622474623781e+00]),
                vectors.LorentzVector([-0.1266910929008660e-02,     0.1042082398896580e-02,    -0.7022482873512663e-03,    -0.2045519905847759e+00]),
                vectors.LorentzVector([ 0.1219812672295137e-03,    -0.5091595601628711e-03,     0.1467817856097107e-02,    -0.1851651194457983e+00]),
                vectors.LorentzVector([ 0.1674475137164211e-02,     0.3402326878938469e-04,     0.2152318927422069e-02,    -0.2904737739132314e+00]),
                vectors.LorentzVector([-0.1462205527979791e-02,     0.6384661230781275e-03,     0.1005399358325202e-02,    -0.2134598759976639e+00]),
                vectors.LorentzVector([-0.1178490912049724e-02,    -0.7861646778555695e-03,     0.8888425631833423e-03,    -0.1948573024707536e+00]),
                vectors.LorentzVector([-0.7634201764963458e-03,    -0.1065933382895628e-02,    -0.6112186838508179e-03,    -0.1758582560513255e+00]),
                vectors.LorentzVector([ 0.6897321690923359e-03,    -0.7887565619660570e-03,    -0.5081255122367074e-03,    -0.1534945899799452e+00]),
                vectors.LorentzVector([ 0.9549065688928160e-03,    -0.7441742800968979e-03,    -0.1853043903801153e-02,    -0.2428870853686338e+00]),
                vectors.LorentzVector([ 0.3923011244565586e-03,     0.2422331014181601e-03,     0.2335445458285408e-02,    -0.2582030700107666e+00]),
                vectors.LorentzVector([-0.9795832381710337e-03,     0.5636957969949398e-04,     0.1226654123452698e-02,    -0.1862106653391436e+00]),
                vectors.LorentzVector([ 0.8181805751421102e-03,     0.5387868493347231e-03,    -0.4237475106015999e-03,    -0.1462625268240969e+00]),
                vectors.LorentzVector([ 0.2786040034765673e-02,     0.8955021698940161e-03,     0.9655352662130770e-04,    -0.3094069455442947e+00]),
                vectors.LorentzVector([ 0.4786066805874474e-03,    -0.6599958883621299e-03,    -0.6424678117632113e-03,    -0.1441327102537679e+00]),
                vectors.LorentzVector([ 0.7146915658674555e-03,    -0.6796215097483205e-03,    -0.2263441657184843e-03,    -0.1422638784860535e+00]),
                vectors.LorentzVector([-0.3953749284791590e-03,     0.2856041461731434e-03,    -0.9683192402034328e-03,    -0.1474968885546181e+00]),
                vectors.LorentzVector([-0.1201028032679855e-03,    -0.1412927108028454e-02,    -0.2765611610570441e-02,    -0.3264872934200733e+00]),
                vectors.LorentzVector([ 0.1100574502856851e-02,    -0.4371110654864841e-03,    -0.1864167037320893e-02,    -0.2424345079179755e+00]),
                vectors.LorentzVector([ 0.1795262230327934e-02,    -0.9113763988578755e-03,     0.3950480353682807e-04,    -0.2248362525381452e+00]),
                vectors.LorentzVector([-0.6409010778289663e-03,    -0.1558254657227798e-02,    -0.2285635914259547e-03,    -0.1972600588983990e+00]),
                vectors.LorentzVector([ 0.6724202149341944e-03,    -0.1821786988428590e-03,     0.8172648592700283e-04,    -0.1221481576719928e+00]),
                vectors.LorentzVector([-0.3264569712789387e-04,     0.4166719452310339e-03,     0.7593802979656003e-03,    -0.1323381913289548e+00]),
                vectors.LorentzVector([ 0.8846940024726268e-03,     0.5197967599592516e-03,    -0.5111351113090440e-03,    -0.1521226890272480e+00]),
                vectors.LorentzVector([ 0.1489822813430809e-02,     0.1859033296997267e-02,    -0.9098200756761861e-03,    -0.2739224230483302e+00]),
                vectors.LorentzVector([-0.2080830866522650e-02,     0.2940388292912083e-02,     0.1273301105043701e-02,    -0.3949308308080461e+00]),
                vectors.LorentzVector([-0.3057596965049211e-03,    -0.4585064151705153e-03,     0.1078164614790584e-03,    -0.1146883391685992e+00]),
        ]),
        analytical_result = 3.21931745143375310e-007j,
        parameter_values = {'m1': 0., 'm2': 0., 'm3': 0., 'm4': 0., 'm5': 0.,
                          'm6': 0., 'm7': 0., 'm8': 0., 'm9': 0., 'm10': 0.,
                          'm11': 0., 'm12': 0., 'm13': 0., 'm14': 0., 'm15': 0.,
                          'm16': 0., 'm17': 0., 'm18': 0., 'm19': 0., 'm20': 0.,
                          'm21': 0., 'm22': 0., 'm23': 0., 'm24': 0., 'm25': 0.,
                          'm26': 0., 'm27': 0., 'm28': 0., 'm29': 0., 'm30': 0.,
                         },
        name = 'Tringigon_P2_euclidean',
    ),
    entry_name = 'Tringigon_P2_euclidean'
)

# Example printout
# ----------------
#hard_coded_topology_collection['DoubleTriangle'].print_topology()
#hard_coded_topology_collection['TriangleBox'].print_topology()
#hard_coded_topology_collection['DoubleBox'].print_topology()
#hard_coded_topology_collection['TriangleBoxBox'].print_topology()

# Example of yaml export and import
# ---------------------------------
#hard_coded_topology_collection['DoubleTriangle'].export_to('ExampleDoubleTriangleExport.yaml')
#test = LoopTopology.import_from('ExampleDoubleTriangleExport.yaml')
#test.print_topology()

# Example of a yaml export and import of a complete TopologyCollection
# --------------------------------------------------------------------
#hard_coded_topology_collection.export_to('ExampleTopologyCollectionExport.yaml')
#test = TopologyCollection.import_from('ExampleTopologyCollectionExport.yaml')
#test['DoubleTriange'].print_topology()

#############################################################################################################
# Define and store hyperparameters
#############################################################################################################

hyperparameters = HyperParameters({

    'General'       :   {
        # can be multiplicative, additive, cutgroups or none
        'deformation_strategy'  :   'additive',
        'topology'              :   'TriangleBoxBox',
        'numerical_threshold'   :   0.,
        # number of digits that should be the same between integrand and rotated version
        'relative_precision'    :   5.,
        'numerical_instability_check': True,
        # statistics will be printed to stderr by default
        'log_to_screen'         :   False,
        'log_file_prefix'       :   'stats/statistics',
        'integration_statistics':   True,
        'statistics_interval'   :   100000,
        'debug'                 :   0
    },

    'Integrator'    :   {
        # The integrator can be vegas or cuhre
        'integrator'        :   'vegas',
        'n_start'           :   int(1.0e6),
        'n_max'             :   int(1.0e9),
        'n_increase'        :   0,
        'n_new'             :   100000,
        'n_min'             :   2,
        'flatness'          :   50.,
        'seed'              :   1,
        'integrated_phase'  :  'imag'
    },

    'Deformation'   :   {
        # positive value: maximum lambda in auto scaling
        # negative value: no auto scaling, lambda is set to abs(lambda)
        'lambda'    :   -0.001,
        # sigma=0 means normal min. sigma large decreases steepness
        'softmin_sigma' : 0,
        # skip hyperboloid checks from lambda scaling
        'skip_hyperboloids' : False,
        'expansion_threshold'   :   0.1,

        'additive'              :   {
            # can be exponential, hyperbolic, or unity
            'mode'  :   'exponential',
            'a_ij'  :   25.0,
        },

        'multiplicative'        :   {
            'M_ij'  :   0.1
        },

        'cutgroups' : {
            'M_ij'  :   0.1,
            'mode'  :   'hyperbolic',
        }
    },

    'Parameterization'   :   {
        # rescale the radius of the ith loop momentum by rescaling^i
        'rescaling' :   1.0,
        # shift the ith loop momentum by e_cm*(i*shift)
        'shift'     :   [0., 0., 0.]

    },

})


def synchronize(root_dir = ''):
    # Synchronise the database of hard-coded topologies to the yaml data base that Rust can easily import
    hard_coded_topology_collection.export_to(os.path.join(root_dir, 'topologies.yaml'))
    hyperparameters.export_to(os.path.join(root_dir, 'hyperparameters.yaml'))

# Main synchronises yaml file to python records
if __name__ == '__main__':
    synchronize()
