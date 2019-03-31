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
        analytical_result = -1.01888646237782899E-063+6.86568168003053015E-063j,
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


# ===========================================================================================
# Decagons with 10 external particles with masses all equal to 500 GeV and no internal masses
# ===========================================================================================

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
        analytical_result = 7.96242898900458028E-015+1.58733080719071658E-005j,
        parameter_values = {'m1': 0.0, 'm2': 0.0, 'm3': 0.0, 'm4': 0.0, 'm5': 0.0,
                            'm6': 0.0, 'm7': 0.0, 'm8': 0.0, 'm9': 0.0, 'm10': 0.0,
                           },
        name = 'Decagon_P2_euclidean_massless',
    ),
    entry_name = 'Decagon_P2_euclidean_massless'
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
