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

    def add_topology(self, topology_to_add):

        assert(isinstance(topology_to_add,LoopTopology))
        if topology_to_add.name is None:
            raise BaseException("Specify the name of the hard-coded topology to add to the collection.")
        self[topology_to_add.name] = topology_to_add


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
                    end_node    = 2,
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
                    end_node    = 2,
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
                vectors.LorentzVector([10.,5.,0.,0.]),
                vectors.LorentzVector([22.,0.,6.,0.]),
                vectors.LorentzVector([-32.,-5.,-6.,0.]),
        ]),
        analytical_result =  -1.95238770805808e-32 - 2.68888267073667e-14j
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
        'deformation_strategy'  :   'generic',
        'topology'              :   'DoubleTriangle',
        'numerical_threshold'   :   0.,
        # number of digits that should be the same between integrand and rotated version
        'relative_precision'    :   5,
        'integration_statistics':   True,
        'statistics_interval'   :   10000,
        'debug'                 :   0
    },

    'Integrator'    :   {
        # The integrator can be vegas or cuhre
        'integrator'        :   'cuhre',
        'n_start'           :   int(1.0e5),
        'n_max'             :   int(1.0e9),
        'n_increase'        :   0,
        'integrated_phase'  :  'real'
    },

    'Deformation'   :   {
        # positive value: maximum lambda in auto scaling
        # negative value: no auto scaling, lambda is set to abs(lambda)
        'lambda'    :   -0.1,

        'additive'              :   {
            # can be exponential or hyperbolic
            'mode'  :   'exponential',
            'a_ij'  :   10.0,
        },

        'multiplicative'        :   {
            'M_ij'  :   0.1
        }

    },

})


def synchronize(root_dir = ''):
    # Synchronise the database of hard-coded topologies to the yaml data base that Rust can easily import
    hard_coded_topology_collection.export_to(os.path.join(root_dir, 'topologies.yaml'))
    hyperparameters.export_to(os.path.join(root_dir, 'hyperparameters.yaml'))

# Main synchronises yaml file to python records
if __name__ == '__main__':
    synchronize()
