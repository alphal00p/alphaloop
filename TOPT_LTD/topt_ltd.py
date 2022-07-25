#!/usr/bin/env python3
import argparse
import itertools
from pprint import pprint, pformat
import multiprocessing
import os
import re
pjoin = os.path.join

root_path = os.path.dirname(os.path.realpath( __file__ ))
import sys
sys.path.insert(0,pjoin(root_path,os.path.pardir))

from LTD.ltd_utils import TopologyGenerator
from cff_generator import cFF_Analyser

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('TOPT_LTD')
logger.setLevel(logging.INFO)

class TOPTLTDException(Exception):
    pass

class TOPTLTD_Analyser(object):

    _ALLOW_NUMERATORS_WHEN_SIMPLIFYING_TOPT_TERMS = False
    _N_CPUS = multiprocessing.cpu_count()
    _MERGE_TERMS = True
    _CAS = 'form' # Other option is 'sympy'
    _FORM_PATH = None
    _SPLIT_ORIENTATIONS = True

    def __init__(self, internal_edges, external_edges, name):
        self.ie = list(internal_edges)
        self.ee = list(external_edges)
        self.name = name
        self.edges = internal_edges+external_edges
        self.iv = sorted(list(set(sum([[a,b] for a,b in self.ie],[]))))
        self.ie_map = { v : [] for v in self.iv }
        for i_e, (a,b) in enumerate(self.ie):
            self.ie_map[a].append( (b,1,i_e) )
            self.ie_map[b].append( (a,-1,i_e) )
        self.ev = {}
        for i_ee, (k, v) in enumerate(self.ee):
            if k in self.iv:
                self.ev[k] = [ (v,-1,i_ee), ]
                raise TOPTLTDException("Define all external edges as incoming.")
            elif v in self.iv:
                self.ev[v] = [ (k,1,i_ee), ]
            else:
                raise TOPTLTDException("External edges not connected to any internal vertex.")
        self.vertices = sorted(list(set(list(self.iv)+[v_infos[0] for all_v_infos in self.ev.values() for v_infos in all_v_infos])))
        for v in self.vertices:
            if v not in self.ev:
                self.ev[v] = []

        self.topo = TopologyGenerator( 
            [ ('q%d'%i, v1, v2) for i, (v1, v2) in enumerate(self.ie) ] + 
            [ ('p%d'%i, v1, v2) for i, (v1, v2) in enumerate(self.ee) ]
        )

        self.loop_momentum_bases = self.topo.loop_momentum_bases()
        self.topo.generate_momentum_flow( loop_momenta = self.loop_momentum_bases[0] )
        self.signature_map = self.topo.get_signature_map()
        # This will be generated when needed.
        self.connected_cluster_to_e_surface_map = None
        self.e_surface_to_connected_cluster_map = None

        # Cache the map from E_surf in the OS energy representation to the node representation
        # when generating the expression
        self.E_surf_to_node_repr_dict = {}

    def generate_orientation(self, topt_ordering):

        #orientation = []
        # More efficient encoding below
        orientation = [None,]*len(self.ie)
        past_vertices = []
        for internal_vertex in topt_ordering:
            for connected_vertex, o, i_e in self.ie_map[internal_vertex]:
                if connected_vertex in past_vertices:
                    #orientation.append( (connected_vertex, internal_vertex) )
                    if orientation[i_e] is not None:
                        raise TOPTLTDException("Edge #%d traversed twice in TOPT ordering: %s"%(i_e, str(topt_ordering)))
                    orientation[i_e] = o
            past_vertices.append(internal_vertex)

        if any(o is None for o in orientation):
            raise TOPTLTDException("Not all edges traversed for TOPT ordering: %s"%(i_e, str(topt_ordering)))

        return tuple(orientation)

    def get_topt_terms(self):

        topt_orientations = {}
        for p in itertools.permutations(self.iv):
            orientation = self.generate_orientation(p)
            if orientation in topt_orientations:
                topt_orientations[orientation].append(p)
            else:
                topt_orientations[orientation] = [p,]
        return topt_orientations

    def get_E_surfaces_for_topt_ordering(self, t):
        
        past_vertices = []
        accumulated_shifts = []
        accumulated_energies = []
        E_surfaces = []
        # Do not add the E-surface embedding the whole graph:
        for internal_vertex in t[:-1]:
            accumulated_shifts.extend([i_ee for _,_,i_ee in self.ev[internal_vertex]])
            # Account for energy-momentum conservation
            while all(ext in accumulated_shifts for ext in range(len(self.ee))):
                for ext in range(len(self.ee)):
                    del accumulated_shifts[accumulated_shifts.index(ext)]
            for connected_vertex, o, i_e in self.ie_map[internal_vertex]:
                if connected_vertex in past_vertices:
                    try:
                        accumulated_energies.remove(i_e)
                    except ValueError:
                        raise TOPTLTDException("Could not find vertex #%d in list even thought it "%i_e+
                        "should already have been encountered when processing TOPT ordering: %s"%(str(t)))
                else:
                    accumulated_energies.append(i_e)
            past_vertices.append(internal_vertex)
            E_surfaces.append( (tuple(sorted(accumulated_energies)),tuple(sorted(accumulated_shifts))) )

            if E_surfaces[-1] not in self.E_surf_to_node_repr_dict:
                # Canonicalise the node representation by always picking the one with the lowest internal vertex id
                if self.iv[0] in past_vertices:
                    self.E_surf_to_node_repr_dict[E_surfaces[-1]] = set(past_vertices)
                else:
                    self.E_surf_to_node_repr_dict[E_surfaces[-1]] = set(self.iv)-set(past_vertices)

            #E_surfaces.append( (tuple(sorted(accumulated_energies)), tuple([]) ) )

        return tuple(E_surfaces)

    def merge_topt_terms_with_numerators_sympy(self, topt_terms, **opts):
        
        import sympy as sp

        if len(topt_terms) == 1:
            return topt_terms

        Es = sp.Matrix(len(self.ie)+1,1, lambda i,j: sp.Symbol('E%d' % (i)))
        ps = sp.Matrix(len(self.ev)+1,1, lambda i,j: sp.Symbol('p%d' % (i)))
        psToSubst = [ ps[i_e] for i_e in range(len(self.ee)-1) ]
        psToSubst.append( -sum(ps[i_e] for i_e in range(len(self.ee)-1)) )

        etas = {}
        for num, denom in topt_terms:
            for E_surf in denom:
                if E_surf not in etas:
                    etas[E_surf] = (sum(Es[i_e] for i_e in E_surf[0]) if len(E_surf[0])>0 else 0) + \
                                (sum(psToSubst[i_e] for i_e in E_surf[1]) if len(E_surf[1])>0 else 0)

        numerator = 0
        for num, denom in topt_terms:
            numerator_piece = 1
            for eta, eta_sp_expr in etas.items():
                if eta not in denom:
                    numerator_piece *= eta_sp_expr
            numerator += numerator_piece
        
        # Now attempt to divide the combined numerator by all possible denominator piece:
        denominator = []
        for eta, eta_sp_expr in etas.items():
            coefs, remainder = sp.reduced(numerator,[eta_sp_expr,])
            coef = coefs[0]
            #coef, remainder = sp.div(numerator, eta_sp_expr, domain ='QQ')
            if remainder == 0:
                numerator = coef
            else:
                denominator.append(eta)

        # Cast the numerator as a decomposition over the E_surfaces in the denominator
        # PROBLEM: E-surfaces in the denominator are not linearly independent so you don't necessarily get a zero remainder
        # denominator = sorted(denominator)
        # coefs, remainder = sp.reduced(numerator,[etas[eta] for eta in denominator])
        # etaVars = sp.Matrix(len(denominator)+1,1, lambda i,j: sp.Symbol('eta_%d' % (i)))
        # numerator = sum( c*etaVars[i_c] for i_c, c in enumerate(coefs) ) + remainder

        num = ( numerator.expand(), tuple([t[0][1][0] for t in topt_terms]) ) 

        return [ [ num, set(denominator) ], ]

    def merge_topt_terms_with_numerators_form(self, topt_terms, **opts):
        
        import form

        if len(topt_terms) == 1:
            return topt_terms

        Es = ['E%d'%i for i in range(len(self.ie))]
        #ps = ['p%d'%i for i in range(len(self.ee))]
        ps = ['p%d'%i for i in range(len(self.ee)-1)]
        ps.append( '(-%s)'%('-'.join(ps[i_e] for i_e in range(len(self.ee)-1))) )
        etas = {}
        for num, denom in topt_terms:
            for E_surf in denom:
                if E_surf not in etas:
                    etas[E_surf] = '+'.join([Es[i_e] for i_e in E_surf[0]]+[ps[i_e] for i_e in E_surf[1]])
        sorted_etas = sorted(list(etas.items()), key=lambda el:el[0])
        numerator = []
        for num, denom in topt_terms:
            numerator_piece = []
            for eta, eta_expr in sorted_etas:
                if eta not in denom:
                    numerator_piece.append('g(%s)'%eta_expr)
            numerator.append('(%s)'%('*'.join(numerator_piece)))
        numerator = '+'.join(numerator)

        # Now attempt to divide the combined numerator by all possible denominator piece:
        denominator = []
        FORM_code = '''#: MaxTermSize 100M
#: WorkSpace 1G
Symbol E0,...,E%(nEs)d;
Symbol p0,...,p%(nPs)d;
Symbol x,y,z,e,markerNum,markerRes;
CFunction f,g,NumTracker,TestRes;

Local n = %(numerator)s;
#do i=1,1
id once g(e?) = e;
if (count(g,1)) redefine i "0";
.sort
#enddo
.sort
Hide n;
Local E = TestRes()*NumTracker(n)*f(n,1);

%(tests)s

.sort
PolyRatFun;
id ifmatch->finalskip NumTracker(y?)*f(z?,1) = NumTracker(z);
id TestRes(?a)*NumTracker(y?)*f(z?,x?) = TestRes(?a,%(last_eta_id)d)*NumTracker(y);
label finalskip;
.sort
Local Num = markerNum*E;
Local Res = markerRes*E;
id markerNum*NumTracker(y?)*TestRes(?a) = y;
id markerRes*NumTracker(y?)*TestRes(?a) = TestRes(?a);
.sort'''%{
            'nEs' : len(self.ie)-1,
            'nPs' : len(self.ev)-1,
            'numerator' : numerator,
            'last_eta_id' : (len(etas)-1),
            'tests' : '\n\n\n'.join([
                '\n'.join([
                    '.sort',
                    'PolyRatFun;',
                    'id ifmatch->skip%d NumTracker(y?)*f(z?,1) = NumTracker(z)*f(z,%s);'%(i_eta,eta_expr),
                    'id TestRes(?a)*NumTracker(y?)*f(z?,x?) = TestRes(?a,%d)*NumTracker(y)*f(y,%s);'%(i_eta-1,eta_expr),
                    'label skip%d;'%i_eta,
                    '.sort',
                    'PolyRatFun f;'
                ])
                for i_eta, (eta, eta_expr) in enumerate(sorted_etas)
            ])
        }
        #print('\n%s\n'%FORM_code)
        with (form.open() if self._FORM_PATH is None else form.open(self._FORM_PATH)) as f:
            try:
                f.write(FORM_code)
                num = f.read('Num')
                surviving_eta_indices = [int(i_str) for i_str in f.read('Res')[8:-1].split(',')]
                numerator = ( num, tuple([t[0][1][0] for t in topt_terms]) ) 
                denominator = [sorted_etas[i][0] for i in surviving_eta_indices]
            except form.formlink.FormError as e:
                print("Error when running FORM.\nError:\n%s\nFORM code:\n%s"%(str(e),FORM_code))
                raise e
            #print(num)
            #print(surviving_eta_indices)
            #stop

        return [ [ numerator, set(denominator) ], ]

    def merge_topt_terms_without_numerators(self, topt_terms, debug=False):
        
        orderings_representatives = { orderings: orderings[0] for (num, orderings), terms in topt_terms}

        while(len(topt_terms) > 1):
            
            all_combinations = list(itertools.combinations( list(range(len(topt_terms))), 2 ))

            # Optimise the sorting of the next step
            # all_combinations = []
            # topt_terms.sort(key=lambda el: el[0][1])
            # for i_a in range(len(topt_terms)):
            #     for i_b in range(len(topt_terms)-1,i_a,-1):
            #         a, b = topt_terms[i_a], topt_terms[i_b]
            #         non_common_e_surfs =  list(a[1].symmetric_difference(b[1]))
            #         if len(non_common_e_surfs) != 2:
            #             continue
            #         # Make sure the same edge does not appear twice in the merged term
            #         if len(set(non_common_e_surfs[0][0]).intersection(set(non_common_e_surfs[1][0])))>0:
            #             continue
            #         all_combinations.append( (i_a, i_b) )

            # # Now sort the combinations by exploring those orderings with the least number of matching first indices.
            # def get_index_of_first_different_element(listA, listB):
            #     for i,(a,b) in enumerate(zip(listA,listB)):
            #         if a!=b:
            #             return i
            #     return len(listA)+1
            # all_combinations.sort(key = lambda el: get_index_of_first_different_element( 
            #     orderings_representatives[topt_terms[el[0]][0][1]], 
            #     orderings_representatives[topt_terms[el[1]][0][1]],
            # ) ,reverse=True)


            possible_cancellations = []
            for i_a, i_b in all_combinations:
                a, b = topt_terms[i_a], topt_terms[i_b]
                if len(a[1])!=len(b[1]):
                    raise TOPTLTDException("All terms generated in merge_topt_terms_with_numerators should have the same number of E-surfaces in the denominator.")
                non_common_e_surfs =  list(a[1].symmetric_difference(b[1]))
                if len(non_common_e_surfs) != 2:
                    continue
                # Make sure the same edge does not appear twice in the merged term
                if len(set(non_common_e_surfs[0][0]).intersection(set(non_common_e_surfs[1][0])))>0:
                    continue
                
                merged_OSEs = tuple(sorted(list(non_common_e_surfs[0][0])+list(non_common_e_surfs[1][0])))
                merged_externals = list(non_common_e_surfs[0][1])+list(non_common_e_surfs[1][1])
                # Account for energy-momentum conservation
                while all(ext in merged_externals for ext in range(len(self.ee))):
                    for ext in range(len(self.ee)):
                        del merged_externals[merged_externals.index(ext)]
                merged_externals = tuple(sorted(merged_externals))
                merged_e_surf_in_num = ( merged_OSEs, merged_externals )
                common_e_surfs = list(a[1].intersection(b[1]))
                if merged_e_surf_in_num in common_e_surfs:
                    possible_cancellations.append( (merged_e_surf_in_num,common_e_surfs,non_common_e_surfs,(i_a, i_b)) )
                    # Add the ordering representative for the merged term
                    orderings_representatives[ tuple(list(a[0][1])+list(b[0][1])) ] = min( list(a[0][1])+list(b[0][1]) )
                    # Uncomment below to always apply the first cancellation found. This will yield a worse structure!!
                    #break

            if len(possible_cancellations) == 0:
                break

            # if len(possible_cancellations) > 1:
            #     for i_pc, ( merged_e_surf_in_num,common_e_surfs,non_common_e_surfs,(i_a, i_b) ) in enumerate(possible_cancellations):
            #         print("Possible cancellation #%d:"%i_pc)
            #         print("   merged_e_surf_in_num              = %s"%str(merged_e_surf_in_num))
            #         print("   common_e_surfs                    = %s"%str(common_e_surfs))
            #         print("   non_common_e_surfs                = %s"%str(non_common_e_surfs))
            #         print("   i_a, i_b, ordering_a, ordering_b  = %d, %d, %s, %s"%(i_a, i_b, str(topt_terms[i_a][0][1]), str(topt_terms[i_b][0][1])))
            #     stop

            # Pick the cancellation that corresponds to merging terms with orderings 
            cancellation_index_to_apply = 0
            current_first_index_differing = len(self.iv)
            for i_cancellation, (merged_e_surf_in_num,common_e_surfs,non_common_e_surfs,(i_a, i_b)) in enumerate(possible_cancellations):
                # Check what is the first index that
                this_first_index_differing = None
                for i_index, (v1, v2) in enumerate(zip(orderings_representatives[topt_terms[i_a][0][1]], orderings_representatives[topt_terms[i_b][0][1]] )):
                    if v1!=v2:
                        this_first_index_differing = i_index
                        break
                if this_first_index_differing is None:
                    this_first_index_differing = -1
                if this_first_index_differing < current_first_index_differing:
                    current_first_index_differing = this_first_index_differing
                    cancellation_index_to_apply = i_cancellation

            merged_e_surf_in_num,common_e_surfs,non_common_e_surfs,(i_a, i_b) = possible_cancellations[cancellation_index_to_apply]
            topt_terms[:] = [ [ (1, tuple(list(topt_terms[i_a][0][1])+list(topt_terms[i_b][0][1])) ) , (set(common_e_surfs) - {merged_e_surf_in_num,}).union(set(non_common_e_surfs)) ] ] + [ t for i_t, t in enumerate(topt_terms) if i_t not in [i_a, i_b] ]
          
        return topt_terms

    def merge_topt_terms_helper(self, args, **opts):

        orientation, topt_terms = args

        # Store the on-shell energy polynomial in the numerator and list of E-surfaces (linear polynomials) in the denominator.
        numerator_and_E_surfaces_per_term = [ [ (1, (tuple(t),) ), set(self.get_E_surfaces_for_topt_ordering(t)) ] for t in topt_terms ]

        # Iteratively merge terms
        if self._MERGE_TERMS:
            if self._ALLOW_NUMERATORS_WHEN_SIMPLIFYING_TOPT_TERMS:
                if self._CAS == 'sympy':
                    numerator_and_E_surfaces_per_term = self.merge_topt_terms_with_numerators_sympy(numerator_and_E_surfaces_per_term,**opts)
                else:
                    numerator_and_E_surfaces_per_term = self.merge_topt_terms_with_numerators_form(numerator_and_E_surfaces_per_term,**opts)
            else:
                numerator_and_E_surfaces_per_term = self.merge_topt_terms_without_numerators(numerator_and_E_surfaces_per_term,**opts)#, debug=(orientation==(-1, 1, -1, 1)))
        return orientation, tuple( ( nE[0], tuple(sorted(list(nE[1]))) ) for nE in numerator_and_E_surfaces_per_term )

    def process_topt_terms(self, topt_terms):

        processed_terms = {}
        logger.info("Processing TOPT terms on %d CPUS and %sconsidering numerators for merging..."%(self._N_CPUS, 'NOT ' if not self._ALLOW_NUMERATORS_WHEN_SIMPLIFYING_TOPT_TERMS else ''))
        
        if self._CAS=='form' and self._MERGE_TERMS and self._ALLOW_NUMERATORS_WHEN_SIMPLIFYING_TOPT_TERMS and not os.path.isfile(os.path.join(root_path,'form.set')):
            with open(os.path.join(root_path,'form.set'),'w') as f:
                f.write(
"""MaxTermSize=10M
WorkSpace=5G
SmallSize=5G
SubTermsInSmall=10M""")

        if self._N_CPUS > 1:
            with multiprocessing.Pool(self._N_CPUS) as p:
                results = p.imap_unordered(self.merge_topt_terms_helper, topt_terms.items())
                for orientation, merged_terms in results:
                    processed_terms[orientation] = merged_terms
                    print("Processing orientation #%d/%d"%(len(processed_terms),len(topt_terms)),end='\r')
        else:
            for o, t in topt_terms.items():
                #if o != (-1, 1, 1, -1, 1, -1, -1):
                #    continue
                orientation, merged_terms = self.merge_topt_terms_helper([o,t])#,debug=True)
                processed_terms[orientation] = merged_terms
                #logger.info("Processing orientation #%d/%d"%(len(processed_terms),len(topt_terms)),end='\r')
        return processed_terms

    def generate_mathematica_cFF_code(self, topt_ltd_terms):

        repl_dict = {}

        repl_dict['name'] = self.name.replace('-','').replace('_','')

        externals_map = {}
        for e_name, sig in self.signature_map.items():
            if not e_name.startswith('p'):
                continue
            if not sum(abs(s) for s in sig[1]) == 1:
                raise TOPTLTDException("Incorrect external momentum signature: %s"%sig)
            try:
                ext_index = sig[1].index(1)
            except ValueError:
                raise TOPTLTDException("Incorrect external momentum signature: %s"%sig)
            externals_map[ext_index] = int(e_name[1:])

        repl_dict['topt_ltd_terms'] = ',\n'.join(
            '<|"Orientation"->%s,"Terms"->{%s}|>'%(
                '{%s}'%(','.join('%d'%i for i in tlt[0])),
                ','.join('%s'%(
                '<|"Num"->(%s),"Orderings"->(%s),"Esurfs"->%s|>'%(
                    re.sub(r'p(\d+)',r'p\1E',re.sub(r'E(\d+)',r'OSE[\1]',str(t[0][0]))),
                    '{%s}'%(','.join(
                        '{%s}'%(','.join( '%d'%o for o in ordering )
                        ) for ordering in t[0][1]
                    )),
                    '{%s}'%(','.join(
                        '<|"OSE"->{%s},"shifts"->{%s}|>'%(
                            ','.join('%d'%i_e for i_e in eta[0]),
                            ','.join('%d'%i_s for i_s in eta[1]),
                        ) for eta in t[1]
                    ))
                )
            ) for t in tlt[1]))
            for tlt in topt_ltd_terms
        )
        repl_dict['topt_OSE_replacement'] = ',\n'.join(
            'OSE[%(id)d] :> Sqrt[(%(lm_str)s).(%(lm_str)s)]'%{
               'id': int(e_name[1:]),
               'lm_str' : ''.join( [ '%sk[[%d]]'%('+' if s > 0 else '-', i_k+1) for i_k, s in enumerate(sig[0]) if s!=0 ]+
                                   [ '%sp[[%d]]'%('+' if s > 0 else '-', externals_map[i_e]+1) for i_e, s in enumerate(sig[1]) if s!=0 ] )
            } for e_name, sig in self.signature_map.items() if e_name.startswith('q')
        )

        repl_dict['shifts_replacement'] = ','.join('pE[%d]:>p%dE'%(i,i) for i in range(len(externals_map)))
        repl_dict['loop_momenta'] = ','.join('k%d'%i for i in range(self.topo.n_loops))
        repl_dict['ext_momenta'] = ','.join('p%d'%i for i in range(len(externals_map)))
        repl_dict['n_internal_edges'] = len(self.ie)
        repl_dict['n_loops'] = self.topo.n_loops
        return \
"""%(name)sTOPTLTDTerms = {%(topt_ltd_terms)s};
ComputeOSEReplacement[k_,p_]:={
%(topt_OSE_replacement)s
};
ComputeTOPTLTD[terms_,NInternalEdges_,numerics_]:=Module[
{ResPerOrientations,OSEReplacement,ShiftsReplacement},
OSEReplacement = ComputeOSEReplacement[{%(loop_momenta)s},{%(ext_momenta)s}]/.numerics;
ShiftsReplacement = {%(shifts_replacement)s}/.numerics;
ResPerOrientations=Table[
Total[Table[
((st["Num"]/.OSEReplacement)/.numerics)/
(Times@@Table[
((Total[Table[OSE[e],{e,eta["OSE"]}]]/.OSEReplacement)+
(Total[Table[pE[s],{s,eta["shifts"]}]]/.ShiftsReplacement))
,{eta, st["Esurfs"]}])
,{st,t["Terms"]}
]]
,{t,terms}
];

( (Complex[0, -1])^(%(n_loops)d) / ((Times@@Table[-2*OSE[ie],{ie,0,NInternalEdges-1}])/.OSEReplacement) )*Total[ResPerOrientations]

];

%(name)scFFNum = ComputeTOPTLTD[%(name)sTOPTLTDTerms,%(n_internal_edges)d,NumericalSample];
Print["%(name)scLTDNum = ", N[%(name)scLTDNum,40]//FullForm];
Print["%(name)scFFNum = ", N[%(name)scFFNum,40]//FullForm];
Print["%(name)sRepresentationsRatio = ", N[%(name)scFFNum/%(name)scLTDNum,40]//FullForm];
"""%repl_dict

    def generate_mathematica_cLTD_code(self):

        externals_map = {}
        for e_name, sig in self.signature_map.items():
            if not e_name.startswith('p'):
                continue
            if not sum(abs(s) for s in sig[1]) == 1:
                raise TOPTLTDException("Incorrect external momentum signature: %s"%sig)
            try:
                ext_index = sig[1].index(1)
            except ValueError:
                raise TOPTLTDException("Incorrect external momentum signature: %s"%sig)
            externals_map[ext_index] = e_name

        repl_dict = {}
        repl_dict['name'] = self.name.replace('-','').replace('_','')
        repl_dict['MM_props'] = '*'.join('\nprop[%s,0](*%s:(%d,%d)*)'%(
                ''.join( [ '%sk%d'%('+' if s > 0 else '-', i_k) for i_k, s in enumerate(sig[0]) if s!=0 ]+
                         [ '%s%s'%('+' if s > 0 else '-', externals_map[i_e]) for i_e, s in enumerate(sig[1]) if s!=0 ] ),
                e_name, self.ie[int(e_name[1:])][0], self.ie[int(e_name[1:])][1]
            ) for e_name, sig in sorted(list(self.signature_map.items()),key=lambda el: el[0]) if e_name.startswith('q'))
        repl_dict['MM_loopmoms'] = ','.join('k%d'%i for i in range(self.topo.n_loops))

        repl_dict['max_rndm'] = self.topo.n_loops*3+(len(externals_map)-1)*4
        repl_dict['num_loop_momenta'] = ',\n'.join('k%d -> {%s}'%(i,','.join('RN[[%d]]'%(i*3+j+1) for j in range(3))) for i in range(self.topo.n_loops))
        repl_dict['ext_momenta'] = ',\n'.join('p%d -> {%s}'%(i,','.join('RN[[%d]]'%(self.topo.n_loops*3+i*3+j+1) for j in range(3))) for i in range(len(externals_map)-1))
        repl_dict['ext_energies'] = ',\n'.join('p%dE -> %s'%(i,'RN[[%d]]'%((self.topo.n_loops+len(externals_map)-1)*3+1+i)) for i in range(len(externals_map)-1))
        repl_dict['energy_repl'] = '{%s}'%(','.join('p%d[0] :> p%dE'%(i,i) for i in range(len(externals_map))))
        repl_dict['mom_conservation'] = 'p%dE->-(%s),p%d->-(%s)'%(
            len(externals_map)-1, '+'.join('p%dE'%i for i in range(len(externals_map)-1)),
            len(externals_map)-1, '+'.join('p%d'%i for i in range(len(externals_map)-1))
        )

        return \
"""(* Loading of the cLTD package *)

cLTDPATH="<PASTE_HERE_YOUR_CLTD_PATH>"(*"/Users/vjhirsch/HEP_programs/cLTD/cLTD.m"*);
FORMPATH="<PASTE_HERE_YOUR_FORM_PATH>"(*"/Users/vjhirsch/HEP_programs/form/bin/form"*);
Get[cLTDPATH];


(* cLTD representation *)

%(name)scLTD = cLTD[%(MM_props)s, loopmom -> {%(MM_loopmoms)s}, EvalAll -> True, "FORMpath" -> FORMPATH];
RNSeed=1;
RN = Table[Prime[i]/Prime[i+1],{i,RNSeed,%(max_rndm)d+RNSeed-1}];
NumericalSample={
%(num_loop_momenta)s,
%(ext_momenta)s,
%(ext_energies)s
};
NumericalSample=Join[NumericalSample,({%(mom_conservation)s}/.NumericalSample)];
%(name)scLTDNum = ((%(name)scLTD[[1]]/.%(energy_repl)s)/.(%(name)scLTD[[2]]/.NumericalSample))/.NumericalSample;
Print["%(name)scLTDNum = ", N[%(name)scLTDNum,40]//FullForm];
"""%repl_dict

    def analyze(self):
        analysis = {
            'topt_ltd_terms' : [],
            'cLTD_mathematica' : None,
            'cFF_mathematica' : None,
        }

        analysis['cLTD_mathematica'] = self.generate_mathematica_cLTD_code()

        # Create all topt orderings
        topt_terms = self.get_topt_terms()

        if not self._SPLIT_ORIENTATIONS:
            logger.info("Merging all orientations within a single list of terms that will be merged. This *will* induce finite differences of numerators then!")
            topt_terms = {(0,) : sum(list(topt_terms.values()),[])}

        logger.info("Original TOPT expression has %d orientations and %d terms, with a maximum of %d terms per orientation."%(
            len(topt_terms.keys()), sum(len(v) for v in topt_terms.values()), max(len(v) for v in topt_terms.values())
        ))

        # Now merge TOPT terms
        processed_topt_terms = self.process_topt_terms(topt_terms)
        logger.info("After merging, the TOPT-LTD expression has %d orientations and %d terms, with a maximum of %d terms per orientation."%(
            len(processed_topt_terms.keys()), sum(len(v) for v in processed_topt_terms.values()), max(len(v) for v in processed_topt_terms.values())
        ))

        analysis['topt_ltd_terms'] = sorted( list(processed_topt_terms.items()), key=lambda el: el[0] )

        logger.info("Maximal number of square roots in E-surfaces in the denominator of any term: %d"%(
            max( max( max(len(d[0]) for d in denom) for num, denom in t) for o,t in analysis['topt_ltd_terms'] )
        ))
        logger.info("Maximal number of E-surfaces in the denominator of any term: %d"%(
            max( max( len(denom) for num, denom in t) for o,t in analysis['topt_ltd_terms'] )
        ))

        analysis['cFF_mathematica'] = self.generate_mathematica_cFF_code(analysis['topt_ltd_terms'])

        return analysis

    def span_graph_from_node(self, seed_node, veto_edges):

        current_graph = { seed_node }
        for i_e, (a,b) in enumerate(self.ie):
            if i_e in veto_edges:
                continue
            if seed_node==a:
                current_graph = current_graph.union(self.span_graph_from_node(b, veto_edges+[i_e,]))
            elif seed_node==b:
                current_graph = current_graph.union(self.span_graph_from_node(a, veto_edges+[i_e,]))
        return current_graph

    def convert_E_surface_into_node_representation(self, E_surface):
        
        if self.e_surface_to_connected_cluster_map is not None and E_surface[0] in self.e_surface_to_connected_cluster_map:
            return self.e_surface_to_connected_cluster_map[ E_surface[0] ]

        DEBUG_THIS_FUNCTION = False

        if (E_surface not in self.E_surf_to_node_repr_dict) or DEBUG_THIS_FUNCTION:

            #raise TOPTLTDException("E-surface '%s' was never encountered."%str(E_surface))

            # We need to generate it and add it to the cache then
            # Pick the first edge being cut and explore all nodes accessible from either side.
            # The sum must cover the whole graph
            node_representation = []
            complete_left_graph, complete_right_graph = {}, {}
            for i_OSE, OSE in enumerate(E_surface[0]):
                seed_node_left, seed_node_right = self.ie[OSE]
                left_graph = self.span_graph_from_node(seed_node_left, list(E_surface[0]))
                right_graph = self.span_graph_from_node(seed_node_right, list(E_surface[0]))
                if left_graph.intersection(right_graph) != set([]):
                        raise TOPTLTDException("The E-surface %s does not separate two disconnected regions."%str(E_surface))
                if min(left_graph) > min(right_graph):
                    left_graph, right_graph = right_graph, left_graph

                if i_OSE == 0:
                    complete_left_graph, complete_right_graph = left_graph, right_graph
                    if left_graph.union(right_graph) == set(self.iv):
                        node_representation.append( tuple(sorted(list(left_graph))) )
                        break
                    else:
                        node_representation.extend( [ tuple(sorted(list(left_graph))), tuple(sorted(list(right_graph))) ] )
                else:
                    node_representation.extend( [ tuple(sorted(list(left_graph))), tuple(sorted(list(right_graph))) ] )

                    complete_left_graph = complete_left_graph.union(left_graph)
                    complete_right_graph = complete_right_graph.union(right_graph)

            if complete_left_graph.union(complete_right_graph) != set(self.iv):
                raise TOPTLTDException("Summing both sides of E-surface %s does not span the whole graph."%str(E_surface))

            node_representation = tuple(sorted(list(set(node_representation))))

            if DEBUG_THIS_FUNCTION and E_surface in self.E_surf_to_node_repr_dict:
                if node_representation != self.E_surf_to_node_repr_dict[E_surface]:
                    raise TOPTLTDException("Inconsistent reconstruction of the node representation of E-surface %s:\n Reconstructed: %s\n Original: %s"%(
                        str(E_surface), str(node_representation), str(self.E_surf_to_node_repr_dict[E_surface])
                    ))
            self.E_surf_to_node_repr_dict[E_surface] = node_representation

        return self.E_surf_to_node_repr_dict[E_surface]

    def generate_connected_cluster_to_e_surface_map(self, cff_analysis, canonicalise = True):

        connected_cluster_to_e_surface_map = {}

        connected_clusters = cff_analysis['connected_clusters']
        
        for cc in connected_clusters:
            OSEs = []
            for node in cc:
                for other_end, _, i_e in self.ie_map[node]:
                    if other_end not in cc:
                        OSEs.append(i_e)
            nodes = set(list(cc))
            if canonicalise:
                # Canonicalize the list of not by taking whichever complement has the smallest node id in it.
                if self.iv[0] not in nodes:
                    nodes = set(self.iv)-nodes
            connected_cluster_to_e_surface_map[ ( tuple( sorted(list(nodes)) ),) ] = tuple(sorted(list(set(list(OSEs)))))


        return connected_cluster_to_e_surface_map

    def is_a_cFF(self, family):

        # First check that the intersection of any two subset of the family is either one of the subset or the empty set
        for i, a in enumerate(family):
            for b in family[i+1:]:
                if a.intersection(b) not in [a,b,set([])]:
                    return False

        # Then check that no set can be written as the union of any number of other sets.
        all_possible_unions = []
        for n_subsets in range(2,len(family)):
            for comb in itertools.combinations(family,n_subsets):
                if set().union(*comb) not in comb:
                    all_possible_unions.append(set().union(*comb))

        for f in family:
            if f in all_possible_unions:
                return False

        return True

    def test_cff_property(self, topt_analysis, cff_analysis, verbosity=1, full_analysis=False):
        
        E_surfaces_violating_connected_cluster = []
        families_violating_cFF = []
        distinct_E_surfaces = set([])
        n_topt_ltd_terms = 0
        first_shown = False
        for o, terms in topt_analysis['topt_ltd_terms']:
            for num, etas in terms:
                n_topt_ltd_terms += 1
                if any(eta[0] not in self.e_surface_to_connected_cluster_map_non_canonicalised for eta in etas):
                    E_surfaces_violating_connected_cluster.extend( [ eta for eta in etas if eta[0] not in self.e_surface_to_connected_cluster_map_non_canonicalised ] )
                    if not first_shown:
                        first_shown = True
                        logger.info("The following E-surface is violating the connected cluster constraints:  %s | %s"%( 
                            str(E_surfaces_violating_connected_cluster[0]), str(self.convert_E_surface_into_node_representation(E_surfaces_violating_connected_cluster[0])) ))
                        logger.info("It appeared in orientation %s with the following terms:\n%s"%(o, '\n'.join( pformat(t) for t in terms )))
                    continue
                else:
                    for eta in etas:
                        distinct_E_surfaces.add(self.e_surface_to_connected_cluster_map_non_canonicalised[eta[0]])
                potential_cFF = tuple([ set(self.e_surface_to_connected_cluster_map_non_canonicalised[eta[0]][0]) for eta in etas ])
                if not self.is_a_cFF(potential_cFF):
                    families_violating_cFF.append(potential_cFF)

        if len(E_surfaces_violating_connected_cluster) == 0:
            logger.info("ALL %d distinct E-surfaces in denominators correspond to a connected cluster."%len(distinct_E_surfaces))
        else:
            logger.info("The following %d E-surfaces%s do *not* correspond to a connected cluster: ( edge representation | node representation )\n%s%s"%(
                len(E_surfaces_violating_connected_cluster), ' (showing first three)' if verbosity <= 1 else '',
                '\n'.join( '%s | %s'%( str(eta), str(self.convert_E_surface_into_node_representation(eta)) ) for eta in 
                (E_surfaces_violating_connected_cluster[:3] if verbosity<=1 else E_surfaces_violating_connected_cluster) ),
                '\n[...]' if verbosity <= 1 and len(E_surfaces_violating_connected_cluster)>3 else ''
                )
            )

        if len(families_violating_cFF) == 0:
            logger.info("ALL E-surfaces in denominators of the %d TOPT LTD terms correspond to a cross-free family."%n_topt_ltd_terms)
        else:
            logger.info("The following %d families%s do *not* correspond to cross-free-family: ( node representation non canonicalised )\n%s%s"%(
                len(families_violating_cFF), ' (showing first three)' if verbosity <= 1 else '',
                '\n'.join( '%s'%( str(tuple([ tuple(sorted(list(es))) for es in family ])) ) for family in 
                (families_violating_cFF[:3] if verbosity<=1 else families_violating_cFF) ),
                '\n[...]' if verbosity <= 1 and len(families_violating_cFF)>3 else ''
                )
            )

        if not full_analysis:
            return

        # First, canonicalise each cFF by making sure the node representation is the one containing the lowest vertex id
        cFF = [ tuple(sorted( list( ( tuple(sorted(list(eta if self.iv[0] in eta else set(self.iv)-eta))), ) for eta in eta_list) )) for eta_list in cff_analysis['cross_free_family'] ] 
        # cFF not canonicalised
        cFF_not_canonicalised = [ tuple(sorted( list( ( tuple(sorted(list(eta))), ) for eta in eta_list) )) for eta_list in cff_analysis['cross_free_family'] ]

        # for acFF in cFF:
        #     print(acFF)
        # stop
        cFF_encountered = []

        etas_violating_cff = []
        for o, terms in topt_analysis['topt_ltd_terms']:
            for num, etas in terms:

                esurf_structure = tuple(sorted([ self.convert_E_surface_into_node_representation(eta) for eta in etas ]))
                try:
                    cFF_index = cFF.index(esurf_structure)
                    cFF_encountered.append(cFF_index)
                except ValueError:
                    etas_violating_cff.append( etas )

        if len(etas_violating_cff) > 1:
            logger.info("There are %d denominator configurations *not* satisfying the cross-free family condition."%len(etas_violating_cff))
            logger.info(("They are ( edge representation | node representation ):\n%s%s" if verbosity > 1 else "The first three are ( edge representation | node representation ):\n%s%s")%(
                '\n'.join( '%s | %s'% (str(etas),str(tuple([ self.convert_E_surface_into_node_representation(eta) for eta in etas]))) for etas in (etas_violating_cff[:3] if verbosity <= 1 else etas_violating_cff)),
                '\n[...]' if verbosity <= 1 else ''
            ))
        else:
            logger.info("ALL TOPT LTD denominator configurations satisfy the cross-free family condition!")

        logger.info("There is a total of %d possible cross-free families, and %d of them were found in the TOPT LTD expression."%(
            len(cFF), len(cFF_encountered)
        ))
        cFF_not_encountered = [ a_cFF for i_cFF, a_cFF in enumerate(cFF_not_canonicalised) if i_cFF not in cFF_encountered ]
        logger.info("The following%s cross-free-families were unaccounted for (in the node representation):\n%s%s"%(
            ' (first at most three)' if verbosity <= 1 else '',
            '\n'.join(
                str(a_cFF) for a_cFF in ( cFF_not_encountered if verbosity > 1 else cFF_not_encountered[:3] )
            ),
            '\n[...]' if verbosity <= 1 else '',
        ))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Compute the cross-free family LTD representation of a graph and analyses thresholds.""")
    requiredNamed = parser.add_argument_group('required named arguments')

    parser.add_argument('--topology', '-t', dest='topology', type=str, default='box', choices=('box','double-box','pentagon','fishnet2x2','triple-box'),
                        help='Specify the topology to run (default: %(default)s).')

    parser.add_argument('--cas', '-cas', dest='cas', type=str, default='form', choices=('form','sympy'),
                        help='Which computer algebra system to use (default: %(default)s).')

    parser.add_argument('--consider_numerator', '-n', dest='consider_numerator', action='store_true', default=False,
                        help='Enable grouping all terms of an orientation and thus consider numerator (default: %(default)s).')

    parser.add_argument('--not_consider_numerator', '-nn', dest='not_consider_numerator', action='store_true', default=False,
                        help='Disable grouping all terms of an orientation and thus consider numerator (default: %(default)s).')

    parser.add_argument('--no_merge_terms', '-nm', dest='merge_terms', action='store_false', default=True,
                        help='Whether to merge terms with same orientation at all (default: %(default)s).')

    parser.add_argument('--verbosity', '-v', dest='verbosity', type=int, default=1,
                        help='verbosity (default: %(default)s).')

    parser.add_argument('--codes', '-c', dest='n_cores', type=int, default=None,
                        help='Number of cores to run on (default: automatically runs on all available).')

    parser.add_argument('--form_path', '-fp', dest='form_path', type=str, default=None,
                        help='Specify FORM path, make sure to have a form.set file there with MaxTermsize=1M at least (default: automatic).')

    parser.add_argument('--no_save_load_results', '-nslr', dest='save_load_results', action='store_false', default=True,
                        help='Disable the saving and loading of results on file  (default: save and loads results into a file named after the topology).')

    parser.add_argument('--clean', '-cl', dest='clean', action='store_true', default=False,
                        help='Remove all results file so that they are not recycled (default: %(default)s).')

    parser.add_argument('--no_output_mathematica', '-nom', dest='output_mathematica', action='store_false', default=True,
                        help='Disable the generation of the mathematica code (default: %(default)s).')

    parser.add_argument('--no_split_contributions_per_orientation', '-nso', dest='split_orientations', action='store_false', default=True,
                        help='Disable  (default: automatically runs on all available).')

    parser.add_argument('--full_cff_analysis', '-f', dest='full_cff_analysis', action='store_true', default=False,
                        help='Enable a full cff analysis (default: %(default)s).')

    args = parser.parse_args()

    TOPTLTD_analyser = None

    if args.topology=='box':
        internal_edges=((1,2),(2,3),(3,4),(4,1))
        external_edges=((101,1),(102,2),(103,3),(104,4))
        TOPTLTD_analyser = TOPTLTD_Analyser(internal_edges, external_edges, args.topology)
        cff_analyzer = cFF_Analyser(internal_edges, external_edges)

    if args.topology=='pentagon':
        internal_edges=((1,2),(2,3),(3,4),(4,5),(5,1))
        external_edges=((101,1),(102,2),(103,3),(104,4),(105,5))
        TOPTLTD_analyser = TOPTLTD_Analyser(internal_edges, external_edges, args.topology)
        cff_analyzer = cFF_Analyser(internal_edges, external_edges)

    if args.topology=='double-box':
        internal_edges=((1,2),(2,3),(3,4),(4,5),(5,6),(6,1),(3,6))
        external_edges=((11,1),(22,2),(44,4),(55,5))
        TOPTLTD_analyser = TOPTLTD_Analyser(internal_edges, external_edges, args.topology)
        cff_analyzer = cFF_Analyser(internal_edges, external_edges)

    if args.topology=='triple-box':
        internal_edges=((1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,1),(8,3),(7,4))
        external_edges=((11,1),(22,2),(33,5),(44,6))
        TOPTLTD_analyser = TOPTLTD_Analyser(internal_edges, external_edges, args.topology)
        cff_analyzer = cFF_Analyser(internal_edges, external_edges)

    if args.topology=='fishnet2x2':
        internal_edges=((1,2),(2,3),(3,9),(9,4),(4,5),(5,6),(6,7),(7,1),(7,8),(8,9),(2,8),(8,5))
        external_edges=((11,1),(22,3),(33,4),(44,6))
        TOPTLTD_analyser = TOPTLTD_Analyser(internal_edges, external_edges, args.topology)
        cff_analyzer = cFF_Analyser(internal_edges, external_edges)

    if args.consider_numerator:
        TOPTLTD_analyser._ALLOW_NUMERATORS_WHEN_SIMPLIFYING_TOPT_TERMS = True
    if args.not_consider_numerator:
        TOPTLTD_analyser._ALLOW_NUMERATORS_WHEN_SIMPLIFYING_TOPT_TERMS = False
    if args.form_path is not None:
        TOPTLTD_analyser._FORM_PATH = args.form_path

    TOPTLTD_analyser._CAS = args.cas

    TOPTLTD_analyser._MERGE_TERMS = args.merge_terms
    if args.n_cores is not None:
        TOPTLTD_analyser._N_CPUS = args.n_cores

    TOPTLTD_analyser._SPLIT_ORIENTATIONS = args.split_orientations

    save_load_filename = pjoin(root_path,'%s_cff_analysis_result.txt'%(args.topology.replace('-','_')))
    if args.clean and os.path.isfile(save_load_filename):
        logger.info("Removing all results file '%s'."%save_load_filename)
        os.remove(save_load_filename)

    if (not args.save_load_results) or (not os.path.isfile(save_load_filename)):
        cff_analysis = cff_analyzer.analyze(args.full_cff_analysis)
    else:
        try:
            logger.info("Recycling cFF analysis results from file '%s'."%save_load_filename)
            cff_analysis = eval(open(save_load_filename,'r').read())
        except Exception as e:
            raise TOPTLTD_analyser("An exception occurred when recycling results from file '%s'. Exception: %s"%(save_load_filename, str(e)))
    
    if args.save_load_results and (not os.path.isfile(save_load_filename)):
        with open(save_load_filename,'w') as f:
            f.write(pformat(cff_analysis))
        logger.info("Raw cFF analysis results written to file '%s'."%save_load_filename)

    # Generate the map from connected cluster to e surface with the clusters non-canonicalised
    TOPTLTD_analyser.connected_cluster_to_e_surface_map_non_canonicalised = TOPTLTD_analyser.generate_connected_cluster_to_e_surface_map(cff_analysis, canonicalise=False)
    TOPTLTD_analyser.e_surface_to_connected_cluster_map_non_canonicalised = { v : k for k,v in TOPTLTD_analyser.connected_cluster_to_e_surface_map_non_canonicalised.items() }
    # Same, but canonicalised this time.
    TOPTLTD_analyser.connected_cluster_to_e_surface_map = TOPTLTD_analyser.generate_connected_cluster_to_e_surface_map(cff_analysis, canonicalise=True)
    TOPTLTD_analyser.e_surface_to_connected_cluster_map = { v : k for k,v in TOPTLTD_analyser.connected_cluster_to_e_surface_map.items() }

    save_load_filename = pjoin(root_path,'%s_topt_analysis_result.txt'%(args.topology.replace('-','_')))
    if args.clean and os.path.isfile(save_load_filename):
        logger.info("Removing all results file '%s'."%save_load_filename)
        os.remove(save_load_filename)
    if (not args.save_load_results) or (not os.path.isfile(save_load_filename)):
        topt_analysis = TOPTLTD_analyser.analyze()
    else:
        try:
            logger.info("Recycling TOPT results from file '%s'."%save_load_filename)
            topt_analysis = eval(open(save_load_filename,'r').read())
        except Exception as e:
            raise TOPTLTDException("An exception occurred when recycling results from file '%s'. Exception: %s"%(save_load_filename, str(e)))
    
    if args.save_load_results and (not os.path.isfile(save_load_filename)):
        with open(save_load_filename,'w') as f:
            f.write(pformat(topt_analysis))
        logger.info("Raw TOPT results written to file '%s'."%save_load_filename)

    if args.verbosity >= 3:
        logger.info("TOPT LTD Terms:\n%s"%pformat(topt_analysis['topt_ltd_terms']))
    if args.output_mathematica:
        mm_filename = pjoin(root_path,'%s_mathematica_comparison.m'%(args.topology.replace('-','_')))
        logger.info("Writing results in a Mathematica notebook named '%s'."%mm_filename)
        with open(mm_filename,'w') as f:
            f.write('%s\n\n(* cFF representation *)\n\n%s'%(
                topt_analysis['cLTD_mathematica']
                ,
                topt_analysis['cFF_mathematica']
            ))

    # Code for debugging features can be placed here
    #print(cff_analyzer.is_completing_a_family([{1, 2},{1, 2, 3, 5, 6}], {3, 5, 6}))

    TOPTLTD_analyser.test_cff_property(topt_analysis, cff_analysis, args.verbosity, args.full_cff_analysis)
