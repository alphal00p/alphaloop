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

class TOPTLTDException(Exception):
    pass

class TOPTLTD_Analyser(object):

    _ALLOW_NUMERATORS_WHEN_SIMPLIFYING_TOPT_TERMS = False
    _N_CPUS = multiprocessing.cpu_count()
    _MERGE_TERMS = True
    _CAS = 'form' # Other option is 'sympy'
    _FORM_PATH = None

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
            #E_surfaces.append( (tuple(sorted(accumulated_energies)), tuple([]) ) )

        return tuple(E_surfaces)

    def merge_topt_terms_with_numerators_sympy(self, topt_terms, **opts):
        
        import sympy as sp

        if len(topt_terms) == 1:
            return topt_terms

        Es = sp.Matrix(len(self.ie)+1,1, lambda i,j: sp.Symbol('E%d' % (i)))
        ps = sp.Matrix(len(self.ev)+1,1, lambda i,j: sp.Symbol('p%d' % (i)))
        etas = {}
        for num, denom in topt_terms:
            for E_surf in denom:
                if E_surf not in etas:
                    etas[E_surf] = (sum(Es[i_e] for i_e in E_surf[0]) if len(E_surf[0])>0 else 0) + \
                                (sum(ps[i_e] for i_e in E_surf[1]) if len(E_surf[1])>0 else 0)

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

        return [ [ numerator.expand(), set(denominator) ], ]

    def merge_topt_terms_with_numerators_form(self, topt_terms, **opts):
        
        import form

        if len(topt_terms) == 1:
            return topt_terms

        Es = ['E%d'%i for i in range(len(self.ie))]
        ps = ['p%d'%i for i in range(len(self.ev))]
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
                numerator = num
                denominator = [sorted_etas[i][0] for i in surviving_eta_indices]
            except form.formlink.FormError as e:
                print("Error when running FORM.\nError:\n%s\nFORM code:\n%s"%(str(e),FORM_code))
                raise e
            #print(num)
            #print(surviving_eta_indices)
            #stop

        return [ [ numerator, set(denominator) ], ]

    def merge_topt_terms_without_numerators(self, topt_terms, debug=False):
        
        while(len(topt_terms) > 1):
            
            for i_a, i_b in itertools.combinations( list(range(len(topt_terms))), 2 ):
                a, b = topt_terms[i_a], topt_terms[i_b]
                if len(a[1])!=len(b[1]):
                    raise TOPTLTDException("All terms generated in merge_topt_terms_with_numerators should have the same number of E-surfaces in the denominator.")
                non_common_e_surfs =  list(a[1].symmetric_difference(b[1]))
                if len(non_common_e_surfs) != 2:
                    continue
                # Make sure the same edge does not appear twice in the merged term
                if len(set(non_common_e_surfs[0][0]).intersection(set(non_common_e_surfs[1][0])))>0 or len(set(non_common_e_surfs[0][1]).intersection(set(non_common_e_surfs[1][1])))>0:
                    continue

                merged_e_surf_in_num = ( tuple(sorted(list(non_common_e_surfs[0][0])+list(non_common_e_surfs[1][0]))),
                                         tuple(sorted(list(non_common_e_surfs[0][1])+list(non_common_e_surfs[1][1])))  )
                
                common_e_surfs = list(a[1].intersection(b[1]))
                if merged_e_surf_in_num in common_e_surfs:
                    topt_terms[:] = [ [ 1, (set(common_e_surfs) - {merged_e_surf_in_num,}).union(set(non_common_e_surfs)) ] ] + [ t for i_t, t in enumerate(topt_terms) if i_t not in [i_a, i_b] ]
                    break
                else:
                    continue
            else:
                break

        return topt_terms

    def merge_topt_terms_helper(self, args):

        orientation, topt_terms = args

        # Store the on-shell energy polynomial in the numerator and list of E-surfaces (linear polynomials) in the denominator.
        numerator_and_E_surfaces_per_term = [ [ 1, set(self.get_E_surfaces_for_topt_ordering(t)) ] for t in topt_terms ]

        # Iteratively merge terms
        if self._MERGE_TERMS:
            if self._ALLOW_NUMERATORS_WHEN_SIMPLIFYING_TOPT_TERMS:
                if self._CAS == 'sympy':
                    numerator_and_E_surfaces_per_term = self.merge_topt_terms_with_numerators_sympy(numerator_and_E_surfaces_per_term)
                else:
                    numerator_and_E_surfaces_per_term = self.merge_topt_terms_with_numerators_form(numerator_and_E_surfaces_per_term)
            else:
                numerator_and_E_surfaces_per_term = self.merge_topt_terms_without_numerators(numerator_and_E_surfaces_per_term)#, debug=(orientation==(-1, 1, -1, 1)))
        return orientation, tuple( ( nE[0], tuple(sorted(list(nE[1]))) ) for nE in numerator_and_E_surfaces_per_term )

    def process_topt_terms(self, topt_terms):

        processed_terms = {}
        print("Processing TOPT terms on %d CPUS and %sconsidering numerators for merging..."%(self._N_CPUS, 'NOT ' if not self._ALLOW_NUMERATORS_WHEN_SIMPLIFYING_TOPT_TERMS else ''))
        
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
                #if o != (-1, 1, -1, 1, -1, 1, -1):
                #    continue
                orientation, merged_terms = self.merge_topt_terms_helper([o,t])
                processed_terms[orientation] = merged_terms
                #print("Processing orientation #%d/%d"%(len(processed_terms),len(topt_terms)),end='\r')
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
                '<|"Num"->(%s),"Esurfs"->%s|>'%(
                    re.sub(r'p(\d+)',r'p\1E',re.sub(r'E(\d+)',r'OSE[\1]',str(t[0]))),
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

( (Complex[0, -1])^(%(n_loops)d) / ((Times@@Table[2*OSE[ie],{ie,0,NInternalEdges-1}])/.OSEReplacement) )*Total[ResPerOrientations]

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
"""cLTDPATH="<PASTE_HERE_YOUR_CLTD_PATH>"(*"/Users/vjhirsch/HEP_programs/cLTD/cLTD.m"*);
FORMPATH="<PASTE_HERE_YOUR_FORM_PATH>"(*"/Users/vjhirsch/HEP_programs/form/bin/form"*);
Get[cLTDPATH];
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
        print("Original TOPT expression has %d orientations and %d terms, with a maximum of %d terms per orientation."%(
            len(topt_terms.keys()), sum(len(v) for v in topt_terms.values()), max(len(v) for v in topt_terms.values())
        ))

        # Now merge TOPT terms
        processed_topt_terms = self.process_topt_terms(topt_terms)
        print("After merging, the TOPT-LTD expression has %d orientations and %d terms, with a maximum of %d terms per orientation."%(
            len(processed_topt_terms.keys()), sum(len(v) for v in processed_topt_terms.values()), max(len(v) for v in processed_topt_terms.values())
        ))

        analysis['topt_ltd_terms'] = sorted( list(processed_topt_terms.items()), key=lambda el: el[0] )

        print("Maximal number of square roots in E-surfaces in the denominator of any term: %d"%(
            max( max( max(len(d[0]) for d in denom) for num, denom in t) for o,t in analysis['topt_ltd_terms'] )
        ))
        print("Maximal number of E-surfaces in the denominator of any term: %d"%(
            max( max( len(denom) for num, denom in t) for o,t in analysis['topt_ltd_terms'] )
        ))

        analysis['cFF_mathematica'] = self.generate_mathematica_cFF_code(analysis['topt_ltd_terms'])

        return analysis

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Compute the cross-free family LTD representation of a graph and analyses thresholds.""")
    requiredNamed = parser.add_argument_group('required named arguments')

    parser.add_argument('--topology', '-t', dest='topology', type=str, default='box', choices=('box','double-box','pentagon'),
                        help='Specify the topology to run (default: %(default)s).')

    parser.add_argument('--cas', '-cas', dest='cas', type=str, default='form', choices=('form','sympy'),
                        help='Which computer algebra system to use (default: %(default)s).')

    parser.add_argument('--consider_numerator', '-n', dest='consider_numerator', action='store_true', default=False,
                        help='Enable grouping all terms of an orientation and thus consider numerator (default: %(default)s).')

    parser.add_argument('--not_consider_numerator', '-nn', dest='not_consider_numerator', action='store_true', default=False,
                        help='Disable grouping all terms of an orientation and thus consider numerator (default: %(default)s).')

    parser.add_argument('--no_merge_terms', '-nm', dest='merge_terms', action='store_false', default=True,
                        help='Whether to merge terms with same orientation at all (default: %(default)s).')

    parser.add_argument('--verbosity', '-v', dest='verbosity', type=int, default=2,
                        help='verbosity (default: %(default)s).')

    parser.add_argument('--codes', '-c', dest='n_cores', type=int, default=None,
                        help='Number of cores to run on (default: automatically runs on all available).')

    parser.add_argument('--form_path', '-fp', dest='form_path', type=str, default=None,
                        help='Specify FORM path, make sure to have a form.set file there with MaxTermsize=1M at least (default: automatic).')

    args = parser.parse_args()

    TOPTLTD_analyser = None

    if args.topology=='box':
        internal_edges=((1,2),(2,3),(3,4),(4,1))
        external_edges=((101,1),(102,2),(103,3),(104,4))
        TOPTLTD_analyser = TOPTLTD_Analyser(internal_edges, external_edges, args.topology)

    if args.topology=='pentagon':
        internal_edges=((1,2),(2,3),(3,4),(4,5),(5,1))
        external_edges=((101,1),(102,2),(103,3),(104,4),(105,5))
        TOPTLTD_analyser = TOPTLTD_Analyser(internal_edges, external_edges, args.topology)

    if args.topology=='double-box':
        internal_edges=((1,2),(2,3),(3,4),(4,5),(5,6),(6,1),(3,6))
        external_edges=((11,1),(22,2),(44,4),(55,5))
        TOPTLTD_analyser = TOPTLTD_Analyser(internal_edges, external_edges, args.topology)
    
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

#    my_f4=cross_f([1,2,3,4,5,6],[[1,2],[2,3],[3,4],[4,5],[5,6],[6,1],[3,6]],[[1,11],[2,22],[4,44],[5,55]])        
#    my_f4.print_family("edges")

    analysis = TOPTLTD_analyser.analyze()
    if args.verbosity >= 3:
        print("TOPT LTD Terms:\n")
        pprint(analysis['topt_ltd_terms'])
    if args.verbosity >= 2:
        mm_filename = pjoin(root_path,'%s_mathematica_comparison.m'%(args.topology.replace('-','_')))
        print("Writing results in a Mathematica notebook named '%s'."%mm_filename)
        with open(mm_filename,'w') as f:
            f.write('\n(* cLTD representation *)\n\n%s\n\n(* cFF representation *)\n\n%s'%(
                analysis['cLTD_mathematica']
                ,
                analysis['cFF_mathematica']
            ))