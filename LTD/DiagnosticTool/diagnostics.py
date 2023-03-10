#!/usr/bin/env python
import numpy
import math
import random

import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath( __file__ )),os.path.pardir))
from ltd_utils import LoopTopology, LoopLine, Propagator, TopologyCollection
import ltd_commons
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

import logging
logging.basicConfig(level=logging.INFO)

class Surface(object):

        """
        Surface contains the label corresponding to the vanishing propagator (self.n_surface), the label corresponding to the vanishing factor (self.marker) and a
        variable (self.exist) establishing if the surface is differenti from the emptyset. Furthermore, it contains a variable establishing if the upper or lower sheet is chosen
        """
        _UPPER=0
        _LOWER=1

        def __init__(self, n, ot_sign, sheet):
            super(Surface, self).__init__()
            super(Surface, self).__init__()
            self.n_surface=n
            self.marker=ot_sign
            self.surface_signs=[]
            self.exist=None
            self.sheet=0
            self.param_variable = None
            self.moms=[]
            self.p0=None
            self.foci=None
            self.n_cut=[]

            if sheet==1:
                self.sheet=self._UPPER
            else:
                self.sheet=self._LOWER

        def is_ellipsoid(self):
            
            assert self.surface_signs[0][0]*self.surface_signs[0][1]!=0,'Zeno screwed up.'
            return ( abs(sum(self.surface_signs[0]) + self.surface_signs[1]) == 3 or
                     (self.surface_signs[1]==0 and abs(sum(self.surface_signs[0]))==2) )
                
        def __str__(self, cut_propagators=None):
            """ Print human-readable information regarding the surface."""

            lines = []

            [(+1,-1),-1]

            lines.append("%-30s : %d"%("Surface ID",self.n_surface))
            if cut_propagators:
                lines.append("%-30s : %s"%("Cut propagators","(%s)"%(','.join('%d'%cp for cp in cut_propagators))))
            lines.append("%-30s : %s"%("Solution type",'+' if self.marker>0 else '-'))
            lines.append("%-30s : %s"%("Surface signature", '(%s %s) %s'%(tuple(
                [{1:'+',-1:'-',0:'0'}[s] for s in list(self.surface_signs[0])+[self.surface_signs[1],]]
            ) ) ))
            lines.append("%-30s : %s"%("Surface momenta", self.moms))
            lines.append("%-30s : %s"%("Surface type", 'Ellipsoid' if self.is_ellipsoid() else 'Hyperboloid'))
            lines.append("%-30s : %s"%("Exists","True" if self.exist else "False"))
            # Only printout hemisphere choice for hyperboloid as it can changed and is not fixed for ellipsoids.
            if not self.is_ellipsoid():
                lines.append("%-30s : %s"%("Hemisphere choice",'+' if self.sheet>0 else '-'))
            lines.append("%-30s : %s"%("Parametrising loop mom. #",'%d'%(self.param_variable+1)))

            return '\n'.join(lines)

        __repr__ = __str__

	def parametric_equation_eval(self, propagators, loopmomenta, n_loops):
            value = 0.
            
            for i in range(0, 2):
                value += self.surface_signs[0][i] * math.sqrt(sum((propagators[self.moms[0][i]].q[j1+1] +
                        sum(propagators[self.moms[0][i]].signature[j] *loopmomenta[j][j1]
                        for j in range(0, n_loops))) ** 2. for j1 in range(0, 3)))

            value +=self.surface_signs[1] * math.sqrt(sum(( propagators[self.moms[1]].q[j1+1] +
                        sum(propagators[self.moms[1]].signature[j] * loopmomenta[j][j1]
                        for j in range(0, n_loops))) ** 2. for j1 in range(0, 3)))

            return value + self.p0


#Momentum info is a class containing, as members, the list of cut propagators (self.cut), the list of propagators (self.prop), the number of loops (n_loops)
#and the target list, which contains the same number of propagators as self.prop; however, the propagator members self.routing, self.q, self.m have been substituted
#with those one would obtain if the choice of loop momenta corresponded to the cut propagators
class Momentum_info(object):

    """
    Momentum info is a class containing, as members, the list of cut propagators (self.cut), the list of propagators (self.prop), the number of loops (n_loops)
    and the target list, which contains the same number of propagators as self.prop; however, the propagator members self.routing, self.q, self.m have been substituted
    with those one would obtain if the choice of loop momenta corresponded to the cut propagators.
    Members:
        -self.cut: contains a list of cut momenta in the form of couplets of numbers (at two loops). The numbers refer to a conventional labelling of the edges.
            Example: for the double triangle, with loop_lines [[0,1],[2],[3,4]], one has
                    self.cut=[0,2] or [1,2] or [0,3] or [0,4] or [1,3] or [1,4] or [2,3] or [2,4]
        -self.prop: contains a list of all propagators of the diagram
        -self.n_loops: number of loops
    """

    def __init__(self, cut, prop, n_loops, signs):
        super(Momentum_info, self).__init__()
        self.cut = cut
        self.prop = prop
        self.n_loops = n_loops
        self.combs = []

    """
    info_prep expresses all the propagator momenta as linear combinations of the cut momenta in self.cut.
    In particular, the output combs[i] contains an object of type propagator whose attributes are
        -the signature of the decomposition of q_i in terms of the cut momenta of self.cut
        -the affine term, depending on the external kinematics, which is what is left out by the previous decomposition
        -the mass of the propagator i
    Example:
        Let q_i=a*q_1+b*q_2+p
        Then combs[i].signature=[a,b], combs[i].q=p
    """


    def info_prep(self):
        mat = []

        for i in range(0, len(self.cut)):
                mat.append(numpy.array(self.prop[self.cut[i]].signature))


        for i in range(0, len(self.prop)):
            #solving the mentioned linear system
            lin_comb = numpy.linalg.solve(numpy.array(mat).T, numpy.array(self.prop[i].signature))  # solving the linear system
            #the external kinematics for each propagator are updated according to this new decomposition
            ext_mom = numpy.array(self.prop[i].q) + sum(
                -lin_comb[j] * numpy.array(self.prop[self.cut[j]].q) for j in range(0, self.n_loops))
            #the target objects are added to a list and then outputted
            finalprop=Propagator(signature=lin_comb,q=ext_mom,m_squared=self.prop[i].m_squared)
            self.combs.append(finalprop)
        return self.combs



class Diagnostic_tool(object):

    """
    The class Diagnostic_tool contains, as members, the list of propagators (self.propagators), the total number of propagators (self.Etot), two lists of surfaces
    the number of loops (self.n_loop), the list of cut propagators (self.cut, as labelled by convention), and the list of propagators in the new basis of cuts as
    produced by the previous class
    """

    def __init__(self, topology):
        super(Diagnostic_tool, self).__init__()
        self.topology=topology

        if topology.n_loops==1:
            self.propagators=topology.loop_lines[0].propagators
            self.Etot=len(topology.loop_lines[0].propagators)
            self.n_loop=topology.n_loops
            self.cut=[[i] for i in range(0,self.Etot)]
            self.sign=[[1] for i in range(0,self.Etot)]

        if topology.n_loops==2:
            self.propagators=[]
            for i_loop_line, loop_line in enumerate(topology.loop_lines):
                for j_prop, prop in enumerate(loop_line.propagators):
                    self.propagators.append(prop)

            #for i in range(0, len(topology.ltd_cut_structure[0])):
            #    for j in range(0, len(topology.loop_lines[i].propagators)):
            #        self.propagators.append(topology.loop_lines[i].propagators[j])


            self.Etot=len(self.propagators)
            self.n_loop=topology.n_loops

            cutty=[]
            signcutty=[]
            for i in range(0, len(topology.ltd_cut_structure)):
                cuttylines=[]
                for j in range(0,len(topology.ltd_cut_structure[i])):
                    if topology.ltd_cut_structure[i][j]!=0:
                        cuttylines.append(j)

                cutty.append([[a,b]
                              for a,b in itertools.product(
                                [i1+sum(len(topology.loop_lines[j].propagators) for j in range(0,cuttylines[0]))
                    for i1 in range(0, len(topology.loop_lines[cuttylines[0]].propagators))],[i2+sum(len(topology.loop_lines[j].propagators) for j in range(0,cuttylines[1])
                                       )
                    for i2 in range(0, len(topology.loop_lines[cuttylines[1]].propagators))])])
                mm=len(topology.loop_lines[cuttylines[0]].propagators)*len(topology.loop_lines[cuttylines[1]].propagators)
                signcutty.append([[topology.ltd_cut_structure[i][cuttylines[0]],topology.ltd_cut_structure[i][cuttylines[1]]]
                                 for l in range(0,mm)
                ])
            self.cut=cutty[0]+cutty[1]+cutty[2]
            self.sign=signcutty[0]+signcutty[1]+signcutty[2]

        self.combs = None
        self.momcheck=None




    def translation_construction(self, surface, loop_mom, n_cut):

        info = Momentum_info(cut=self.cut[n_cut], prop=self.propagators, n_loops=self.n_loop, signs=self.sign[n_cut])
        self.combs = info.info_prep()

        mom1 = 0
        mom2 = 0
        mom3 = [0, 0, 0, 0]
        tot_mom = [0, 0, 0, 0]
        chosen_loop = None
        signs = []
        othersign = 0
        moms=[]
        othermom=0
        surface.p0=self.combs[surface.n_surface].q[0]
        surface.n_cut=n_cut
        #Running though the coefficients of the decomposition of the chosen (n_surface) propagator in terms of cut momenta
        for i in range(0, self.n_loop):

            #the first non zero coefficient of this decomposition is chosen; this coefficient corresponds to a cut propagator of momentum q
            if self.combs[surface.n_surface].signature[i] != 0:

                #running through the loop momentum routing of q=a k_1+b k_2 + p
                for i1 in range(0, self.n_loop):
                    #the first non zero loop momentum contained in q is chosen, say k_1
                    if self.propagators[self.cut[n_cut][i]].signature[i1] != 0:
                        #its label and sign are saved
                        chosen_loop = i1
                        signs.append(self.sign[n_cut][i]*self.combs[surface.n_surface].signature[i])

                        #everything in q which is not k_1 is saved in mom1
                        mom1 = numpy.array(self.propagators[self.cut[n_cut][i]].q) + sum(self.propagators[self.cut[n_cut][i]].signature[i2] * numpy.array(loop_mom[i2]) for i2 in range(i1 + 1, self.n_loop))
                        moms.append(self.cut[n_cut][i])
                        #running through all the remaining coefficients of the decomposition of q in the cut momenta
                        for j in range(i + 1, self.n_loop):
                            #the non-zero coefficients are chosen, say q' has a non-zero coefficient
                            if self.combs[surface.n_surface].signature[j] != 0:
                                #and those whose momentum routing comprises k_1 (in this chosen example) are selected
                                if self.propagators[self.cut[n_cut][j]].signature[i1] != 0:
                                    #everything in q' which is not k_1 is saved in mom2
                                    mom2 = numpy.array(self.propagators[self.cut[n_cut][j]].q) + sum(self.propagators[self.cut[n_cut][j]].signature[i3] * numpy.array(loop_mom[i3]) for i3 in range(0, i1) + range(i1 + 1, self.n_loop))
                                    signs.append(self.sign[n_cut][j]*self.combs[surface.n_surface].signature[j])
                                    moms.append(self.cut[n_cut][j])
                                else:
                                    #if this fails, we save it in mom3 and its sign in other_sign
                                    mom3 = numpy.array(self.propagators[self.cut[n_cut][j]].q) + sum(self.propagators[self.cut[n_cut][j]].signature[i3] * numpy.array(loop_mom[i3]) for i3 in range(0, self.n_loop))
                                    othersign = self.sign[n_cut][j]*self.combs[surface.n_surface].signature[j]
                                    othermom=self.cut[n_cut][j]
                        #if the previous cycle has failed, the only possibility is that the propagator's momentum itself contains k_1 as loop momentum
                        if self.propagators[surface.n_surface].signature[i1] != 0:
                            #again everything is saved
                            mom2 = numpy.array(self.propagators[surface.n_surface].q) + sum(self.propagators[surface.n_surface].signature[i3] * numpy.array(loop_mom[i3]) for i3 in range(0, i1) + range(i1 + 1, self.n_loop))
                            signs.append(surface.marker)
                            moms.append(surface.n_surface)
                        #if this if clause has failed, it means that the previous one has succeded; again we save everything
                        else:
                            mom3 = numpy.array(self.propagators[surface.n_surface].q) + sum(self.propagators[surface.n_surface].signature[i3] * numpy.array(loop_mom[i3]) for i3 in range(0, self.n_loop))
                            othersign = surface.marker
                            othermom=surface.n_surface
                        break
                break


        
       # given the previous information, we can write |k_1+mom1|+|k_1+mom2|=p^0-|k_2+mom3| (again this holds in the specific chosen example)
        #thus we can calculate the translation necessary to make this equation into |k_1+h|+|k_1-h|=p^0-|k_2+mom3|
        surface.moms=[moms,othermom]
        tot_mom = -(numpy.array(mom1) + numpy.array(mom2)) / 2.
        surface.param_variable=chosen_loop
        surface.surface_signs=[signs, othersign]
     #   print([tot_mom, mom3, mom1, mom2])
        self.momcheck=[tot_mom, mom3, mom1, mom2]
        return [tot_mom, mom3, mom1, mom2]


    #The following function constructs the rotation matrix necessary to allign the simmetry axis of the ellipsoid with the z axis
    def rot_matrix_construction(self, surface, loop_mom, n_cut):

        #The vector corresponding to the simmetry axis of the ellipsoid is maj_axis, i.e. the vector h in the previous example
        maj_axis = numpy.array(self.translation_construction(surface, loop_mom,n_cut)[0]) + numpy.array(
            self.translation_construction(surface, loop_mom, n_cut)[2])

        #h is normalized to unity if it is non-zero. If it is, the identity matrix is returned
        mom_norm = (sum(maj_axis[j] ** 2. for j in range(1, 4))) ** (0.5)
        if mom_norm != 0:
            mom_unit = [(1 / mom_norm) * maj_axis[j] for j in range(1, 4)]
        else:
            self.rot_matrix = numpy.identity(3)
            return

        z_axis = [0, 0, 1.]
        #angle between the z axis and h
        angle = math.acos(mom_unit[2])

        #rotation axis obtained as the normalized cross product e_z x h
        rot_axis = numpy.cross(z_axis, mom_unit)
        rot_axis = 1. / (sum(rot_axis[i] ** 2 for i in range(0, 3)) ** 0.5) * numpy.array(rot_axis)

        #matrices necessary to the decomposition of a rotation matrix into identity, simmetric matrix and skew-simmetric matrix
        A_mat = numpy.matrix(
            [[0, -rot_axis[2], rot_axis[1]], [rot_axis[2], 0, -rot_axis[0]], [-rot_axis[1], rot_axis[0], 0]])
        aa_mat = numpy.matrix([[rot_axis[i] * rot_axis[j] for j in range(0, 3)] for i in range(0, 3)])
        #construction of the rotation matrix
        self.rot_matrix = math.cos(angle) * (numpy.identity(3) - aa_mat) + aa_mat + math.sin(angle) * A_mat




    def Get_surface_point(self, surface, kx, ky, a, c, surface_type):
        if surface.sheet==1:
            g=-1
        else:
            g=1

        #parametrization of the z component of the chosen loop momentum
        if surface_type == 1:
            kz = g * (c * (1. - (1. / a) * (kx**2.+ky**2.))) ** 0.5
        else:

            kz = g * (c * (1. + (1. / a) *  (kx**2.+ky**2.) )) ** 0.5

        return kz


    #Determines the existence of a surface (and changes its attribute self.exist) if one of randomly generate n_points produces a non null result
    def determine_existence(self, n_points, surf, n_cut):
        for k in range(0, n_points):
            loopmomenta = [[random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)],
                           [random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)]]

            ll1 = self.get_parametrization(u=random.uniform(0, 1), v=random.uniform(0, 1), loop_momenta=loopmomenta,
                                           surface=surf, n_cut=n_cut)

            if ll1 != None:
                surf.exist=1
                return None
        surf.exist=0

    def get_parametrization(self, u, v, surface, loop_momenta, n_cut):
        
        #The function takes an input vector of loop momenta and uses the necessary ones in between them to calculate the z component of the parametrizing loop variable
        #IN THE AFFINELY MAPPED space, and then maps back the obtained point of the surface to the original space
        #surface type and signs




        p_tr = numpy.array(self.translation_construction(surface, loop_momenta, n_cut))

        surface_type = surface.surface_signs[0][0] * surface.surface_signs[0][1]
        sursign=surface.surface_signs[0][0] + surface.surface_signs[0][1]


        p0=(-self.combs[surface.n_surface].q[0] - surface.surface_signs[1] * sum(
            p_tr[1][l] ** 2 for l in range(1, 4)) ** 0.5) / 2.
        p0_sq = (p0) ** 2.
        p_tr3d = numpy.array([p_tr[0][1], p_tr[0][2], p_tr[0][3]]) + numpy.array(
            [p_tr[2][1],
             p_tr[2][2],
             p_tr[2][3]])

        mod_p_sq = sum(p_tr3d[i] ** 2. for i in range(0, 3))  # /4.

        #parameters of the surface
        ra=p0_sq - mod_p_sq - self.propagators[surface.n_surface].m_squared ** 2.
        rc=(p0_sq / (p0_sq - mod_p_sq)) * ra
        a = abs(ra)
        c = abs(rc)

        surface_type = surface.surface_signs[0][0] * surface.surface_signs[0][1]
        sursign=surface.surface_signs[0][0] + surface.surface_signs[0][1]

        #General conditions to satisfy
        if ra>0 and surface_type==-1:
            return None

        if ra<0 and surface_type==1:
            return None

        if sursign==2 and p0<0:
            return None

        if sursign==-2 and p0>0:
            return None

        if p0*surface.surface_signs[0][0]<0 and surface_type==-1 and surface.sheet==0:
            surface.sheet=1


        if p0*surface.surface_signs[0][0]>0 and surface_type==-1 and surface.sheet==1:
            surface.sheet=0


        if surface_type==-1:
            kx = u*2.-1.
            ky = v*2.-1.
        else:
            kx = (u*2.-1.) * math.sqrt(a)
            ky = (v*2. -1.) * math.sqrt(a-kx**2)

        if (kx**2.+ky**2.<a and surface_type==1) or surface_type==-1:
            kz=self.Get_surface_point(surface=surface, kx=kx, ky=ky, a=a, c=c, surface_type=surface_type)
        elif kx ** 2. + ky ** 2. > a and surface_type == 1:
            return None

        vec = [kx, ky, kz]

        self.rot_matrix_construction(surface, loop_momenta, n_cut)
        
        vec=numpy.array(numpy.dot(self.rot_matrix, vec)) + numpy.array([p_tr[0][1], p_tr[0][2], p_tr[0][3]])
        
        # Stupid fix for when self.rot_matrix is the identity then the return type of the operation above 
        # does not have the same shape... I'll let zeno fix this correctly :S
        if isinstance(vec[0], float):
            vec = numpy.array([vec])

        foci1 = numpy.array(numpy.dot(self.rot_matrix, [0., 0.,math.sqrt(mod_p_sq)])) + numpy.array(
            [p_tr[0][1], p_tr[0][2], p_tr[0][3]])
        foci2=numpy.array(numpy.dot(self.rot_matrix, [0.,0.,-math.sqrt(mod_p_sq)])) + numpy.array([p_tr[0][1], p_tr[0][2], p_tr[0][3]])
        surface.foci=[foci1,foci2]
            
        if surface.param_variable == 0:
            vec = [vec[0], [loop_momenta[1][1], loop_momenta[1][2], loop_momenta[1][3]]]
        else:
            vec = [[loop_momenta[0][1], loop_momenta[0][2], loop_momenta[0][3]], vec[0]]

        return vec

    def all_surfaces(self, n_points):
        all_surf=[]
        l=0
        #print("LAAA")
        for i in range(0, len(self.cut)):
            info = Momentum_info(cut=self.cut[i], prop=self.propagators, n_loops=self.n_loop, signs=self.sign[i])
            self.combs = info.info_prep()
            for j in range(0, self.Etot):
                if self.combs[j].q[0] != 0:
                    #print("LOOOO")
                    surp = Surface(j, 1, 1)
                    surm = Surface(j, -1, -1)

                    p_tr = numpy.array(self.translation_construction(surp, [[0,0,0,0],[0,0,0,0]], i))
                    spl=surp.surface_signs[0][0]+surp.surface_signs[0][1]+surp.surface_signs[1]

                    if abs(spl)>1:
                        #print("HERE")
                        check=0
                        surp.exist=0
                        onsh=self.combs[surp.n_surface].q[0]**2-sum(self.combs[surp.n_surface].q[i1]**2 for i1 in range(1,4))

                        if onsh>0 and (-self.combs[surp.n_surface].q[0])*spl>0:
                            #surp.exist ==1
                            check=1
                        if n_points!=0 and self.n_loop>1:
                            self.determine_existence(n_points, surp, i)
                        else:
                            surp.exist=check
                        if surp.exist == check:
                            if check==1:
                                all_surf.append([self.cut[i], surp])
                        else:
                            print(surp.exist)
                            print(check)
                            print(surp.moms)
                            print(surp.surface_signs)
                            print("ZENOFUCKERY: Numerical existence of ellipsoid does not coincide with analytical")
                            #sys.exit(1)

                    else:
                        self.determine_existence(n_points, surp, i)
                        if surp.exist == 1:
                            all_surf.append([self.cut[i],surp])
                            #print(surp.surface_signs)
                            #print(spl)

                    p_tr = numpy.array(self.translation_construction(surm, [[0,0,0,0],[0,0,0,0]], i))
                    spl=surm.surface_signs[0][0]+surm.surface_signs[0][1]+surm.surface_signs[1]

                    if abs(spl)>1:
                        check=0
                        surm.exist=0
                        onsh=self.combs[surm.n_surface].q[0]**2-sum(self.combs[surm.n_surface].q[i1]**2 for i1 in range(1,4))

                        if onsh>0 and (-self.combs[surm.n_surface].q[0])*spl>0:
                            check=1
                        if n_points!=0:
                            self.determine_existence(n_points, surm, i)
                        else:
                            surm.exist=check
                        if surm.exist==check:
                            if check==1:
                                all_surf.append([self.cut[i],surm])
                        else:
                            print(surm.exist)
                            print(check)
                            print(surm.moms)
                            print(surm.surface_signs)
                            print("ERRRROOOOOORRRR: Numerical existence of ellipsoid does not coincide with analytical")
                            sys.exit(1)
                    else:
                        self.determine_existence(n_points, surm, i)
                        if surm.exist == 1:
                            all_surf.append([self.cut[i],surm])

        
        return all_surf

    def generate_surface_points(self, n_points, surface, n_cut):
        t=0
        point_list=[]
        if surface.exist==1:
            while t!=n_points:
                loopmomenta = [[random.uniform(0, 0.1), random.uniform(-1,1), random.uniform(-1,1), random.uniform(-1,1)],
                               [random.uniform(0, 1), random.uniform(-1,1), random.uniform(-1,1), random.uniform(-1,1)]]

                ll1 = diag.get_parametrization(u=random.uniform(-1, 1), v=random.uniform(-1, 1), loop_momenta=loopmomenta,
                                               surface=surface, n_cut=n_cut)

                if ll1 != None:
                    t+=1
                    point_list.append(ll1)
            return point_list

        else:
            return None

    def generate_limit(self, length, direction, surface, n_cut):
        n_points=int(length/0.001)
        lll=self.generate_surface_points(n_points=1, surface=surface, n_cut=n_cut)
        if lll==None:
            return None
        limit_list=[]
        for i in range(0,n_points):
            limit_list.append(numpy.array(lll)+(0.001*i-(length/2.))*numpy.array(direction))
        return limit_list

    def check_similarity(self, all_surfaces):
        non_redundant_sur=[]
        for i in range(0,len(all_surfaces)):
            l=0

            for j in range(0, len(non_redundant_sur)):
                if all_surfaces[i][1].moms[0][0] == non_redundant_sur[j][0][1].moms[0][0] and all_surfaces[i][1].moms[0][1] == \
                        non_redundant_sur[j][0][1].moms[0][1] and all_surfaces[i][1].moms[1] == non_redundant_sur[j][0][1].moms[1]:
                    if (all_surfaces[i][1].surface_signs[0][0] == non_redundant_sur[j][0][1].surface_signs[0][0] and
                        all_surfaces[i][1].surface_signs[0][1] == non_redundant_sur[j][0][1].surface_signs[0][1] and
                        all_surfaces[i][1].surface_signs[1] == non_redundant_sur[j][0][1].surface_signs[1]) or (
                            all_surfaces[i][1].surface_signs[0][0] == -non_redundant_sur[j][0][1].surface_signs[0][0] and
                            all_surfaces[i][1].surface_signs[0][1] == -non_redundant_sur[j][0][1].surface_signs[0][1] and
                            all_surfaces[i][1].surface_signs[1] == -non_redundant_sur[j][0][1].surface_signs[1]):
                        l += 1
                        non_redundant_sur[j].append(all_surfaces[i])

                if all_surfaces[i][1].moms[0][1] == non_redundant_sur[j][0][1].moms[0][0] and all_surfaces[i][1].moms[0][0] == \
                        non_redundant_sur[j][0][1].moms[0][1] and all_surfaces[i][1].moms[1] == non_redundant_sur[j][0][1].moms[1]:
                    if (all_surfaces[i][1].surface_signs[0][1] == non_redundant_sur[j][0][1].surface_signs[0][0] and
                        all_surfaces[i][1].surface_signs[0][0] == non_redundant_sur[j][0][1].surface_signs[0][1] and
                        all_surfaces[i][1].surface_signs[1] == non_redundant_sur[j][0][1].surface_signs[1]) or (
                            all_surfaces[i][1].surface_signs[0][1] == -non_redundant_sur[j][0][1].surface_signs[0][0] and
                            all_surfaces[i][1].surface_signs[0][0] == -non_redundant_sur[j][0][1].surface_signs[0][1] and
                            all_surfaces[i][1].surface_signs[1] == -non_redundant_sur[j][0][1].surface_signs[1]):
                        l += 1
                        non_redundant_sur[j].append(all_surfaces[i])
                if all_surfaces[i][1].moms[0][0] == non_redundant_sur[j][0][1].moms[1] and all_surfaces[i][1].moms[0][1] == \
                        non_redundant_sur[j][0][1].moms[0][1] and all_surfaces[i][1].moms[1] == non_redundant_sur[j][0][1].moms[0][0]:
                    if (all_surfaces[i][1].surface_signs[0][0] == non_redundant_sur[j][0][1].surface_signs[1] and
                        all_surfaces[i][1].surface_signs[0][1] == non_redundant_sur[j][0][1].surface_signs[0][1] and
                        all_surfaces[i][1].surface_signs[1] == non_redundant_sur[j][0][1].surface_signs[0][0]) or (
                            all_surfaces[i][1].surface_signs[0][0] == -non_redundant_sur[j][0][1].surface_signs[1] and
                            all_surfaces[i][1].surface_signs[0][1] == -non_redundant_sur[j][0][1].surface_signs[0][1] and
                            all_surfaces[i][1].surface_signs[1] == -non_redundant_sur[j][0][1].surface_signs[0][0]):
                        l += 1
                        non_redundant_sur[j].append(all_surfaces[i])
                if all_surfaces[i][1].moms[0][0] == non_redundant_sur[j][0][1].moms[0][0] and all_surfaces[i][1].moms[0][1] == \
                        non_redundant_sur[j][0][1].moms[1] and all_surfaces[i][1].moms[1] == non_redundant_sur[j][0][1].moms[0][1]:
                    if (all_surfaces[i][1].surface_signs[0][0] == non_redundant_sur[j][0][1].surface_signs[0][0] and
                        all_surfaces[i][1].surface_signs[0][1] == non_redundant_sur[j][0][1].surface_signs[1] and
                        all_surfaces[i][1].surface_signs[1] == non_redundant_sur[j][0][1].surface_signs[0][1]) or (
                            all_surfaces[i][1].surface_signs[0][0] == -non_redundant_sur[j][0][1].surface_signs[0][0] and
                            all_surfaces[i][1].surface_signs[0][1] == -non_redundant_sur[j][0][1].surface_signs[1] and
                            all_surfaces[i][1].surface_signs[1] == -non_redundant_sur[j][0][1].surface_signs[0][1]):
                        l += 1
                        non_redundant_sur[j].append(all_surfaces[i])


            if l==0:
                newset=[]
                newset.append(all_surfaces[i])
                non_redundant_sur.append(newset)

        return non_redundant_sur


def print_all_surfaces(classified_surfaces, show_group_members=False):

    res = ["Surface groups content:"]
    for i_group, surfaces_group in enumerate(classified_surfaces):
        res.append('')
        res.append('='*40)
        res.append('')        
        res.append('Defining surface of the surfaces group #%d:'%(i_group+1))
        res.append(surfaces_group[0][1].__str__(cut_propagators = surfaces_group[0][0]))
        if show_group_members:
            res.append('Identical to the following %d other surfaces:'%(len(surfaces_group[1:])))
            for i_member, surface in enumerate(surfaces_group[1:]):
                res.append('-'*40)
                res.append('Identical surface #%d:'%(i_member+1))
                res.append(surface[1].__str__(cut_propagators=surface[0]))
        else:
            res.append('Identical to %d other surfaces.'%(len(surfaces_group[1:])))
    
    res.append('')
    res.append('='*40)
    res.append('')    
    res.append('A total of %d groups of identical surfaces were found (%d elliptic, %d hyperbolic) for a total of %d surfaces.'%(
            len(classified_surfaces),
            len([1 for s in classified_surfaces if s[0][1].is_ellipsoid()]),
            len([1 for s in classified_surfaces if not s[0][1].is_ellipsoid()]),
            sum(len(surface_group) for surface_group in classified_surfaces)
        )
    )

    return '\n'.join(res)

if __name__ == '__main__':
   

    if len(sys.argv)>=2 and not sys.argv[1].startswith('--'):
        topology_name = sys.argv[1]
    else:
        topology_name = "DoubleTriangle"

    options = {
        'n_trial_points' : 2000,
        'seed' : 1,
        'show_group_members' : False
    }
    for arg in sys.argv:
        if arg.startswith('--'):
            try:
                opt, value = arg[2:].split('=')
            except:
                opt, value = arg[2:], None
            
            if opt == 'seed':
                options['seed'] = int(value)
            elif opt == 'n_points':
                options['n_trial_points'] = int(value)
            elif opt == 'show_group_members':
                options['show_group_members'] = (value is None or value.upper() in ['T','TRUE','ON'])
            else:
                logging.warning("Option '%s' is not reckognized and will be ignored."%opt)

    random.seed(options['seed'])

    double_triangle =ltd_commons.hard_coded_topology_collection[topology_name]

    diag=Diagnostic_tool(double_triangle)
   
   
    # Determining existence of surfaces!
    logging.info("Analysing surfaces...")
    #sur=Surface(n=0,ot_sign=-1,sheet=1)
    #diag.determine_existence(100,sur,7)

    # Make sure that N_points is used as follows:
    # Ellipsoids: still run existence check but if exact result says does not exist but numerical random tests says it does, then crash.
    # Hyperboloids: simply do the rnadom check existence test, but maybe eventually improve sampling heuristics.

    all_surfaces=diag.all_surfaces(options['n_trial_points'])
    logging.info('='*40)    


    classified_surfaces = diag.check_similarity(all_surfaces)

    logging.info(print_all_surfaces(classified_surfaces,show_group_members=options['show_group_members']))

    """
    #surcheck=Surface(n=3,ot_sign=1,sheet=1)
    surcheck=classified_surfaces[0][0][1]
    points=diag.generate_surface_points(n_points=10000, surface=surcheck, n_cut=0)
    #print(points)
    fig = plt.figure()
    fig2 = plt.figure()
    ax = Axes3D(fig)
    ax2 = Axes3D(fig2)

    for i in range(0,10000):
        xs = points[i][0][0]
        ys = points[i][0][1]
        zs = points[i][0][2]
        ax.scatter(xs, ys, zs)

    for i in range(0,10000):
        xs = points[i][1][0]
        ys = points[i][1][1]
        zs = points[i][1][2]
        ax2.scatter(xs, ys, zs)


    plt.show()
    """


