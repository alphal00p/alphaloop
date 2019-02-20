import numpy
import math
import random



#The object propagator contains members corresponding to loop momentum routing (self.routing), external kinematics (self.q) and mass (self.m)
class Propagator(object):
        def __init__(self, momentum_routing, external_kinematics, mass):
            super(Propagator, self).__init__()
            self.routing=momentum_routing
            self.p=external_kinematics
            self.m=mass

        def on_shellness(self):
            return self.p[0]**2-sum(self.p[i]**2 for i in range(1,4))


#Surface contains the label corresponding to the vanishing propagator (self.n_surface), the label corresponding to the vanishing factor (self.marker) and a
#variable (self.exist) establishing if the surface is differenti from the emptyset. Furthermore, it contains a variable establishing if the upper or lower sheet is chosen
class Surface(object):

        _UPPER=0
        _LOWER=1

        def __init__(self, n, ot_sign, sheet):
            super(Surface, self).__init__()
            self.n_surface=n
            self.marker=ot_sign
            self.exist=None
            self.sheet=0

            if sheet==1:
                self.sheet=self._UPPER
            else:
                self.sheet=self._LOWER

        #this function changes the self.exist variable whenever the set is non-empty

        def Is_surface(self, prop):
            if prop[self.n_surface].p[0] < 0 and prop[self.n_surface].on_shellness() > 0:
                self.exist=1
            if prop[self.n_surface].on_shellness() < 0:
                self.exist=-1


#Momentum info is a class containing, as members, the list of cut propagators (self.cut), the list of propagators (self.prop), the number of loops (n_loops)
#and the target list, which contains the same number of propagators as self.prop; however, the propagator members self.routing, self.q, self.m have been substituted
#with those one would obtain if the choice of loop momenta corresponded to the cut propagators
class Momentum_info(object):
    def __init__(self, cut, prop, n_loops):
        super(Momentum_info, self).__init__()
        self.cut = cut
        self.prop = prop
        self.n_loops = n_loops
        self.combs = []

    def info_prep(self):

        #mat is a matrix which contains, as rows, the routing of the cut propagator. for each propagator corresponding to routing (a,b) and total momentum q the system mat*(x,y)=(a,b)
        #then (x,y) are the coefficients of the decomposition of q in the basis of the cuts

        mat = [self.prop[i].routing for i in self.cut]

        for i in range(0, len(self.prop)):
            #solving the mentioned linear system
            lin_comb = numpy.linalg.solve(numpy.array(mat).T, numpy.array(self.prop[i].routing))  # solving the linear system
            #the external kinematics for each propagator are updated according to this new decomposition
            ext_mom = numpy.array(self.prop[i].p) + sum(
                -lin_comb[j] * numpy.array(self.prop[self.cut[j]].p) for j in range(0, self.n_loops))
            #the target objects are added to a list and then outputted
            finalprop=Propagator(lin_comb,ext_mom,self.prop[i].m)
            self.combs.append(finalprop)

        return self.combs

#The class Diagnostic_tool contains, as members, the list of propagators (self.propagators), the total number of propagators (self.Etot), two lists of surfaces
#the number of loops (self.n_loop), the list of cut propagators (self.cut, as labelled by convention), and the list of propagators in the new basis of cuts as
#produced by the previous class

class Diagnostic_tool(object):
    def __init__(self, propagators, Etot, n_loop, cut):
        super(Diagnostic_tool, self).__init__()
        self.propagators = propagators
        self.Etot = Etot
        self.n_loop = n_loop
        self.rot_matrix = []
        self.cut = cut
        info = Momentum_info(cut=cut, prop=propagators, n_loops=n_loop)
        self.combs = info.info_prep()



    def translation_construction(self, surface, loop_mom):

        mom1 = 0
        mom2 = 0
        mom3 = [0, 0, 0, 0]
        tot_mom = [0, 0, 0, 0]
        chosen_loop = 0
        signs = []
        othersign = 0

        #Running though the coefficients of the decomposition of the chosen (n_surface) propagator in terms of cut momenta
        for i in range(0, self.n_loop):

            #the first non zero coefficient of this decomposition is chosen; this coefficient corresponds to a cut propagator of momentum q
            if self.combs[surface.n_surface].routing[i] != 0:

                #running through the loop momentum routing of q=a k_1+b k_2 + p
                for i1 in range(0, self.n_loop):
                    #the first non zero loop momentum contained in q is chosen, say k_1
                    if self.propagators[self.cut[i]].routing[i1] != 0:
                        #its label and sign are saved
                        chosen_loop = i1
                        signs.append(self.combs[surface.n_surface].routing[i])
                        #everything in q which is not k_1 is saved in mom1
                        mom1 = numpy.array(self.propagators[self.cut[i]].p) + sum(self.propagators[self.cut[i]].routing[i2] * numpy.array(loop_mom[i2]) for i2 in range(i1 + 1, self.n_loop))

                        #running through all the remaining coefficients of the decomposition of q in the cut momenta
                        for j in range(i + 1, self.n_loop):
                            #the non-zero coefficients are chosen, say q' has a non-zero coefficient
                            if self.combs[surface.n_surface].routing[j] != 0:
                                #and those whose momentum routing comprises k_1 (in this chosen example) are selected
                                if self.propagators[self.cut[j]].routing[i1] != 0:
                                    #everything in q' which is not k_1 is saved in mom2
                                    mom2 = numpy.array(self.propagators[self.cut[j]].p) + sum(self.propagators[self.cut[j]].routing[i3] * numpy.array(loop_mom[i3]) for i3 in range(0, i1 - 1) + range(i1 + 1, self.n_loop))
                                    signs.append(self.combs[surface.n_surface].routing[j])
                                else:
                                    #if this fails, we save it in mom3 and its sign in other_sign
                                    mom3 = numpy.array(self.propagators[self.cut[j]].p) + sum(self.propagators[self.cut[j]].routing[i3] * numpy.array(loop_mom[i3]) for i3 in range(0, i1 - 1) + range(i1 + 1, self.n_loop))
                                    othersign = self.combs[surface.n_surface].routing[j]
                        #if the previous cycle has failed, the only possibility is that the propagator's momentum itself contains k_1 as loop momentum
                        if self.propagators[surface.n_surface].routing[i1] != 0:
                            #again everything is saved
                            mom2 = numpy.array(self.propagators[surface.n_surface].p) + sum(self.propagators[surface.n_surface].routing[i3] * numpy.array(loop_mom[i3]) for i3 in range(0, i1 - 1) + range(i1 + 1, self.n_loop))
                            signs.append(surface.marker)
                        #if this if clause has failed, it means that the previous one has succeded; again we save everything
                        else:
                            mom3 = numpy.array(self.propagators[surface.n_surface].p) + sum(self.propagators[surface.n_surface].routing[i3] * numpy.array(loop_mom[i3]) for i3 in range(0, i1 - 1) + range(i1 + 1, self.n_loop))
                            othersign = surface.marker
                        break
                break


        #given the previous information, we can write |k_1+mom1|+|k_1+mom2|=p^0-|k_2+mom3| (again this holds in the specific chosen example)
        #thus we can calculate the translation necessary to make this equation into |k_1+h|+|k_1-h|=p^0-|k_2+mom3|
        tot_mom = -(numpy.array(mom1) + numpy.array(mom2)) / 2.
        return [tot_mom, chosen_loop, mom3, [signs, othersign], mom1]


    #The following function constructs the rotation matrix necessary to allign the simmetry axis of the ellipsoid with the z axis
    def rot_matrix_construction(self, surface, loop_mom):

        #The vector corresponding to the simmetry axis of the ellipsoid is maj_axis, i.e. the vector h in the previous example
        maj_axis = numpy.array(self.translation_construction(surface, loop_mom)[0]) + numpy.array(
            self.translation_construction(surface, loop_mom)[4])

        #h is normalized to unity if it is non-zero. If it is, the identity matrix is returned
        mom_norm = (sum(maj_axis[j] ** 2 for j in range(1, 4))) ** (0.5)
        if mom_norm != 0:
            mom_unit = [(1 / mom_norm) * maj_axis[j] for j in range(1, 4)]
        else:
            self.rot_matrix = numpy.identity(3)
            return

        z_axis = [0, 0, 1]

        #angle between the z axis and h
        angle = math.acos(mom_unit[2])

        #rotation axis obtained as the normalized cross product e_z x h
        rot_axis = numpy.cross(z_axis, mom_unit)
        rot_axis = 1 / (sum(rot_axis[i] ** 2 for i in range(0, 3)) ** 0.5) * numpy.array(rot_axis)

        #matrices necessary to the decomposition of a rotation matrix into identity, simmetric matrix and skew-simmetric matrix
        A_mat = numpy.matrix(
            [[0, -rot_axis[2], rot_axis[1]], [rot_axis[2], 0, -rot_axis[0]], [-rot_axis[1], rot_axis[0], 0]])
        aa_mat = numpy.matrix([[rot_axis[i] * rot_axis[j] for j in range(0, 3)] for i in range(0, 3)])

        #construction of the rotation matrix
        self.rot_matrix = math.cos(angle) * (numpy.identity(3) - aa_mat) + aa_mat + math.sin(angle) * A_mat


    #The following function takes an input vector of loop momenta and uses the necessary ones in between them to calculate the z component of the parametrizing loop variable
    #IN THE AFFINELY MAPPED space, and then maps back the obtained point of the surface to the original space
    def Get_surface_point(self, surface, loop_mom):
        #constructing the rotation matrix
        self.rot_matrix_construction(surface, loop_mom)
        #constructing the translation vector
        p_tr = numpy.array(self.translation_construction(surface, loop_mom))

        #determining if the parametrized surface is an ellipsoid or an hyperboloid
        surface_type = p_tr[3][0][0] * p_tr[3][0][1]

        # determining the parameters of the surface
        p0_sq = ((self.combs[surface.n_surface].p[0] - p_tr[3][1] * sum(
            p_tr[2][l] ** 2 for l in range(1, 4)) ** 0.5) / 2.) ** 2

        p_tr3d = numpy.array([p_tr[0][1], p_tr[0][2], p_tr[0][3]]) + numpy.array(
            [self.translation_construction(surface, loop_mom)[4][1],self.translation_construction(surface, loop_mom)[4][2],self.translation_construction(surface, loop_mom)[4][3]])

        mod_p_sq = sum(p_tr3d[i] ** 2 for i in range(0, 3))  # /4.

        a = abs(p0_sq - mod_p_sq - self.propagators[surface.n_surface].m ** 2)
        c = abs((p0_sq / (p0_sq - mod_p_sq)) * a)

        #hyperboloids exist only if the on_shellness of the momentum is less than the mass square
        # ellipsoids exist only if the major axis has positive length
        if self.combs[surface.n_surface].p[0] - p_tr[3][1] * sum(
                p_tr[2][l] ** 2 for l in range(1, 4)) ** 0.5 < 0 and surface_type == 1:
            print("the surface does not exist for these values of loop momenta")
            return

        if p0_sq - mod_p_sq - self.propagators[surface.n_surface].m ** 2>0 and surface_type==-1:
            print("the surface does not exist for these values of the external momenta")
            return
        #parametrization of the z component of the chosen loop momentum
        if surface_type == 1:

            kz = (-1) ** (surface.sheet) * (c * (1 - (1 / a) * sum((loop_mom[p_tr[1]][j]) ** 2 for j in range(1,3)))) ** 0.5
        else:

            kz = (-1) ** (surface.sheet) * (c * (1 + (1 / a) * sum((loop_mom[p_tr[1]][j]) ** 2 for j in range(1, 3)))) ** 0.5
        # THE INPUT LOOP MOMENTA ARE PLUGGED IN THE ROTATED SURFACE
        vec = [loop_mom[p_tr[1]][1], loop_mom[p_tr[1]][2], kz]
        #transforming back the vector on the surface
        vec = numpy.array(numpy.dot(self.rot_matrix, vec)) + numpy.array([p_tr[0][1], p_tr[0][2], p_tr[0][3]])

        if p_tr[1] == 0:
            vec = [vec[0], [loop_mom[1][1], loop_mom[1][2], loop_mom[1][3]]]
        else:
            vec = [[loop_mom[0][1], loop_mom[0][2], loop_mom[0][3]], vec[0]]

        return vec




class trial_function(object):
    def __init__(self):
        super(trial_function, self).__init__()

    def Eval(self, loop_value, p1, p2, p0, sign):
        mod_loop1 = sum((loop_value[0][i] + p1[i]) ** 2 for i in range(0, 3))
        mod_loop2 = sum((loop_value[0][i] + p2[i]) ** 2 for i in range(0, 3))
        return math.sqrt(mod_loop1 + 0 ** 2) + ((-1) ** sign) * math.sqrt(mod_loop2 + 0 ** 2) - p0

    def Eval2(self, loop_value, p1, p2, p3, p0):
        mod_loop1 = sum((loop_value[0][i] + p1[i]) ** 2 for i in range(0, 3))
        mod_loop2 = sum((loop_value[1][i] + p2[i]) ** 2 for i in range(0, 3))
        mod_loop3 = sum((loop_value[1][i] + loop_value[0][i] + p3[i]) ** 2 for i in range(0, 3))
        return math.sqrt(mod_loop3) + math.sqrt(mod_loop1 + 0 ** 2) + math.sqrt(mod_loop2 + 0 ** 2) - p0


prop1=Propagator([1, 0], [0.25, 0, 0, 1], 0)
prop2=Propagator([0, 1], [0, 1, 0, 3], 0)
prop3=Propagator([1, 1], [0, 6, 0.5, 0], 0)
prop4=Propagator([0, 1], [3, 0, 2, 0], 0)
prop5=Propagator([1, 0], [0.5, 1, 1, 0], 0)


Diag = Diagnostic_tool(propagators=[prop1,prop2,prop3,prop4,prop5], Etot=5, n_loop=2, cut=[0, 2])

sur=Surface(n=4, ot_sign=-1, sheet=1)
trial_f = trial_function()
vec_collection=[]
for i in range(0,100):
    x=random.uniform(-3,3)
    y=random.uniform(-3,3)
    vec = Diag.Get_surface_point(sur, [[0, x, y, 1.1], [0, 0.3, 1.23, 0]])
    vec_collection.append(trial_f.Eval(vec, [0, 0, 1], [1, 1, 0], 0.25,1))

print(vec_collection)


prop11=Propagator([1, 0], [0, 0, 0, 1], 0)
prop21=Propagator([0, 1], [0, 1, 0, 3], 0)
prop31=Propagator([1, 1], [7.5, 1, 0.5, 0], 0)
prop41=Propagator([0, 1], [5, 0, 2, 0], 0)
prop51=Propagator([1, 0], [0.5, 1, 1, 0], 0)


Diag = Diagnostic_tool(propagators=[prop11,prop21,prop31,prop41,prop51], Etot=5, n_loop=2, cut=[0, 1])


sur2=Surface(n=2, ot_sign=1, sheet=1)
vec_collection2=[]
for i in range(0,100):
    x=random.uniform(-1,1)
    y=random.uniform(-1,1)
    vec = Diag.Get_surface_point(sur2, [[0, x, y, 1.1], [0, 0.3, 1.23, 0]])
    vec_collection2.append(trial_f.Eval2(vec,[0,0,1],[1,0,3],[1, 0.5, 0], 7.5))

print(vec_collection2)
