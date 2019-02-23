import numpy
import math
import random



class Propagator(object):
        """
        The object propagator contains members corresponding to loop momentum routing (self.routing), external kinematics (self.q) and mass (self.m)
        """
        def __init__(self, momentum_routing, external_kinematics, mass):
            super(Propagator, self).__init__()
            self.routing=momentum_routing
            self.p=external_kinematics
            self.m=mass

        def on_shellness(self):
            return self.p[0]**2-sum(self.p[i]**2 for i in range(1,4))



class Surface(object):

        """
        Surface contains the label corresponding to the vanishing propagator (self.n_surface), the label corresponding to the vanishing factor (self.marker) and a
        variable (self.exist) establishing if the surface is differenti from the emptyset. Furthermore, it contains a variable establishing if the upper or lower sheet is chosen
        """
        _UPPER=0
        _LOWER=1

        def __init__(self, n, ot_sign, sheet):
            super(Surface, self).__init__()
            self.n_surface=n
            self.marker=ot_sign
            self.surface_signs=[]
            self.exist=None
            self.sheet=0
            self.param_variable = None

            if sheet==1:
                self.sheet=self._UPPER
            else:
                self.sheet=self._LOWER





#Momentum info is a class containing, as members, the list of cut propagators (self.cut), the list of propagators (self.prop), the number of loops (n_loops)
#and the target list, which contains the same number of propagators as self.prop; however, the propagator members self.routing, self.q, self.m have been substituted
#with those one would obtain if the choice of loop momenta corresponded to the cut propagators
class Momentum_info(object):

    """
    Momentum info is a class containing, as members, the list of cut propagators (self.cut), the list of propagators (self.prop), the number of loops (n_loops)
    and the target list, which contains the same number of propagators as self.prop; however, the propagator members self.routing, self.q, self.m have been substituted
    with those one would obtain if the choice of loop momenta corresponded to the cut propagators
    """

    def __init__(self, cut, prop, n_loops, signs):
        super(Momentum_info, self).__init__()
        self.cut = cut
        self.prop = prop
        self.n_loops = n_loops
        self.combs = []
        self.signs=signs

    def info_prep(self):

        """
        mat is a matrix which contains, as rows, the routing of the cut propagator. for each propagator corresponding to routing (a,b) and total momentum q the system mat*(x,y)=(a,b)
        then (x,y) are the coefficients of the decomposition of q in the basis of the cuts
        """

        mat = []

        for i in range(0, len(self.cut)):
                mat.append(numpy.array(self.prop[self.cut[i]].routing))


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



class Diagnostic_tool(object):

    """
    The class Diagnostic_tool contains, as members, the list of propagators (self.propagators), the total number of propagators (self.Etot), two lists of surfaces
    the number of loops (self.n_loop), the list of cut propagators (self.cut, as labelled by convention), and the list of propagators in the new basis of cuts as
    produced by the previous class
    """

    def __init__(self, propagators, Etot, n_loop, cut, sign):
        super(Diagnostic_tool, self).__init__()
        self.propagators = propagators
        self.Etot = Etot
        self.n_loop = n_loop
        self.rot_matrix = []
        self.cut = cut
        info = Momentum_info(cut=cut, prop=propagators, n_loops=n_loop, signs=sign)
        self.combs = info.info_prep()
        self.sign=sign



    def translation_construction(self, surface, loop_mom):

        mom1 = 0
        mom2 = 0
        mom3 = [0, 0, 0, 0]
        tot_mom = [0, 0, 0, 0]
        chosen_loop = None
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
                        signs.append(self.sign[i]*self.combs[surface.n_surface].routing[i])
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
                                    signs.append(self.sign[j]*self.combs[surface.n_surface].routing[j])
                                else:
                                    #if this fails, we save it in mom3 and its sign in other_sign
                                    mom3 = numpy.array(self.propagators[self.cut[j]].p) + sum(self.propagators[self.cut[j]].routing[i3] * numpy.array(loop_mom[i3]) for i3 in range(0, i1 - 1) + range(i1 + 1, self.n_loop))
                                    othersign = self.sign[j]*self.combs[surface.n_surface].routing[j]
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


        """
        given the previous information, we can write |k_1+mom1|+|k_1+mom2|=p^0-|k_2+mom3| (again this holds in the specific chosen example)
        thus we can calculate the translation necessary to make this equation into |k_1+h|+|k_1-h|=p^0-|k_2+mom3|
        """
        tot_mom = -(numpy.array(mom1) + numpy.array(mom2)) / 2.
        surface.param_variable=chosen_loop
        surface.surface_signs=[signs, othersign]
        return [tot_mom, mom3, mom1, mom2]


    #The following function constructs the rotation matrix necessary to allign the simmetry axis of the ellipsoid with the z axis
    def rot_matrix_construction(self, surface, loop_mom):

        #The vector corresponding to the simmetry axis of the ellipsoid is maj_axis, i.e. the vector h in the previous example
        maj_axis = numpy.array(self.translation_construction(surface, loop_mom)[0]) + numpy.array(
            self.translation_construction(surface, loop_mom)[2])

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




    def Get_surface_point(self, surface, kx, ky, a, c, surface_type):

        #parametrization of the z component of the chosen loop momentum
        if surface_type == 1:

            kz = (-1) ** (surface.sheet) * (c * (1 - (1 / a) * (kx**2+ky**2))) ** 0.5
        else:

            kz = (-1) ** (surface.sheet) * (c * (1 + (1 / a) *  (kx**2+ky**2) )) ** 0.5

        return kz


    #Determines the existence of a surface (and changes its attribute self.exist) if one of randomly generate n_points produces a non null result
    def determine_existence(self, n_points, surf):
        for k in range(0, n_points):
            loopmomenta = [[random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)],
                           [random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)]]

            ll1 = diag.get_parametrization(u=random.uniform(0, 1), v=random.uniform(0, 1), loop_momenta=loopmomenta,
                                           surface=surf)

            if ll1 != None:
                surf.exist=1
                return None
        surf.exist=None


    def get_parametrization(self, u, v, surface, loop_momenta):
        """
        The function takes an input vector of loop momenta and uses the necessary ones in between them to calculate the z component of the parametrizing loop variable
        IN THE AFFINELY MAPPED space, and then maps back the obtained point of the surface to the original space
        """
        p_tr = numpy.array(self.translation_construction(surface, loop_momenta))
        p0=(-self.combs[surface.n_surface].p[0] - surface.surface_signs[1] * sum(    #CHANGED SIGN TO THE ENERGY COMPONENT
            p_tr[1][l] ** 2 for l in range(1, 4)) ** 0.5) / 2.
        p0_sq = (p0) ** 2

        p_tr3d = numpy.array([p_tr[0][1], p_tr[0][2], p_tr[0][3]]) + numpy.array(
            [p_tr[2][1],
             p_tr[2][2],
             p_tr[2][3]])

        mod_p_sq = sum(p_tr3d[i] ** 2 for i in range(0, 3))  # /4.

        #parameters of the surface
        ra=p0_sq - mod_p_sq - self.propagators[surface.n_surface].m ** 2
        rc=(p0_sq / (p0_sq - mod_p_sq)) * ra
        a = abs(ra)
        c = abs(rc)

        #surface type ans signs
        surface_type = surface.surface_signs[0][0] * surface.surface_signs[0][1]
        sursign=surface.surface_signs[0][0] + surface.surface_signs[0][1]

        #General conditions to satisfy

        if ra>0 and surface_type==-1:
            return None

        if ra<0 and surface_type==1:
            return None

        if sursign==2 and self.combs[surface.n_surface].p[0]>0:
            return None

        if sursign==-2 and self.combs[surface.n_surface].p[0]<0:
            return None

        if p0<0 and surface_type==-1 and surface.sheet==0:
            surface.sheet=1

        if p0>0 and surface_type==-1 and surface.sheet==1:
            surface.sheet=0

        kx=u*a
        ky=v*a
        kz=self.Get_surface_point(surface=surface, kx=kx, ky=ky, a=a, c=c, surface_type=surface_type)
        vec = [kx, ky, kz]

        self.rot_matrix_construction(surface, loop_momenta)

        vec=numpy.array(numpy.dot(self.rot_matrix, vec)) + numpy.array([p_tr[0][1], p_tr[0][2], p_tr[0][3]])


        if surface.param_variable == 0:
            vec = [vec[0], [loop_momenta[1][1], loop_momenta[1][2], loop_momenta[1][3]]]
        else:
            vec = [[loop_momenta[0][1], loop_momenta[0][2], loop_momenta[0][3]], vec[0]]

        return vec



class Evaluator(object):
    """
    The class takes a set of propagators, a list of cut propagators and the signs corresponding to the residue which is going to be picked up
    The class is initialized by saving these informations and immediately decomposing each propagator in the cut momentum basis, by using
    the Momentum_info class
    """
    def __init__(self, propagators, cut, signs):
        super(Evaluator,self).__init__()
        self.propagators=propagators
        self.cut=cut
        info=Momentum_info(cut=cut, n_loops=len(cut), prop=propagators, signs=signs)
        self.combs=info.info_prep()
        self.signs=signs

    def evaluate_propagator(self, n_prop, loop_momenta):
        """
        The evaluate_propagator functions takes the label corresponding to the ropagator one wants to evaluate and the loop momenta and returns
        the propagator as evaluated on the hyperboloid sheets identified by the cut momenta
        """
        #If the chosen propagator is a cut propagator, then the half-propagator is returned
        if n_prop in self.cut:
            return 2.*math.sqrt(
                sum((
                    sum(
                        self.propagators[n_prop].routing[i]*loop_momenta[i][j] for i in range(0, len(self.cut)))+self.propagators[n_prop].p[j])**2 for j in range(1,4))
            )
        #if the chosen propagator is not a cut propagator, then a different evaluation is enacted
        else:
            #assume e.g. that k_1+p and k_1+k_2+p' are cut and k_2+p'' is the chosen momenta, with residue signs [-1,1]
            #then self.combs[n_prop].p[0]=(p''-p'+p)[0]
            energy=self.combs[n_prop].p[0]+sum(
                #and the following is +(+1)*sqrt((\vec3d{k_1+k_2+p'})**2)-(-1)*sqrt((\vec3d{k_1+p})**2)
                self.combs[n_prop].routing[i]*self.signs[i]*math.sqrt(
                    sum((
                        sum(
                            self.propagators[self.cut[i]].routing[k]*loop_momenta[k][j] for k in range(0,len(self.cut))
                        ) + self.propagators[self.cut[i]].p[j])**2 for j in range(1,4))
                ) for i in range(0,len(self.cut)))
            #thus energy=(p''-p'+p)[0]+|\vec3d{k_1+k_2+p'}|+|\vec3d{k_1+p}|
            square_mom=math.sqrt(sum((sum(self.propagators[n_prop].routing[i]*loop_momenta[i][j] for i in range(0, len(self.cut)))+self.propagators[n_prop].p[j])**2 for j in range(1,4)))
            #while square_mom=|\vec3d{k_2+p''}|
            return energy**2-square_mom**2



class trial_function(object):
    """
    Trial function to try the results of the algorithm
    """

    def __init__(self):
        super(trial_function, self).__init__()

    def Eval(self, loop_value, p1, p2, p0, k):
        mod_loop1 = sum((loop_value[k][i] + p1[i]) ** 2 for i in range(0, 3))
        mod_loop2 = sum((loop_value[k][i] + p2[i]) ** 2 for i in range(0, 3))
        return math.sqrt(mod_loop1 + 0 ** 2) + math.sqrt(mod_loop2 + 0 ** 2) - p0

    def Eval2(self, loop_value, p1, p2, p3, p0):
        mod_loop1 = sum((loop_value[0][i] + p1[i]) ** 2 for i in range(0, 3))
        mod_loop2 = sum((loop_value[1][i] + p2[i]) ** 2 for i in range(0, 3))
        mod_loop3 = sum((loop_value[1][i] + loop_value[0][i] + p3[i]) ** 2 for i in range(0, 3))
        return math.sqrt(mod_loop3) + math.sqrt(mod_loop1 + 0 ** 2) + math.sqrt(mod_loop2 + 0 ** 2) - p0

    def Eval3(self, loop_value, p1, p2, p3, p0):
        mod_loop1 = sum((loop_value[0][i] + p1[i]) ** 2 for i in range(0, 3))
        mod_loop2 = sum((loop_value[1][i] + p2[i]) ** 2 for i in range(0, 3))
        mod_loop3 = sum((loop_value[1][i] + loop_value[0][i] + p3[i]) ** 2 for i in range(0, 3))
        return -math.sqrt(mod_loop3) + math.sqrt(mod_loop1 + 0 ** 2) + math.sqrt(mod_loop2 + 0 ** 2) - p0

    def Eval4(self, loop_value, p1, p2, p3, p0):
        mod_loop1 = sum((loop_value[0][i] + p1[i]) ** 2 for i in range(0, 3))
        mod_loop2 = sum((loop_value[1][i] + p2[i]) ** 2 for i in range(0, 3))
        mod_loop3 = sum((loop_value[1][i] + loop_value[0][i] + p3[i]) ** 2 for i in range(0, 3))
        return math.sqrt(mod_loop3) + math.sqrt(mod_loop1 + 0 ** 2) - math.sqrt(mod_loop2 + 0 ** 2) - p0


if __name__ == '__main__':

    prop11=Propagator([1, 0], [0, 0, 0, 0], 0)
    prop21=Propagator([1, 0], [1, 0.25, 0, 0.5], 0)
    prop31=Propagator([1, 1], [0, 0, 0, 0], 0)
    prop41=Propagator([0, 1], [0, 0, 0, 0], 0)
    prop51=Propagator([0, 1], [-1, -0.25, 0, -0.5], 0)

    trial_f=trial_function()

    double_triangle=[prop11,prop21,prop31,prop41,prop51]

    cuts=[[0,2],[0,3],[0,4],[1,2],[1,3],[1,4],[2,4],[2,3]]
    signs=[[-1,1],[1,1],[1,1],[-1,1],[1,1],[1,1],[1,1],[1,1]]

#Determining existence of surfaces!

    for i in range(0,8):
        ifo = Momentum_info(cut=cuts[i],prop=double_triangle, n_loops=2, signs=signs[i])
        ifo.info_prep()


        diag = Diagnostic_tool(propagators=double_triangle,Etot=5,n_loop=2,cut=cuts[i], sign=signs[i])
        for j in range(0,5):
            if ifo.combs[j].p[0]!=0:
                surp=Surface(j,1,1)
                surm=Surface(j,-1,1)
                diag.determine_existence(100,surp)
                if surp.exist==1:
                    print(j)
                    print(cuts[i])
                    print("+")


                print("-------")
                diag.determine_existence(100,surm)

                if surm.exist==1:
                    print(j)
                    print(cuts[i])
                    print("-")

                print("-------")



#Determining if I get the right points on the existent surface

    diag1=Diagnostic_tool(propagators=double_triangle,Etot=5,n_loop=2,cut=[0,2],sign=[-1,1])
    check1=[]

    for i in range(0,100):
        loopm=[[random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)],
          [random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)]]
        surr1=Surface(n=1,ot_sign=-1,sheet=1)
        vec=diag1.get_parametrization(u=random.uniform(0,1),v=random.uniform(0,1),surface=surr1, loop_momenta=loopm)
    #print(vec)
        if vec!=None:
            check1.append(trial_f.Eval(vec,[0,0,0],[0.25,0,0.5],1,0))



    check2=[]


    for i in range(0,100):
        surr1=Surface(n=4,ot_sign=1,sheet=1)
        loopm = [[random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)],
             [random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)]]
        vec=diag1.get_parametrization(u=random.uniform(0,1),v=random.uniform(0,1),surface=surr1, loop_momenta=loopm)

        if vec!=None:
            check2.append(trial_f.Eval2(vec,[0,0,0],[-0.25,0,-0.5],[0,0,0],1))


    check3=[]

    diag2=Diagnostic_tool(propagators=double_triangle,Etot=5,n_loop=2,cut=[0,3],sign=[1,1])


    for i in range(0,100):
        surr1=Surface(n=4,ot_sign=1,sheet=1)
        loopm = [[random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)],
             [random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)]]
        vec=diag2.get_parametrization(u=random.uniform(0,1),v=random.uniform(0,1),surface=surr1, loop_momenta=loopm)

        if vec!=None:
            check3.append(trial_f.Eval(vec,[0,0,0],[-0.25,0,-0.5],1,1))



    check4=[]

    diag3=Diagnostic_tool(propagators=double_triangle,Etot=5,n_loop=2,cut=[0,2],sign=[-1,1])


    for i in range(0,100):
        surr1=Surface(n=4,ot_sign=-1,sheet=1)
        loopm = [[random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)],
             [random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)]]
        vec=diag3.get_parametrization(u=random.uniform(0,1),v=random.uniform(0,1),surface=surr1, loop_momenta=loopm)

        if vec!=None:
            check4.append(trial_f.Eval4(vec,[0,0,0],[-0.25,0,-0.5],[0,0,0],1))



    check5=[]

    diag4=Diagnostic_tool(propagators=double_triangle,Etot=5,n_loop=2,cut=[1,3],sign=[1,1])


    for i in range(0,100):
        surr1=Surface(n=2,ot_sign=1,sheet=1)
        loopm = [[random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)],
             [random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)]]
        vec=diag4.get_parametrization(u=random.uniform(0,1),v=random.uniform(0,1),surface=surr1, loop_momenta=loopm)

        if vec!=None:
            check5.append(trial_f.Eval2(vec,[0.25,0,0.5],[0,0,0],[0,0,0],1))



    check6=[]

    diag5=Diagnostic_tool(propagators=double_triangle,Etot=5,n_loop=2,cut=[1,3],sign=[1,1])


    for i in range(0,100):
        surr1=Surface(n=2,ot_sign=-1,sheet=1)
        x1=random.uniform(0, 1)
        x2=random.uniform(0, 1)
        x3=random.uniform(0, 1)
        x4=random.uniform(0, 1)
        loopm = [[random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), x1],
             [random.uniform(0, 1), x2, x3, x4]]
        u=random.uniform(0, 1)
        v=random.uniform(0, 1)
        vec=diag5.get_parametrization(u=u,v=v,surface=surr1, loop_momenta=loopm)

        if vec!=None:
            check6.append(trial_f.Eval3(vec,[0.25,0,0.5],[0,0,0],[0,0,0],1))
        #print(trial_f.Eval3(vec,[0.25,0,0.5],[0,0,0],[0,0,0],1))
            if trial_f.Eval3(vec,[0.25,0,0.5],[0,0,0],[0,0,0],1)>0.01:
                surr1m = Surface(n=2, ot_sign=-1, sheet=-1)
                vecw = diag5.get_parametrization(u=u, v=v, surface=surr1m, loop_momenta=loopm)
                print(trial_f.Eval3(vecw,[0.25,0,0.5],[0,0,0],[0,0,0],1))





    print("111111111")
    print(check1)
    print("22222222222")
    print(check2)
    print("33333333")
    print(check3)
    print("444444444")
    print(check4)
    print("5555555555")
    print(check5)
    print("6666666")
    print(check6)


#Eval_propagators


    prop111=Propagator([1, 0], [1, 0, 0, 0], 0)
    prop211=Propagator([1, 0], [1, 1, 1, 1], 0)
    prop311=Propagator([1, 1], [2, 0, 0, 0], 0)
    prop411=Propagator([0, 1], [1, 0, 0, 0], 0)
    prop511=Propagator([0, 1], [-1, -0.25, 0, -0.5], 0)

    trial_f=trial_function()

    double_triangle1=[prop111,prop211,prop311,prop411,prop511]


    evals=Evaluator(propagators=double_triangle1, cut=[0,2], signs=[-1,1])


    print(evals.evaluate_propagator(n_prop=3,loop_momenta=[[0,0,0,0],[0,0,0,0]]))
