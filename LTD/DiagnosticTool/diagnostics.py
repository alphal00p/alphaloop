import numpy
import math


#cut is a list which contains the number of the propagators which are cut
#propagators are given as lsits [[loop momenta flow],[ext momenta], mass]
#Etot is the total number of edges


class Momentum_info(object):
        def __init__(self, cut, prop, n_loops):
            super(Momentum_info, self).__init__()
            self.cut=cut
            self.prop=prop
            self.n_loops=n_loops
            self.combs=[]

        def info_prep(self):   #info_prep expresses all the non-cut propagators as linear combinations of the cut propagators

            mat=[self.prop[i][0] for i in self.cut]  #mat contains as row the momentum routing of the cut propagators
            print(self.prop)
            for i in range(0,len(self.prop)):

                    lin_comb=numpy.linalg.solve(numpy.array(mat).T, numpy.array(self.prop[i][0]))  #solving the linear system
                    ext_mom=numpy.array(self.prop[i][1])+sum(-lin_comb[j]*numpy.array(self.prop[j][1]) for j in self.cut)
                    self.combs.append([lin_comb,ext_mom,self.prop[i][2]])

            return self.combs


mominfo=Momentum_info([0,1], [[[1,1],[0.5,0.12,0.32,0.05], 0.25],[[0,1],[0.2,0.1,0.5,0], 0.25],[[1,1],[0.5,0.12,0.32,0.05], 0.25],[[0,1],[0.2,0.1,0.5,0], 0.25],[[1,0],[0,0,0,0], 0.25]], 2)

print(mominfo.info_prep())



class Diagnostic_tool(object):
         def __init__(self, propagators, Etot, n_loop, cut):
            super(Diagnostic_tool, self).__init__()
            self.propagators=propagators
            self.Etot=Etot
            self.hyperboloids=[]
            self.ellipsoids=[]
            self.n_loop=n_loop
            self.rot_matrix=[]
            self.cut=cut
            info=Momentum_info(cut, propagators, n_loop)
            self.combs=info.info_prep()

         def OnShell(self,vec):  #outputs a vector having information on the sign of the energy component of the input vec and its onshell value
            onshellness=vec[0]**2-sum(vec[i]**2 for i in range(1,4))

            if vec[0]==0 or onshellness==0:
                return [0,0]
            elif onshellness>0:
                if vec[0]>0:
                    return [1,1]
                else:
                    return [1,-1]
            else:
                if vec[0]>0:
                    return [-1,1]
                else:
                    return [-1,-1]


         def List_Surfaces(self): #catalogue of surfaces according to their onshell value and energy
            for i in range(0,len(self.combs)):
                if i not in self.cut:
                    print(self.OnShell(self.combs[i][1]))
                    if self.OnShell(self.combs[i][1])[1]!=0:
                        if self.OnShell(self.combs[i][1])[0]==-1:
                            self.hyperboloids.append(i)
                        elif self.OnShell(self.combs[i][1])[0]==1 and self.OnShell(self.combs[i][1])[1]!=0:
                            self.ellipsoids.append(i)

         def print_surfaces(self):
             print(self.hyperboloids)
             print(self.ellipsoids)



         def translation_construction(self,n_surface, loop_mom):  #constructs the translation part of the affine map necessary to map the ellipsoid/hyperboloid to a surface with simmetry axis coinciding with z

             mom1=0
             mom2=0
             tot_mom=0
             chosen_loop=0

             for i in range(0,self.n_loop):

                 if self.combs[n_surface][0][i]!=0: #takes the first available cut momenta q

                    for i1 in range(0,self.n_loop):
			#takes the first available loop momentum on which q depends and
                        if self.propagators[self.cut[i]][0][i1]!=0:  choses it to parametrize the surface (in particular, its z component)
                            chosen_loop=i1
				#constructs the momentum term corresponding to the chosen cut
                            mom1=numpy.array(self.propagators[self.cut[i]][1])+sum(self.propagators[self.cut[i]][0][i2]*numpy.array(loop_mom[i2]) for i2 in range(i1+1,self.n_loop))  
				#constructs the momentum term corresponding to the only other square root (in the surface equation) containing the chosen loop momentum
                            for j in range(i+1,self.n_loop):
                                if self.combs[n_surface][0][j]!=0:
                                    if self.propagators[self.cut[j]][0][i1] != 0:
                                        mom2 = numpy.array(self.propagators[self.cut[j]][1]) + sum(self.propagators[self.cut[j]][0][i3] * numpy.array(loop_mom[i3]) for i3 in range(0, i1-1)+range(i1+1,self.n_loop))
                            if self.propagators[n_surface][0][i1]!=0:
                                mom2=numpy.array(self.propagators[n_surface][1])+sum(self.propagators[n_surface][0][i3] * numpy.array(loop_mom[i3]) for i3 in range(0, i1-1)+range(i1+1,self.n_loop))
                            break
                    break
             #print(mom1)
             #print(mom2)
		#construct the translation vector
             tot_mom=(numpy.array(mom1)-numpy.array(mom2))/2

             return [tot_mom,chosen_loop]
       


		#construct the rotation necessary to rotate the translation vector and make it coincide with the z axis
         def rot_matrix_construction(self, n_surface, loop_mom):
             maj_axis=self.translation_construction(n_surface,loop_mom)[0]
		#vector corresponding to the major axis of the ellipsoid
             mom_norm=1/(sum(maj_axis[j]**2 for j in range(1,4)))**(0.5)
             mom_unit=[mom_norm*maj_axis[j] for j in range(1,4)]
		#z axis
             z_axis=[0,0,1]
		#rotation angle
             angle=math.acos(mom_unit[0])
		#rotation axis
             rot_axis=numpy.cross(x_axis,mom_unit)
		#skew-simmetric and hermitian matrix decomposition of rotation matrix
             A_mat=numpy.matrix([[0,-rot_axis[2],rot_axis[1]],[rot_axis[2],0,-rot_axis[0]],[-rot_axis[1], rot_axis[0], 0]])
             aa_mat=numpy.matrix([[rot_axis[i]*rot_axis[j] for j in range(0,3)]for i in range(0,3)])
		#rotation matrix
             self.rot_matrix=math.cos(angle)*(numpy.identity(3)-aa_mat)+aa_mat+math.sin(angle)*A_mat

		#given the matrix and the translation vector, provides a point on the surface
         def Get_surface_point(self,n_surface, loop_mom, surface_sign):
             p0_sq=(self.combs[n_surface][1][0]/2.)**2
             p_tr=self.translation_construction(n_surface,loop_mom)
		#3d translation vector
             p_tr3d=[p_tr[0][1],p_tr[0][2],p_tr[0][3]]
             mod_p_sq=sum(p_tr[0][i]**2 for i in range(1,4))
		#major and minor axis length of the ellipsoid
             a=p0_sq-mod_p_sq-self.propagators[n_surface][2]**2
             c=(p0_sq/(p0_sq-mod_p_sq))*a
		#parametrization of the ellipsoid vector
             kz=(-1)**(surface_sign)*(c*(1-(1/a)*sum((loop_mom[p_tr[1]][j])**2 for j in range(1,2))))**0.5
             vec=[loop_mom[p_tr[1]][1],loop_mom[p_tr[1]][2],kz]
		#mapping back the vector
             vec=numpy.array(numpy.dot(numpy.matrix(self.rot_matrix).T, vec))-numpy.array(p_tr3d)
             return vec


Diag=Diagnostic_tool([[[1,0],[1,0,0,1], 0.25],[[0,1],[2,1,0,0], 0.25],[[1,1],[3,0,1,0], 0.25],[[0,1],[4,1,1,0], 0.25],[[1,0],[5,0,1,1], 0.25]],5,2, [0,1])
print(Diag.OnShell([-2,0.6,0.6,0.6]))
"""
Diag.List_Surfaces()

Diag.print_surfaces()
"""
print(Diag.translation_construction(2,[[5,5,5,5],[0.2,0.2,0.2,0.2]]))

"""
Diag.rot_matrix_construction(4,[[5,5,5,5],[0.2,0.2,0.2,0.2]])
print(Diag.Get_surface_point(4,[[5,0.1,0.1,0.1],[0.2,0.2,0.2,0.2]],1))
"""
