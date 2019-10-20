#!/usr/bin/env python
import os
import sys
pjoin = os.path.join
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(root_path,os.pardir))
sys.path.insert(0, pjoin(root_path,os.pardir,os.pardir))

import vectors; import math; sqrt = math.sqrt

import ltd_utils
from ltd_utils import TopologyGenerator, TopologyCollection
from analytical_expressions  import ladder_phi, analytic_four_point_ladder, analytic_three_point_ladder, analytic_two_point_ladder

euclidean_2_to_2_points = [
    vectors.LorentzVector([7.,0.,0.,4.*sqrt(7.)])*(1./sqrt(63.)),
    vectors.LorentzVector([-7.,0.,0.,5.*sqrt(7.)])*(1./sqrt(63.)),
    vectors.LorentzVector([-25.,-sqrt(78.),-4.*sqrt(39.),-4.*sqrt(7)])*(1./sqrt(63.)),
    vectors.LorentzVector([25.,sqrt(78.),4.*sqrt(39.),-5.*sqrt(7)])*(1./sqrt(63.))
]

physical_2_to_2_points = [
    vectors.LorentzVector([29.,0.,0.,sqrt(721.)])*(1./sqrt(120.)),
    vectors.LorentzVector([31.,0.,0.,-sqrt(721.)])*(1./sqrt(120.)),
    vectors.LorentzVector([-29.,4.*sqrt(366./103.),6.*sqrt(854./103.),-43.*sqrt(7./103.)])*(1./sqrt(120.)),
    vectors.LorentzVector([-31.,-4.*sqrt(366./103.),-6.*sqrt(854./103.),43.*sqrt(7./103.)])*(1./sqrt(120.))
]


def test_PS_point(PS):
    print "M1=",PS[0].square()
    print "M2=",PS[1].square()
    print "M3=",PS[2].square()
    print "M4=",PS[3].square()
    print "s=",(PS[0]+PS[1]).square()
    print "t=",(PS[0]+PS[2]).square()
    print "Sum=",(PS[0]+PS[1]+PS[2]+PS[3])

def generate_two_to_two_at_fixed_angle_positive_s(
    M1sq, M2sq, M3sq, M4sq, s, costheta, cosphi):
    """ Generates a point with given masses, positive s and angles theta and phi between p1 and p3."""
    sintheta = sqrt(1.-costheta**2)
    sinphi = sqrt(1.-cosphi**2)
    return [
  vectors.LorentzVector([
      (M1sq - M2sq + s)/(2.*sqrt(s)), 
      0, 
      0, 
      sqrt(M1sq**2 + (M2sq - s)**2 - 2*M1sq*(M2sq + s))/(2.*sqrt(s)), \

  ]),
  vectors.LorentzVector([
      (-M1sq + M2sq + s)/(2.*sqrt(s)), 
      0, 
      0, 
      -sqrt(M1sq**2 + (M2sq - s)**2 - 2*M1sq*(M2sq + s))/(2.*sqrt(s)), 
  ]),
  vectors.LorentzVector([
      -(M3sq - M4sq + s)/(2.*sqrt(s)), 
      (cosphi*sqrt(M3sq**2 + (M4sq - s)**2 - 2*M3sq*(M4sq + \
s))*sintheta)/(2.*sqrt(s)), 
      (sqrt(M3sq**2 + (M4sq - s)**2 - 2*M3sq*(M4sq + \
s))*sinphi*sintheta)/(2.*sqrt(s)), 
      (costheta*sqrt(M3sq**2 + (M4sq - s)**2 - 2*M3sq*(M4sq + \
s)))/(2.*sqrt(s)), 
  ]),
  vectors.LorentzVector([
      -(-M3sq + M4sq + s)/(2.*sqrt(s)), 
      -(cosphi*sqrt(M3sq**2 + (M4sq - s)**2 - 2*M3sq*(M4sq + \
s))*sintheta)/(2.*sqrt(s)), 
      -(sqrt(M3sq**2 + (M4sq - s)**2 - 2*M3sq*(M4sq + \
s))*sinphi*sintheta)/(2.*sqrt(s)), 
      -(costheta*sqrt(M3sq**2 + (M4sq - s)**2 - 2*M3sq*(M4sq + \
s)))/(2.*sqrt(s)), 
  ]),
    ]

def generate_two_to_two_at_fixed_t_positive_s(
    M1sq, M2sq, M3sq, M4sq, s, t, cosphi):
    """ Generates a point with given masses, positive s, t and angle phi between p1 and p3."""
    sinphi = sqrt(1.-cosphi**2)
    E3 = (M3sq - M4sq + ((M1sq - M2sq + s)/(2.*sqrt(s)) + (-M1sq + M2sq + \
s)/(2.*sqrt(s)))**2)/(2.*((M1sq - M2sq + s)/(2.*sqrt(s)) + (-M1sq + \
M2sq + s)/(2.*sqrt(s))))
    E4 = (-M3sq + M4sq + ((M1sq - M2sq + s)/(2.*sqrt(s)) + (-M1sq + M2sq + \
s)/(2.*sqrt(s)))**2)/(2.*((M1sq - M2sq + s)/(2.*sqrt(s)) + (-M1sq + \
M2sq + s)/(2.*sqrt(s))))
    p3r = -sqrt(((M1sq**2 + (M2sq - s)**2 - 2*M1sq*(M2sq + s))*(M4sq**2 + \
(-M3sq + ((M1sq - M2sq + s)/(2.*sqrt(s)) + (-M1sq + M2sq + \
s)/(2.*sqrt(s)))**2)**2 - 2*M4sq*(M3sq + ((M1sq - M2sq + \
s)/(2.*sqrt(s)) + (-M1sq + M2sq + \
s)/(2.*sqrt(s)))**2)))/s)/(2.*sqrt(((M1sq**2 + (M2sq - s)**2 - \
2*M1sq*(M2sq + s))*((M1sq - M2sq + s)/(2.*sqrt(s)) + (-M1sq + M2sq + \
s)/(2.*sqrt(s)))**2)/s))
    costheta = (-2*sqrt(M1sq**2 + (M2sq - s)**2 - 2*M1sq*(M2sq + s))*((M1sq - M2sq + \
s)/(2.*sqrt(s)) + (-M1sq + M2sq + s)/(2.*sqrt(s)))*(((M1sq - M2sq + \
s)**2*(-M1sq + M2sq + s))/(8.*s**1.5) + ((-M1sq + M2sq + s)*(-M3sq + \
(M1sq**2 + (M2sq - s)**2 - 2*M1sq*(M2sq + s))/(4.*s) + \
t))/(2.*sqrt(s)) + ((M1sq - M2sq + s)*(-M4sq + (-M1sq + M2sq + \
s)**2/(4.*s) + (M1sq**2 + (M2sq - s)**2 - 2*M1sq*(M2sq + s))/(4.*s) + \
t))/(2.*sqrt(s))))/(sqrt(s)*sqrt(((M1sq**2 + (M2sq - s)**2 - \
2*M1sq*(M2sq + s))*((M1sq - M2sq + s)/(2.*sqrt(s)) + (-M1sq + M2sq + \
s)/(2.*sqrt(s)))**2)/s)*sqrt(((M1sq**2 + (M2sq - s)**2 - 2*M1sq*(M2sq \
+ s))*(M4sq**2 + (-M3sq + ((M1sq - M2sq + s)/(2.*sqrt(s)) + (-M1sq + \
M2sq + s)/(2.*sqrt(s)))**2)**2 - 2*M4sq*(M3sq + ((M1sq - M2sq + \
s)/(2.*sqrt(s)) + (-M1sq + M2sq + s)/(2.*sqrt(s)))**2)))/s))
    print(1.-costheta**2)
    sintheta = sqrt(1.-costheta**2)
    return [
  vectors.LorentzVector([
      (M1sq - M2sq + s)/(2.*sqrt(s)), 
      0, 
      0, 
      sqrt(M1sq**2 + (M2sq - s)**2 - 2*M1sq*(M2sq + s))/(2.*sqrt(s)), \

  ]),
  vectors.LorentzVector([
      (-M1sq + M2sq + s)/(2.*sqrt(s)), 
      0, 
      0, 
      -sqrt(M1sq**2 + (M2sq - s)**2 - 2*M1sq*(M2sq + s))/(2.*sqrt(s)), 
  ]),
  vectors.LorentzVector([
      -E3, 
      -(cosphi*sqrt(1 - costheta**2)*p3r), 
      -(sqrt(1 - costheta**2)*p3r*sinphi), 
      -(costheta*p3r), 
  ]),
  vectors.LorentzVector([
      -E4, 
      cosphi*sqrt(1 - costheta**2)*p3r, 
      sqrt(1 - costheta**2)*p3r*sinphi, 
      costheta*p3r, 
  ]),
    ]

def generate_two_to_two_at_fixed_angle_negative_s(
    M1sq, M2sq, M3sq, M4sq, minus_s, costheta, cosphi):
    """ Generates a point with given masses, positive s and angles theta and phi between p1 and p3."""
    s = -minus_s
    sintheta = sqrt(1.-costheta**2)
    sinphi = sqrt(1.-cosphi**2)
    sectheta = 1./sintheta
    cos2theta = 2.*costheta**2-1.
    cos4theta = 2.*(2.*costheta**2-1.)**2-1.
    tantheta = sintheta/costheta

    return [
  vectors.LorentzVector([
      sqrt(M1sq**2 + 2*M1sq*(-M2sq + s) + (M2sq + s)**2)/(2.*sqrt(s)), 
      0, 
      0, 
      (-M1sq + M2sq + s)/(2.*sqrt(s)), 
  ]),
  vectors.LorentzVector([
      -sqrt(M1sq**2 + 2*M1sq*(-M2sq + s) + (M2sq + s)**2)/(2.*sqrt(s)), 
      0, 
      0, 
      (M1sq - M2sq + s)/(2.*sqrt(s)), 
  ]),
  vectors.LorentzVector([
      -sqrt(M3sq**2 - 2*M3sq*M4sq + 2*cos2theta*M3sq*s + (M4sq + s)**2)/(2.*sqrt(costheta**2*s)), 
      -(cosphi*(-M3sq + M4sq + s)*tantheta)/(2.*sqrt(s)), 
      -((-M3sq + M4sq + s)*sinphi*tantheta)/(2.*sqrt(s)), 
      -(-M3sq + M4sq + s)/(2.*sqrt(s)), 
  ]),
  vectors.LorentzVector([
      sqrt(M3sq**2 - 2*M3sq*M4sq + 2*cos2theta*M3sq*s + (M4sq + s)**2)/(2.*sqrt(costheta**2*s)), 
      (cosphi*(-M3sq + M4sq + s)*tantheta)/(2.*sqrt(s)), 
      ((-M3sq + M4sq + s)*sinphi*tantheta)/(2.*sqrt(s)), 
      -(M3sq - M4sq + s)/(2.*sqrt(s)), 
  ]),
    ]

def generate_two_to_two_at_fixed_t_negative_s(
    M1sq, M2sq, M3sq, M4sq, minus_s, t, cosphi):
    """ Generates a point with given masses, negative s, t and angle phi between p1 and p3.""" 
    s = -minus_s 
    sinphi = sqrt(1.-cosphi**2)

    E1=sqrt(M1sq**2 + 2*M1sq*(-M2sq + s) + (M2sq + s)**2)/(2.*sqrt(s))
    p1z=(-M1sq + M2sq + s)/(2.*sqrt(s))
    p2z=(M1sq - M2sq + s)/(2.*sqrt(s))

    E3=-(M4sq*p1z + M3sq*p2z + p1z**2*p2z + p1z*p2z**2 + E1**2*(p1z + p2z) - \
p1z*t - p2z*t)/(2.*E1*(p1z + p2z))
    E4=(M4sq*p1z + M3sq*p2z + p1z**2*p2z + p1z*p2z**2 + E1**2*(p1z + p2z) - \
p1z*t - p2z*t)/(2.*E1*(p1z + p2z))
    p3r=-sqrt(E1**4*(p1z + p2z)**2 + (M4sq*p1z + M3sq*p2z + (p1z + \
p2z)*(p1z*p2z - t))**2 + 2*E1**2*(p1z + p2z)*(M4sq*p1z - M3sq*(2*p1z \
+ p2z) + (p1z + p2z)*(p1z*p2z - t)))/(2.*sqrt(E1**2*(p1z + p2z)**2))
    p4z=-p1z - p2z
    costheta=(sqrt(E1**2*(p1z + p2z)**2)*(-M3sq + M4sq + (p1z + p2z)**2))/((p1z + \
p2z)*sqrt(E1**4*(p1z + p2z)**2 + (M4sq*p1z + M3sq*p2z + (p1z + \
p2z)*(p1z*p2z - t))**2 + 2*E1**2*(p1z + p2z)*(M4sq*p1z - M3sq*(2*p1z \
+ p2z) + (p1z + p2z)*(p1z*p2z - t))))

    return [
  vectors.LorentzVector([
      E1, 
      0, 
      0, 
      p1z, 
  ]),
  vectors.LorentzVector([
      -E1, 
      0, 
      0, 
      p2z, 
  ]),
  vectors.LorentzVector([
      E3, 
      cosphi*sqrt(1 - costheta**2)*p3r, 
      sqrt(1 - costheta**2)*p3r*sinphi, 
      costheta*p3r, 
  ]),
  vectors.LorentzVector([
      E4, 
      -(cosphi*sqrt(1 - costheta**2)*p3r), 
      -(sqrt(1 - costheta**2)*p3r*sinphi), 
      -(costheta*p3r) + p4z, 
  ])
    ]

if __name__ == '__main__':
    
#    print "Hard-coded benchmark points:"
#    test_PS_point(euclidean_2_to_2_points)
#    test_PS_point(physical_2_to_2_points)
#    print ""
#    print "Generated points wiht positive s:"
#    test_PS_point(generate_two_to_two_at_fixed_angle_positive_s(-1.,-2.,-3.,-4.,5,0.1,0.2))
#    test_PS_point(generate_two_to_two_at_fixed_t_positive_s(1., 2., 3., 4., 30., -5., 0.2))
#    print "Generated points wiht negative s:"
#    test_PS_point(generate_two_to_two_at_fixed_angle_negative_s(-1.,-2.,-3.,-4.,-9.,2./3.,1./3.))
#    test_PS_point(generate_two_to_two_at_fixed_t_negative_s(1., 2., 3., 4., -15.,-5., 0.2))
    
    # Now generating topologies for each of the points to be scanned.
    
    box_file_name = 'box_scan.yaml'
    double_box_file_name = 'double_box_scan.yaml'
    triple_box_file_name = 'triple_box_scan_yaml'
    n_points = 21
    n_analytic_points = 210
    M1sq, M2sq, M3sq, M4sq, s, cosphi = -5., -1., -1., -1., -1., 1./3.
    

#    costheta_angles = [-1.+ 2.*(i/float(n_points)) for i in range(n_points+1)]
    # Trim the endpoints
#    costheta_angles = costheta_angles[1:-1]
    
    t_values = [-6. + 5.*i for i in range(n_points+1)]
    t_values_analytic_line = [-6. + 0.5*i for i in range(n_analytic_points+1)]

    box = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 4), ('p2', 4, 3), ('p3', 3, 2), ('p4', 2, 1),
    ])
    doublebox = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 6)
    ])
    triplebox = TopologyGenerator([
        ('q1', 101, 1), ('q2', 102, 2), ('q3', 103, 3), ('q4', 104, 4),
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 9), ('p6', 9, 8), ('p7', 8, 6),
        ('p8', 8, 4), ('p9', 4, 3), ('p10', 3, 9)
    ])

    box_topologies = TopologyCollection()
    doublebox_topologies = TopologyCollection()
    triplebox_topologies = TopologyCollection()
   
    analytic_box = []
    analytic_doublebox = []
    analytic_triplebox = []

#    for i_cos, costheta in enumerate(costheta_angles):
    for i, t in enumerate(t_values):
#        PS = generate_two_to_two_at_fixed_angle_negative_s(M1sq, M2sq, M3sq, M4sq, s, costheta, cosphi)
        PS = generate_two_to_two_at_fixed_t_negative_s(M1sq, M2sq, M3sq, M4sq, s, t, cosphi)
#        PS = generate_two_to_two_at_fixed_t_negative_s(-5., -1., -1., -1., -1.,106., 1./3.)
        external_momenta = { 'q1': PS[0], 'q2': PS[1] , 'q3': PS[2], 'q4': -PS[0]-PS[1]-PS[2] }
        loop_masses = {} # Massless internal lines

        for topology_name, master_topo, loop_momenta_names, topologies_collection, analytic_container in [
            ('box', box, ('p4',), box_topologies, analytic_box),
            ('doublebox', doublebox, ('p4','p2'), doublebox_topologies, analytic_doublebox),
            ('triplebox', triplebox, ('p4','p2','p6'), triplebox_topologies, analytic_triplebox),
        ]:
            #print(topology_name, t)
            #test_PS_point(PS)
            #print(PS[0])
            #print(PS[1])
            #print(PS[2])
            #print(PS[3])
            n_loops = len(loop_momenta_names)
            analytic_result = analytic_four_point_ladder( 
                    PS[0].square(), PS[1].square(), PS[2].square(), PS[3].square(),
                    (PS[0]+PS[1]).square(), (PS[1]+PS[2]).square(), n_loops)
#            analytic_container.append((costheta,analytic_result))
            analytic_container.append((t,analytic_result))            
            #print('aa',analytic_result)
            low_level_topology = master_topo.create_loop_topology(
                "scan_%d"%(i+1), 
                ext_mom=external_momenta, 
                mass_map=loop_masses, 
                loop_momenta_names=loop_momenta_names, 
                analytic_result = analytic_result
            )
            topologies_collection.add_topology(low_level_topology,entry_name = "scan_%d"%(i+1))
#            print('ana',low_level_topology.analytic_result) 
   
    analytic_line_box = []
    analytic_line_doublebox = []
    analytic_line_triplebox = []
    for i, t in enumerate(t_values_analytic_line):
        u = M1sq+M2sq+M3sq+M4sq-s-t
        if u==0.: continue
        analytic_line_box.append((t, analytic_four_point_ladder( M1sq, M2sq, M3sq, M4sq, s, u ,1)))
        analytic_line_doublebox.append((t, analytic_four_point_ladder( M1sq, M2sq, M3sq, M4sq, s, u,2)))
        analytic_line_triplebox.append((t, analytic_four_point_ladder( M1sq, M2sq, M3sq, M4sq, s, u,3)))


    # Now write out the files
    box_topologies.export_to(os.path.join(root_path,'box','topologies.yaml'))
    doublebox_topologies.export_to(os.path.join(root_path,'doublebox','topologies.yaml'))
    triplebox_topologies.export_to(os.path.join(root_path,'triplebox','topologies.yaml'))

    box_analytic = open(os.path.join(root_path,'box','analytic_result.dat'),'w')
    box_analytic.write('\n'.join('%.16e %.16e'%(x,float(abs(y))) for x, y in analytic_box))
    box_analytic.close()

    doublebox_analytic = open(os.path.join(root_path,'doublebox','analytic_result.dat'),'w')
    doublebox_analytic.write('\n'.join('%.16e %.16e'%(x,float(abs(y))) for x, y in analytic_doublebox))
    doublebox_analytic.close()

    triplebox_analytic = open(os.path.join(root_path,'triplebox','analytic_result.dat'),'w')
    triplebox_analytic.write('\n'.join('%.16e %.16e'%(x,float(abs(y))) for x, y in analytic_triplebox))
    triplebox_analytic.close()

    box_line_analytic = open(os.path.join(root_path,'box','analytic_line_result.dat'),'w')
    box_line_analytic.write('\n'.join('%.16e %.16e'%(x,float(abs(y))) for x, y in analytic_line_box))
    box_line_analytic.close()

    doublebox_line_analytic = open(os.path.join(root_path,'doublebox','analytic_line_result.dat'),'w')
    doublebox_line_analytic.write('\n'.join('%.16e %.16e'%(x,float(abs(y))) for x, y in analytic_line_doublebox))
    doublebox_line_analytic.close()

    triplebox_line_analytic = open(os.path.join(root_path,'triplebox','analytic_line_result.dat'),'w')
    triplebox_line_analytic.write('\n'.join('%.16e %.16e'%(x,float(abs(y))) for x, y in analytic_line_triplebox))
    triplebox_line_analytic.close()
