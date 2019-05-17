import sys
sys.path.insert(0, "../LTD")
import vectors,topologies

mu_sq=1e3*1j

mytop=topologies.hard_coded_topology_collection['manual_Box_massless_1000']
#mytop=topologies.hard_coded_topology_collection['manual_Box_no_ellipse']

qs=[prop.q for prop in mytop.loop_lines[0].propagators]
ms=[prop.m_squared for prop in mytop.loop_lines[0].propagators]

ps=mytop.external_kinematics
    
s= (qs[1]-qs[3]).square()
t= (qs[0]-qs[2]).square()

m1=ps[0].square()
m3=ps[2].square()

def collinear_x(k,p):
        keta = k[0]*p[0] + k[1]*p[1] + k[2]*p[2] + k[3]*p[3]
        peta = p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]
        return  keta/peta
    

def numerator(loop_momenta):
    ct_numerator = 1.0;

    k=vectors.LorentzVector([complex(r,i)for r,i in loop_momenta[0]])
    props=[(k+q).square()- m for q,m in zip(qs,ms)]

    x4 = collinear_x(k+qs[2],ps[3])
    x4b= 1.0 - x4

    
    ct_numerator -= (props[0] * props[1])/((x4b * t + x4 * m1) * (x4 * s + x4b * m3))\
                   * mu_sq* (mu_sq - props[2] - props[3])/(props[2]-mu_sq)/(props[3]- mu_sq)

    #print(ct_numerator)
    return (ct_numerator.real, ct_numerator.imag)
