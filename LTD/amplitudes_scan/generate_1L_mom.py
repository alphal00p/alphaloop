#!/usr/bin/python 

import numpy as np
import itertools
import sys

#Define scalar products
def sp2(p):
    return sp(p,p)

def sp(p,q):
    return p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3]
def v_sq(p):
    return p[1]*p[1]+p[2]*p[2]+p[3]*p[3]

#boost the 4-vector p by the velocity v
def boost(p,v):
    vsq = v[0]**2+v[1]**2+v[2]**2
    if vsq == 0:
        return p
    
    gamma = 1.0/np.sqrt(1-vsq)
    n = v / np.sqrt(vsq)
    

    Lambda = np.array([[gamma, -gamma*v[0], -gamma*v[1], -gamma*v[2]],
            [-gamma*v[0], 1.0+(gamma-1.0)*n[0]**2,(gamma-1.0)*n[0]*n[1],(gamma-1.0)*n[0]*n[2]],
            [-gamma*v[1], (gamma-1.0)*n[1]*n[0],1.0+(gamma-1.0)*n[1]**2,(gamma-1.0)*n[1]*n[2]],
            [-gamma*v[2], (gamma-1.0)*n[2]*n[0],(gamma-1.0)*n[2]*n[1],1.0+(gamma-1.0)*n[2]**2]])
                
    return Lambda.dot(p)

#Start from one momentum and end up with two with mass m1 and m2
#p has to be have positive mass 
def split_mom(p,msq1,msq2):
    #p->p1+p2
    s = sp2(np.array(p))
    p1 = np.random.rand(4)*2.0-1.0
    p2 = -p1
    
    v2_norm = (pow(msq1,2)+pow(s-msq2,2)-2.0*msq1*(s+msq2))/(4.0*s)

    #Redefine p1
    v2 = v_sq(p1)
    p1=p1*np.sqrt(v2_norm/v2)
    p1[0] = np.sign(p[0]) * np.sqrt(v2_norm+msq1)
   
    #Redefine p2
    v2 = v_sq(p2)
    p2=p2*np.sqrt(v2_norm/v2)
    p2[0] = np.sign(p[0]) * np.sqrt(v2_norm+msq2)
  
    v_boost = -np.array(p[1:])/np.sqrt(sp2(p)+v_sq(p))
    return (boost(p1,v_boost),boost(p2,v_boost))
 
#def fix_sp_after_split(p1,p2,q,sp1q):
    
#def generate_moms(external_masses, virtual_masses,final_state_number, seed=0):
def generate_moms(external_masses, final_state_number, seed=0):
    np.random.seed(seed);
    m2=external_masses
    sij = np.sort(np.random.rand(final_state_number-1))
    sij[final_state_number-2] = m2[1+final_state_number]

    #sij=virtual_masses
    externals = []

    #INCOMING
    (p1,p2) = split_mom([1.0,0.0,0.0,0.0],m2[0],m2[1])
    externals += [p1,p2]

    #OUTGOING
    q=np.array([1.0,0.0,0.0,0.0])
    for i in range(2,1+final_state_number):
        (p,q) = split_mom(q,m2[i],sij[i-2])
        externals += [-p]
    externals += [-q]
    
    return externals

if __name__ == "__main__":
    format_type = ""
    
    #Read the input parameters
    if len(sys.argv) == 1 or len(sys.argv) > 4:
        print("Need at least the number of final states.")
        print("The second integer will be used as the seed for the random numbers")
        sys.exit(1)  
    elif len(sys.argv) == 2:
        N=int(sys.argv[1])
        seed = 0
    elif len(sys.argv) == 3:
        N=int(sys.argv[1])
        seed = int(sys.argv[2])
    elif len(sys.argv) == 4:
        N=int(sys.argv[1])
        seed = int(sys.argv[2])
        format_type = sys.argv[3];
    
    np.random.seed(seed)
    
    #print("=> Start generating PS points for a 2->{} process\n".format(N))

    #mass fraction compare to the energy of the c.o.m.
    m2 = [0.0]*(2+N)
    sij = np.sort(np.random.rand(N-1))
    sij[N-2] = m2[1+N]

    externals = generate_moms(m2,3,seed)
    externals[2:] = list(itertools.permutations(externals[2:]))[np.random.randint(5)]
    #print("Sum incoming: {}".format(sum(p for p in externals[:2])))
    #print("Sum outgoing: {}".format(sum(p for p in externals[2:])))

    ss = "msq = ["
    for (i,mom) in enumerate(externals):
        ss += " {},".format(sp2(mom))
    ss = ss[:-1] + " ]"
    print(ss)
    
    for i in range(1+N):
        for j in range(i+1,1+N):
            sij=sp2(externals[i]+externals[j]) 
            print("s({},{})={}".format(i+1,j+1,sij))
    print("Generated PS:")
    for (i,mom) in enumerate(externals):
        if format_type == "fortran":
            if i>=2:
                print("{}d0 {}d0 {}d0 {}d0".format(-mom[0],-mom[1],-mom[2],-mom[3]))
            else:
                print("{}d0 {}d0 {}d0 {}d0".format(mom[0],mom[1],mom[2],mom[3]))
        elif format_type == "mathematica":
            if i==0:
                print("{");
            print("{{{}}}".format((mom).tolist()).replace("[","").replace("]",""))

            if i==len(externals)-1:
                print("}")
            else:
                print(",")
        elif format_type == "rust":
            print("\tvectors.LorentzVector({}),".format(mom.tolist()))
        else:
            print("\t%s" %str(mom.tolist()))
