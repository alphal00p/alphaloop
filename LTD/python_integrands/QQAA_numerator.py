import integrand
import numpy as np
from gamma_chains import gamma_chain

compute_gamma_chain = gamma_chain.compute_chain

class Diagram(integrand.Diagram):

    def __init__(self, *args, **opts):
        super(Diagram, self).__init__(*args, **opts)

    def __call__(self, PS_point, loop_momenta, invariants, spin_index):
        
        # shortcuts
        invs    = invariants
        ps      = PS_point
        l       = loop_momenta[0]
        params  = self.parameters 
        prop    = self.prop
        si      = spin_index  # spin_index helps us to sum over all spins/polarizations

        def compute_polarization(p): # only for massless case
            # spherical coordinates of the spatial part of p
            rho=np.linalg.norm(p.space())
            theta=np.arccos(p[3]/rho)
            phi = np.arctan2(p[2],p[1])
            phase = np.exp(1j*phi)

            # define gamma matrices in Dirac representation, only for the use of computing polarization vectors and possibly change of basis between Dirac and Weyl representations
            gamma0 = np.array(
                      [[1., 0., 0., 0.],
                       [0., 1., 0., 0.],
                       [0., 0., -1., 0.],
                       [0., 0., 0., -1.]])
            gamma5 = np.array(
                      [[0., 0., 1., 0.],
                       [0., 0., 0., 1.],
                       [1., 0., 0., 0.],
                       [0., 1., 0., 0.]])

            # spinors as documented in Britto 2011, rewritten in Dirac representation, note the four-momentum is contravariant whereas in Britto they are covariant
            u_plus = np.sqrt(rho/2.)*np.array([-np.sqrt(1-np.cos(theta))/phase, np.sqrt(1+np.cos(theta)), np.sqrt(1-np.cos(theta))/phase, -np.sqrt(1+np.cos(theta))])
            u_minus = -np.sqrt(rho/2.)*np.array([np.sqrt(1+np.cos(theta)), np.sqrt(1-np.cos(theta))*phase, np.sqrt(1+np.cos(theta)), np.sqrt(1-np.cos(theta))*phase]) 
            u_bar_plus = np.matmul(np.conj(u_plus),gamma0)
            u_bar_minus = np.matmul(np.conj(u_minus),gamma0)

            # photon polarization vectors
            rotation_matrix = np.matmul(
                np.array([[1., 0., 0., 0.],
                          [0., np.cos(phi), -np.sin(phi), 0.],
                          [0., np.sin(phi), np.cos(phi), 0.],
                          [0., 0., 0., 1.]]),
                np.array([[1., 0., 0., 0.],
                          [0., np.cos(theta), 0., np.sin(theta)],
                          [0., 0., 1., 0.],
                          [0., -np.sin(theta), 0., np.cos(theta)]])
                )
            a1 = np.matmul(rotation_matrix, np.array([0.,1.,0.,0.]))
            a2 = np.matmul(rotation_matrix, np.array([0.,0.,1.,0.]))

            return [u_plus, u_minus, u_bar_plus, u_bar_minus, a1, a2]

        u_p1 = compute_polarization(ps[0])[0:2]
        vbar_p2 = compute_polarization(ps[1])[2:4]
        a_mu = compute_polarization(-ps[3])[4:6]
        a_nu = compute_polarization(-ps[2])[4:6]

        # For the call to gamma_chains, it is necessary to cast LorentzVectors into lists
        vecs = map(list, [a_nu[si[2]], a_mu[si[3]], ps[0]+ps[3], ps[0]-l, -ps[1]-l, ps[0]+ps[3]-l])
        u = list(u_p1[si[0]])
        vbar = list(vbar_p2[si[1]])

        # Tree diagram    
        if self.identifier == 'D0':
            return (
                -complex(0.,1.)*params['alpha_ew']*4*np.pi
                *(1./invs['s14'])
                *compute_gamma_chain(vbar,u,[1, 3, 2], vecs)
            )

        # Define a common propagator factor
        box = prop(l, 0.)*prop(ps[0]-l, 0.)*prop(ps[1]+l, 0.)*prop(ps[0]+ps[3]-l, 0.)

        # Loop diagrams
                
        if self.identifier == 'D1N':
            return (
                -params['alpha_ew']*4*np.pi*params['alpha_s']*4*np.pi*params['C_F']
                *(1./invs['s14'])
                *prop(ps[1]+l, 0.)*prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.)
                *compute_gamma_chain(vbar,u, [-1, 5, 1, 6, -1, 3, 2], vecs)
                /box
            )

        elif self.identifier == 'D2N':
            return (
                -params['alpha_ew']*4*np.pi*params['alpha_s']*4*np.pi*params['C_F']
                *(1./invs['s14'])
                *prop(ps[0]-l, 0.)*prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.)
                *compute_gamma_chain(vbar, u, [1, 3, -1, 6, 2, 4, -1], vecs)
                /box
            )

        elif self.identifier == 'D3N':
            return (
                -params['alpha_ew']*4*np.pi*params['alpha_s']*4*np.pi*params['C_F']
                *prop(ps[1]+l, 0.)*prop(ps[0]-l, 0.)*prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.) 
                *compute_gamma_chain(vbar, u, [-1, 5, 1, 6, 2, 4, -1], vecs)
                /box
            )

        elif self.identifier == 'D4N':
            return (
                -params['alpha_ew']*4*np.pi*params['alpha_s']*4*np.pi*params['C_F']
                *(1./invs['s14'])**2
                *prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.)
                *compute_gamma_chain(vbar, u, [1, 3, -1, 6, -1, 3, 2], vecs)
                /box
            )

        else:
            raise NotImplementedError

class Counterterm(integrand.Counterterm):

    def __init__(self, *args, **opts):
        super(Counterterm, self).__init__(*args, **opts)

    def __call__(self, PS_point, loop_momenta, invariants, spin_index):

        # shortcuts
        invs    = invariants
        ps      = PS_point
        l       = loop_momenta[0]
        params  = self.parameters 
        prop    = self.prop


        # Instantialize diagrams for use
        D0 = Diagram('D0', self.parameters)
        D1 = Diagram('D1N', self.parameters)
        D2 = Diagram('D2N', self.parameters)
        D4 = Diagram('D4N', self.parameters)

        # Define a common propagator factor
        box = prop(l, 0.)*prop(ps[0]-l, 0.)*prop(ps[1]+l, 0.)*prop(ps[0]+ps[3]-l, 0.)

        # For the call to gamma_chains, it is necessary to cast LorentzVectors into lists
        if self.identifier == 'UV1N':
            return (
                D1(PS_point, loop_momenta, invariants, spin_index)
                *(prop(ps[1]+l, 0.)*prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.))**(-1)
                *(prop(l, invs['m_UV2']))**3
            )

        elif self.identifier == 'UV2N':
            return (
                D2(PS_point, loop_momenta, invariants, spin_index)
                *(prop(ps[0]-l, 0.)*prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.))**(-1)
                *(prop(l, invs['m_UV2']))**3
            )

        elif self.identifier == 'UV4N':
            return (
                D4(PS_point, loop_momenta, invariants, spin_index)
                *(prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.))**(-1)
                *(prop(l-ps[0]/2-ps[3]/2, invs['m_UV2']))**2
            )

        elif self.identifier == 'SoftN':
            return (
                D0(PS_point, loop_momenta, invariants, spin_index)
                *2*complex(0.,1.)*params['alpha_s']*4*np.pi*params['C_F']
                *invs['s12']
                *prop(l, 0.)*prop(ps[1]+l, 0.)*prop(ps[0]-l,0.)
                /box
            )

        elif self.identifier == 'Coll1N':
            return (
                D0(PS_point, loop_momenta, invariants, spin_index)
                *(-2)*complex(0.,1.)*params['alpha_s']*4*np.pi*params['C_F']
                *(prop(ps[0]-l, 0.)*prop(l, 0.)-(prop(l, invs['m_UV2']))**2)
                /box
            )

        elif self.identifier == 'Coll2N':
            return (
                D0(PS_point, loop_momenta, invariants, spin_index)
                *(-2)*complex(0.,1.)*params['alpha_s']*4*np.pi*params['C_F']
                *(prop(ps[1]+l, 0.)*prop(l, 0.)-(prop(l, invs['m_UV2']))**2) 
                /box
            )
        else:
            raise NotImplementedError

class Amplitude(integrand.Amplitude):

    def __init__ (self, identifier, n_loops=1,**opts):
        super(Amplitude, self).__init__(identifier, n_loops=n_loops, **opts) 

        if self.identifier != 'qq_aa':
            raise BaseException("AmplitudeQQAA only implements type identifier 'qq_aa', not '%s'."%self.identifier)
        if self.n_loops not in [0, 1]:
            raise BaseException("AmplitudeQQAA only available at tree and one-loop for now")

        # Now create contributing diagrams
        self.diagrams = []
        if self.n_loops==1:
            self.diagrams.append(Diagram('D1N', self.parameters))
            self.diagrams.append(Diagram('D2N', self.parameters))
            self.diagrams.append(Diagram('D3N', self.parameters))
            self.diagrams.append(Diagram('D4N', self.parameters))
        elif self.n_loops==2:
            # Example of where you'd implement the two-loop amplitude
            raise NotImplementedError
        else:
            raise NotImplementedError

        # And now the counterterms
        self.counterterms = []
        if self.n_loops==1:
            self.counterterms.append(Counterterm('UV1N', self.parameters))
            self.counterterms.append(Counterterm('UV2N', self.parameters))
            self.counterterms.append(Counterterm('UV4N', self.parameters))
            self.counterterms.append(Counterterm('SoftN', self.parameters))
            self.counterterms.append(Counterterm('Coll1N', self.parameters))
            self.counterterms.append(Counterterm('Coll2N', self.parameters))
        elif self.n_loops==2:
            # Example of where you'd implement the two-loop amplitude
            raise NotImplementedError
        else:
            raise NotImplementedError

    def compute_invariants(self, PS_point, loop_momenta):

        invariants = {}
        for i, pi in enumerate(PS_point):
            for j, pj in enumerate(PS_point):
                invariants['s%d%d'%(i+1,j+1)] = (pi+pj).square()
        
        invariants['m_UV2'] = 1000.
        
        for i, pi in enumerate(PS_point):
            invariants['sk%d'%(i+1)] = (loop_momenta[0]+pi).square()

        return invariants

    def __call__(self, PS_point, loop_momenta, spin_index=[0,0,0,0]):
        
        res = 0.
        
        invariants = self.compute_invariants(PS_point, loop_momenta)

        for diag in self.diagrams:
            res += diag(PS_point, loop_momenta, invariants, spin_index)

        for CT in self.counterterms:
            res -= CT(PS_point, loop_momenta, invariants, spin_index)

        return res 

