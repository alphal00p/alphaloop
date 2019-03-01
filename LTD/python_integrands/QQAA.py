import integrand
import numpy as np
from gamma_chains import gamma_chain

compute_gamma_chain = gamma_chain.compute_chain

class Diagram(integrand.Diagram):

    def __init__(self, *args, **opts):
        super(Diagram, self).__init__(*args, **opts)

    def __call__(self, PS_point, loop_momenta, invariants):
        
        # shortcuts
        invs    = invariants
        ps      = PS_point
        l       = loop_momenta[0]
        params  = self.parameters 
        prop    = self.prop

        # spinor wave function as documented in Britto 2011
        def polar(v):
            rho=np.linalg.norm(v)
            theta=np.arccos(v[2]/rho)
            if np.sin(theta) > 10**(-3): # avoid getting 0 
                phase = v[1]/(rho*np.sin(theta))+complex(0., 1.)*v[0]/(rho*np.sin(theta))
            else:
                phase = 1
            return [rho, np.cos(theta), phase]

        p2pol = polar(ps[1].space())
        p1pol = polar(ps[0].space())
        vbar = np.asarray([-np.sqrt(p2pol[0]*(1+p2pol[1])), complex(0., 1.)*np.sqrt(p2pol[0]*(1-p2pol[1]))*p2pol[2], 0, 0])
        u    = np.asarray([0, 0, -np.sqrt(p1pol[0]*(1+p1pol[1])), -complex(0., 1.)*np.sqrt(p1pol[0]*(1-p1pol[1]))/p1pol[2]]) 
        # Convert the eigenfunctions to Dirac representation
        gamma5 = np.asarray(
                  [[-1, 0, 0, 0],
                  [0, -1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])
        gamma0 = np.asarray(
                  [[0, 0, 1, 0],
                  [0, 0, 0, 1],
                  [1, 0, 0, 0],
                  [0, 1, 0, 0]])
        conversion_u = 1/np.sqrt(2)*(np.identity(4)-np.matmul(gamma5, gamma0))
        conversion_vbar = 1/np.sqrt(2)*(np.identity(4)+np.matmul(gamma5, gamma0))
        u = list(np.matmul(conversion_u, u))
        vbar = list(np.matmul(vbar, conversion_vbar))

        # photon polarization vectors
        def perp(vt, ref1 = [0., 0., 1.], ref2 = [0., 1., 0.]):
            if np.cross(vt, ref1).any():
                res = np.cross(vt, ref1)/np.linalg.norm(np.cross(vt, ref1))
            else :
                res = np.cross(vt, ref2)/np.linalg.norm(np.cross(vt, ref2))
            return res

        A3 = [0] + list(perp(ps[2].space()))
        A4 = [0] + list(perp(ps[3].space()))

        # For the call to gamma_chains, it is necessary to cast LorentzVectors into lists
        vecs = [A3, A4, list(ps[0]+ps[3]), list(ps[0]-l), list(-ps[1]-l), list(ps[0]+ps[3]-l)]

        # Tree diagram    
        if self.identifier == 'D0':
            return (
                -complex(0.,1.)*params['alpha_ew']*4*np.pi
                *(1./invs['s14'])
                *compute_gamma_chain(vbar,u,[1, 3, 2], vecs)
            )

        # Loop diagrams
                
        elif self.identifier == 'D1':
            return (
                -params['alpha_ew']*4*np.pi*params['alpha_s']*4*np.pi*params['C_F']
                *(1./invs['s14'])
                *prop(ps[1]+l, 0.)*prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.)
                *compute_gamma_chain(vbar,u, [-1, 5, 1, 6, -1, 3, 2], vecs)
            )

        elif self.identifier == 'D2':
            return (
                -params['alpha_ew']*4*np.pi*params['alpha_s']*4*np.pi*params['C_F']
                *(1./invs['s14'])
                *prop(ps[0]-l, 0.)*prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.)
                *compute_gamma_chain(vbar, u, [1, 3, -1, 6, 2, 4, -1], vecs)
            )

        elif self.identifier == 'D3':
            return (
                -params['alpha_ew']*4*np.pi*params['alpha_s']*4*np.pi*params['C_F']
                *prop(ps[1]+l, 0.)*prop(ps[0]-l, 0.)*prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.) 
                *compute_gamma_chain(vbar, u, [-1, 5, 1, 6, 2, 4, -1], vecs)
            )

        elif self.identifier == 'D4':
            return (
                -params['alpha_ew']*4*np.pi*params['alpha_s']*4*np.pi*params['C_F']
                *(1./invs['s14'])**2
                *prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.)
                *compute_gamma_chain(vbar, u, [1, 3, -1, 6, -1, 3, 2], vecs)
            )

        else:
            raise NotImplementedError

class Counterterm(integrand.Counterterm):

    def __init__(self, *args, **opts):
        super(Counterterm, self).__init__(*args, **opts)

    def __call__(self, PS_point, loop_momenta, invariants):

        # shortcuts
        invs    = invariants
        ps      = PS_point
        l       = loop_momenta[0]
        params  = self.parameters 
        prop    = self.prop


        # Instantialize diagrams for use
        D0 = Diagram('D0', self.parameters)
        D1 = Diagram('D1', self.parameters)
        D2 = Diagram('D2', self.parameters)
        D4 = Diagram('D4', self.parameters)

        # For the call to gamma_chains, it is necessary to cast LorentzVectors into lists
        if self.identifier == 'UV1':
            return (
                D1(PS_point, loop_momenta, invariants)
                *(prop(ps[1]+l, 0.)*prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.))**(-1)
                *(prop(l, params['m_IR2']))**3
            )

        elif self.identifier == 'UV2':
            return (
                D2(PS_point, loop_momenta, invariants)
                *(prop(ps[0]-l, 0.)*prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.))**(-1)
                *(prop(l, params['m_IR2']))**3
            )

        elif self.identifier == 'UV4':
            return (
                D4(PS_point, loop_momenta, invariants)
                *(prop(l, 0.)*prop(ps[0]+ps[3]-l, 0.))**(-1)
                *(prop(l-ps[0]/2-ps[3]/2, params['m_IR2']))**2
            )

        elif self.identifier == 'Soft':
            return (
                D0(PS_point, loop_momenta, invariants)
                *2*complex(0.,1.)*params['alpha_s']*4*np.pi*params['C_F']
                *invs['s12']
                *prop(l, 0.)*prop(ps[1]+l, 0.)*prop(ps[0]-l,0.)
            )

        elif self.identifier == 'Coll1':
            return (
                D0(PS_point, loop_momenta, invariants)
                *(-2)*complex(0.,1.)*params['alpha_s']*4*np.pi*params['C_F']
                *(prop(ps[0]-l, 0.)*prop(l, 0.)-(prop(l, params['m_IR2']))**2)
            )

        elif self.identifier == 'Coll2':
            return (
                D0(PS_point, loop_momenta, invariants)
                *(-2)*complex(0.,1.)*params['alpha_s']*4*np.pi*params['C_F']
                *(prop(ps[1]+l, 0.)*prop(l, 0.)-(prop(l, params['m_IR2']))**2) 
            )
        else:
            raise NotImplementedError

class Amplitude(integrand.Amplitude):

    def __init__ (self, identifier, n_loops=1, **opts):
        super(Amplitude, self).__init__(identifier, n_loops=n_loops, **opts) 

        if self.identifier != 'qq_aa':
            raise BaseException("AmplitudeQQAA only implements type identifier 'qq_aa', not '%s'."%self.identifier)
        if self.n_loops != 1:
            raise BaseException("AmplitudeQQAA only available at one-loop for now")

        # Now create contributing diagrams
        self.diagrams = []
        if self.n_loops==1:
            self.diagrams.append(Diagram('D1', self.parameters))
            self.diagrams.append(Diagram('D2', self.parameters))
            self.diagrams.append(Diagram('D3', self.parameters))
            self.diagrams.append(Diagram('D4', self.parameters))
        elif self.n_loops==2:
            # Example of where you'd implement the two-loop amplitude
            raise NotImplementedError
        else:
            raise NotImplementedError

        # And now the counterterms
        self.counterterms = []
        if self.n_loops==1:
            self.counterterms.append(Counterterm('UV1', self.parameters))
            self.counterterms.append(Counterterm('UV2', self.parameters))
            self.counterterms.append(Counterterm('UV4', self.parameters))
            self.counterterms.append(Counterterm('Soft', self.parameters))
            self.counterterms.append(Counterterm('Coll1', self.parameters))
            self.counterterms.append(Counterterm('Coll2', self.parameters))
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
        
        invariants['internal_masses'] = 0.
        
        for i, pi in enumerate(PS_point):
            invariants['sk%d'%(i+1)] = (loop_momenta[0]+pi).square()

        return invariants

    def __call__(self, PS_point, loop_momenta):
        
        res = 0.

        invariants = self.compute_invariants(PS_point, loop_momenta)

        for diag in self.diagrams:
            res += diag(PS_point, loop_momenta, invariants)
           # print (diag.identifier, "%.1f"%np.log(np.absolute(diag(PS_point, loop_momenta, invariants))))


        for CT in self.counterterms:
            res -= CT(PS_point, loop_momenta, invariants)
           # print (CT.identifier, "%.1f"%np.log(np.absolute(CT(PS_point, loop_momenta, invariants))))

        return res

