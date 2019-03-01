import vectors

class Parametrisation(object):

    def __init__(self, n_loops=1, **opts):
        self.n_loops = n_loops

    def __call__(self, random_variables, *args, **opts):
        raise BaseException("This function should be implemented in a daughter class.")

class CartesianParametrisation(Parametrisation):

    def __init__(self, *args, **opts):
        super(CartesianParametrisation, self).__init__(*args, **opts)

    def __call__(self, random_variables):
        
        loop_momenta = []
        jacobian     = 1.
        # Map unit hypercube to infinite hypercube
        for i, _ in enumerate(range(self.n_loops)):
            # Energy set to zero for now as it will be set by LTD cut conditions
            loop_momenta.append(vectors.LorentzVector([0.,]+
                    [(-1./x+1./(1.-x)) for x in random_variables[i*3:(i+1)*3]])) 
            jacobian /= 1./x**2 + 1./(1.-x)**2

        return loop_momenta, jacobian
