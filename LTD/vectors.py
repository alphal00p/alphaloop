import numpy as np
import math

class InvalidOperation(Exception):
    pass

class Vector(np.ndarray):

    def __new__(cls, *args, **opts):

        if args and isinstance(args[0], Vector):
            foo = args[0].get_copy()
        else:
            foo = np.asanyarray(*args, **opts).view(cls)
        return foo

    def eps(self):

        try: return np.finfo(self.dtype).eps
        except: return 0

    def huge(self):

        if np.issubdtype(self.dtype, np.inexact):
            return np.finfo(self.dtype).max
        elif np.issubdtype(self.dtype, np.integer):
            return np.iinfo(self.dtype).max
        else:
            raise ValueError

    def __eq__(self, other):

        eps = max(self.eps(), other.eps())
        # numpy.allclose uses abs(self-other) which would compute the norm
        # for LorentzVector, thus view self and other as numpy.ndarray's
        return np.allclose(
            self.view(type=np.ndarray), other.view(type=np.ndarray),
            math.sqrt(eps), 0.
        )

    def __ne__(self, other):

        return not self.__eq__(other)

    def __hash__(self):

        return tuple(x for x in self).__hash__()

    def get_copy(self):

        # The vector instantiated by get_copy() should be modified
        # without changing the previous instance, irrespectively of the
        # (presumably few) layers that compose entries of the vector
        # return copy.deepcopy(self)
        return copy.copy(self)

    def square(self):

        return self.dot(self)

    def __abs__(self):

        return math.sqrt(self.square())

    def normalize(self):

        self.__idiv__(abs(self))
        return

    def project_onto(self, v):

        return (self.dot(v) / v.square()) * v

    def component_orthogonal_to(self, v):

        return self - self.project_onto(v)

    # Specific to 3D vectors
    def cross(self, v):

        assert len(self) == 3
        assert len(v) == 3
        return Vector([
            self[1] * v[2] - self[2] * v[1],
            self[2] * v[0] - self[0] * v[2],
            self[0] * v[1] - self[1] * v[0]
         ])


class LorentzVector(Vector):

    def __new__(cls, *args, **opts):

        if len(args) == 0:
            return super(LorentzVector, cls).__new__(cls, [0., 0., 0., 0.], **opts)
        return super(LorentzVector, cls).__new__(cls, *args, **opts)

    def space(self):
        """Return the spatial part of this LorentzVector."""

        return self[1:].view(type=Vector)

    def dot(self, v):
        """C ompute the Lorentz scalar product."""
        ## The implementation below allows for a check but it should be done upstream and
        ## significantly slows down the code here.
        # pos = self[0]*v[0]
        # neg = self.space().dot(v.space())
        # if pos+neg != 0 and abs(2*(pos-neg)/(pos+neg)) < 100.*self.eps(): return 0
        # return pos - neg
        return self[0]*v[0] - self[1]*v[1] - self[2]*v[2] - self[3]*v[3]

    def square_almost_zero(self):
        """Check if the square of this LorentzVector is zero within numerical accuracy."""

        return (self.square() / self.view(Vector).square()) ** 2 < self.eps()

    def rho2(self):
        """Compute the radius squared."""

        return self.space().square()

    def rho(self):
        """Compute the radius."""

        return abs(self.space())

    def set_square(self, square, negative=False):
        """Change the time component of this LorentzVector
        in such a way that self.square() = square.
        If negative is True, set the time component to be negative,
        else assume it is positive.
        """

        # Note: square = self[0]**2 - self.rho2(),
        # so if (self.rho2() + square) is negative, self[0] is imaginary.
        # Letting math.sqrt fail if data is not complex on purpose in this case.
        self[0] = math.sqrt(self.rho2() + square)
        if negative: self[0] *= -1
        return self

    def rotoboost(self, p, q):
        """Apply the Lorentz transformation that sends p in q to this vector."""

        # NOTE: when applying the same Lorentz transformation to many vectors,
        #       this function goes many times through the same checks.

        # Compute squares
        p2 = p.square()
        q2 = q.square()
        # Numerical tolerances
        p_eps = math.sqrt(p.eps())
        q_eps = math.sqrt(q.eps())
        # Check if both Lorentz squares are small compared to the euclidean squares,
        # in which case the alternative formula should be used
        if p.square_almost_zero() and q.square_almost_zero():
            # Use alternative formula
            if p == self:
                for i in range(len(self)):
                    self[i] = q[i]
            else:
                print("Error in vectors.rotoboost: missing formula")
                print("Boosting %s (%.9e)" % (str(self), self.square()))
                print("p = %s (%.9e)" % (str(p), p2))
                print("q = %s (%.9e)" % (str(q), q2))
                print("Eq. (4.14) of arXiv:0706.0017v2, p. 26 not implemented")
                raise NotImplementedError
            return self
        else:
            # Check that the two invariants are close,
            # else the transformation is invalid
            if abs(p2-q2)/(abs(p2)+abs(q2)) > (p_eps+q_eps):
                print("Error in vectors.rotoboost: nonzero, unequal squares")
                print("p = %s (%.9e)" % (str(p), p2))
                print("q = %s (%.9e)" % (str(q), q2))
                raise InvalidOperation
            # Compute scalar products
            pq = p + q
            pq2 = pq.square()
            p_s = self.dot(p)
            pq_s = self.dot(pq)
            # Assemble vector
            self.__iadd__(2 * ((p_s/q2) * q - (pq_s/pq2) * pq))
            return self

    def pt(self, axis=3):
        """Compute transverse momentum."""

        return math.sqrt(
            sum(self[i]**2 for i in range(1, len(self)) if i != axis) )

    def pseudoRap(self):
        """Compute pseudorapidity."""

        pt = self.pt()
        if pt < self.eps() and abs(self[3]) < self.eps():
            return self.huge()*(self[3]/abs(self[3]))
        th = math.atan2(pt, self[3])
        return -math.log(math.tan(th/2.))

    def rap(self):
        """Compute rapidity in the lab frame. (needs checking)"""

        if self.pt() < self.eps() and abs(self[3]) < self.eps():
            return self.huge()*(self[3]/abs(self[3]))

        return .5*math.log((self[0]+self[3])/(self[0]-self[3]))

    def getdelphi(self, p2):
        """Compute the phi-angle separation with p2."""

        pt1 = self.pt()
        pt2 = p2.pt()
        if pt1 == 0. or pt2 == 0.:
            return self.huge()
        tmp = self[1]*p2[1] + self[2]*p2[2]
        tmp /= (pt1*pt2)
        if abs(tmp) > (1.0+math.sqrt(self.eps())):
            print("Cosine larger than 1. in phase-space cuts.")
            raise ValueError
        if abs(tmp) > 1.0:
            return math.acos(tmp/abs(tmp))
        return math.acos(tmp)

    def deltaR(self, p2):
        """Compute the deltaR separation with momentum p2."""

        delta_eta = self.pseudoRap() - p2.pseudoRap()
        delta_phi = self.getdelphi(p2)
        return math.sqrt(delta_eta**2 + delta_phi**2)

    def boostVector(self):

        if self == LorentzVector():
            return Vector([0.] * 3)
        if self[0] <= 0. or self.square() < 0.:
            print("Attempting to compute a boost vector from")
            print("%s (%.9e)" % (str(self), self.square()))
            raise InvalidOperation
        return self.space()/self[0]

    def cosTheta(self):

        ptot = self.rho()
        assert (ptot > 0.)
        return self[3] / ptot

    def phi(self):

        return math.atan2(self[2], self[1])
    
    def boost(self, boost_vector, gamma=-1.):
        """Transport self into the rest frame of the boost_vector in argument.
        This means that the following command, for any vector p=(E, px, py, pz)
            p.boost(-p.boostVector())
        transforms p to (M,0,0,0).
        """

        b2 = boost_vector.square()
        if gamma < 0.:
            gamma = 1.0 / math.sqrt(1.0 - b2)

        bp = self.space().dot(boost_vector)
        gamma2 = (gamma-1.0) / b2 if b2 > 0 else 0.
        factor = gamma2*bp + gamma*self[0]
        self_space = self.space()
        self_space += factor*boost_vector
        self[0] = gamma*(self[0] + bp)

