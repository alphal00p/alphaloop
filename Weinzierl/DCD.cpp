#include "DCD.h"

/*==========================================================
  =                    Real 4-vector  
  ==========================================================*/

//Constructor
DIdeform::R4vector::R4vector()
{
  std::vector<my_real>(4, 0.).swap(v);
}

DIdeform::R4vector::R4vector(const std::vector<my_real> &w)
{
  if (w.size() != 4)
  {
    printf("R4vector needs to be initialized with a 4d vector\n");
    exit(1);
  }
  v = w;
}

//Vector functions
my_real &DIdeform::R4vector::operator[](const int mu)
{
  return this->v.at(mu);
}

const my_real &DIdeform::R4vector::operator()(const int mu) const
{
  return this->v.at(mu);
}

DIdeform::R4vector DIdeform::R4vector::dual() const
{
  std::vector<my_real> dual_v(4);
  dual_v[0] = (*this)(0);
  for (int mu = 1; mu < 4; mu++)
    dual_v[mu] = -(*this)(mu);
  return dual_v;
}

//Operators
DIdeform::R4vector DIdeform::R4vector::operator+(const DIdeform::R4vector &rhs)
{
  std::vector<my_real> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = (*this)(mu) + rhs(mu);
  return result;
}

DIdeform::R4vector DIdeform::R4vector::operator-()
{
  std::vector<my_real> negative(4, 0.);
  for (int mu = 0; mu < 4; mu++)
    negative[mu] = -(*this)(mu);
  return negative;
}

DIdeform::R4vector DIdeform::R4vector::operator-(const DIdeform::R4vector &rhs)
{
  std::vector<my_real> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = (*this)(mu)-rhs(mu);
  return result;
}

DIdeform::R4vector DIdeform::operator*(const my_real a, const DIdeform::R4vector &w)
{
  std::vector<my_real> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = a * w(mu);
  return result;
}

std::ostream &DIdeform::operator<<(std::ostream &os, DIdeform::R4vector rhs)
{
  os << "[";
  for (int mu = 0; mu < 4; mu++)
  {
    if (rhs(mu) > 0.0)
      os << "+";
    os << rhs(mu);
    if (mu < 3)
      os << ", ";
  }
  os << "]";
  return os;
}

//Scalar product with metric (+,-,-,-)
my_real DIdeform::R4vector::operator*(const DIdeform::R4vector &rhs)
{
  my_real result = (*this)(0) * rhs(0);
  for (int mu = 1; mu < 4; mu++)
    result -= (*this)(mu)*rhs(mu);
  return result;
}
/*--------------- End Real 4-vector --------------------*/

/*==========================================================
  =                  Complex 4-vector
  ==========================================================*/

//Constructor
DIdeform::C4vector::C4vector()
{
  std::vector<my_comp>(4, 0.).swap(v);
}

DIdeform::C4vector::C4vector(const std::vector<my_comp> &w)
{
  if (w.size() != 4)
  {
    printf("R4vector needs to be initialized with a 4d vector\n");
    exit(1);
  }
  v = w;
}

//Vector functions
my_comp &DIdeform::C4vector::operator[](const int mu)
{
  return this->v.at(mu);
}

const my_comp &DIdeform::C4vector::operator()(const int mu) const
{
  return this->v.at(mu);
}

DIdeform::C4vector DIdeform::C4vector::dual() const
{
  std::vector<my_comp> dual_v(4);
  dual_v[0] = (*this)(0);
  for (int mu = 1; mu < 4; mu++)
    dual_v[mu] = -(*this)(mu);
  return dual_v;
}

//Operators
DIdeform::C4vector DIdeform::C4vector::operator+(const DIdeform::C4vector &rhs)
{
  std::vector<my_comp> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = (*this)(mu) + rhs(mu);
  return result;
}

DIdeform::C4vector DIdeform::C4vector::operator-()
{
  std::vector<my_comp> negative(4, 0.);
  my_comp zero(0., 0.);
  for (int mu = 0; mu < 4; mu++)
    negative[mu] = zero - (*this)(mu);
  return negative;
}

DIdeform::C4vector DIdeform::C4vector::operator-(const DIdeform::C4vector &rhs)
{
  std::vector<my_comp> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = (*this)(mu)-rhs(mu);
  return result;
}

DIdeform::C4vector DIdeform::operator*(const my_comp a, const DIdeform::C4vector &w)
{
  std::vector<my_comp> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = a * w(mu);
  return result;
}

std::ostream &DIdeform::operator<<(std::ostream &os, DIdeform::C4vector rhs)
{
  os << "[";
  for (int mu = 0; mu < 4; mu++)
  {
    os << rhs(mu);
    if (mu < 3)
      os << ", ";
  }
  os << "]";
  return os;
}

//Scalar product with metric (+,-,-,-)
my_comp DIdeform::C4vector::operator*(const DIdeform::C4vector &rhs)
{
  my_comp result = (*this)(0) * rhs(0);
  for (int mu = 1; mu < 4; mu++)
    result -= (*this)(mu)*rhs(mu);
  return result;
}

/*--------------- End Complex 4-vector --------------------*/

/*==========================================================
  =                   Hypercube mapping
  ==========================================================*/
my_real alpha = 1.0e+02;

DIdeform::R4vector DIdeform::hypcub_mapping(std::vector<my_real> x)
{
  DIdeform::R4vector momentum;
  for (int mu = 0; mu < 4; mu++)
      momentum[mu] = alpha * log(x[mu] / (1 - x[mu]));
      //momentum[mu] = alpha * (1. / (1. - x[mu]) - 1. / x[mu]);

  return momentum;
}

//my_real DIdeform::hypcub_jacobian(DIdeform::R4vector& x){
my_real DIdeform::hypcub_jacobian(std::vector<my_real> x)
{
  my_real jacobian = 1.;
  for (int mu = 0; mu < 4; mu++)
    jacobian = jacobian * (alpha / (x[mu] * (1 - x[mu])));
    //jacobian = jacobian * alpha * (1. / pow(x[mu], 2) + 1. / pow(1. - x[mu], 2));
  return jacobian;
}

/*--------------- End Hypercube mapping --------------------*/

/***********************************************************************
   START: Class ContourDeform 

   Definition of all the function needed from the class ContourDefrom
   
   It is divided in:

    - Constructors
    - Public Functions

    - Hepler Functions
    - Helper Functions - Derivatives
    - Exterior Regions
    - Exterior Regions - Derivatives
    - Interior Regions 
    - Interior Regions - Derivatives
    - Deformation Vector
    - Scaling Parameter
    - Scaling Parameter - Derivatives



***********************************************************************/

/*==========================================================
  =                   Constructors  
  ==========================================================*/
DIdeform::ContourDeform::ContourDeform(std::vector<DIdeform::R4vector> &Qs)
{
  qi = Qs;
  pi = Qs;

  //Identify topology
  legs = qi.size();
  for (int i = 0; i < legs; i++)
  {
    pi[i] = qi[i] - qi[(i - 1 + legs) % legs];
    mi.push_back(0.0); // Set everything massless;
    if (pi[i](0) < 0.)
      A = i;
  }

  if (pi[0](0) > 0.)
  {
    std::printf("p0 has to be an incoming momenta! i.e p0(0) < 0.\n");
    exit(1);
  }

  //Compute P+ and P-
  set_PpPm(Qs);
  //Compute Misq
  set_global_var();
  //std::printf("mu_P: %.5e\n",mu_P);
  //Resize grad varibales
  std::vector<std::vector<my_real>>(4, std::vector<my_real>(4, 0.)).swap(gradk_int);
  std::vector<std::vector<my_real>>(4, std::vector<my_real>(4, 0.)).swap(gradk_ext);
  std::vector<std::vector<my_real>>(4, std::vector<my_real>(4, 0.)).swap(gradk0);
  std::vector<my_real>(4, 0.).swap(gradlambda);
}

bool DIdeform::ContourDeform::get_PpPm(short int plus_or_minus, std::vector<DIdeform::R4vector> &qs)
{
  int size = qs.size();
  if (size == 1)
    return true;

  //Remove vectors that lies in the forward/backward light-cone of another
  DIdeform::R4vector qij;
  std::vector<DIdeform::R4vector>::iterator pos;
  for (int i = 0; i < size; i++)
  {
    for (int j = i + 1; j < size; j++)
    {
      qij = qs[i] - qs[j];
      if (qij * qij >= 0.0)
      {
        if (plus_or_minus * qij(0) < 0.0)
        {
          qs.erase(qs.begin() + j);
          j--;
        }
        else
        {
          qs.erase(qs.begin() + i);
          j=i;
        }
        if (--size == 1)
          return true;
      }
    }
  }

  //Find the pair with the smallest space-like separation
  int pos_i = 0, pos_j = 0;
  my_real sp;
  for (int i = 0; i < size; i++)
  {
    for (int j = i + 1; j < size; j++)
    {
      qij = qs[i] - qs[j];
      if (i == 0 && j == 1)
      {
        sp = -qij * qij;
        pos_i = i;
        pos_j = j;
      }
      else if (sp > -qij * qij)
      {
        sp = -qij * qij;
        pos_i = i;
        pos_j = j;
      }
    }
  }
  //Update qs
  if (plus_or_minus == +1)
    qs[pos_i] = zp(qs[pos_i] + qs[pos_j], qs[pos_i] - qs[pos_j]);
  if (plus_or_minus == -1)
    qs[pos_i] = zm(qs[pos_i] + qs[pos_j], qs[pos_i] - qs[pos_j]);
  qs.erase(qs.begin() + pos_j);

  //Repeat
  get_PpPm(plus_or_minus, qs);
}
void DIdeform::ContourDeform::set_PpPm(std::vector<DIdeform::R4vector> &Qs)
{
  // //Note that this definition is valid only for massless momenta
  // Pp = zp(qi[A] + qi[0], qi[A] - qi[0]);
  // Pm = zm(qi[A - 1] + qi[legs - 1], qi[A - 1] - qi[legs - 1]);

  std::vector<DIdeform::R4vector> qs;
  qs = Qs;
  get_PpPm(+1, qs);
  Pp = qs[0];

  qs = Qs;
  get_PpPm(-1, qs);
  Pm = qs[0];
  
  //Check the validity of Pp and Pm
  my_real eps = std::numeric_limits<my_real>::epsilon();
  DIdeform::R4vector vv;
  for (int i = 0; i < 4; i++)
  {
    vv = Qs[i] - Pp;
    if (vv * vv + eps < 0.0 || vv(0) < 0.0)
    {
      std::printf("The build of P+ has failed. One of the Qs is not in its forward light-cone\n");
      std::printf("\tValue: %.10e, v(0): %.10e\n", vv * vv, vv(0));
      std::cout << "\t" << Pp << std::endl;
      exit(1);
    }
    vv = Qs[i] - Pm;
    if (vv * vv < 0.0 || vv(0) > 0.0)
    {
      std::printf("The build of P- has failed. One of the Qs is not in its forward light-cone\n");
      std::printf("\tValue: %.10e, v(0): %.10e\n", vv * vv, vv(0));
      std::cout << "\t" << Pm << std::endl;
      exit(1);
    }
  }
}

void DIdeform::ContourDeform::set_global_var()
{
  //Center of mass energy squared
  Ecmsq = (pi[0] + pi[A]) * (pi[0] + pi[A]);
  mu_P = (Pm - Pp) * (Pm - Pp);

  //Set M1sq and M2sq
  M1sq = pow(M1f, 2) * Ecmsq;
  M2sq = pow(M2f, 2) * std::max(mu_P, Ecmsq);
  M3sq = pow(M3f, 2) * std::max(mu_P, Ecmsq);

  //Set orthogonal vectors ei's
  std::vector<DIdeform::R4vector>(4).swap(ei);
  for (int mu = 0; mu < 4; mu++)
    ei[mu][mu] = E_soft * sqrt(Ecmsq);
}

/*==========================================================
  =                   Public Functions  
  ==========================================================*/

void DIdeform::ContourDeform::loop_momentum(DIdeform::R4vector &loopmomentum)
{
  l = loopmomentum;
  
  //Set l-qi
  std::vector<DIdeform::R4vector>(legs).swap(lqi);
  for (int i = 0; i < legs; i++)
    lqi[i] = l - qi[i];

  //Set K+, K- and K_center
  kp = l - Pp;
  km = l - Pm;
  k_center = 0.5 * (kp + km);
}

void DIdeform::ContourDeform::deform(DIdeform::C4vector &k, my_comp &jacobian)
{
  set_k0();
  set_lambda();

  DIdeform::C4vector ell;
  my_comp ii(0.0, 1.0);
  for (int mu = 0; mu < 4; mu++)
    ell[mu] = l(mu) + ii * lambda * k0(mu);

  std::vector<std::vector<my_comp>> grad(4, std::vector<my_comp>(4, 0.));
  for (int mu = 0; mu < 4; mu++)
  {
    grad[mu][mu] += 1.;
    for (int nu = 0; nu < 4; nu++)
      grad[mu][nu] += ii * (lambda * gradk0[mu][nu] + gradlambda[nu] * k0[mu]);
  }

  k = ell;
  jacobian = Determinant(grad);
}

/*==========================================================
  =                 Helper Functions 
  ==========================================================*/

my_real DIdeform::ContourDeform::hp(DIdeform::R4vector &k, my_real msq, my_real Msq)
{
  my_real Ek = k(0);                 //energy part of k
  my_real vksq = pow(Ek, 2) - k * k; //vector part of k squared
  my_real val = pow(+Ek - sqrt(vksq + msq), 2);

  return val / (val + Msq);
}

my_real DIdeform::ContourDeform::hm(DIdeform::R4vector &k, my_real msq, my_real Msq)
{
  my_real Ek = k(0);                 //energy part of k
  my_real vksq = pow(Ek, 2) - k * k; //vector part of k squared
  my_real val = pow(-Ek - sqrt(vksq + msq), 2);

  return val / (val + Msq);
}
my_real DIdeform::ContourDeform::h0(DIdeform::R4vector &k, my_real msq, my_real Msq)
{
  my_real Ek = k(0);                 //energy part of k
  my_real vksq = pow(Ek, 2) - k * k; //vector part of k squared
  my_real val = pow(abs(Ek) - sqrt(vksq + msq), 2);

  return val / (val + Msq);
}

my_real DIdeform::ContourDeform::ht(my_real t, my_real Msq)
{
  return t < 0 ? 0 : t / (t + Msq);
}

/*==========================================================
  =        Helper Functions - Derivatives 
  ==========================================================*/

DIdeform::R4vector DIdeform::ContourDeform::gradloghp(DIdeform::R4vector &k, my_real msq, my_real Msq)
{
  my_real Ek = k(0);                 //energy part of k
  my_real vksq = pow(Ek, 2) - k * k; //vector part of k squared
  my_real w = +Ek - sqrt(vksq + msq);

  DIdeform::R4vector gradw = -(1. / sqrt(vksq + msq)) * k;
  gradw[0] = +1.;

  DIdeform::R4vector grad = 2.0 * Msq / w / (pow(w, 2) + Msq) * gradw;
  return grad;
}

DIdeform::R4vector DIdeform::ContourDeform::gradloghm(DIdeform::R4vector &k, my_real msq, my_real Msq)
{
  my_real Ek = k(0);                 //energy part of k
  my_real vksq = pow(Ek, 2) - k * k; //vector part of k squared
  my_real w = -Ek - sqrt(vksq + msq);

  DIdeform::R4vector gradw = -(1. / sqrt(vksq + msq)) * k;
  gradw[0] = -1.;

  DIdeform::R4vector grad = 2.0 * Msq / w / (pow(w, 2) + Msq) * gradw;
  return grad;
}

DIdeform::R4vector DIdeform::ContourDeform::gradlogh0(DIdeform::R4vector &k, my_real msq, my_real Msq)
{
  my_real Ek = k(0);                 //energy part of k
  my_real vksq = pow(Ek, 2) - k * k; //vector part of k squared
  my_real w = abs(Ek) - sqrt(vksq + msq);

  DIdeform::R4vector gradw = -(1. / sqrt(vksq + msq)) * k;
  gradw[0] = Ek > 0 ? +1. : -1.;

  DIdeform::R4vector grad = 2.0 * Msq / w / (pow(w, 2) + Msq) * gradw;
  return grad;
}

my_real DIdeform::ContourDeform::dloght(my_real t, my_real Msq)
{
  return t < 0 ? 0 : Msq / t / (t + Msq);
}

/*==========================================================
  =                 Exterior Region 
  ==========================================================*/

DIdeform::R4vector DIdeform::ContourDeform::zp(DIdeform::R4vector x, DIdeform::R4vector y)
{
  my_real Ey = y(0);                        //Energy part of y
  my_real vymod = sqrt(pow(Ey, 2) - y * y); //Modulus of the vector part of y

  if (vymod == 0.0)
    return 0.5 * x;
  else
    return 0.5 * (x + (1. / vymod) * ((y * y) * g0mu - Ey * y));
}

DIdeform::R4vector DIdeform::ContourDeform::zm(DIdeform::R4vector x, DIdeform::R4vector y)
{
  my_real Ey = y(0);                        //Energy part of y
  my_real vymod = sqrt(pow(Ey, 2) - y * y); //Modulus of the vector part of y

  if (vymod == 0.0)
    return 0.5 * x;
  else
    return 0.5 * (x - (1. / vymod) * ((y * y) * g0mu - Ey * y));
}

void DIdeform::ContourDeform::set_k_ext()
{
  //Set c+ and c-
  cp = 1.;
  cm = 1.;
  for (int i = 0; i < legs; i++)
  {
    cp = cp * hm(lqi[i], pow(mi[i], 2), M3sq);
    cm = cm * hp(lqi[i], pow(mi[i], 2), M3sq);
  }
  //k_ext
  //k_ext=(cp*kp+cm*km).dual();
  k_ext = (cp * kp + cm * km).dual();
}

/*==========================================================
  =          Exterior Region - Derivatives
  ==========================================================*/

void DIdeform::ContourDeform::set_gradk_ext()
{
  //Set gradc+ and gradc-
  DIdeform::R4vector gradcp, gradcm;

  for (int i = 0; i < legs; i++)
  {
    gradcp = gradcp + gradloghm(lqi[i], pow(mi[i], 2), M3sq);
    gradcm = gradcm + gradloghp(lqi[i], pow(mi[i], 2), M3sq);
  }
  //Using c+ and c- set from set_k_ext
  gradcp = cp * gradcp;
  gradcm = cm * gradcm;

  //gradk_ext
  for (int mu = 0; mu < 4; mu++)
  {
    gradk_ext[mu][mu] = gmunu[mu] * (cp + cm);
    for (int nu = 0; nu < 4; nu++)
    {
      if (mu == nu)
      {
        gradk_ext[mu][nu] += gmunu[mu] * gradcp(nu) * kp(mu);
        gradk_ext[mu][nu] += gmunu[mu] * gradcm(nu) * km(mu);
      }
      else
      {
        gradk_ext[mu][nu] = gmunu[mu] * gradcp(nu) * kp(mu);
        gradk_ext[mu][nu] += gmunu[mu] * gradcm(nu) * km(mu);
      }
    }
  }
}

/*==========================================================
  =                 Interior Region 
  ==========================================================*/

//Interior region
my_real DIdeform::ContourDeform::g(DIdeform::R4vector &k, my_real Msq)
{
  DIdeform::R4vector dualk = k.dual();
  return gamma1 * Msq / (k * dualk + Msq);
}

my_real DIdeform::ContourDeform::dij(int i, int j)
{
  DIdeform::R4vector qij = qi[i] - qi[j];
  my_real eps = std::numeric_limits<my_real>::epsilon();
  my_real msq = pow(mi[j], 2);

  if (msq <= eps && i == j)
  { //Case 1
    return 1;
  }
  else if (msq <= eps && //Case 2
           qij(0) <= 0.0 &&
           abs(qij * qij) <= eps)
  {
    return hp(lqi[j], msq, M1sq);
  }
  else if (msq <= eps && //Case 3
           qij(0) > 0.0 &&
           abs(qij * qij) <= eps)
  {
    return hm(lqi[j], msq, M1sq);
  }
  else //Otherwise
    return std::max(h0(lqi[j], msq, M1sq), ht(-2 * (lqi[j] * lqi[i]), M1sq));
}

my_real DIdeform::ContourDeform::dijk(int i, int j, int k)
{
  DIdeform::R4vector qij = qi[i] - qi[j];
  my_real eps = std::numeric_limits<my_real>::epsilon();
  my_real msq = pow(mi[k], 2);

  DIdeform::R4vector lvij = l - vij(i, j);
  my_real zij = qij * qij - pow(mi[i] + mi[j], 2);

  return ht(zij, M1sq) * std::max(h0(lqi[k], msq, M1sq), ht(-2 * (lqi[k] * lvij), M1sq));
}

my_real DIdeform::ContourDeform::ci(int i)
{
  my_real c_i = g(k_center, M2sq);
  //std::cout<< c_i << " "<< M2sq << std::endl;
  for (int j = 0; j < legs; j++)
    c_i = c_i * dij(i, j);

  return c_i;
}

my_real DIdeform::ContourDeform::cij(int i, int j)
{
  my_real c_ij = g(k_center, M2sq);

  for (int k = 0; k < legs; k++)
    c_ij = c_ij * dijk(i, j, k);

  return c_ij;
}

DIdeform::R4vector DIdeform::ContourDeform::vij(int i, int j)
{
  DIdeform::R4vector qij = qi[i] - qi[j];

  return 0.5 * (qi[i] + qi[j] - ((mi[i] - mi[j]) / sqrt(qij * qij)) * qij);
}

void DIdeform::ContourDeform::set_k_int()
{
  /* k_int = -c_i lqi - c_ij lvij + k_soft
     lvij = l - vij
     k_soft to be implemented
  */
  DIdeform::R4vector k1; // Zero vector
  //Massless case
  for (int i = 0; i < legs; i++)
  {
    k1 = k1 + ci(i) * lqi[i];
    //std::cout << k1 << std::endl;
  }

  k_int = -k1;

  //Massive case
  //.... k2 and ksoft
}

/*==========================================================
  =          Interior Region - Derivatives
  ==========================================================*/

DIdeform::R4vector DIdeform::ContourDeform::gradlogg(DIdeform::R4vector &k, my_real Msq)
{
  DIdeform::R4vector dualk = k.dual();
  DIdeform::R4vector grad = (-2.0 / (k * dualk + Msq)) * k;
  return grad;
}

DIdeform::R4vector DIdeform::ContourDeform::gradlogdij(int i, int j)
{
  DIdeform::R4vector qij = qi[i] - qi[j];
  my_real eps = std::numeric_limits<my_real>::epsilon();
  my_real msq = pow(mi[j], 2);

  DIdeform::R4vector grad; // At the moment is a zero vector

  if (msq <= eps && i == j) //Case 1
  {
  }
  else if (msq <= eps && //Case 2
           qij(0) <= 0.0 &&
           abs(qij * qij) <= eps)
  {
    grad = gradloghp(lqi[j], msq, M1sq);
  }
  else if (msq <= eps && //Case 3
           qij(0) > 0.0 &&
           abs(qij * qij) <= eps)
  {
    grad = gradloghm(lqi[j], msq, M1sq);
  }
  else
  { //Otherwise
    if (h0(lqi[j], msq, M1sq) > ht(-2 * (lqi[j] * lqi[i]), M1sq))
      grad = gradlogh0(lqi[j], msq, M1sq);
    else
      grad = dloght(-2 * (lqi[j] * lqi[i]), M1sq) *
             (-2) * (lqi[j].dual() + lqi[i].dual());
  }

  return grad;
}

DIdeform::R4vector DIdeform::ContourDeform::gradlogdijk(int i, int j, int k)
{
  DIdeform::R4vector qij = qi[i] - qi[j];
  my_real eps = std::numeric_limits<my_real>::epsilon();
  my_real msq = pow(mi[k], 2);

  DIdeform::R4vector lvij = l - vij(i, j);
  my_real zij = qij * qij - pow(mi[i] + mi[j], 2);

  DIdeform::R4vector grad; // At the moment is a zero vector
  if (h0(lqi[k], msq, M1sq) > ht(-2 * (lqi[k] * lvij), M1sq))
    grad = gradlogh0(lqi[k], msq, M1sq);
  else
    grad = dloght(-2 * (lqi[k] * lvij), M1sq) * (-2) * (lqi[k].dual() + lvij.dual());

  return grad;
}

DIdeform::R4vector DIdeform::ContourDeform::gradlogci(int i)
{
  DIdeform::R4vector grad = gradlogg(k_center, M2sq);

  for (int j = 0; j < legs; j++)
    grad = grad + gradlogdij(i, j);

  return grad;
}

DIdeform::R4vector DIdeform::ContourDeform::gradlogcij(int i, int j)
{
  DIdeform::R4vector grad = gradlogg(k_center, M2sq);

  for (int k = 0; k < legs; k++)
    grad = grad + gradlogdijk(i, j, k);

  return grad;
}

void DIdeform::ContourDeform::set_gradk_int()
{
  /* MASSLESS CASE */
  DIdeform::R4vector gradci;
  my_real c_i;
  //Reset gradk_int
  std::vector<std::vector<my_real>>(4, std::vector<my_real>(4, 0.)).swap(gradk_int);

  for (int i = 0; i < legs; i++)
  {
    c_i = ci(i);
    gradci = c_i * gradlogci(i);

    for (int mu = 0; mu < 4; mu++)
    {
      gradk_int[mu][mu] += -c_i;
      for (int nu = 0; nu < 4; nu++)
      {
        gradk_int[mu][nu] += -gradci[nu] * lqi[i](mu);
      }
    }
  }
}

/*==========================================================
  =                   Deformation Vector
  ==========================================================*/

void DIdeform::ContourDeform::set_k0()
{
  set_k_ext();
  set_k_int();

  set_gradk_ext();
  set_gradk_int();

  k0 = k_ext + k_int;

  for (int mu = 0; mu < 4; mu++)
    for (int nu = 0; nu < 4; nu++)
      gradk0[mu][nu] = gradk_ext[mu][nu] + gradk_int[mu][nu];
}

/*==========================================================
  =                 Scaling Parameter 
  ==========================================================*/

my_real DIdeform::ContourDeform::Xi(int i)
{
  return pow((k0 * lqi[i]) / (k0 * k0), 2);
}

my_real DIdeform::ContourDeform::Yi(int i)
{
  return (lqi[i] * lqi[i] - pow(mi[i], 2)) / (k0 * k0);
}

my_real DIdeform::ContourDeform::lambda_i(int i)
{
  my_real lambdai = lambda_max;
  my_real xi = Xi(i);
  my_real yi = Yi(i);

  if (2 * xi < yi)
    lambdai = sqrt(yi / 4.);
  else if (yi < 0)
    lambdai = sqrt(xi - yi / 2.);
  else if (yi < 2 * xi)
    lambdai = sqrt(xi - yi / 4.);

  return lambdai;
}

my_real DIdeform::ContourDeform::lambda_coll()
{
  /*MASSLESS CASE*/
  my_real C = 0.;
  for (int i = 0; i < legs; i++)
    C += ci(i);

  return 1 / (4. * C);
}

void DIdeform::ContourDeform::set_lambda()
{
  //Start by assuming lambda_max is the smallest
  short int which = 5;
  lambda = lambda_max;
  gradlambda = {0., 0., 0., 0.};

  my_real lambda_tmp;
  DIdeform::R4vector grad;

  //Find smallest lambda_i
  for (int i = 0; i < legs; i++)
  {
    lambda_tmp = lambda_i(i);
    if (lambda_tmp < lambda)
    {
      which = i;
      lambda = lambda_tmp;
      grad = gradlambda_i(i);
      for (int mu = 0; mu < 4; mu++)
      {
        gradlambda[mu] = grad(mu);
      }
    }
  }

  //Compare with lambda_coll
  lambda_tmp = lambda_coll();
  if (lambda_tmp < lambda)
  {
    which = 4;
    lambda = lambda_tmp;
    grad = gradlambda_coll();
    for (int mu = 0; mu < 4; mu++)
    {
      gradlambda[mu] = grad(mu);
    }
  }
  //std::cout << "which: " << which << std::endl;
}

/*==========================================================
  =          Scaling Parameter - Derivatives
  ==========================================================*/

DIdeform::R4vector DIdeform::ContourDeform::gradXi(int i)
{
  DIdeform::R4vector grad;
  my_real k0sq = k0 * k0;
  my_real k0lqi = k0 * lqi[i];

  for (int mu = 0; mu < 4; mu++)
  {
    grad[mu] += gmunu[mu] * (2.0 / (k0lqi)) * k0(mu);
    for (int nu = 0; nu < 4; nu++)
    {
      grad[mu] += (gmunu[nu] * (-4.0 / (k0sq)*k0(nu) + 2.0 / (k0lqi)*lqi[i](nu)) * gradk0[nu][mu]);
    }
  }
  grad = Xi(i) * grad;

  return grad;
}

DIdeform::R4vector DIdeform::ContourDeform::gradYi(int i)
{
  my_real w = (lqi[i] * lqi[i] - pow(mi[i], 2));
  my_real k0sq = k0 * k0;

  DIdeform::R4vector grad;
  for (int mu = 0; mu < 4; mu++)
  {
    grad[mu] += gmunu[mu] * (2.0 / w) * lqi[i](mu);
    for (int nu = 0; nu < 4; nu++)
    {
      grad[mu] += gmunu[nu] * (-2.0 / k0sq) * k0(nu) * gradk0[nu][mu];
    }
  }
  grad = Yi(i) * grad;

  return grad;
}

DIdeform::R4vector DIdeform::ContourDeform::gradlambda_i(int i)
{
  my_real xi = Xi(i);
  my_real yi = Yi(i);

  DIdeform::R4vector grad;

  if (2 * xi < yi)
    grad = 0.25 * gradYi(i);
  else if (yi < 0)
    grad = gradXi(i) - 0.5 * gradYi(i);
  else if (yi < 2 * xi)
    grad = gradXi(i) - 0.25 * gradYi(i);

  grad = 0.5 / lambda_i(i) * grad;

  return grad;
}

DIdeform::R4vector DIdeform::ContourDeform::gradlambda_coll()
{
  /*MASSLESS CASE*/
  DIdeform::R4vector gradC;
  for (int i = 0; i < legs; i++)
    gradC = gradC + ci(i) * gradlogci(i);

  return -4.0 * pow(lambda_coll(), 2) * gradC;
}
/*--------------- END - Class ContourDefrom   --------------------*/

// Function Determinant(M,d).
// This function calculates the determinant of any input matrix using LU decomposition.
// Code by W. Gong and D.E. Soper
my_comp DIdeform::Determinant(const std::vector<std::vector<my_comp>> &bb)
{
  // Define matrix related variables.
  int dimension = bb.at(0).size();
  int maxsize = dimension - 1;
  static const my_comp one(1.0, 0.0);
  my_comp determinant;
  int indx[dimension], d;
  // Initialize the determinant.
  determinant = one;
  // Inintialize the matrix to be decomposed with the transferred matrix b.
  std::vector<std::vector<my_comp>> aa = bb;

  // Define parameters used in decomposition.
  int i, imax, j, k, flag = 1;
  my_comp dumc, sum;
  my_real aamax, dumr;
  my_real vv[maxsize];
  // Initialize the parity parameter.
  d = 1;
  // Get the implicit scaling information.
  for (i = 0; i < dimension; i++)
  {
    aamax = 0;
    for (j = 0; j < dimension; j++)
    {
      if (abs(aa[i][j]) > aamax)
        aamax = abs(aa[i][j]);
    }
    // Set a flag to check if the determinant is zero.
    if (aamax == 0)
      flag = 0;
    // Save the scaling.
    vv[i] = 1.0 / aamax;
  }
  if (flag == 1)
  {
    for (j = 0; j < dimension; j++)
    {
      for (i = 0; i < j; i++)
      {
        sum = aa[i][j];
        for (k = 0; k < i; k++)
          sum = sum - aa[i][k] * aa[k][j];
        aa[i][j] = sum;
      }
      //Initialize for the search for largest pivot element.
      aamax = 0;
      for (i = j; i < dimension; i++)
      {
        sum = aa[i][j];
        for (k = 0; k < j; k++)
          sum = sum - aa[i][k] * aa[k][j];
        aa[i][j] = sum;
        // Figure of merit for the pivot.
        dumr = vv[i] * abs(sum);
        // Is it better than the best so far?
        if (dumr >= aamax)
        {
          imax = i;
          aamax = dumr;
        }
      } // End for(i = j; i < dimension; i++)
      // See if we need to interchange rows.
      if (j != imax)
      {
        for (k = 0; k < dimension; k++)
        {
          dumc = aa[imax][k];
          aa[imax][k] = aa[j][k];
          aa[j][k] = dumc;
        }
        // Change the parity of d.
        d = -d;
        // Interchange the scale factor.
        vv[imax] = vv[j];
      } // End if(j != imax)
      indx[j] = imax;
      if (j != dimension - 1)
      {
        dumc = my_real(1.0) / aa[j][j];
        for (i = j + 1; i < dimension; i++)
          aa[i][j] = aa[i][j] * dumc;
      }
    } // End for(j = 0; j < dimension; j++)
  }   // END if(flag == 1)
  // Calculate the determinant using the decomposed matrix.
  if (flag == 0)
    determinant = 0.0;
  else
  { // Multiply the diagonal elements.
    for (int diagonal = 0; diagonal < dimension; diagonal++)
    {
      determinant = determinant * aa[diagonal][diagonal];
    }
    determinant = (my_real)d * determinant;
  } // End if(flag == 0) ... else
  return determinant;
}
// End of function Determinant(M,d).

#ifdef _DCD_MAIN_

int main()
{
  // DIdeform::R4vector p1({0.5, 0.5, 0.0, 1.0});
  //  DIdeform::R4vector p2({0.5, -0.5, 0.0, -2.0});
  //  DIdeform::R4vector p3({0.5, -10.0, 0.5, -4.0});
  DIdeform::R4vector p1({-4.7213384875835902e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, -4.6142211817746778e+02});
  DIdeform::R4vector p2({-5.0290194983056193e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, 4.6142211817746778e+02});
  DIdeform::R4vector p3({4.5162178137689028e+02, -3.3466005003229913e+02, -1.3177690342268895e+00, -4.4307423883444024e+01});
  DIdeform::R4vector p4({5.2341401721203044e+02, 3.3466005003229913e+02, 1.3177690342268895e+00, 4.4307423883444024e+01});

  std::vector<DIdeform::R4vector> Qs(4);
  //Qs[0] = Qs[0];
  //Qs[1] = -p2;
  //Qs[2] = p1 - p3;
  //Qs[3] = p1;

  Qs[0] = p1;
  Qs[1] = Qs[0] + p2;
  Qs[2] = Qs[1] + p3;
  Qs[3] = 0.0 * Qs[3];

  DIdeform::R4vector shift(p1);
  for (int i = 0; i < 4; i++)
  {
    std::cout << "Q[" << i << "]" << Qs[i] << std::endl;
  }

  //DIdeform::R4vector l({0.5, 0.5, 0.5, 0.5});
  DIdeform::R4vector l({-96518.69521385277, -974059.78677852161, -486540.90951631521, -324033.96315526153});
  //DIdeform::R4vector l({0.371187,0.367224,-0.0345374,-0.0480211});
  //DIdeform::R4vector l({1652e+10, 4500e+10, 2315e+10, 4521e+10});
  DIdeform::ContourDeform contour(Qs);
  contour.loop_momentum(l);

  DIdeform::C4vector ell;
  my_comp jacobian;
  contour.deform(ell, jacobian);

  std::cout << "Pp:\t" << contour.Pp << std::endl;
  std::cout << "Pm:\t" << contour.Pm << std::endl;
  
  std::cout << "l:\t" << l << std::endl;
  std::cout << "ell:\t" << ell << std::endl;
  std::cout << "J:\t" << jacobian << std::endl;

  my_real eps = sqrt(std::numeric_limits<my_real>::epsilon());
  std::vector<std::vector<my_comp>> grad(4, std::vector<my_comp>(4, 0.));
  std::vector<DIdeform::R4vector> ev(4);
  for (int mu = 0; mu < 4; mu++)
    ev[mu][mu] = 1.0;

  DIdeform::R4vector dl;
  DIdeform::C4vector dell;

  //Numeric Jacobian
  my_real ep;
  for (int mu = 0; mu < 4; mu++)
  {
    ep = std::max(l[mu] * eps, eps);
    dl = l + ep * ev[mu];

    contour.loop_momentum(dl);
    contour.deform(dell, jacobian);

    dell = 1. / ep * (dell - ell);
    for (int nu = 0; nu < 4; nu++)
    {
      grad[mu][nu] = dell(nu);
    }
  }

std::cout << "NJ:\t" << DIdeform::Determinant(grad) << std::endl;

/*
  std::cout << "g(kc): " << contour.g(contour.k_center, contour.M2sq) << std::endl;
  std::cout << "cp: " << contour.cp << " " << contour.Pp << std::endl;
  std::cout << "cm: " << contour.cm << " " << contour.Pm << std::endl;
  for (int i = 0; i < 4; i++)
  {
    std::cout << "c" << i << ": " << contour.ci(i) << std::endl;
    std::cout << "x" << i << ": " << contour.Xi(i) << std::endl;
    std::cout << "y" << i << ": " << contour.Yi(i) << std::endl;
    std::cout << std::endl;
  }
  my_comp denominator = 1.;
  DIdeform::C4vector prop_mom;
  for (int i = 0; i < 4; i++)
  {
    for (int mu = 0; mu < 4; mu++)
      prop_mom[mu] = l(mu) - Qs[i](mu);
    std::cout << Qs[i] << std::endl;
    std::cout << prop_mom * prop_mom << std::endl;
  }
   
  //Numerator
  my_comp r = pow(ell(3), 2);
  std::cout << jacobian * r / denominator << std::endl;

  my_real x, dx;
  for (int i = 2; i < 3; i++)
  {
    contour.loop_momentum(l);
    contour.set_k0();
    int j=0;
    std::cout << "A[d"<<j << i << "]\t" << contour.dij(j,i) * contour.gradlogdij(j,i) << std::endl;

    std::cout << "N\t[";
    x = contour.dij(j,i);
    for (int mu = 0; mu < 4; mu++)
    {
      ep = l[mu] * eps;
      dl = l + ep * ev[mu];

      contour.loop_momentum(dl);
      contour.set_k0();
      dx = contour.dij(j,i);

      dx = 1. / ep * (dx - x);
      std::cout << dx << ", ";
    }
    std::cout << std::endl;
  }
  //k_int
  contour.loop_momentum(l);
  contour.set_k0();
  contour.set_lambda();
  DIdeform::R4vector v, dv;
  std::cout << "A:\t[k_int]\n";
  for (int nu = 0; nu < 4; nu++)
  {
    for (int mu = 0; mu < 4; mu++)
      std::cout << contour.gradk_int[mu][nu] << ", ";
    std::cout << std::endl;
  }
  std::cout << "\nN:\t[\n";
  v = contour.k_int;
  for (int mu = 0; mu < 4; mu++)
  {
    ep = l[mu] * eps;
    dl = l + ep * ev[mu];

    contour.loop_momentum(dl);
    contour.set_k0();
    dv = contour.k_int;

    dv = 1. / ep * (dv - v);
    for (int nu = 0; nu < 4; nu++)
    {
      std::cout << dv[nu] << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  //k_ext
  contour.loop_momentum(l);
  contour.set_k0();
  contour.set_lambda();
  std::cout << "A:\t[k_ext]\n";
  for (int nu = 0; nu < 4; nu++)
  {
    for (int mu = 0; mu < 4; mu++)
      std::cout << contour.gradk_ext[mu][nu] << ", ";
    std::cout << std::endl;
  }
  std::cout << "\nN:\t[\n";
  v = contour.k_ext;
  for (int mu = 0; mu < 4; mu++)
  {
    ep = l[mu] * eps;
    dl = l + ep * ev[mu];

    contour.loop_momentum(dl);
    contour.set_k0();
    dv = contour.k_ext;

    dv = 1. / ep * (dv - v);
    for (int nu = 0; nu < 4; nu++)
    {
      std::cout << dv[nu] << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  //k0
  contour.loop_momentum(l);
  contour.set_k0();
  contour.set_lambda();
  std::cout << "A:\t[k0]\n";
  for (int nu = 0; nu < 4; nu++)
  {
    for (int mu = 0; mu < 4; mu++)
      std::cout << contour.gradk0[mu][nu] << ", ";
    std::cout << std::endl;
  }
  std::cout << "\nN:\t[\n";
  v = contour.k0;
  for (int mu = 0; mu < 4; mu++)
  {
    ep = l[mu] * eps;
    dl = l + ep * ev[mu];

    contour.loop_momentum(dl);
    contour.set_k0();
    dv = contour.k0;

    dv = 1. / ep * (dv - v);
    for (int nu = 0; nu < 4; nu++)
    {
      std::cout << dv[nu] << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
*/
  return 0;
}
#endif