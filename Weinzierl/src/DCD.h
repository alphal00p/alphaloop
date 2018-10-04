/*
Direct Contour Deformation from Sebastian Becker and Stefan Weinzierl
 */

#ifndef _DCD_H_
#define _DCD_H_

#include <stdio.h>
#include <iostream>
#include <vector>
#include <complex>
#include <limits>
#include <cmath>

/* Define type  */
typedef double                my_real;
typedef std::complex<my_real> my_comp;

namespace DIdeform
{

  class R4vector
  {
  public:
    std::vector<my_real> v;

    R4vector();
    R4vector(const std::vector<my_real> &w);

    my_real &operator[](const int mu);
    const my_real &operator()(const int mu) const;
    R4vector dual() const;
    my_real square() const;

    //Operator
    R4vector operator+(const R4vector &) const;
    R4vector operator-();
    R4vector operator-(const R4vector &) const;
    my_real operator*(const R4vector &) const;
    R4vector operator=(const std::vector<my_real> &);
  };

  R4vector operator*(const my_real, const R4vector &);
  std::ostream &operator<<(std::ostream &os, DIdeform::R4vector rhs);

  class C4vector
  {
  public:
    std::vector<my_comp> v;

    C4vector();
    C4vector(const std::vector<my_comp> &w);

    my_comp &operator[](const int mu);
    const my_comp &operator()(const int mu) const;
    C4vector dual() const;
    my_comp square() const;

    //Operator
    C4vector operator+(const C4vector &) const;
    C4vector operator-();
    C4vector operator-(const C4vector &) const;
    my_comp operator*(const C4vector &) const;
    C4vector operator=(const std::vector<my_comp> &);
  };

  C4vector operator*(const my_comp, const C4vector &);
  std::ostream &operator<<(std::ostream &os, DIdeform::C4vector rhs);
  
  //Mixed operator C4vector R4vector
  C4vector operator*(const my_comp, const R4vector &);
  
  C4vector operator+(const R4vector &, const C4vector &);
  C4vector operator+(const C4vector &, const R4vector &);
  
  C4vector operator-(const R4vector &, const C4vector &);
  C4vector operator-(const C4vector &, const R4vector &);
    
  my_comp operator*(const R4vector &, const C4vector &);
  my_comp operator*(const C4vector &, const R4vector &);

  //Map the hypercube to the momentum space
  R4vector hypcub_mapping(std::vector<my_real> x);
  my_real hypcub_jacobian(std::vector<my_real> x);

  my_comp Determinant(const std::vector<std::vector<my_comp>> &bb);

  /*The contour deformation depends only on the q_i momenta
    We thake p_n to be an incoming momenta and q_n to be zero.
    Note that the external OUTGOING momenta can be recovered by
      p_i= q_i - q_(i-1)
   */
  class ContourDeform
  {
  public:
    //Metric tensor
    std::vector<my_real> gmunu = {1., -1., -1., -1.};

    //Topology information
    int A;    //Position of second incoming momentum
    int legs; // Number of external legs

    //Global variables
    my_real M1f = 0.035, M2f = 0.7, M3f = 0.035, M4f = 0.035;
    my_real gamma1 = 0.7;
    my_real gamma2 = 0.008;
    my_real E_soft = 0.03;
    my_real Ecmsq, mu_Psq, M1sq, M2sq, M3sq, M4sq;
    my_comp mu_UVsq = my_comp(0.0, -1e5); // should scale with problem
    short int which_hypercube_map = 0;
    short int channel_id = -1;

    //Auxiliary vectors
    R4vector Pp, kp; //zp(qa+q0, qa-q0)
    R4vector Pm, km; //zm(q(a-1)+q(n-1), q(a-1)-q(n-1))
    R4vector g0mu = R4vector({1., 0., 0., 0.});

    //Coefficents
    my_real cp, cm; // prod(i=1,n, hpm(ki,mi))

    //Helper functions
    my_real hp(R4vector &k, my_real mass, my_real Msq);
    my_real hm(R4vector &k, my_real mass, my_real Msq);
    my_real h0(R4vector &k, my_real mass, my_real Msq);
    my_real ht(my_real t, my_real Msq);

    R4vector gradloghp(R4vector &k, my_real mass, my_real Msq);
    R4vector gradloghm(R4vector &k, my_real mass, my_real Msq);
    R4vector gradlogh0(R4vector &k, my_real mass, my_real Msq);
    my_real dloght(my_real t, my_real Msq);

    //Exterior region
    R4vector zp(R4vector x, R4vector y);
    R4vector zm(R4vector x, R4vector y);
    R4vector k_ext;
    std::vector<std::vector<my_real>> gradk_ext;
    void set_k_ext();
    void set_gradk_ext();

    //Interior region
    my_real g(R4vector &k, my_real Msq);
    my_real dij(int i, int j);
    my_real dijk(int i, int j, int k);
    my_real ci(int i);
    my_real cij(int i, int j);
    R4vector vij(int i, int j);

    R4vector gradlogg(R4vector &k, my_real Msq);
    R4vector gradlogdij(int i, int j);
    R4vector gradlogdijk(int i, int j, int k);
    R4vector gradlogci(int i);
    R4vector gradlogcij(int i, int j);

    R4vector k_center;
    R4vector k_int;
    std::vector<std::vector<my_real>> gradk_int;
    void set_k_int();
    void set_gradk_int();

    //Soft region
    std::vector<R4vector> ei;
    R4vector k_soft;

    //Deformation vector
    R4vector k0;
    std::vector<std::vector<my_real>> gradk0;
    void set_k0();

    //Scaling parameter
    my_real lambda;
    my_real lambda_max = 1.;
    std::vector<my_real> gradlambda;

    my_real Xi(int i);
    my_real Yi(int i);
    my_real lambda_i(int i);
    my_real lambda_coll();
    my_real lambda_UV();
    //my_real lambda_uv(int i);
    void set_lambda();

    R4vector gradXi(int i);
    R4vector gradYi(int i);
    R4vector gradlambda_i(int i);
    R4vector gradlambda_coll();
    R4vector gradlambda_UV();

    //Conformal Mapping
    R4vector weinzierl_mapping(std::vector<my_real> x, my_real& jacobian);

    //Derivatives

    //Constructor Function
    void set_PpPm(std::vector<R4vector> &Qs);
    bool test_PpPm(std::vector<R4vector> &Qs);
    void get_PpPm(short int plus_or_minus, std::vector<R4vector> &qs);
    void set_global_var();
    //void Find_Topology();

  public:
    R4vector l; //loop momentum
    R4vector UV_offset;

    std::vector<R4vector> qi;  // qi are given as input
    std::vector<R4vector> pi;  // pi are the external momenta
    std::vector<R4vector> lqi; // l-qi is the propagator momentum

    bool masslessQ = false;  // If true al masses are dropped
    std::vector<my_real> mi; // masses squared

    //Constructor
    ContourDeform(std::vector<R4vector> &);
    
    //Callable Functions
    R4vector hypcub_mapping(std::vector<my_real> x, my_real& jacobian);
    void loop_momentum(R4vector &loopmomentum);
    void deform(C4vector &k, my_comp &jacobian);
  };
}; // namespace DIdeform

#endif
