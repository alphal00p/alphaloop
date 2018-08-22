#include "cuba.h"
#include "DCD.h"
#include "Integrand.h"
#include <fstream>

using namespace std;

extern my_real alpha;

DIdeform::R4vector l;
DIdeform::C4vector prop_mom;
DIdeform::C4vector ell;

my_comp hypercube_jacobian;
my_comp deformation_jacobian;
my_comp Njacobian;
my_comp jacobian;

 //Defined in cuba.c 
extern int which_integrand;
extern std::vector<DIdeform::R4vector> Qs;
extern DIdeform::ContourDeform * deformer;
//ofstream file("points.txt", ios::out);

void PS_points(DIdeform::R4vector &, DIdeform::R4vector &, DIdeform::R4vector &, DIdeform::R4vector &, int seed = 0);
my_comp numeric_jacobian(DIdeform::ContourDeform &deformer, DIdeform::R4vector l);

int Integrand(const int *ndim, const cubareal xx[],
              const int *ncomp, cubareal ff[], void *userdata)
{

#define k0 xx[0]
#define k1 xx[1]
#define k2 xx[2]
#define k3 xx[3]
#define f_real ff[0]
#define f_imag ff[1]


  alpha = std::sqrt(deformer->mu_P);
  
  l = DIdeform::hypcub_mapping({k0, k1, k2, k3});
  hypercube_jacobian = DIdeform::hypcub_jacobian({k0, k1, k2, k3});

  deformer->loop_momentum(l);
  deformer->deform(ell, deformation_jacobian);

  my_comp function;
  my_real factor = 1.;

  //Check agreement with Numeric Jacobian
  /*Njacobian = numeric_jacobian(contour, l);
    if (std::abs((Njacobian - deformation_jacobian) / deformation_jacobian) > 0.05)
    {
      char s_point[100];
      std::sprintf(s_point, "{%+.5e, %+.5e, %+.5e, %+.5e}",
                   l(0), l(1), l(2), l(3));

      file << s_point
           << "\t" << deformation_jacobian
           << "\t" << Njacobian << "\n";
      file.flush();
    }
    */
  jacobian = hypercube_jacobian * deformation_jacobian;
  function = jacobian;

  switch (which_integrand)
  {
  /* === BOX 6d === */
  case 0:
    function *= cuba_integrand::box1L_6d(ell, Qs);
    break;

  /* === BOX OFFSHELL === */
  case 1:
    function *= cuba_integrand::box1L_offshell(ell, Qs);
    break;

  /* === BOX SUBTRACTED === */
  case 2:
    function *= cuba_integrand::box1L_subtracted(ell, Qs);
    break;
  /* === BOX SUBTRACTED === */
  case 3:
    function *= cuba_integrand::box1L_subtracted(ell, Qs);
    break;
  };
  
  f_real = function.real();
  f_imag = function.imag();

  return 0;
  }

 //Numeric Jacobian
  my_comp numeric_jacobian(DIdeform::ContourDeform &deformer, DIdeform::R4vector l)
  {
    my_comp Njacobian;

    my_real eps = std::sqrt(std::numeric_limits<my_real>::epsilon());
    std::vector<std::vector<my_comp>> grad(4, std::vector<my_comp>(4, 0.));
    std::vector<DIdeform::R4vector> ev(4);
    for (int mu = 0; mu < 4; mu++)
      ev[mu][mu] = 1.0;

    DIdeform::R4vector dl;
    DIdeform::C4vector ell;
    DIdeform::C4vector dell;

    //Initial position
    deformer.deform(ell, Njacobian);

    my_real ep;
    for (int mu = 0; mu < 4; mu++)
    {
      ep = std::max(l[mu] * eps, eps);
      dl = l + ep * ev[mu];

      deformer.loop_momentum(dl);
      deformer.deform(dell, Njacobian);

      dell = 1. / ep * (dell - ell);
      for (int nu = 0; nu < 4; nu++)
      {
        grad[mu][nu] = dell(nu);
      }
    }
    Njacobian = DIdeform::Determinant(grad);
    return Njacobian;
  }
