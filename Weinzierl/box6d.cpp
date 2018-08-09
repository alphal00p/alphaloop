#include "cuba.h"
#include "DCD.h"
#include <fstream>

#ifndef ROTATION_ANGLE
#define ROTATION_ANGLE M_PI/2.
#endif

using namespace std;

DIdeform::R4vector p1({0.5, 0.5, 0.0, 0.0});
DIdeform::R4vector p2({0.5,-0.5, 0.0, 0.0});
DIdeform::R4vector p3({0.5, 0.5 * std::cos(ROTATION_ANGLE ), 0.5 * std::sin(ROTATION_ANGLE), 0.0});


ofstream file("points.txt",ios::out);

int Integrand(const int *ndim, const cubareal xx[],
	      const int *ncomp, cubareal ff[], void *userdata) {

#define k0 xx[0]
#define k1 xx[1]
#define k2 xx[2]
#define k3 xx[3]
#define f_real ff[0]
#define f_imag ff[1]

  DIdeform::R4vector l;
  DIdeform::C4vector prop_mom;
  DIdeform::C4vector ell;
  my_comp hypercube_jacobian;
  my_comp deformation_jacobian;
  my_comp Njacobian;
  my_comp jacobian;
  std::vector<DIdeform::R4vector> Qs(4);

  Qs[1] = -p2;
  Qs[2] = p1 - p3;
  Qs[3] = p1;
  Qs[0] = 0.0 * Qs[0];
  
  DIdeform::R4vector shift(Qs[3]);
  for(int i=0; i<4; i++)
    Qs[i] = Qs[i];// - shift;

  //Variable map
  l = DIdeform::hypcub_mapping({k0,k1,k2,k3});
  hypercube_jacobian = DIdeform::hypcub_jacobian({k0, k1, k2, k3});

  DIdeform::ContourDeform contour(Qs);
  contour.lambda_max = 1.0;
  
  contour.loop_momentum(l);
  contour.deform(ell,deformation_jacobian);
  deformation_jacobian = abs(deformation_jacobian);

  my_comp f;
  my_real factor = 1.;
  
  //Numeric Jacobian
  my_real eps = sqrt(std::numeric_limits<my_real>::epsilon());
  std::vector<std::vector<my_comp>> grad(4, std::vector<my_comp>(4, 0.));
  std::vector<DIdeform::R4vector> ev(4);
  for (int mu = 0; mu < 4; mu++)
    ev[mu][mu] = 1.0;

  DIdeform::R4vector dl;
  DIdeform::C4vector dell;

  my_real ep;
  for (int mu = 0; mu < 4; mu++)
  {
    ep = std::max(l[mu] * eps, eps);
    dl = l + ep * ev[mu];

    contour.loop_momentum(dl);
    contour.deform(dell, Njacobian);
    
    dell = 1. / ep * (dell - ell);
    for (int nu = 0; nu < 4; nu++)
    {
      grad[mu][nu] = dell(nu);
    }
  }
  Njacobian = abs(DIdeform::Determinant(grad));
  
  if (abs((Njacobian - deformation_jacobian) / deformation_jacobian) > 0.05)
  {
    file << "{"
         << l(0) << ","
         << l(1) << ","
         << l(2) << ","
         << l(3) << "}"
         << "\t" << deformation_jacobian 
         << "\t" << Njacobian
         << std::endl;
    file.flush();
         //exit(1);
  }
    
 // if (abs(Njacobian - jacobian) > .1)
 //   cout << "BOOM" << endl;
  //deformation_jacobian = Njacobian;

  if (!true)
  {
    jacobian = hypercube_jacobian * deformation_jacobian;
    
    //Denominator
    my_comp denominator = 1.;
    for (int i = 0; i < 4; i++)
    {
      for (int mu = 0; mu < 4; mu++)
        prop_mom[mu] = ell(mu) - Qs[i](mu);
      denominator = denominator * (prop_mom * prop_mom);
  }
    //Numerator
    my_comp r = pow(ell(3), 2);

    factor = (4.0 * M_PI) / 2.0;
    f = jacobian * r / denominator;
  }
  else
  {
    jacobian = hypercube_jacobian * deformation_jacobian;
    f = jacobian * 1. / (pow((ell * ell.dual()) / pow(1, 2), 3) + 50.);
  }

  f_real = f.real() * factor;
  f_imag = f.imag() * factor;

  return 0;
}
