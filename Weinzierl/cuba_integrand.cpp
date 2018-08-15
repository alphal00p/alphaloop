#include "cuba.h"
#include "DCD.h"
#include <fstream>

#ifndef ROTATION_ANGLE
#define ROTATION_ANGLE M_PI/8.
#endif

#define WHITCH_INTEGRAND 2 //0:box_6d, 1:box_offshell, 2:box_subtracted, 3:test_function
#define PS_SEED 22

using namespace std;

my_comp ii(0.0, 1.0);
DIdeform::R4vector p1, p2, p3, p4;
extern my_real alpha;

ofstream file("points.txt", ios::out);
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

    DIdeform::R4vector l;
    DIdeform::C4vector prop_mom;
    DIdeform::C4vector ell;

    my_comp hypercube_jacobian;
    my_comp deformation_jacobian;
    my_comp Njacobian;
    my_comp jacobian;

    PS_points(p1,p2,p3,p4,PS_SEED);
    std::vector<DIdeform::R4vector> Qs(4);

    Qs[0] = p1;
    Qs[1] = Qs[0] + p2;
    Qs[2] = Qs[1] + p3;
    Qs[3] = 0.0 * Qs[3];

    DIdeform::R4vector shift(Qs[0]);
    for (int i = 0; i < 4; i++)
      Qs[i] = Qs[i] - 0.0 * shift;

    DIdeform::ContourDeform contour(Qs);
    contour.lambda_max = 1.0;

 // cout << k0 << endl;
 //  std::printf("{%f, %f, %f, %f}\n",k0,k1,k2,k3);
 //  exit(1); 
    //Variable map
    alpha = std::sqrt(contour.mu_P);
    l = DIdeform::hypcub_mapping({k0, k1, k2, k3});
    hypercube_jacobian = DIdeform::hypcub_jacobian({k0, k1, k2, k3});
    
    contour.loop_momentum(l);
    contour.deform(ell, deformation_jacobian);

    my_comp f;
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

    switch (WHITCH_INTEGRAND)
    {
    /* === BOX 6d === */
    case 0: 
    {
    //Denominator
    my_comp denominator = 1.;
    for (int i = 0; i < 4; i++)
    {
      for (int mu = 0; mu < 4; mu++)
        prop_mom[mu] = ell(mu) - Qs[i](mu);
      denominator = denominator * (prop_mom * prop_mom);
    }
    //Numerator
    my_comp r = std::pow(ell(3), 2);
    factor = (4.0 * M_PI) / 2.0;
    f = factor * jacobian * r / denominator;
  }
  break;
  
  /* === BOX OFFSHELL === */
  case 1: 
  {
    //Denominator
    my_comp denominator = 1.;
    for (int i = 0; i < 4; i++)
    {
      for (int mu = 0; mu < 4; mu++)
        prop_mom[mu] = ell(mu) - Qs[i](mu);
      denominator = denominator * (prop_mom * prop_mom);
    }
    f = /*1e10 */ (1.0/ii/std::pow(M_PI,2)) * jacobian  / denominator;
  }
  break;
  
  /* === BOX SUBTRACTED === */
  case 2: 
  {
    //Denominator
    my_comp denominator = 1.;
    for (int i = 0; i < 4; i++)
    {
      for (int mu = 0; mu < 4; mu++)
        prop_mom[mu] = ell(mu) - Qs[i](mu);
      denominator = denominator * (prop_mom * prop_mom);
    }
    
    //Regulator
    my_comp F=1;
    my_real sij;
    for (int i = 0; i < 4; i++)
    {
      for (int mu = 0; mu < 4; mu++)
        prop_mom[mu] = ell(mu) - Qs[(i+0)%4](mu);
      sij = (Qs[i % 2] - Qs[i % 2 + 2]) * (Qs[i % 2] - Qs[i % 2 + 2]);
      F -= (prop_mom * prop_mom) / sij;
    }

    f = /*1e10 */ (1.0/ii/std::pow(M_PI,2)) * F * jacobian / denominator;
  }
  break;

    /* === TEST INTEGRAL === */
  case 3:
  {
    f = jacobian * 1. / (std::pow((ell * ell.dual()) / std::pow(1, 2), 3) + 50.);
  }
  break;
  }
  
  f_real = f.real();
  f_imag = f.imag();

  return 0;
}

//Phasespace points
void PS_points(DIdeform::R4vector &p1, DIdeform::R4vector &p2, DIdeform::R4vector &p3, DIdeform::R4vector &p4, int seed)
{
  switch (seed)
  {
  case 1666:
  {
    //seed 666
    p1 = DIdeform::R4vector({-4.7213384875835902e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, -4.6142211817746778e+02});
    p2 = DIdeform::R4vector({-5.0290194983056193e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, 4.6142211817746778e+02});
    p3 = DIdeform::R4vector({4.5162178137689028e+02, -3.3466005003229913e+02, -1.3177690342268895e+00, -4.4307423883444024e+01});
    p4 = DIdeform::R4vector({5.2341401721203044e+02, 3.3466005003229913e+02, 1.3177690342268895e+00, 4.4307423883444024e+01});
  }
  break;

  case 12:
  {
    //seeds 2 massive
    p1 = DIdeform::R4vector({-4.7809952420694083e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, -4.6752449673455959e+02});
    p2 = DIdeform::R4vector({-5.0850678957797919e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, 4.6752449673455959e+02});
    p3 = DIdeform::R4vector({4.5782801395958194e+02, 1.3758384614384497e+02, 8.1217573038820291e+01, -3.0672606911725950e+02});
    p4 = DIdeform::R4vector({5.2877829982533808e+02, -1.3758384614384497e+02, -8.1217573038820291e+01, 3.0672606911725950e+02});
  }
  break;

  case 22:
  {
    //seeds 2 massless
    p1 = DIdeform::R4vector({-4.8678217388064405e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, -4.8678217388064405e+02});
    p2 = DIdeform::R4vector({-4.8678217388064405e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, 4.8678217388064405e+02});
    p3 = DIdeform::R4vector({4.8678217388064405e+02, 1.9365322696179936e+02, 1.1431607376733305e+02, -4.3172577844468481e+02});
    p4 = DIdeform::R4vector({4.8678217388064405e+02, -1.9365322696179936e+02, -1.1431607376733305e+02, 4.3172577844468481e+02});
  }
  break;

  case 133:
  {
    //seed 33 massive
    p1 = DIdeform::R4vector({-4.3403531061888162e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, -4.2235843884552497e+02});
    p2 = DIdeform::R4vector({-4.6731857534665693e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, 4.2235843884552497e+02});
    p3 = DIdeform::R4vector({4.1184646746703140e+02, -3.7233508163465089e+01, 2.1500619135292814e+02, 1.7889526632871485e+02});
    p4 = DIdeform::R4vector({4.8950741849850709e+02, 3.7233508163465089e+01, -2.1500619135292814e+02, -1.7889526632871485e+02});
  }
  break;
  case 233:
  {
    //seed 33 massless
    p1 = DIdeform::R4vector({-3.9756551766525951e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, -3.9756551766525951e+02});
    p2 = DIdeform::R4vector({-3.9756551766525951e+02, -0.0000000000000000e+00, -0.0000000000000000e+00, 3.9756551766525951e+02});
    p3 = DIdeform::R4vector({3.9756551766525951e+02, -5.2461217332126608e+01, 3.0293912899098285e+02, 2.5205960731275806e+02});
    p4 = DIdeform::R4vector({3.9756551766525951e+02, 5.2461217332126608e+01, -3.0293912899098285e+02, -2.5205960731275806e+02});
  }
  break;

  default:
  {
    //Angular seeds
    p1 = DIdeform::R4vector({-0.5, -0.5, 0.0, 0.0});
    p2 = DIdeform::R4vector({-0.5, 0.5, 0.0, 0.0});
    p3 = DIdeform::R4vector({0.5, -0.5 * std::cos(ROTATION_ANGLE), -0.5 * std::sin(ROTATION_ANGLE), 0.0});
    p4 = DIdeform::R4vector({0.5, 0.5 * std::cos(ROTATION_ANGLE), 0.5 * std::sin(ROTATION_ANGLE), 0.0});
  }
  break;
  }
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
