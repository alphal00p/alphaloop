#include "DI_PathDef.h"
#include "cuba.h"

#ifndef ROTATION_ANGLE
#define ROTATION_ANGLE M_PI/4.
#endif


/***************** Hypercube mapping *****************/

//Cuba variable mapping
//Map from [0,1] to (-inf,inf)
cubareal VarChange1(cubareal x){ return 1.*log(x/(1-x)); } 
cubareal Jacobian1(cubareal x) { return 1./(x*(1-x)); }
//Map from [0,1] to  [0,inf)
cubareal VarChange2(cubareal x){ return x/(1.0-x); }
cubareal Jacobian2(cubareal x) { return pow(1-x,-2); }

/******************************************************/



int Integrand(const int *ndim, const cubareal xx[],
	      const int *ncomp, cubareal ff[], void *userdata) {

#define k0 xx[0]
#define k1 xx[1]
#define k2 xx[2]
#define k3 xx[3]
#define f_real ff[0]
#define f_imag ff[1]
  
  //Set Variables
  PathDeformation deformed(ROTATION_ANGLE);
  deformed.set_loop_mom(VarChange1(k0)
			,VarChange1(k1)
			,VarChange1(k2)
			,VarChange1(k3));
  //Deform path
  deformed.deform_loop_path();
  
  //Jacobian product
  complex<cubareal> JJ=Jacobian1(k0)*Jacobian1(k1)*Jacobian1(k2)*Jacobian1(k3);
  JJ=JJ*deformed.jacobian;

  //Denominator
  deformed.compute_box_propagators();
  auto denominator = ( deformed.prop1 
		       * deformed.prop2 
		       * deformed.prop3
		       * deformed.prop4);
    
  //Numerator
  complex<cubareal> r=pow(deformed.ell[3],2);
  
  //Integral
  auto integrand=JJ*r/denominator;
  cubareal factor=1.0*(4.0*M_PI)/2.0;
  f_real=integrand.real()*factor;
  f_imag=integrand.imag()*factor;
  
  return 0;
}

