/*
  Path deformation for integration in direct space of 1-loop integrals
*/
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
using namespace std;

#include <vector>
#include <math.h>
#include <complex>


#define dimensions 4

/*Typedef for vectors */
typedef double my_real;
typedef complex<double> my_complex;

typedef vector<complex<my_real> > cvector;
typedef vector<my_real> rvector;

#ifndef ROTATION_ANGLE
#define ROTATION_ANGLE M_PI/4.
#endif

/********************** Class **********************/

class PathDeformation
{
 private:
  
  /*deformation global variables*/
  rvector v     = rvector(dimensions,0.);
  rvector kappa = rvector(dimensions,0.);
  my_real M1=0.,M2=0.,M3=0.;
  my_real x=0., xbar=0.;

  /*Check if variables are initialized*/
  bool momentaQ=false;
  
public:
  
  /*Tolopogy variable*/
  //external momenta
  rvector p1 = rvector(dimensions, 0.);
  rvector p2 = rvector(dimensions,0.);
  rvector p3 = rvector(dimensions,0.);
  rvector p4 = rvector(dimensions,0.);
  //Loop momentum
  rvector k = rvector(dimensions,0.);
  //deformed loop momentum
  cvector ell = cvector(dimensions,0.);
  complex<my_real> jacobian;
  
  //Propagator momenta (l-qi)^2
  rvector q1,q2,q3,q4;
  vector<rvector> Q;
  
  //Propagators
  complex<my_real> prop1;
  complex<my_real> prop2;
  complex<my_real> prop3;
  complex<my_real> prop4;
  void compute_box_propagators();
  
  //Constructor
  PathDeformation(my_real angle);
  bool init_momenta(my_real angle);
  
  // Check if ell is inside the negative light cone starting from Q[j];
  bool ellinCminus(int j);
  // Check if ell is inside the positive light cone starting from Q[j];
  bool ellinCplus(int j);

  /* Deformation Functions */ 
  inline my_real hm(rvector q);  //h- function
  inline my_real hp(rvector q);  //h+ function

  inline my_real g(rvector& q);  //g(l) function
  inline my_real gp(rvector& q); //g+ function
  inline my_real gm(rvector& q); //g- function

  /* Derivatives */ 
  rvector gradloghm(const rvector k); // Gradient of log(h-(k));
  rvector gradloghp(const rvector k); // Gradient of log(h+(k));
    
  rvector gradlogg(const rvector k);   // Gradient of log(g(k));
  rvector gradloggp(const rvector k);  // Gradient of log(g+(k));
  rvector gradloggm(const rvector k);  // Gradient of log(g-(k));

  // l is deformed to ell with the corresponding jacobian
  bool loop_setQ=false, deformedQ=false;
  void set_loop_mom(my_real l0,my_real l1,my_real l2,my_real l3);
  void deform_loop_path();

  
};
