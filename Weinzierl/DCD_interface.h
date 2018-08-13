#ifndef _DCD_INTERFACE_H_
#define _DCD_INTERFACE_H_

#include "DCD.h"


#define _OP_M1_FACTOR 0
#define _OP_M2_FACTOR 1
#define _OP_M3_FACTOR 2
#define _OP_GAMMA1_FACTOR 3
#define _OP_GAMMA2_FACTOR 4

/***********************************************************
 * The error are coded as follows:
 * 
 *    return 1   -> wrong dimension for the input
 *    return 99  -> No function found for the corresponding option
 *    return 101 -> deformer has not been created.
 * 
 * *********************************************************/

extern "C"
{
  //set factors using option ids (python function)
  int set_factor_int(short int op_id , int v[], int d);
  int set_factor_double(short int op_id, double v[],int d);
  
  //set factors function (internal)
  int set_M(double v[], int d, double & M, bool & external_M);
  int set_gamma(double v[], int d, double & gamma, bool & external_gamma);
  
  //First one needs to give the Qs
  int append_Q(double v[], int d);
  
  //Optionally one can set Pp and Pm
  int set_Pp(double v[], int d);
  int set_Pm(double v[], int d);
  
  // Set up the deforme with the above settings
  int init();

  // Apply deformation for one given loop momentum point
  int deform_loop_momentum(double l[], int d);
  
  //Get the results of the deformation
  double *get_jacobian();
  my_real *get_deformed_loop_momentum();
  void clear();
}

#endif