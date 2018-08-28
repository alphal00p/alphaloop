#ifndef _DCD_INTERFACE_H_
#define _DCD_INTERFACE_H_

#include "DCD.h"


#define _OP_M1_FACTOR 11
#define _OP_M2_FACTOR 12
#define _OP_M3_FACTOR 13
#define _OP_M4_FACTOR 14

#define _OP_GAMMA1_FACTOR 21
#define _OP_GAMMA2_FACTOR 22

#define _OP_MAPPING 31
#define _OP_CHANNEL_ID 41

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
  int update_M(double v[], int d, my_real & deformer_M);
  int set_gamma(double v[], int d, double & gamma, bool & external_gamma);
  int update_gamma(double v[], int d, my_real & deformer_gamma);

  int set_mapping(int v[], int d, int &mapping, bool &external_mapping);
  int update_mapping(int v[], int d, short int& deformer_mapping);
  int set_channel_id(int v[], int d, int &channel, bool &external_channel_id);
  int update_channel_id(int v[], int d, short int& channel_id);

  //First one needs to give the Qs
  int append_Q(double v[], int d);
  
  //Optionally one can set Pp and Pm
  int set_Pp(double v[], int d);
  int set_Pm(double v[], int d);
  
  // Set up the deforme with the above settings
  int init();

  double* hypcub_mapping(double x[], int d, double* jacobian);

  // Apply deformation for one given loop momentum point
  int deform_loop_momentum(double l[], int d);
  
  //Get the results of the deformation
  double *get_jacobian();
  my_real *get_deformed_loop_momentum();
  void clear();
}

#endif