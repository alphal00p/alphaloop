#ifndef _INTEGRAND_INTERFACE_H_
#define _INTEGRAND_INTERFACE_H_

#include "DCD.h"

#define _OP_INTEGRAND_ID 1
#define _OP_CHANNEL_ID 2
#define _OP_UVSQ 3


extern "C"
{
  int set_factor_int(short int op_id, int value);
  int set_factor_complex(short int op_id, double real, double imag);
  int append_Q(double v[]);
  void reset_Q();
  int evaluate(double* loop_momentum_real, double* loop_momentum_imag, double& f_real, double& f_imag);
}

#endif