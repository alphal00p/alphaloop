#include "DCD.h"
#include "integrand_interface.h"

extern "C"
{

int integrand_id;
int channel_id;
my_comp mu_UVsq, s12, s23;
std::vector<DIdeform::R4vector> Qs;

int set_factor_int(short int op_id, int value)
{
    switch (op_id) 
    {
        case _OP_INTEGRAND_ID:
          integrand_id = value;
          return 0;
        case _OP_CHANNEL_ID:
        {
          channel_id = value;
          return 0;
        }
        default:
          return 99;
    }
}

int set_factor_complex(short int op_id, double real, double imag)
{
    switch (op_id) 
    {
        case _OP_UVSQ:
          mu_UVsq.real(real);
          mu_UVsq.imag(imag);
          return 0;
        case _OP_S12:
          s12.real(real);
          s12.imag(imag);
          return 0;
        case _OP_S23:
          s23.real(real);
          s23.imag(imag);
          return 0;
        default:
          return 99;
    }
}

//First one needs to give the Qs
int append_Q(double v[])
{
  DIdeform::R4vector q;
  for (int i = 0; i < 4; i++)
  {
      q[i] = v[i];
  }
  Qs.push_back(q);

  return 0;
}

void reset_Q()
{
  Qs.clear();
}

int evaluate(double* loop_momentum_real, double* loop_momentum_imag, double& f_real, double& f_imag)
{

  DIdeform::C4vector ell;
  for (int i = 0; i < 4; i++) {
    ell[i].real(loop_momentum_real[i]);
    ell[i].imag(loop_momentum_imag[i]);
  }

  my_comp function;
  switch (integrand_id)
  {
  /* === BOX 6d === */
  case 0:
    function = box1L_6d(ell, Qs);
    break;

  /* === BOX OFF-SHELL === */
  case 1:
    function = box1L_offshell(ell, Qs);
    break;

  /* === BOX ALL ON-SHELL SUBTRACTED === */
  case 2:
    function = box1L_subtracted(ell, Qs);
    break;
   /* === BOX ALL ON-SHELL SUBTRACTED === */
 
  case 3:
    function = box1L_subtracted_ch(ell, Qs, channel_id);
    break;
  /* === BOX ONE OFF-SHELL SUBTRACTED === */
 
  case 4:
    function = box1L_one_offshell_subtracted(ell, Qs, mu_UVsq, s12, s23);
    break;
  default:
    return 1;
  };

  f_real = function.real();
  f_imag = function.imag();

  return 0;
}

}