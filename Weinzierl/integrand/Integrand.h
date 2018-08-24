#ifndef _INTEGRAND_H_
#define _INTEGRAND_H_

#include "DCD.h"

namespace cuba_integrand
{
//box1L
my_comp box1L_6d(DIdeform::C4vector & ell, std::vector<DIdeform::R4vector> & Qs);
my_comp box1L_offshell(DIdeform::C4vector & ell, std::vector<DIdeform::R4vector> & Qs);
my_comp box1L_subtracted(DIdeform::C4vector & ell, std::vector<DIdeform::R4vector> & Qs);
my_comp box1L_subtracted_ch(DIdeform::C4vector &ell, std::vector<DIdeform::R4vector> &Qs, int ch_id);
my_comp box1L_one_offshell_subtracted(DIdeform::C4vector & ell, std::vector<DIdeform::R4vector> & Qs, my_comp& mu_UVsq);

};//namespace Integrand

#endif