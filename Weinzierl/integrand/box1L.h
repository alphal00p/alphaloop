#ifndef _BOX1L_H_
#define _BOX1L_H_

#include "DCD.h"

//box1L
my_comp box1L_6d(DIdeform::C4vector & ell, std::vector<DIdeform::R4vector> & Qs);
my_comp box1L_offshell(DIdeform::C4vector & ell, std::vector<DIdeform::R4vector> & Qs, int ch_id);
my_comp box1L_subtracted(DIdeform::C4vector & ell, std::vector<DIdeform::R4vector> & Qs);
my_comp box1L_subtracted_ch(DIdeform::C4vector &ell, std::vector<DIdeform::R4vector> &Qs, int ch_id);
my_comp box1L_one_offshell_subtracted(DIdeform::C4vector & ell, std::vector<DIdeform::R4vector> & Qs, my_comp& mu_UVsq, my_comp& s12,  my_comp& s23);

#endif
