#include "Integrand.h"

using namespace DIdeform;

namespace cpp_integrand
{
extern my_real s12,s23;
//Imaginary I
my_comp ii(0.0, 1.0);

//Collinear auxiliary vectors
DIdeform::R4vector eta1({1, 1, 0, 0});
DIdeform::R4vector eta2({1, -1, 0, 0});

/* One loop box if 6-dimension.
 * Two dimensions are integrated out analyticly in the rest frame 
 * where
 * 
 * p1 = 0.5 * (-1,-1, 0,  0 )
 * p2 = 0.5 * (-1,-1, 0,  0 )
 * p1 = 0.5 * ( 1, cos(x), sin(x),  0 )
 * p1 = 0.5 * ( 1,-cos(x),-sin(x),  0 )
 * */
my_comp box1L_6d(C4vector &ell, std::vector<R4vector> &Qs)
{
    //    C4vector prop_mom;
    my_real factor = (4.0 * M_PI) / 2.0;

    //Denominator
    my_comp denominator = 1.;
    for (int i = 0; i < 4; i++)
        denominator *= (ell - Qs[i]).square();
    //Numerator
    my_comp r = std::pow(ell(3), 2);
    return factor * r / denominator;
}

/* One loop box with off-shell external momenta */
my_comp box1L_offshell(C4vector &ell, std::vector<R4vector> &Qs)
{
    C4vector prop_mom;
    my_comp factor = 1.0 / ii / std::pow(M_PI, 2);

    //Denominator
    my_comp denominator = 1.;
    for (int i = 0; i < 4; i++)
        denominator *= (ell - Qs[i]).square();

    return factor / denominator;
}

/* One loop box with on-shell external momenta 
 * The poles have been removed by the subtraction term
 * (1-A13/s - A24/t)
 * */
my_comp box1L_subtracted(C4vector &ell, std::vector<R4vector> &Qs)
{
    C4vector prop_mom;
    my_comp factor = 1.0 / ii / std::pow(M_PI, 2);

    //Denominator
    my_comp denominator = 1.;
    for (int i = 0; i < 4; i++)
        denominator *= (ell - Qs[i]).square();

    //Regulator
    my_comp F = 1;
    my_real sij;
    for (int i = 0; i < 4; i++)
    {
        sij = (Qs[i % 2] - Qs[i % 2 + 2]).square();
        F -= (ell - Qs[i]).square() / sij;
    }

    return factor * F / denominator;
}

/* One loop box with on-shell external momenta 
 * The poles have been removed by the subtraction term
 * (1-A13/s - A24/t)
 * Channel options
 */
 my_comp box1L_subtracted_ch(C4vector &ell, std::vector<R4vector> &Qs, int ch_id)
{
    my_comp factor = 1.0 / ii / std::pow(M_PI, 2);

    //Denominator
    my_comp denominator = 1.;
    for (int i = 0; i < 4; i++)
        denominator *= (ell - Qs[i]).square();

    //Regulator
    my_comp F = 1;
    my_real sij;
    for (int i = 0; i < 4; i++)
    {
        sij = (Qs[i % 2] - Qs[i % 2 + 2]).square();
        F -= (ell - Qs[i]).square() / sij;
    }

    //Channelling a la Weinzierl
    my_comp MC_factor=0.0;
    my_comp tmp = 0;
    for (int i = 0; i < 3; i++)
    {
        tmp = std::abs((ell - Qs[i]).square()) * std::abs((ell - Qs[i + 1]).square());
        MC_factor += tmp;

        if (i == ch_id) {
            factor /= tmp;
            factor /= tmp;
        }
    }

    return factor * F / denominator / MC_factor / MC_factor;
}

/* One loop box with one off-shell external momenta 
 * The poles have been removed by the subtraction term
 * I = (1-A13/s - A24/t)/(A1 A2 A3 A4) 
 *     - UV_term_2/(A2 A3 s t xbar_2) 
 *     - UV_term_3/(A3 A4 s t x_3)
 * */

inline my_comp box1L_collinear_term(bool x_xbar,
                                    const DIdeform::R4vector &q1,
                                    const DIdeform::R4vector &q2,
                                    const DIdeform::C4vector &ell,
                                    const my_comp &mu2)
{
    //True: compute with x otherwise with 1-x
    DIdeform::R4vector p=q2-q1;
    //TODO: global s and t, ideally sij and sij_inv
    my_comp cFactor=1.0/s12/s23;
    
    //Divide by x or 1-x
    DIdeform::R4vector eta;
    eta = eta1 * p != 0 ? eta1 : eta2; 

    my_comp x = (eta * ell)/(eta * p);
    x = x_xbar ? x : 1.0 - x;
    cFactor *= 1.0 / x;

    //Divide by A1 A2
    my_comp A1, A2;
    A1 = (ell - q1).square();
    cFactor *= 1.0 / A1;

    A2 = (ell - q2).square();
    cFactor *= 1.0 / A2;

    //UV factor
    cFactor *= mu2 / (mu2 + A1);

    return cFactor;
    }

my_comp box1L_one_offshell_subtracted(C4vector &ell, std::vector<R4vector> &Qs, my_comp& mu_UVsq)
{
    C4vector prop_mom;
    my_comp factor = 1.0 / ii / std::pow(M_PI, 2);
    
    //Denominator
    my_comp denominator = 1.;
    for (int i = 0; i < 4; i++)
        denominator *= (ell - Qs[i]).square();

    //Soft Regulator
    my_comp soft_F = 1;
    my_real sij;
    for (int i = 2; i < 4; i++)
    { //As for the massless box but without A0
        sij = (Qs[i % 2] - Qs[i % 2 + 2]).square();
        soft_F -= (ell - Qs[i]).square() / sij;
    }
    
    //Collinear
    DIdeform::R4vector p;
    my_comp coll_F = 0;
    coll_F -= box1L_collinear_term(true , Qs[3], Qs[0], ell, mu_UVsq);
    coll_F -= box1L_collinear_term(false, Qs[1], Qs[2], ell, mu_UVsq);
    return factor * (soft_F / denominator - coll_F);
    }


}; //namespace Integrand
