#include "DCD.h"

using namespace DIdeform;

namespace cuba_integrand
{
//Imaginary I
my_comp ii(0.0, 1.0);

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
    C4vector prop_mom;
    my_real factor = (4.0 * M_PI) / 2.0;

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
    {
        for (int mu = 0; mu < 4; mu++)
            prop_mom[mu] = ell(mu) - Qs[i](mu);
        denominator = denominator * (prop_mom * prop_mom);
    }

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
    {
        for (int mu = 0; mu < 4; mu++)
            prop_mom[mu] = ell(mu) - Qs[i](mu);
        denominator = denominator * (prop_mom * prop_mom);
    }

    //Regulator
    my_comp F = 1;
    my_real sij;
    for (int i = 0; i < 4; i++)
    {
        for (int mu = 0; mu < 4; mu++)
            prop_mom[mu] = ell(mu) - Qs[(i + 0) % 4](mu);
        sij = (Qs[i % 2] - Qs[i % 2 + 2]) * (Qs[i % 2] - Qs[i % 2 + 2]);
        F -= (prop_mom * prop_mom) / sij;
    }

    return factor * F / denominator;
}


/* One loop box with on-shell external momenta 
 * The poles have been removed by the subtraction term
 * (1-A13/s - A24/t)
 * Channel options
 * */
my_comp box1L_subtracted_ch(C4vector &ell, std::vector<R4vector> &Qs, int ch_id)
{
    C4vector prop_mom;
    my_comp factor = 1.0 / ii / std::pow(M_PI, 2);

    //Denominator
    my_comp denominator = 1.;
    for (int i = 0; i < 4; i++)
    {
        for (int mu = 0; mu < 4; mu++)
            prop_mom[mu] = ell(mu) - Qs[i](mu);
        denominator = denominator * (prop_mom * prop_mom);
    }

    //Regulator
    my_comp F = 1;
    my_real sij;
    for (int i = 0; i < 4; i++)
    {
        for (int mu = 0; mu < 4; mu++)
            prop_mom[mu] = ell(mu) - Qs[(i + 0) % 4](mu);
        sij = (Qs[i % 2] - Qs[i % 2 + 2]) * (Qs[i % 2] - Qs[i % 2 + 2]);
        F -= (prop_mom * prop_mom) / sij;
    }
    
    //Channelling 
    my_comp MC_factor=0.0;
    for (int i = 0; i < 4; i++)
    {
        for (int mu = 0; mu < 4; mu++)
            prop_mom[mu] = ell(mu) - Qs[i](mu);
        MC_factor += std::pow(prop_mom * prop_mom, 4);
        if (i == ch_id)
            factor *= pow(prop_mom * prop_mom,4);
    }

    return factor * F / denominator / MC_factor;
}
}; //namespace Integrand
