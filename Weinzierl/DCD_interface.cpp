#include "DCD.h"
#include "DCD_interface.h"

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
    //Deformer
    DIdeform::ContourDeform *deformer = NULL;

    //Deformer arguments
    std::vector<DIdeform::R4vector> Qs;
    DIdeform::R4vector Pp, Pm, loop_momentum;
    double M1_factor, M2_factor, M3_factor, M4_factor;
    double gamma1, gamma2;
    
    //Additional arguments flags
    bool external_Pp = false;
    bool external_Pm = false;
    bool external_M1 = false;
    bool external_M2 = false;
    bool external_M3 = false;
    bool external_M4 = false;
    bool external_gamma1 = false;
    bool external_gamma2 = false;

    //Deformer outputs
    DIdeform::C4vector deformed_loop_momentum;
    my_comp deformation_jacobian;
    double pyLoop[8];
    double jacobian[2];
    
    
    //set factors using option ids (python function)
    int set_factor_int(short int op_id , int v[], int d){
        return 99;
    }
    
    int set_factor_double(short int op_id, double v[],int d){
        switch (op_id)
        {
        case _OP_M1_FACTOR:
            return set_M(v, d, M1_factor,external_M1);
       
        case _OP_M2_FACTOR:
            return set_M(v, d, M2_factor,external_M2);
       
        case _OP_M3_FACTOR:
            return set_M(v, d, M3_factor,external_M3);
       
        case _OP_M4_FACTOR:
            return set_M(v, d, M4_factor,external_M4);
       
        case _OP_GAMMA1_FACTOR:
            return set_gamma(v, d, gamma1, external_gamma1);
       
        case _OP_GAMMA2_FACTOR:
            return set_gamma(v, d, gamma2, external_gamma2);
       
        default:
            return 99;
        };
        return 0;
    }

    //Gives the option to initialize with give P+ and P-
    int set_Pp(double v[], int d)
    {
        //Check dimension
        if (d != 4)
            return 1;

        for (int i = 0; i < 4; i++)
            Pp[i] = v[i];

        //If deformer exists set new P+, otherwise do it later
        if (deformer == NULL)
            external_Pp = true;
        else
            deformer->set_global_var();

        return 0;
    }

    int set_Pm(double v[], int d)
    {
        //Check dimension
        if (d != 4)
            return 1;

        for (int i = 0; i < 4; i++)
            Pm[i] = v[i];

        //If deformer exists set new P-, otherwise do it later
        if (deformer == NULL)
            external_Pm = true;
        else
            deformer->set_global_var();

        return 0;
    }

    int set_M(double v[], int d, double &M, bool &external_M)
    {
        //Check if there is only 1 input
        if (d != 1)
            return 1;
        M = v[0];
        
        //If deformer exists set new M, otherwise do it later
        if (deformer == NULL)
            external_M = true;
        else
            deformer->set_global_var();
        return 0;
    }

    int set_gamma(double v[], int d, double &gamma, bool &external_gamma)
    {
        //Check if there are 1 inputs
        if (d != 1)
            return 1;
        gamma = v[0];

        //If deformer exists set new gamma, otherwise do it later
        if (deformer == NULL)
            external_gamma = true;
        else
            deformer->set_global_var();
        return 0;
    }

    //First one needs to give the Qs
    int append_Q(double v[], int d)
    {
        //Check dimension
        if (d != 4)
            return 1;

        DIdeform::R4vector q;
        for (int i = 0; i < 4; i++)
        {
            q[i] = v[i];
            v[i] = 10.;
        }
        Qs.push_back(q);

        return 0;
    }
    /*
     * This function needs the Qs vectors from the propagators
     * and the P+ and P- vector the the external deformation=
     */
    int init()
    {
        //Check if deformer exists
        if (deformer != NULL)
            return 101;

        //If Qs has not been initialized return 1
        if (Qs.size() < 3)
            return 1;

        //create the deformer
        deformer = new DIdeform::ContourDeform(Qs);

        //overwrite P+- when necessary
        if (external_Pp || external_Pm ||
            external_M1 || external_M2 ||
            external_M3 || external_M4 ||
            external_gamma1 || external_gamma2)
        {
            if (external_Pp)
                deformer->Pp = Pp;
            if (external_Pm)
                deformer->Pm = Pm;
            if (external_M1)
                deformer->M1f = M1_factor;
            if (external_M2)
                deformer->M2f = M2_factor;
            if (external_M3)
                deformer->M3f = M3_factor;
            if (external_M4)
                deformer->M4f = M4_factor;
            if (external_gamma1)
                deformer->gamma1 = gamma1;
            if (external_gamma2)
                deformer->gamma2 = gamma2;
            deformer->set_global_var();
        }

        //std::printf("M1_factor:\t%f   |   ", deformer->M1f);
        //std::printf("M2_factor:\t%f   |   ", deformer->M2f);
        //std::printf("M3_factor:\t%f\n", deformer->M3f);
        //std::printf("gamma1:   \t%f   |   ", deformer->gamma1);
        //std::printf("gamma2:   \t%f   |\n", deformer->gamma2);
        return 0;
    }

    //Create the deformation
    int deform_loop_momentum(double l[], int d)
    {
        //Check if deformer exists
        if (deformer == NULL)
            return 101;

        //Check dimension
        if (d != 4)
            return 1;

        for (int i = 0; i < 4; i++)
        {
            loop_momentum[i] = l[i];
        }
        deformer->loop_momentum(loop_momentum);
        deformer->deform(deformed_loop_momentum,
                         deformation_jacobian);

        return 0;
    }

    //Get the deformation
    double *get_jacobian()
    {
        jacobian[0] = deformation_jacobian.real();
        jacobian[1] = deformation_jacobian.imag();
        return jacobian;
    }

    //Get the deformed loop
    my_real *get_deformed_loop_momentum()
    {
        for (int i = 0; i < 4; i++)
        {
            pyLoop[i] = deformed_loop_momentum[i].real();
            pyLoop[i + 4] = deformed_loop_momentum[i].imag();
        }
        return pyLoop;
    }

    /* 
     *  This function needs to be called when the computation is over
     *  or when one wants to compute the deformation for a new set of
     *  Qs.
     */
    void clear()
    {
        //Free momory allocated to class element
        if (deformer != NULL)
            delete deformer;
        deformer = NULL;

        //Reset variable for new call
        external_Pp = false;
        external_Pm = false;
        external_M1 = false;
        external_M2 = false;
        external_M3 = false;
        external_M4 = false;
        external_gamma1 = false;
        external_gamma2 = false;
        
        //Clear Qs vector
        Qs.clear();
    }
}
