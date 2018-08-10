#include "DCD.h"

extern "C"
{
    //Deformer
    DIdeform::ContourDeform *deformer;

    //Deformer arguments
    std::vector<DIdeform::R4vector> Qs;
    DIdeform::R4vector Pp, Pm, loop_momentum;
    bool external_Pp = false;
    bool external_Pm = false;

    //Deformer outputs
    DIdeform::C4vector deformed_loop_momentum;
    my_comp deformation_jacobian;
    double pyLoop[8];
    double jacobian[2];
    
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


    //Gives the option to initialize with give P+ and P-
    int set_Pp(double v[], int d)
    {
        //Check dimension
        if (d != 4)
            return 1;

        for (int i = 0; i < 4; i++)
            Pp[i] = v[i];
        external_Pp = true;

        return 0;
    }

    int set_Pm(double v[], int d)
    {
        //Check dimension
        if (d != 4)
            return 1;

        for (int i = 0; i < 4; i++)
            Pm[i] = v[i];
        external_Pm = true;

        return 0;
    }

    int test()
    {
        return Qs.size();
    }

    /*
     * This function needs the Qs vectors from the propagators
     * and the P+ and P- vector the the external deformation=
     */
    int init()
    {
        //If Qs has not been initialized return 1
        if (Qs.size() < 3)
            return 1;
        //create the deformer
        deformer = new DIdeform::ContourDeform(Qs);

        //overwrite P+- when necessary
        if (external_Pp && external_Pm)
        {
            deformer->Pp = Pp;
            deformer->Pm = Pm;
            deformer->set_global_var();
        }
        else
            return 2;

        return 0;
    }

    //Create the deformation
    int deform_loop_momentum(double l[], int d)
    {
        //Check dimension
        if (d != 4)
            return 1;

        for (int i = 0; i < 4; i++){
            loop_momentum[i] = l[i];
        }
        deformer->loop_momentum(loop_momentum);
        deformer->deform(deformed_loop_momentum,
                         deformation_jacobian);

        return 0;
    }

    //Get the deformation
    double* get_jacobian()
    {
        jacobian[0] = deformation_jacobian.real();
        jacobian[1] = deformation_jacobian.imag();
        return jacobian;
    }

    //Get the deformed loop
    my_real* get_deformed_loop_momentum()
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
        delete deformer;
        //Reset variable for new call
        external_Pp = false;
        external_Pm = false;
        //Clear Qs vector
        std::vector<DIdeform::R4vector>().swap(Qs);
    }
}
