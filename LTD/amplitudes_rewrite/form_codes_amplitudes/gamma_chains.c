%(numbertype)s compute_chain%(mode)s(int chain_length, %(numbertype)s *vecs[])
{
    // vecs is pointer to complex momenter
    %(numbertype)s res = 0.0;
    %(numbertype)s vbar[4];
    %(numbertype)s vbartmp[4]; 
    vbar[0]= vecs[0][0];
    vbar[1]= vecs[0][1];
    vbar[2]= vecs[0][2];
    vbar[3]= vecs[0][3];    

    int i;

    // multiplication from right with slashed momentum
    for (i = 1; i < chain_length - 1; i++)
    {
        // there must be a better way
        vbartmp[0] = vbar[3] * (vecs[i][1] + I * vecs[i][2]) + vbar[2] * (vecs[i][0] + vecs[i][3]);
        vbartmp[1] = vbar[2] * (vecs[i][1] - I * vecs[i][2]) + vbar[3] * (vecs[i][0] - vecs[i][3]);
        vbartmp[2] = vbar[1] * (-vecs[i][1] - I * vecs[i][2]) + vbar[0] * (vecs[i][0] - vecs[i][3]);
        vbartmp[3] = vbar[0] * (-vecs[i][1] + I * vecs[i][2]) + vbar[1] * (vecs[i][0] + vecs[i][3]);

        vbar[0] = vbartmp[0];
        vbar[1] = vbartmp[1];
        vbar[2] = vbartmp[2];
        vbar[3] = vbartmp[3];
    }
    // finalproduct
    res = vbar[0] * vecs[chain_length - 1][0] + vbar[1] * vecs[chain_length - 1][1] + vbar[2] * vecs[chain_length - 1][2] + vbar[3] * vecs[chain_length - 1][3];
    return res;
}