Auto S cMi1L1, cmi1L2, cmi2L2, alarm;
S eulergamma, log4pi, pi;
CF uvid;

#procedure TruncateExpansion()
    argument rat;
        id ep^{`MAXPOLE'+`SELECTEDEPSILONORDER'+1} = 0;
    endargument;
#endprocedure

#procedure IntegrateUV1L()
    if ((count(vxs,1) == 1) && match(vxs(k1?vector_,-k1?)));
* reduce the numerator
        repeat id g(k1?,k1?)*uvprop(k1?,n1?) = uvprop(k1, n1-1) + mUV2^2 * uvprop(k1, n1);
        id uvprop(k1?,n?) = uvprop(n);

* 1-loop IBP
        id uvprop(n1?{<1}) = 0;
        repeat id uvprop(n1?{>1}) = uvprop(-1 + n1)*rat((2 + (4-2*ep) - 2*n1), 2* (-1 + n1)) / mUV2^2;
        id uvprop(1) = uvid(1,1) * rat(1,ep);
    endif;
    id vxs(k1?,-k1?) = 1;
#endprocedure

#procedure IntegrateUV2L()
    
* two vxs -> 2 loops
    if (count(vxs,1) != 2) goto end;
    id vxs(k1?,k2?,k3?)*vxs(k4?,k5?,k6?) = map(kk1,k1,1)*map(kk2,k2,1)*map(kk3,k3,1)*map(kk1,-k1,-1)*map(kk2,-k2,-1)*map(kk3,-k3,-1);

    id map(kk1,k1?,nn1?)*map(kk2,k2?,nn2?)*map(kk3,k3?,nn3?)*uvprop(k1?,n1?)*uvprop(k2?,n2?)*uvprop(k3?,n3?) = 
        uvid(2, 1, n1, n2, n3)*map(kk1,k1,nn1)*map(kk2,k2,nn2)*map(kk3,k3,nn3);

* reduce the numerator
    repeat id g(k1?,k1?)*map(kk1,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = (uvid(2,1,n1-1,n2,n3) + mUV2^2 * uvid(2,1,n1,n2,n3))*map(kk1,k1,n);
    repeat id g(k1?,k1?)*map(kk2,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = (uvid(2,1,n1,n2-1,n3) + mUV2^2 * uvid(2,1,n1,n2,n3))*map(kk2,k1,n);
    repeat id g(k1?,k1?)*map(kk3,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = (uvid(2,1,n1,n2,n3-1) + mUV2^2 * uvid(2,1,n1,n2,n3))*map(kk3,k1,n);

    repeat id g(k1?,k2?)*map(kk1,k1?,nn1?)*map(kk2,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*map(kk1,k1,nn1)*map(kk2,k2,nn2)*(
        uvid(2,1,n1,n2,n3-1) - uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2-1,n3) - mUV2^2 * uvid(2,1,n1,n2,n3)
    );

    repeat id g(k1?,k2?)*map(kk1,k1?,nn1?)*map(kk3,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*map(kk1,k1,nn1)*map(kk3,k2,nn2)*(
        uvid(2,1,n1,n2-1,n3) - uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2,n3-1) - mUV2^2 * uvid(2,1,n1,n2,n3)
    );

    repeat id g(k1?,k2?)*map(kk2,k1?,nn1?)*map(kk3,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*map(kk2,k1,nn1)*map(kk3,k2,nn2)*(
        uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2-1,n3) - uvid(2,1,n1,n2,n3-1) - mUV2^2 * uvid(2,1,n1,n2,n3)
    );

    id map(?a) = 1;

    #define nloop "1"
    #do loop=1,1
* zero sectors
        id uvid(2,1,n1?{<1},n2?{<1},n3?{<1}) = 0;
        id uvid(2,1,n1?{<1},n2?{<1},n3?) = 0;
        id uvid(2,1,n1?{<1},n2?,n3?{<1}) = 0;
        id uvid(2,1,n1?,n2?{<1},n3?{<1}) = 0;

* sector mappings   
        id uvid(2, 1, n1?{>0}, n2?{<1}, n3?{>0}) = uvid(2, 1, n2, n3, n1);
        id uvid(2, 1, n1?{>0}, n2?{>0}, n3?{<1}) = uvid(2, 1, n3, n2, n1);

* actual reduction
        id ifmatch->end uvid(2, 1, n1?{<0}, n2?{>0}, n3?{>1}) = 
            uvid(2, 1, 1 + n1, -1 + n2, n3)*rat(1,1) + 
            uvid(2, 1, 1 + n1, n2, -1 + n3)*rat((-3 + (4-2*ep) - 3*n1),2*(-1 + n3)) + 
            uvid(2, 1, 2 + n1, -1 + n2, -1 + n3)*rat((-1 - n1),2*(-1 + n3)) + 
            uvid(2, 1, 2 + n1, n2, -2 + n3)*rat((1 + n1),2*(-1 + n3)) - 
            uvid(2, 1, 2 + n1, n2, -1 + n3)*mUV2^2*rat(3*(1 + n1),2*(-1 + n3));

        id ifmatch->end uvid(2, 1, n1?{<0}, n2?{>1}, n3?{>0}) = 
            uvid(2, 1, 1 + n1, -1 + n2, n3)*rat((-3 + (4-2*ep) - 3*n1),2*(-1 + n2)) +
            uvid(2, 1, 1 + n1, n2, -1 + n3)*rat(1,1) +
            uvid(2, 1, 2 + n1, -2 + n2, n3)*rat((1 + n1),2*(-1 + n2)) + 
            uvid(2, 1, 2 + n1, -1 + n2, -1 + n3)*rat((-1 - n1),2*(-1 + n2)) - 
            uvid(2, 1, 2 + n1, -1 + n2, n3)*mUV2^2*rat(3*(1 + n1),2*(-1 + n2));

        id ifmatch->end uvid(2, 1, n1?{<1}, n2?{>0}, n3?{>1})  = 
            uvid(2, 1, n1, n2, -1 + n3)*mUV2^-2*rat(2 + (4-2*ep) - n1 - 2*n3,2*(-1 + n3)) + 
            uvid(2, 1, 1 + n1, -1 + n2, -1 + n3)*mUV2^-2*rat(n1,2*(-1 + n3)) - 
            uvid(2, 1, 1 + n1, n2, -2 + n3)*mUV2^2*rat(n1,2*(-1 + n3)) - 
            uvid(2, 1, 1 + n1, n2, -1 + n3)*rat(n1,2*(-1 + n3));

        id ifmatch->end uvid(2, 1, n1?{<1}, n2?{>1}, n3?{>0})  = 
            uvid(2, 1, n1, -1 + n2, n3)*mUV2^-2*rat(2 + (4-2*ep) - n1 - 2*n2,2*(-1 + n2)) - 
            uvid(2, 1, 1 + n1, -2 + n2, n3)*mUV2^-2*rat(n1,2*(-1 + n2)) + 
            uvid(2, 1, 1 + n1, -1 + n2, -1 + n3)*mUV2^-2*rat(n1,2*(-1 + n2)) - 
            uvid(2, 1, 1 + n1, -1 + n2, n3)*rat(n1,2*(-1 + n2));

        id ifmatch->end uvid(2, 1, n1?{<0}, 1, 1)  = 
            uvid(2, 1, 1 + n1, 0, 1)*rat(n1,-2 + (4-2*ep) - n1) - 
            uvid(2, 1, 1 + n1, 1, 0)*rat(n1,-2 + (4-2*ep) - n1) + 
            uvid(2, 1, 1 + n1, 1, 1)*mUV2^2*rat(-3 + (4-2*ep) - 2*n1,-2 + (4-2*ep) - n1) + 
            uvid(2, 1, 1 + n1, 2, 0)*mUV2^2*rat(2,-2 + (4-2*ep) - n1) - 
            uvid(2, 1, 2 + n1, 0, 1)*mUV2^2*rat(1 + n1,2 - (4-2*ep) + n1) - 
            uvid(2, 1, 2 + n1, 1, 0)*mUV2^2*rat(1 + n1,-2 + (4-2*ep) - n1) - 
            uvid(2, 1, 2 + n1, 1, 1)*mUV2^4*rat(3*(1 + n1),-2 + (4-2*ep) - n1);

        id ifmatch->end uvid(2, 1, n1?{>0}, n2?{>0}, n3?{>1})  = 
            uvid(2, 1, -1 + n1, 1 + n2, -1 + n3)*mUV2^-2*rat(n2,3*(-1 + n3)) + 
            uvid(2, 1, n1, n2, -1 + n3)*mUV2^-2*rat((3 + (4-2*ep) - 3*n3),3*(-1 + n3)) - 
            uvid(2, 1, n1, 1 + n2, -2 + n3)*mUV2^-2*rat(n2,3*(-1 + n3)) + 
            uvid(2, 1, 1 + n1, -1 + n2, -1 + n3)*mUV2^-2*rat(n1,3*(-1 + n3)) - 
            uvid(2, 1, 1 + n1, n2, -2 + n3)*mUV2^-2*rat(n1,3*(-1 + n3));          

        id ifmatch->end uvid(2, 1, n1?{>0}, n2?{>1}, n3?{>0}) = 
            uvid(2, 1, -1 + n1, n2, n3)*mUV2^-2*rat(1,3) + 
            uvid(2, 1, n1, -1 + n2, n3)*mUV2^-2*rat(3 + (4-2*ep) - 3*n2,3*(-1 + n2)) - 
            uvid(2, 1, n1, n2, -1 + n3)*mUV2^-2*rat(1,3) - 
            uvid(2, 1, 1 + n1, -2 + n2, n3)*mUV2^-2*rat(2*n1,3*(-1 + n2)) + 
            uvid(2, 1, 1 + n1, -1 + n2, -1 + n3)*mUV2^-2*rat(2*n1,3*(-1 + n2)); 
        
        id ifmatch->end uvid(2, 1, n1?{>1}, n2?{>0}, n3?{>0}) = 
            uvid(2, 1, -2 + n1, 1 + n2, n3)*mUV2^-2*rat(-2*n2,3*(-1 + n1)) + 
            uvid(2, 1, -1 + n1, n2, n3)*mUV2^-2*rat((3 + (4-2*ep) - 3*n1),3*(-1 + n1)) + 
            uvid(2, 1, -1 + n1, 1 + n2, -1 + n3)*mUV2^-2*rat(2*n2,3*(-1 + n1)) + 
            uvid(2, 1, n1, -1 + n2, n3)*mUV2^-2*rat(1,3) - 
            uvid(2, 1, n1, n2, -1 + n3)*mUV2^-2*rat(1,3);

        id ifmatch->end uvid(2, 1, n1?{>0}, n2?{>0}, n3?{>1})  = 
            uvid(2, 1, n1, 1 + n2, -1 + n3)*rat(n2,-1 + n3) + 
            uvid(2, 1, 1 + n1, -1 + n2, n3)*rat(1,1) - 
            uvid(2, 1, 1 + n1, 1 + n2, -2 + n3)*rat(n2,-1 + n3) + 
            uvid(2, 1, 2 + n1, -1 + n2, -1 + n3)*rat((-1 - n1),-1 + n3) + 
            uvid(2, 1, 2 + n1, n2, -2 + n3)*rat((1 + n1),-1 + n3); 

        label end;

        id uvid(2,1,0,1,1) = tp2P011;
        id uvid(2,1,1,1,1) = tp2P111;

        if (match(uvid(2, 1, ?a)));
            redefine loop "0";
        endif;
    
        .sort:reduction-`nloop++';
    #enddo
* MIs are uvid(2,1,0,2,2) and uvid(2,1,2,2,1)
    id tp2P011 = uvid(2,1)*rat(1,ep^2)*mUV2^4*rat(4,(4-2*ep)^2 - 4*(4-2*ep) + 4);
    id tp2P111 = mUV2^2*rat(9, (4-2*ep)^2 - 5*(4-2*ep) + 6)*(uvid(2,2) - uvid(2,1)*rat(1,ep^2)*rat(-1,3));
#endprocedure

#procedure IntegrateUV()
* match the topology, working down in loop count
    #do i = {2L,1L}
        #call IntegrateUV`i'()
    #enddo

    if (count(vxs,1,uvprop,1));
        Print "Unsubstituted UV graph: %t";
        exit "Critical error";
    endif;
#endprocedure

#procedure Masters()
    id uvid(1,1) = rat(ep,1) * mUV2^2 * -(rat(1,ep) + 1 + (1 + 1/12*pi^2)*rat(ep,1) + (1 + 1/12*pi^2 - 1/3*z3)*rat(ep^2,1));
    id uvid(2,1) = rat(ep^2,1) * (rat(1,ep^2) + pi^2/6 - (-2*z3)/3 * rat(ep, 1));
    id uvid(2,2) = 2 * ( -8453430797319808/21639331090712481 + 28375860251062315/42167665933747557 * rat(ep, 1));
#endprocedure

#procedure SubstituteMasters()
    Multiply replace_(D, 4 - 2 * ep);
    id ep^n1? = rat(ep^n1,1);

    B+ uvid;
    .sort:masters-1;
    PolyRatFun rat(expand,ep,{`MAXPOLE'+`SELECTEDEPSILONORDER'});
    Keep brackets;

* Normalize with alphaLoop convention
    id uvid(n?,n1?) = uvid(n,n1) * (i_*16*pi^2)^-n;

* add log(mUV^2/mu^2)-dependence
    id uvid(n?,n1?) = uvid(n,n1) * (1 - logmUVmu * rat(ep, 1) + 1/2 * logmUVmu^2 * rat(ep^2, 1) - 1/6 * logmUVmu^3 * rat(ep^3, 1))^n;

    .sort
    #call TruncateExpansion()

    B+ uvid;
    .sort:normalization;
    Keep brackets;

    #call Masters()
    .sort:substitution;

    #call TruncateExpansion()
    .sort:truncation;

* fill in constants, they have not been substituted before to have exact cancellation
*    id pi = 30246273033735921/9627687726852338;
   id z3 = 15408859284534249/12818743641862727;

#endprocedure