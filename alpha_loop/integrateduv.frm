Auto S cMi1L1, cmi1L2, cmi2L2, alarm;
S eulergamma, log4pi, pi;
CF uvid;

* TODO: set the expansion depth to n+selected, where n is the power of the deepest pole of the entire diagram
#define MAXPOLE "3"

#procedure TruncateExpansion()
    argument rat;
        id ep^{`MAXPOLE'+`SELECTEDEPSILONORDER'+1} = 0;
    endargument;
#endprocedure

#procedure IntegrateUV1L()
    if (match(vxs(k1?vector_,-k1?)));
* reduce the numerator
        repeat id g(k1?,k1?)*uvprop(k1?,n1?) = uvprop(k1, n1-1) + mUV2^2 * uvprop(k1, n1);
        id uvprop(k1?,n?) = uvprop(n);

* 1-loop IBP
        id uvprop(n1?{<1}) = 0;
        repeat id uvprop(n1?{>1}) = uvprop(-1 + n1)*rat((2 + (4-2*ep) - 2*n1), 2* (-1 + n1)) / mUV2^2;
        id uvprop(1) = uvid(1,1) * rat(1, ep) * mUV2^2 * rat(2, (4-2*ep) - 2);
    endif;
    id vxs(k1?,-k1?) = 1;
#endprocedure

#procedure IntegrateUV2L()
    
* two vxs -> 2 loops
    if (count(vxs,1) != 2) goto end;
    id vxs(k1?,k2?,k3?)*vxs(k4?,k5?,k6?) = map(kk1,k1,1)*map(kk2,k2,1)*map(kk3,k3,1)*map(kk1,-k1,-1)*map(kk2,-k2,-1)*map(kk3,-k3,-1);

*    id ifnomatch->test vxs(k1?,k2?,k3?)*vxs(-k1?,-k2?,-k3?) = map(kk1,k1,1)*map(kk2,k2,1)*map(kk3,k3,1)*map(kk1,-k1,-1)*map(kk2,-k2,-1)*map(kk3,-k3,-1);
    id map(kk1,k1?,nn1?)*map(kk2,k2?,nn2?)*map(kk3,k3?,nn3?)*uvprop(k1?,n1?)*uvprop(k2?,n2?)*uvprop(k3?,n3?) = 
        uvid(2, 1, n1, n2, n3)*map(kk1,k1,nn1)*map(kk2,k2,nn2)*map(kk3,k3,nn3);

* reduce the numerator
    repeat id g(k1?,k1?)*map(kk1,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = (uvid(2,1,n1-1,n2,n3) + mUV2^2 * uvid(2,1,n1,n2,n3))*map(kk1,k1,n);
    repeat id g(k1?,k1?)*map(kk2,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = (uvid(2,1,n1,n2-1,n3) + mUV2^2 * uvid(2,1,n1,n2,n3))*map(kk2,k1,n);
    repeat id g(k1?,k1?)*map(kk3,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = (uvid(2,1,n1,n2,n3-1) + mUV2^2 * uvid(2,1,n1,n2,n3))*map(kk3,k1,n);

    repeat id g(k1?,k2?)*map(kk1,k1?,nn1?)*map(kk2,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*map(kk1,k1,nn1)*map(kk2,k2,nn2)*(
        uvid(2,1,n1,n2,n3-1) - uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2-1,n3) + mUV2^2 * uvid(2,1,n1,n2,n3)
    );

    repeat id g(k1?,k2?)*map(kk1,k1?,nn1?)*map(kk3,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*map(kk1,k1,nn1)*map(kk3,k2,nn2)*(
        uvid(2,1,n1,n2-1,n3) - uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2,n3-1) + mUV2^2 * uvid(2,1,n1,n2,n3)
    );

    repeat id g(k1?,k2?)*map(kk2,k1?,nn1?)*map(kk3,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*map(kk2,k1,nn1)*map(kk3,k2,nn2)*(
        uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2-1,n3) - uvid(2,1,n1,n2,n3-1) + mUV2^2 * uvid(2,1,n1,n2,n3)
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

#procedure Mastermi1L1()
    id uvid(1,1) = rat(ep,1) * (
        cMi1L1EpsM1logmUV0*logmUV^0*rat(1,ep^1)
        +cMi1L1Eps0logmUV0*logmUV^0*rat(ep^0,1)
        +cMi1L1Eps0logmUV1*logmUV^1*rat(ep^0,1)
        +cMi1L1Eps1logmUV0*logmUV^0*rat(ep^1,1)
        +cMi1L1Eps1logmUV1*logmUV^1*rat(ep^1,1)
        +cMi1L1Eps1logmUV2*logmUV^2*rat(ep^1,1)
        +cMi1L1Eps2logmUV0*logmUV^0*rat(ep^2,1)
        +cMi1L1Eps2logmUV1*logmUV^1*rat(ep^2,1)
        +cMi1L1Eps2logmUV2*logmUV^2*rat(ep^2,1)
        +cMi1L1Eps2logmUV3*logmUV^3*rat(ep^2,1)
        +cMi1L1Eps3logmUV0*logmUV^0*rat(ep^3,1)
        +cMi1L1Eps3logmUV1*logmUV^1*rat(ep^3,1)
        +cMi1L1Eps3logmUV2*logmUV^2*rat(ep^3,1)
        +cMi1L1Eps3logmUV3*logmUV^3*rat(ep^3,1)
        +cMi1L1Eps3logmUV4*logmUV^4*rat(ep^3,1)
        +alarmMi1L1*rat(ep^4,1)
    );

    id cMi1L1EpsM1logmUV0= i_/(16*pi^2);
    id cMi1L1Eps0logmUV0 = i_*(-eulergamma + log4pi)/(16*pi^2);
    id cMi1L1Eps0logmUV1 = -i_/(16*pi^2);
    id cMi1L1Eps1logmUV0 = i_*(6*eulergamma^2 + pi^2 - 12*eulergamma*log4pi + 6*log4pi^2)/(192*pi^2);
    id cMi1L1Eps1logmUV1 = i_*(eulergamma - log4pi)/(16*pi^2);
    id cMi1L1Eps1logmUV2 = i_/(32*pi^2);
#endprocedure

#procedure Mastermi1L2()
    id uvid(2,1) = rat(ep^2,1) * (
	 cMi1L2EpsM2logmUV0*logmUV^0*rat(1,ep^2)
	 +cMi1L2EpsM1logmUV0*logmUV^0*rat(1,ep^1)
	 +cMi1L2EpsM1logmUV1*logmUV^1*rat(1,ep^1)
	 +cMi1L2Eps0logmUV0*logmUV^0*rat(ep^0,1)
	 +cMi1L2Eps0logmUV1*logmUV^1*rat(ep^0,1)
	 +cMi1L2Eps0logmUV2*logmUV^2*rat(ep^0,1)
	 +cMi1L2Eps1logmUV0*logmUV^0*rat(ep^1,1)
	 +cMi1L2Eps1logmUV1*logmUV^1*rat(ep^1,1)
	 +cMi1L2Eps1logmUV2*logmUV^2*rat(ep^1,1)
	 +cMi1L2Eps1logmUV3*logmUV^3*rat(ep^1,1)
	 +cMi1L2Eps2logmUV0*logmUV^0*rat(ep^2,1)
	 +cMi1L2Eps2logmUV1*logmUV^1*rat(ep^2,1)
	 +cMi1L2Eps2logmUV2*logmUV^2*rat(ep^2,1)
	 +cMi1L2Eps2logmUV3*logmUV^3*rat(ep^2,1)
	 +cMi1L2Eps2logmUV4*logmUV^4*rat(ep^2,1)
	 +alarmMi1L2*rat(ep^3,1)
	);

    id cMi1L2EpsM2logmUV0=-1/(256*pi^4);
    id cMi1L2EpsM1logmUV0=-(-eulergamma + log4pi)/(128*pi^4);
    id cMi1L2EpsM1logmUV1=1/(128*pi^4);
    id cMi1L2Eps0logmUV0=-(12*eulergamma^2 + pi^2 - 24*eulergamma*log4pi + 12*log4pi^2)/(1536*pi^4);
    id cMi1L2Eps0logmUV1=-(eulergamma - log4pi)/(64*pi^4);
    id cMi1L2Eps0logmUV2=-1/(128*pi^4);
#endprocedure

#procedure Mastermi2L2()
    id uvid(2,2) = 2 * (
	 cMi2L2Eps0logmUV0*logmUV^0*rat(ep^0,1)
	 +cMi2L2Eps1logmUV0*logmUV^0*rat(ep^1,1)
	 +cMi2L2Eps1logmUV1*logmUV^1*rat(ep^1,1)
	 +cMi2L2Eps2logmUV0*logmUV^0*rat(ep^2,1)
	 +cMi2L2Eps2logmUV1*logmUV^1*rat(ep^2,1)
	 +cMi2L2Eps2logmUV2*logmUV^2*rat(ep^2,1)
	 +alarmMi2L2*rat(ep^3,1)
	);

* TODO: use analytic result instead
    id cMi2L2Eps0logmUV0 = rat(335456055722,21413414437668207);
    id cMi2L2Eps1logmUV0 = rat(821288315359,23993180133650287);
    id cMi2L2Eps1logmUV1 = rat(-703172337945,22443059883227614);
    id cMi2L2Eps2logmUV0 = rat(1862483703443,19882032248356283);
    id cMi2L2Eps2logmUV1 = rat(-2629675244357,38411767616149931);
    id cMi2L2Eps2logmUV2 = rat(703172337945,22443059883227614);
#endprocedure

#procedure SubstituteMasters()
    Multiply replace_(D, 4 - 2 * ep);
    id ep^n1? = rat(ep^n1,1);

    B+ uvid;
    .sort:masters-1;
    PolyRatFun rat(expand,ep,{`MAXPOLE'+`SELECTEDEPSILONORDER'});
    Keep brackets;

* normalize with 1/(4 pi e^-gamma)^ep
    id uvid(n?,n1?) = uvid(n,n1) * (1
        +(eulergamma - log4pi)*rat(ep,1) + 1/2*(eulergamma^2 - 2*eulergamma*log4pi + log4pi^2)*rat(ep^2, 1)
        +1/6*(eulergamma^3 - 3*eulergamma^2*log4pi + 3*log4pi*log4pi^2 - log4pi^3)*rat(ep^3, 1))^n;

* add mu^2-dependence
    id uvid(n?,n1?) = uvid(n,n1) * (1
        + logmu * rat(ep, 1) + 1/2 * logmu^2 * rat(ep^2, 1) + 1/6 * logmu^3 * rat(ep^3, 1))^n;

    .sort
    #call TruncateExpansion()

    B+ uvid;
    .sort:normalization;
    Keep brackets;

    #call Mastermi1L1()
    #call Mastermi1L2()
    #call Mastermi2L2()
    .sort:substitution;

    #call TruncateExpansion()
    .sort:truncation;
#endprocedure