Auto S cMi1L1;
S mapu;
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
        repeat id k1?.k1?*uvprop(k1?,n1?) = uvprop(k1, n1-1) + mUV^2 * uvprop(k1, n1);
        id uvprop(k1?,n?) = uvprop(n);

* 1-loop IBP
        id uvprop(n1?{<1}) = 0;
        repeat id uvprop(n1?{>1}) = uvprop(-1 + n1)*rat((2 + (4-2*ep) - 2*n1), 2* (-1 + n1)) / mUV^2;
        id uvprop(1) = uvid(1,1) * rat(1, ep) * mUV^2 * rat(2, (4-2*ep) - 2);
    endif;
    id vxs(k1?,-k1?) = 1;
#endprocedure

#procedure IntegrateUV2L()
    id ifnomatch->no2L vxs(k1?,k2?,k3?)*vxs(-k1?,-k2?,-k3?) = map(kk1,k1,1)*map(kk2,k2,1)*map(kk3,k3,1)*map(kk1,-k1,-1)*map(kk2,-k2,-1)*map(kk3,-k3,-1);
    id map(kk1,k1?)*map(kk2,k2?)*map(kk3,k3?)*uvprop(k1?,n1?)*uvprop(k2?,n2?)*uvprop(k3?,n3?) = 
        uvid(2, 1, n1, n2, n3)*map(kk1,k1)*map(kk2,k2)*map(kk3,k3);

* reduce the numerator
    repeat id k1?.k1?*map(kk1,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = uvid(2,1,n1-1,n2,n3) + mUV^2 * uvid(2,1,n1,n2,n3);
    repeat id k1?.k1?*map(kk2,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = uvid(2,1,n1,n2-1,n3) + mUV^2 * uvid(2,1,n1,n2,n3);
    repeat id k1?.k1?*map(kk3,k1?,n?)*uvid(2,1,n1?,n2?,n3?) = uvid(2,1,n1,n2,n3-1) + mUV^2 * uvid(2,1,n1,n2,n3);

    repeat id k1?.k2?*map(kk1,k1?,nn1?)*map(kk2,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*(
        uvid(2,1,n1,n2,n3-1) - uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2-1,n3) - mUV^2 * uvid(2,1,n1,n2,n3)
    );

    repeat id k1?.k2?*map(kk1,k1?,nn1?)*map(kk3,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*(
        uvid(2,1,n1,n2-1,n3) - uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2,n3-1) - mUV^2 * uvid(2,1,n1,n2,n3)
    );

    repeat id k1?.k2?*map(kk2,k1?,nn1?)*map(kk3,k2?,nn2?)*uvid(2,1,n1?,n2?,n3?) = 1/2*nn1*nn2*(
        uvid(2,1,n1-1,n2,n3) - uvid(2,1,n1,n2-1,n3) - uvid(2,1,n1,n2,n3-1) - mUV^2 * uvid(2,1,n1,n2,n3)
    );

    id map(?a) = 1;

*   TODO: 2-loop IBP on uvid(2,1,n1,n2,n3)
    exit "Critical error: 2-loop UV subgraph not implemented";

    label no2L;
#endprocedure

#procedure IntegrateUV()
* match the topology, working down in loop count
    #do i = {2L,1L}
        #call IntegrateUV`i'()
    #enddo

    if (count(vxs,1,uvprop,1));
        Print "Unsubstited UV graph: %t";
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
*    id cMi1L1Eps2logmUV0 = i_*rat(267137567007244,17222977639926303);
*    id cMi1L1Eps2logmUV1 = i_*rat(-238365146153033,13782143439795685);
*    id cMi1L1Eps2logmUV2 = i_*rat(31591818988785,5106723491205427);
*    id cMi1L1Eps2logmUV3 = i_*rat(-26474513963629,25084126035084908);
*    id cMi1L1Eps3logmUV0 = i_*rat(99723921272160,7862278376418437);
*    id cMi1L1Eps3logmUV1 = i_*rat(-267137567007244,17222977639926303);
*    id cMi1L1Eps3logmUV2 = i_*rat(238365146153033,27564286879591370);
*    id cMi1L1Eps3logmUV3 = i_*rat(-10530606329595,5106723491205427);
*    id cMi1L1Eps3logmUV4 = i_*rat(3349661396909,12694975820195403);
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

    .sort
* at this stage all poles should be substituted
    argument rat;
        id ep^{`SELECTEDEPSILONORDER'+1} = 0;
    endargument;
#endprocedure