*---- defnition of eps-expanded mi1L1: lMUV=log(mUV^2)
	id mi1L1 = (
	 cMi1L1EpsM1lMUV0*rat(lMUV^0,1)*rat(1,eps^1) 
	 +cMi1L1Eps0lMUV0*rat(lMUV^0,1)*rat(eps^0,1) 
	 +cMi1L1Eps0lMUV1*rat(lMUV^1,1)*rat(eps^0,1) 
	 +cMi1L1Eps1lMUV0*rat(lMUV^0,1)*rat(eps^1,1) 
	 +cMi1L1Eps1lMUV1*rat(lMUV^1,1)*rat(eps^1,1) 
	 +cMi1L1Eps1lMUV2*rat(lMUV^2,1)*rat(eps^1,1) 
	 +cMi1L1Eps2lMUV0*rat(lMUV^0,1)*rat(eps^2,1) 
	 +cMi1L1Eps2lMUV1*rat(lMUV^1,1)*rat(eps^2,1) 
	 +cMi1L1Eps2lMUV2*rat(lMUV^2,1)*rat(eps^2,1) 
	 +cMi1L1Eps2lMUV3*rat(lMUV^3,1)*rat(eps^2,1) 
	 +cMi1L1Eps3lMUV0*rat(lMUV^0,1)*rat(eps^3,1) 
	 +cMi1L1Eps3lMUV1*rat(lMUV^1,1)*rat(eps^3,1) 
	 +cMi1L1Eps3lMUV2*rat(lMUV^2,1)*rat(eps^3,1) 
	 +cMi1L1Eps3lMUV3*rat(lMUV^3,1)*rat(eps^3,1) 
	 +cMi1L1Eps3lMUV4*rat(lMUV^4,1)*rat(eps^3,1) 
	 +cMi1L1Eps4lMUV0*rat(lMUV^0,1)*rat(eps^4,1) 
	 +cMi1L1Eps4lMUV1*rat(lMUV^1,1)*rat(eps^4,1) 
	 +cMi1L1Eps4lMUV2*rat(lMUV^2,1)*rat(eps^4,1) 
	 +cMi1L1Eps4lMUV3*rat(lMUV^3,1)*rat(eps^4,1) 
	 +cMi1L1Eps4lMUV4*rat(lMUV^4,1)*rat(eps^4,1) 
	 +cMi1L1Eps4lMUV5*rat(lMUV^5,1)*rat(eps^4,1) 
	 +alarmMi1L1Eps5*rat(eps^5,1)
	);
*
*
*---- BEGIN: defintion of the numeric coefficients for mi1L1: order 100 digits accurat 
*
	 id cMi1L1EpsM1lMUV0 = rat(132049792606502*i_,20852467428353107);
	 id cMi1L1Eps0lMUV0 = rat(930733834984243*i_,75225176822851668);
	 id cMi1L1Eps0lMUV1 = rat(-132049792606502*i_,20852467428353107);
	 id cMi1L1Eps1lMUV0 = rat(238365146153033*i_,13782143439795685);
	 id cMi1L1Eps1lMUV1 = rat(-1057101110939383*i_,85438623805262522);
	 id cMi1L1Eps1lMUV2 = rat(66024896303251*i_,20852467428353107);
	 id cMi1L1Eps2lMUV0 = rat(267137567007244*i_,17222977639926303);
	 id cMi1L1Eps2lMUV1 = rat(-238365146153033*i_,13782143439795685);
	 id cMi1L1Eps2lMUV2 = rat(31591818988785*i_,5106723491205427);
	 id cMi1L1Eps2lMUV3 = rat(-26474513963629*i_,25084126035084908);
	 id cMi1L1Eps3lMUV0 = rat(99723921272160*i_,7862278376418437);
	 id cMi1L1Eps3lMUV1 = rat(-267137567007244*i_,17222977639926303);
	 id cMi1L1Eps3lMUV2 = rat(238365146153033*i_,27564286879591370);
	 id cMi1L1Eps3lMUV3 = rat(-10530606329595*i_,5106723491205427);
	 id cMi1L1Eps3lMUV4 = rat(3349661396909*i_,12694975820195403);
	 id cMi1L1Eps4lMUV0 = rat(168637982940181*i_,23208841390718796);
	 id cMi1L1Eps4lMUV1 = rat(-99723921272160*i_,7862278376418437);
	 id cMi1L1Eps4lMUV2 = rat(133568783503622*i_,17222977639926303);
	 id cMi1L1Eps4lMUV3 = rat(-147392469819533*i_,51132915875907201);
	 id cMi1L1Eps4lMUV4 = rat(10530606329595*i_,20426893964821708);
	 id cMi1L1Eps4lMUV5 = rat(-1275309116435*i_,24166649219167214);
*
*---- END: defintion of the coefficients for mi1L1
*
* 
*---- BEGIN: defintion of the analytc coefficients for mi1L1
*
*	cMi1L1EpsM1lMUV0=imI/(16*Pi^2)
*
*
*	cMi1L1Eps0lMUV0=(imI*(-EulerGamma + Log[4*Pi]))/(16*Pi^2)
*
*
*	cMi1L1Eps0lMUV1=-imI/(16*Pi^2)
*
*
*	cMi1L1Eps1lMUV0=(imI*(6*EulerGamma^2 + Pi^2 - 12*EulerGamma*Log[4*Pi] + 6*Log[4*Pi]^2))/(192*Pi^2)
*
*
*	cMi1L1Eps1lMUV1=(imI*(EulerGamma + Log[1/(4*Pi)]))/(16*Pi^2)
*
*
*	cMi1L1Eps1lMUV2=imI/(32*Pi^2)
*
*
*	cMi1L1Eps2lMUV0=(imI*(-2*EulerGamma^3 + 6*EulerGamma^2*Log[4*Pi] + Pi^2*Log[4*Pi] - EulerGamma*(Pi^2 + 6*Log[4*Pi]^2) + 2*(8*Log[2]^3 + 12*Log[2]^2*Log[Pi] + Log[64]*Log[Pi]^2 + Log[Pi]^3 - 2*Zeta[3])
*	))/(192*Pi^2)
*
*
*	cMi1L1Eps2lMUV1=-(imI*(6*EulerGamma^2 + Pi^2 - 12*EulerGamma*Log[4*Pi] + 6*Log[4*Pi]^2))/(192*Pi^2)
*
*
*	cMi1L1Eps2lMUV2=(imI*(-EulerGamma + Log[4*Pi]))/(32*Pi^2)
*
*
*	cMi1L1Eps2lMUV3=-imI/(96*Pi^2)
*
*
*	cMi1L1Eps3lMUV0=(imI*(20*EulerGamma^4 + 3*Pi^4 - 80*EulerGamma^3*Log[4*Pi] + 20*Pi^2*Log[4*Pi]^2 + 20*EulerGamma^2*(Pi^2 + 6*Log[4*Pi]^2) - 40*EulerGamma*(Pi^2*Log[4*Pi] + 2*(8*Log[2]^3 + 12*Log[2]^2*
*	Log[Pi] + Log[64]*Log[Pi]^2 + Log[Pi]^3 - 2*Zeta[3])) + 20*Log[4*Pi]*(8*Log[2]^3 + 12*Log[2]^2*Log[Pi] + Log[64]*Log[Pi]^2 + Log[Pi]^3 - 8*Zeta[3])))/(7680*Pi^2)
*
*
*	cMi1L1Eps3lMUV1=(imI*(2*EulerGamma^3 - 6*EulerGamma^2*Log[4*Pi] - Pi^2*Log[4*Pi] + EulerGamma*(Pi^2 + 6*Log[4*Pi]^2) - 2*(8*Log[2]^3 + 12*Log[2]^2*Log[Pi] + Log[64]*Log[Pi]^2 + Log[Pi]^3 - 2*Zeta[3]))
*	)/(192*Pi^2)
*
*
*	cMi1L1Eps3lMUV2=(imI*(6*EulerGamma^2 + Pi^2 - 12*EulerGamma*Log[4*Pi] + 6*Log[4*Pi]^2))/(384*Pi^2)
*
*
*	cMi1L1Eps3lMUV3=(imI*(EulerGamma + Log[1/(4*Pi)]))/(96*Pi^2)
*
*
*	cMi1L1Eps3lMUV4=imI/(384*Pi^2)
*
*
*	cMi1L1Eps4lMUV0=(imI*(-12*EulerGamma^5 + 60*EulerGamma^4*Log[4*Pi] + 9*Pi^4*Log[4*Pi] - 20*EulerGamma^3*(Pi^2 + 6*Log[4*Pi]^2) - 3*EulerGamma*(3*Pi^4 + 20*Pi^2*Log[4*Pi]^2 + 20*Log[4*Pi]*(8*Log[2]^3 +
*	 12*Log[2]^2*Log[Pi] + Log[64]*Log[Pi]^2 + Log[Pi]^3 - 8*Zeta[3])) + 60*EulerGamma^2*(Pi^2*Log[4*Pi] + 2*(8*Log[2]^3 + 12*Log[2]^2*Log[Pi] + Log[64]*Log[Pi]^2 + Log[Pi]^3 - 2*Zeta[3])) + 20*Pi^2*(8*Lo
*	g[2]^3 + 12*Log[2]^2*Log[Pi] + Log[64]*Log[Pi]^2 + Log[Pi]^3 - 2*Zeta[3]) + 12*(32*Log[2]^5 + 80*Log[2]^4*Log[Pi] + 80*Log[2]^3*Log[Pi]^2 + Log[Pi]^5 + 40*Log[2]^2*(Log[Pi]^3 - 2*Zeta[3]) - 20*Log[Pi]
*	^2*Zeta[3] + 10*Log[2]*(Log[Pi]^4 - 8*Log[Pi]*Zeta[3]) - 24*Zeta[5])))/(23040*Pi^2)
*
*
*	cMi1L1Eps4lMUV1=-(imI*(20*EulerGamma^4 + 3*Pi^4 - 80*EulerGamma^3*Log[4*Pi] + 20*Pi^2*Log[4*Pi]^2 + 20*EulerGamma^2*(Pi^2 + 6*Log[4*Pi]^2) - 40*EulerGamma*(Pi^2*Log[4*Pi] + 2*(8*Log[2]^3 + 12*Log[2]^2
*	*Log[Pi] + Log[64]*Log[Pi]^2 + Log[Pi]^3 - 2*Zeta[3])) + 20*Log[4*Pi]*(8*Log[2]^3 + 12*Log[2]^2*Log[Pi] + Log[64]*Log[Pi]^2 + Log[Pi]^3 - 8*Zeta[3])))/(7680*Pi^2)
*
*
*	cMi1L1Eps4lMUV2=(imI*(-2*EulerGamma^3 + 6*EulerGamma^2*Log[4*Pi] + Pi^2*Log[4*Pi] - EulerGamma*(Pi^2 + 6*Log[4*Pi]^2) + 2*(8*Log[2]^3 + 12*Log[2]^2*Log[Pi] + Log[64]*Log[Pi]^2 + Log[Pi]^3 - 2*Zeta[3])
*	))/(384*Pi^2)
*
*
*	cMi1L1Eps4lMUV3=-(imI*(6*EulerGamma^2 + Pi^2 - 12*EulerGamma*Log[4*Pi] + 6*Log[4*Pi]^2))/(1152*Pi^2)
*
*
*	cMi1L1Eps4lMUV4=(imI*(-EulerGamma + Log[4*Pi]))/(384*Pi^2)
*
*
*	cMi1L1Eps4lMUV5=-imI/(1920*Pi^2)
*
*
*---- END: defintion of the analytc coefficients for mi1L1
