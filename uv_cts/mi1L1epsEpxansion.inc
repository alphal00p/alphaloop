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
	 id cMi1L1EpsM1lMUV0 = rat(3753287685595852886*i_,592695434564981160875);
	 id cMi1L1EpsM0lMUV0 = rat(762137634968669473*i_,61598639910666565099);
	 id cMi1L1EpsM0lMUV1 = rat(-3753287685595852886*i_,592695434564981160875);
	 id cMi1L1Eps1lMUV0 = rat(923351794526504479*i_,53387699850157443397);
	 id cMi1L1Eps1lMUV1 = rat(-762137634968669473*i_,61598639910666565099);
	 id cMi1L1Eps1lMUV2 = rat(576014374954367863*i_,181921088324491682682);
	 id cMi1L1Eps2lMUV0 = rat(1362742903630679197*i_,87859191131897149395);
	 id cMi1L1Eps2lMUV1 = rat(-923351794526504479*i_,53387699850157443397);
	 id cMi1L1Eps2lMUV2 = rat(2061011697536681054*i_,333156405308478518067);
	 id cMi1L1Eps2lMUV3 = rat(-1158158248837043577*i_,1097334119988482864726);
	 id cMi1L1Eps3lMUV0 = rat(673070460489262038*i_,53065174932988313567);
	 id cMi1L1Eps3lMUV1 = rat(-1362742903630679197*i_,87859191131897149395);
	 id cMi1L1Eps3lMUV2 = rat(890366431681052425*i_,102961008129343655431);
	 id cMi1L1Eps3lMUV3 = rat(-941049777501783509*i_,456353685129811648265);
	 id cMi1L1Eps3lMUV4 = rat(108385788986963215*i_,410774346240489478193);
	 id cMi1L1Eps4lMUV0 = rat(956623842732035481*i_,131655577524449151731);
	 id cMi1L1Eps4lMUV1 = rat(-673070460489262038*i_,53065174932988313567);
	 id cMi1L1Eps4lMUV2 = rat(2100765274346925365*i_,270882405287681646017);
	 id cMi1L1Eps4lMUV3 = rat(-912356673578020461*i_,316511807529973429019);
	 id cMi1L1Eps4lMUV4 = rat(44728035633278509*i_,86761845665812257671);
	 id cMi1L1Eps4lMUV5 = rat(-21677157797392643*i_,410774346240489478193);
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
