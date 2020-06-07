*---- defnition of eps-expanded mi1L2: lMUV=log(mUV^2)
	id mi1L2 = (
	 cMi1L2EpsM2lMUV0*rat(lMUV^0,1)*rat(1,eps^2) 
	 +cMi1L2EpsM1lMUV0*rat(lMUV^0,1)*rat(1,eps^1) 
	 +cMi1L2EpsM1lMUV1*rat(lMUV^1,1)*rat(1,eps^1) 
	 +cMi1L2Eps0lMUV0*rat(lMUV^0,1)*rat(eps^0,1) 
	 +cMi1L2Eps0lMUV1*rat(lMUV^1,1)*rat(eps^0,1) 
	 +cMi1L2Eps0lMUV2*rat(lMUV^2,1)*rat(eps^0,1) 
	 +cMi1L2Eps1lMUV0*rat(lMUV^0,1)*rat(eps^1,1) 
	 +cMi1L2Eps1lMUV1*rat(lMUV^1,1)*rat(eps^1,1) 
	 +cMi1L2Eps1lMUV2*rat(lMUV^2,1)*rat(eps^1,1) 
	 +cMi1L2Eps1lMUV3*rat(lMUV^3,1)*rat(eps^1,1) 
	 +cMi1L2Eps2lMUV0*rat(lMUV^0,1)*rat(eps^2,1) 
	 +cMi1L2Eps2lMUV1*rat(lMUV^1,1)*rat(eps^2,1) 
	 +cMi1L2Eps2lMUV2*rat(lMUV^2,1)*rat(eps^2,1) 
	 +cMi1L2Eps2lMUV3*rat(lMUV^3,1)*rat(eps^2,1) 
	 +cMi1L2Eps2lMUV4*rat(lMUV^4,1)*rat(eps^2,1) 
	 +alarmMi1L2Eps3*rat(eps^3,1)
	);
*
*
*---- BEGIN: defintion of the numeric coefficients for mi1L2: order 100 digits accurat 
*
	 id cMi1L2EpsM2lMUV0 = rat(1053830132972*i_^2,26279074651403377);
	 id cMi1L2EpsM1lMUV0 = rat(3339596090290*i_^2,21311861810762421);
	 id cMi1L2EpsM1lMUV1 = rat(-2249444293247*i_^2,28046889463912230);
	 id cMi1L2Eps0lMUV0 = rat(24451607791211*i_^2,65707409260937486);
	 id cMi1L2Eps0lMUV1 = rat(-7908143437641*i_^2,25233180235885567);
	 id cMi1L2Eps0lMUV2 = rat(2249444293247*i_^2,28046889463912230);
	 id cMi1L2Eps1lMUV0 = rat(16300708140685*i_^2,26105439591783829);
	 id cMi1L2Eps1lMUV1 = rat(-25321691305119*i_^2,34022767500029765);
	 id cMi1L2Eps1lMUV2 = rat(7908143437641*i_^2,25233180235885567);
	 id cMi1L2Eps1lMUV3 = rat(-2249444293247*i_^2,42070334195868345);
	 id cMi1L2Eps2lMUV0 = rat(35128989601063*i_^2,41642764367418656);
	 id cMi1L2Eps2lMUV1 = rat(-21652387683479*i_^2,17338053470147957);
	 id cMi1L2Eps2lMUV2 = rat(25321691305119*i_^2,34022767500029765);
	 id cMi1L2Eps2lMUV3 = rat(-5272095625094*i_^2,25233180235885567);
	 id cMi1L2Eps2lMUV4 = rat(655292079547*i_^2,24511259838894524);
*
*---- END: defintion of the coefficients for mi1L2
*
* 
*---- BEGIN: defintion of the analytc coefficients for mi1L2
*
*	cMi1L2EpsM2lMUV0=imI^2/(256*Pi^4)
*
*
*	cMi1L2EpsM1lMUV0=(imI^2*(-EulerGamma + Log[4*Pi]))/(128*Pi^4)
*
*
*	cMi1L2EpsM1lMUV1=-imI^2/(128*Pi^4)
*
*
*	cMi1L2Eps0lMUV0=(imI^2*(12*EulerGamma^2 + Pi^2 - 24*EulerGamma*Log[4*Pi] + 12*Log[4*Pi]^2))/(1536*Pi^4)
*
*
*	cMi1L2Eps0lMUV1=(imI^2*(EulerGamma + Log[1/(4*Pi)]))/(64*Pi^4)
*
*
*	cMi1L2Eps0lMUV2=imI^2/(128*Pi^4)
*
*
*	cMi1L2Eps1lMUV0=(imI^2*(-4*EulerGamma^3 + 32*Log[2]^3 + 48*Log[2]^2*Log[Pi] + 24*Log[2]*Log[Pi]^2 + 4*Log[Pi]^3 + 12*EulerGamma^2*Log[4*Pi] + Pi^2*Log[4*Pi] - EulerGamma*(Pi^2 + 12*Log[4*Pi]^2) - 2*Ze
*	ta[3]))/(768*Pi^4)
*
*
*	cMi1L2Eps1lMUV1=-(imI^2*(12*EulerGamma^2 + Pi^2 - 24*EulerGamma*Log[4*Pi] + 12*Log[4*Pi]^2))/(768*Pi^4)
*
*
*	cMi1L2Eps1lMUV2=(imI^2*(-EulerGamma + Log[4*Pi]))/(64*Pi^4)
*
*
*	cMi1L2Eps1lMUV3=-imI^2/(192*Pi^4)
*
*
*	cMi1L2Eps2lMUV0=(imI^2*(240*EulerGamma^4 + 7*Pi^4 - 960*EulerGamma^3*Log[4*Pi] + 120*Pi^2*Log[4*Pi]^2 + 120*EulerGamma^2*(Pi^2 + 12*Log[4*Pi]^2) + 240*Log[4*Pi]*(8*Log[2]^3 + 12*Log[2]^2*Log[Pi] + Log
*	[64]*Log[Pi]^2 + Log[Pi]^3 - 2*Zeta[3]) - 240*EulerGamma*(32*Log[2]^3 + 48*Log[2]^2*Log[Pi] + 24*Log[2]*Log[Pi]^2 + 4*Log[Pi]^3 + Pi^2*Log[4*Pi] - 2*Zeta[3])))/(92160*Pi^4)
*
*
*	cMi1L2Eps2lMUV1=(imI^2*(4*EulerGamma^3 - 32*Log[2]^3 - 48*Log[2]^2*Log[Pi] - 24*Log[2]*Log[Pi]^2 - 4*Log[Pi]^3 - 12*EulerGamma^2*Log[4*Pi] - Pi^2*Log[4*Pi] + EulerGamma*(Pi^2 + 12*Log[4*Pi]^2) + 2*Zet
*	a[3]))/(384*Pi^4)
*
*
*	cMi1L2Eps2lMUV2=(imI^2*(12*EulerGamma^2 + Pi^2 - 24*EulerGamma*Log[4*Pi] + 12*Log[4*Pi]^2))/(768*Pi^4)
*
*
*	cMi1L2Eps2lMUV3=(imI^2*(EulerGamma + Log[1/(4*Pi)]))/(96*Pi^4)
*
*
*	cMi1L2Eps2lMUV4=imI^2/(384*Pi^4)
*
*
*---- END: defintion of the analytc coefficients for mi1L2
