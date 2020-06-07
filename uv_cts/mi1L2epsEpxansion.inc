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
	 id cMi1L2EpsM2lMUV0 = rat(115804891537979393*i_^2,2887795000833488180494);
	 id cMi1L2EpsM1lMUV0 = rat(164377612643638857*i_^2,1048987024397929152763);
	 id cMi1L2EpsM1lMUV1 = rat(-115692154855406038*i_^2,1442491858461459690459);
	 id cMi1L2Eps0lMUV0 = rat(321683385408215851*i_^2,864441391255224539716);
	 id cMi1L2Eps0lMUV1 = rat(-142554264587898449*i_^2,454859914985757447942);
	 id cMi1L2Eps0lMUV2 = rat(115692154855406038*i_^2,1442491858461459690459);
	 id cMi1L2Eps1lMUV0 = rat(340104754187323081*i_^2,544674749016305542867);
	 id cMi1L2Eps1lMUV1 = rat(-903841203288474388*i_^2,1214420424998004940155);
	 id cMi1L2Eps1lMUV2 = rat(142554264587898449*i_^2,454859914985757447942);
	 id cMi1L2Eps1lMUV3 = rat(-77165682131128477*i_^2,1443194679439101890353);
	 id cMi1L2Eps2lMUV0 = rat(292305378538579976*i_^2,346505952492378668711);
	 id cMi1L2Eps2lMUV1 = rat(-154450769515073901*i_^2,123675769134740586643);
	 id cMi1L2Eps2lMUV2 = rat(260474432472042686*i_^2,349979033742780400439);
	 id cMi1L2Eps2lMUV3 = rat(-372401921398758530*i_^2,1782381328236515114463);
	 id cMi1L2Eps2lMUV4 = rat(192970573669107870*i_^2,7218081860128436051447);
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
