#include "DI_PathDef.h"
#include "VectorProducts.h"

cvector zeroCvector(dimensions,0.);
rvector zeroRvector(dimensions,0.);

/*Constructor*/

PathDeformation::PathDeformation(my_real const angle){
  init_momenta(angle);
  loop_setQ=false;
  deformedQ=false;
}

/* path point to be deformed*/
void PathDeformation::set_loop_mom(my_real l0,my_real l1,my_real l2,my_real l3){
  k={l0,l1,l2,l3};
  loop_setQ=true;
  deformedQ=false;
}


/*Setup the external momenta and global variables 
  independent from the loop momentum*/
bool PathDeformation::init_momenta(my_real angle){
  if(momentaQ) return true;
  
  p1.at(0)= 0.5; p1.at(1)= 0.5;
  p2.at(0)= 0.5; p2.at(1)=-0.5;
  
  p3.at(0)= 0.5; p3.at(1)= 0.5*cos(angle); p3.at(2)= 0.5*sin(angle);
  p4.at(0)= 0.5; p4.at(1)=-0.5*cos(angle); p4.at(2)=-0.5*sin(angle);

  q1=zeroRvector;  q2=q1-p2;
  q3=p1-p3;        q4=p1;
  Q={q1,q2,q3,q4};

  //Set the constant values M1 and M2
  my_real sp12=SP(p1,p2);
  M1=0.0025*sp12;
  M2=sp12;
  M3=sp12;

  //Set the vector v from eq(22)
  rvector q_a1=q4-q2;
  rvector q_a1bar=q2-q1;
  rvector q_a2=q1-q1;
  rvector q_a2bar=q1-q2;
  my_real a1=SP(q_a1,p2)/sp12;
  my_real a1bar=SP(q_a1bar,p2)/sp12;
  my_real a2=SP(q_a2,p1)/sp12;
  my_real a2bar=SP(q_a2bar,p1)/sp12;
  
  
  rvector v1=a1*q1 + a1bar*q4;
  rvector v2=a2*q2 + a2bar*q1;
  v=.5*(v1+v2);

  return true;
}

void PathDeformation::compute_box_propagators(){
  if(!deformedQ){
    printf("Before computing the propagators it's necessary to deform the path!\n");
    exit(1);
  }
  vector<cvector> Qc(4,zeroCvector);
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      Qc[i][j]=Q[i][j];
    }
  }  
  prop1=square(ell-Qc[0]);
  prop2=square(ell-Qc[1]);
  prop3=square(ell-Qc[2]);
  prop4=square(ell-Qc[3]);
}

/*
 *  arXiv:0812.2686v1 
 *  case for A=1, where the two incoming momenta are adjacent 
 */
 
// Check if ell is inside the negative light cone starting from Q[j];
bool PathDeformation::ellinCminus(int j){
  return (square(k - Q[j]) > 0.0 && k[0] < Q[j][0]);
}
// Check if ell is inside the positive light cone starting from Q[j];
bool PathDeformation::ellinCplus(int j){
  return (square(k - Q[j]) > 0.0 && k[0] > Q[j][0]);
}  

/*************************
/* Deformation Functions *
*************************/
//h- function
inline my_real PathDeformation::hm(rvector q){
  my_real absvk=0.;
  for(int i=1; i < dimensions; i++) absvk+=pow(q.at(i),2);
  absvk=sqrt(absvk);
  my_real Ek=q.at(0);

  if(Ek < -absvk) return 0;
  my_real val=pow(absvk+Ek,2);

  return val/(val+M1);
  
}

//h+ function
inline my_real PathDeformation::hp(rvector q){
  my_real absvk=0.;
  for(int i=1; i < dimensions; i++) absvk+=pow(q.at(i),2);
  absvk=sqrt(absvk);
  my_real Ek=q.at(0);

  if(Ek > absvk) return 0;
  my_real val=pow(absvk-Ek,2);

  return val/(val+M1);
  
}

//g(l) function
inline my_real PathDeformation::g(rvector& q){
  my_real gamma1=.7;

  my_real l0=q.at(0);
  my_real v0=v.at(0);
  
  rvector vect=q-v;
  vect.at(0)=0.;
  
  my_real val1=pow(l0-v0,2);
  my_real val2=-SP(vect,vect);
    
  return gamma1*M2/(val1+val2+M2);
}

//g+ function
inline my_real PathDeformation::gp(rvector& q){
  my_real gamma2=1;
  
  my_real l0=q.at(0);
  my_real v0=v.at(0);

  rvector vect=q-v;
  vect.at(0)=0.;
  
  my_real E=l0-v0;
  my_real w=sqrt(-SP(vect,vect)+M2);
  my_real val=pow(1+E/w,2);

  return gamma2/(1+val);
}

//g- function
inline my_real PathDeformation::gm(rvector& q){
  my_real gamma2=1;
  
  my_real l0=q.at(0);
  my_real v0=v.at(0);

  rvector vect=q-v;
  vect.at(0)=0.;
  
  my_real E=l0-v0;
  my_real w=sqrt(-SP(vect,vect)+M2);
  my_real val=pow(1-E/w,2);

  return gamma2/(1+val);
}


/*************************
/*       Derivative      *
*************************/

// Gradient of log(hm(k));
rvector PathDeformation::gradloghm(const rvector k)
{
  rvector grad;
  my_real Ek = k.at(0);
  my_real absvk = sqrt( k.at(1)*k.at(1) + k.at(2)*k.at(2) + k.at(3)*k.at(3));
  if(Ek < -absvk) return zeroRvector;
  
  my_real w = absvk + Ek;
  rvector gradw=(1/absvk)*k;
  gradw.at(0)=1.0;
  grad = (2.0*M1/w/(pow(w,2) + M1))*gradw;
  return grad;

}

// Gradient of log(hp(k));
rvector PathDeformation::gradloghp(const rvector k)
{
  rvector grad;
  my_real Ek = k.at(0);
  my_real absvk = sqrt( k.at(1)*k.at(1) + k.at(2)*k.at(2) + k.at(3)*k.at(3));
  if( absvk < Ek) return zeroRvector;
  else
    {
      my_real w = absvk - Ek;
      rvector gradw=(1/absvk)*k;
      gradw.at(0)=-1.0;
      grad = (2.0*M1/w/(pow(w,2) + M1))*gradw;
      return grad;
    }
}

// Gradient of log(g(k));
// Note that grad propto the contravariant components of k since we differentiate 
// the euclidian square of k.
rvector PathDeformation::gradlogg(const rvector k)
{
  //v=0
  rvector grad;
  my_real kEsq =k.at(0)*k.at(0)+k.at(1)*k.at(1)+k.at(2)*k.at(2)+k.at(3)*k.at(3);
  grad = (-2.0/( kEsq + M2 ))*k;
  return grad;
}

// Gradient of log(g+(k));
rvector PathDeformation::gradloggp(const rvector k)
{
  //v=0
  rvector grad=zeroRvector;
  my_real omega = sqrt( k.at(1)*k.at(1) + k.at(2)*k.at(2) + k.at(3)*k.at(3) + M2 );
  my_real E = k.at(0);
  my_real ratio = E/omega;
  my_real factor = 2.0*(1 + ratio)/(1.0 + pow((1.0 + ratio), 2))/omega;
  grad.at(0) = - factor;
  for(int mu = 1; mu<=3; mu++) grad.at(mu) = factor*ratio/omega * k.at(mu);
  return grad;
}

// Gradient of log(g-(k));
rvector PathDeformation::gradloggm(const rvector k)
{
  //v=0
  rvector grad=zeroRvector;
  my_real omega = sqrt( k.at(1)*k.at(1) + k.at(2)*k.at(2) + k.at(3)*k.at(3) + M2 );
  my_real E = k.at(0);
  my_real ratio = E/omega;
  my_real factor = 2.0*(1 - ratio)/(1.0 + pow((1.0 - ratio), 2))/omega;
  grad.at(0) =  factor;
  for(int mu = 1; mu<=3; mu++) grad.at(mu) = - factor*ratio/omega * k.at(mu);
  return grad;
}


// l is deformed to ell with the corresponding jacobian
void PathDeformation::deform_loop_path(){

  /*==================================================
    =        COMPUTE KAPPA_0 AND ITS GRADIENT        =
    ==================================================*/
  if(!loop_setQ){
    printf("You need to select a point to be deformed using PathDeformation::set_loop_mom(my_real l0,my_real l1,my_real l2,my_real l3).");
    exit(0);
  }
  
  vector< my_real> cj(4);
  vector<rvector> gradlogcj(4,zeroRvector);
  vector<rvector> gradcj(4,zeroRvector);
  rvector tmp_cv;
  
  
  rvector kappa=zeroRvector;
  vector<rvector> gradkappa(dimensions,zeroRvector);
  
  //C_J functions
  cj[2]=hm(k-q2)*hp(k-q4)*hm(k-q1)*hp(k-q1)*g(k);  //c3
  cj[0]=hm(k-q3)*hp(k-q3)*g(k); //c1
  cj[1]=hp(k-q3)*hp(k-q1)*g(k); //c2
  cj[3]=hm(k-q3)*hm(k-q1)*g(k); //cN

  //logC_3 derivative
  tmp_cv = gradloghm(k-q2);
  tmp_cv += gradloghp(k-q4);
  tmp_cv += gradloghm(k-q1);
  tmp_cv += gradloghp(k-q1);
  tmp_cv += gradlogg(k);
  gradlogcj[2]=tmp_cv;
  //logC_1 derivative
  tmp_cv = gradloghm(k-q3);
  tmp_cv += gradloghp(k-q3);
  tmp_cv += gradlogg(k);
  gradlogcj[0]=tmp_cv;
  //logC_2 derivative
  tmp_cv = gradloghp(k-q3);
  tmp_cv += gradloghp(k-q1);
  tmp_cv += gradlogg(k);
  gradlogcj[1]=tmp_cv;
  //logC_4 derivative
  tmp_cv = gradloghm(k-q3);
  tmp_cv += gradloghm(k-q1);
  tmp_cv += gradlogg(k);
  gradlogcj[3] = tmp_cv;

  //Check regions
  bool askQ;
  askQ=(! ellinCminus(1)) && (! ellinCplus(3)) &&
    (! ellinCminus(0)) && (! ellinCplus(0)) ;
  if(!askQ){
    cj[2]=0.0;
    gradlogcj[2]=zeroRvector;    
  }

  askQ=(! ellinCminus(3-1)) && (! ellinCplus(3));
  if(!askQ){
    cj[0]=0.0;
    gradlogcj[0]=zeroRvector;    
  }

  askQ=(! ellinCplus(3-1)) && (! ellinCplus(1-1));
  if(!askQ){
    cj[1]=0.0;
    gradlogcj[1]=zeroRvector;    
  }

  askQ=(! ellinCminus(3-1)) && (! ellinCminus(1-1));
  if(!askQ){
    cj[2]=0.0;
    gradlogcj[2]=zeroRvector;    
  }

  for(int j=0; j<4; j++){
    tmp_cv=cj.at(j)*gradlogcj.at(j);
    gradcj.at(j)=tmp_cv;
  }
  
  //Update the value of kappa and gradkappa from the c_j functions
  for(int j=0; j<4; j++){
    kappa+=(-cj[j])*(k-Q[j]);  // kappa update
    for(int mu=0; mu<4; mu++){ // gradkappa update
      gradkappa[mu][mu]+=(-1)*cj[j];
      for(int nu=0; nu<4; nu++){
	gradkappa[mu][nu]+=(-1.)*(k[mu]-Q[j][mu])*gradcj[j][nu];
      }
    }
  }

  // Start to compute the contribution from the c+ c- functions
  my_real cp, cm;
  rvector gradlogcp=zeroRvector;
  rvector gradlogcm=zeroRvector;

  //Set the variable x and xbar from eq(16)
  my_real sp12=SP(p1,p2);
  x    = SP(k-q2,p2)/sp12;
  xbar = SP(k-q1,p1)/sp12;
  //C+ and C-
  cp = (x+xbar) > 0 && (! ellinCminus(3))? hm(k-q4)*(x+xbar)*gm(k) : 0;
  cm = (x+xbar) < 0 && (! ellinCplus(1))? -hp(k-q2)*(x+xbar)*gp(k) : 0;
  //gradcplus grandcminus
  gradlogcp  = gradloghm(k - q4);
  gradlogcp += (1.0/(x + xbar)/sp12)*(p1+p2);
  gradlogcp += gradloggm(k);

  gradlogcm  = gradloghp(k - q2);
  gradlogcm += (1.0/(x + xbar)/sp12)*(p1+p2);
  gradlogcm += gradloggp(k);

  //Update kappa
  kappa+=(cp-cm)*(p1+p2);
  //Update gradkappa
  for(int mu=0; mu<4; mu++)
    for(int nu=0; nu<4; nu++){
      gradkappa[mu][nu]+=+(p1[mu]+p2[mu])*cp*gradlogcp[nu];
      gradkappa[mu][nu]+=-(p1[mu]+p2[mu])*cm*gradlogcm[nu];
    }

  /*==================================================
    =        COMPUTE LAMBDA AND ITS GRADIENT        =
    ==================================================*/
  // So far kappa and gradkappa are really kappa_0 and its gradient.
  // Now we need the lambda_i.
  my_real lambda_i; 
  int which_lambda = -1, nverts= 4;

  // Set our maximum possible value for lambda. One could change the value.
  my_real lambda = 1; // to .1
  
  // Compute sum of all the c_j and their gradients
  my_real cjsum=0; 
  rvector gradcjsum=zeroRvector;
  for(int j=0; j<4; j++) {
    cjsum+=cj[j];
    gradcjsum+=gradcj[j];
  }

  if( 0.25/cjsum < lambda){//smallest between lambda_0 and lambda
    lambda = 0.25/cjsum;
    which_lambda = nverts; // lambda_0 is the smallest so far.
  }
  
  my_real kappa2, kqi2, dotkqikappa; //kqi alias for k-Q[i]
  for(int i = 0; i < nverts; i++){
    kappa2 = SP(kappa,kappa); 
    dotkqikappa = SP(kappa,k-Q[i]);
    kqi2 = SP(k-Q[i],k-Q[i]);
    
    if( kappa2*kqi2 > 2.0*pow(dotkqikappa,2) ){
      // Region 1, 2.0*pow(dotkqikappa,2). < kappa2*kqi2 
      lambda_i = 0.5*sqrt(kqi2/kappa2);
    }else if( kappa2*kqi2 > 0.0 ){
      // Region 2, 0.0 < kappa2*kqi2 < 2.0*pow(dotkqikappa,2).
      lambda_i = sqrt(4.0*pow(dotkqikappa,2) - kappa2*kqi2);
      lambda_i *= 0.5/abs(kappa2);
    }else{
      // Region 3, kappa2*kqi2 < 0.0.
      lambda_i = sqrt(4.0*pow(dotkqikappa,2) - 2.0*kappa2*kqi2);
      lambda_i *= 0.5/abs(kappa2);    
    }
    
    // We now have lambda_i, select smallest compared to lambda
    if(lambda_i < lambda){ 
      lambda = lambda_i;
      which_lambda = i; //lamda_i is the smallest so far
    }
  }

  // We have lambda. Now calculate its gradient.
  rvector gradlambda;

  if(which_lambda < 0)
    //Trivial case
    gradlambda = zeroRvector;
  else if(which_lambda == nverts){
    //Grad corresponding to lambda_0
    gradlambda = (- 0.25/pow(cjsum,2))*gradcjsum;
  }else if(which_lambda < nverts){
    //Grad corresponding to lambda_i, i==which_lambda
    kappa2 = SP(kappa,kappa);
    dotkqikappa = SP(kappa,k-Q[which_lambda]);
    kqi2 = SP(k-Q[which_lambda],k-Q[which_lambda]);

    // We need the gradient of kappa2.
    rvector gradkappa2=zeroRvector;
    for(int nu = 0; nu < 4; nu++){
      my_real temp = (2.0*kappa[0]*gradkappa[0][nu]);
      for(int mu = 1; mu < 4; mu++)
	temp += ((-2.0)*kappa[mu]*gradkappa[mu][nu]);
      gradkappa2.at(nu) = temp;
    }

    // We need the gradient of kqi2, aka (k-Q[i])^2
    rvector gradkqi2=zeroRvector;
    gradkqi2[0] = 2.0*(k[0] - Q[which_lambda][0]);
    for(int nu = 1; nu < 4; nu++)
      gradkqi2[nu] = - 2.0*(k[nu] - Q[which_lambda][nu]);
    
    // We need the gradient of dotkqikappa, for which we need
    // also kappac, the covariant components of kappa.
    rvector graddotkqikappa=zeroRvector;
    rvector kappacov={kappa[0], -kappa[1], -kappa[2], -kappa[3]};
    for(int nu = 0; nu < 4; nu++){
      my_real temp = kappacov[nu];
      temp += ((k[0] - Q[which_lambda][0])*gradkappa[0][nu]);
      for(int mu = 1; mu < 4; mu++)
	temp = ((-1.)*(k[mu] - Q[which_lambda][mu])*gradkappa[mu][nu]);
      graddotkqikappa[nu] = temp;
    }

    // Now we have the ingredients, so we calculate gradlambda.
    rvector gradlambdasq=zeroRvector, gradnumerator=zeroRvector;
    my_real numerator;
    if( kappa2*kqi2 > 2.0*pow(dotkqikappa,2) ){
      // Region 1.
      numerator = kappa2*kqi2;
      gradnumerator = kappa2*gradkqi2 + kqi2*gradkappa2;
    }else if( kappa2*kqi2 > 0.0 ){
      // Region 2, 0.0 < kappa2*kqi2 < 2.0*pow(dotkqikappa,2).
      numerator = 4.0*pow(dotkqikappa,2) - kqi2*kappa2;
      gradnumerator = 8.0*dotkqikappa*graddotkqikappa 
	- kappa2*gradkqi2 - kqi2*gradkappa2;
    }else{
      // Region 3, kappa2*kqi2 < 0.0.
      numerator = 4.0*pow(dotkqikappa,2) - 2.0*kqi2*kappa2;
      gradnumerator = 8.0*dotkqikappa*graddotkqikappa
	- 2.0*kappa2*gradkqi2 - 2.0*kqi2*gradkappa2;
    }
    
    gradlambdasq = gradnumerator - (2.0*numerator/kappa2)*gradkappa2;
    gradlambdasq = (1.0/4.0/pow(kappa2,2))*gradlambdasq;
    gradlambda = (1.0/2.0/lambda)*gradlambdasq;
  }
  
  /*==================================================
    =      COMPUTE FINAL KAPPA AND ITS GRADIENT      =
    ==================================================*/
  // We now have lambda and gradlambda. 
  // Use this to get the new kappa and gradkappa.
  for(int nu = 0; nu < 4 ; nu++)
    for(int mu = 0; mu < 4;mu++)
      gradkappa[mu][nu] = lambda*gradkappa[mu][nu] + kappa[mu]*gradlambda[nu];
  
  kappa = lambda*kappa;

  /*==================================================
    =                    OUTPUT                      =
    ==================================================*/
  my_complex ii(0.,1.);
  
  //New Momentum along the modified path
  for(int mu = 0; mu <= 3; mu++)
    ell[mu]=k[mu]+ii*kappa[mu];
  
  //Compure the Jacobian
  vector<cvector> grad(dimensions,zeroCvector);
  for(int mu = 0; mu <= 3; mu++){
    grad[mu][mu]+=1.;
    for(int nu = 0; nu <= 3; nu++)
      grad[mu][nu]+=ii*gradkappa[mu][nu];
  }
  jacobian = Determinant(grad);

  deformedQ=true;
  return void();
}



// Function Determinant(M,d).
// This function calculates the determinant of any input matrix using LU decomposition.
// Code by W. Gong and D.E. Soper
my_complex Determinant(const vector<cvector > & matrix)
{
  // Define matrix related variables.
  int dimension=matrix.at(0).size();
  int maxsize=dimension-1;
  static const my_complex one(1.0, 0.0);
  my_complex determinant;
  my_complex aa[dimension][dimension];
  int indx[dimension], d;
  // Initialize the determinant.
  determinant = one;
  // Inintialize the matrix to be decomposed with the transferred matrix b.
  for(int i=0; i < dimension; i++)
    {
      for(int j=0; j < dimension; j++)
	{
	  aa[i][j] = matrix[i][j];
	}
    }
  // Define parameters used in decomposition.
  int i, imax, j, k, flag=1;
  my_complex dumc,sum;
  my_real aamax, dumr;
  my_real vv[maxsize];
  // Initialize the parity parameter.
  d = 1;
  // Get the implicit scaling information.
  for(i = 0; i < dimension; i++)
    {
      aamax=0;
      for(j = 0; j < dimension; j++)
	{
	  if(abs(aa[i][j]) > aamax) aamax=abs(aa[i][j]);
	}
      // Set a flag to check if the determinant is zero.
      if(aamax == 0) flag=0;
      // Save the scaling.
      vv[i]=1.0/aamax;
    }
  if(flag == 1)
    {
      for(j = 0; j < dimension; j++)
	{
	  for(i = 0; i < j; i++)
	    {
	      sum = aa[i][j];
	      for(k = 0; k < i; k++) sum = sum - aa[i][k]*aa[k][j];
	      aa[i][j] = sum;
	    }
	  //Initialize for the search for largest pivot element.
	  aamax=0;
	  for(i = j; i < dimension; i++)
	    {
	      sum = aa[i][j];
	      for(k = 0; k < j; k++) sum = sum - aa[i][k]*aa[k][j];
	      aa[i][j]=sum;
	      // Figure of merit for the pivot.
	      dumr = vv[i] * abs(sum);
	      // Is it better than the best so far?
	      if(dumr >= aamax)
		{
		  imax=i;
		  aamax=dumr;
		}
	    }  // End for(i = j; i < dimension; i++)
	  // See if we need to interchange rows.
	  if(j != imax)
	    {
	      for(k = 0; k < dimension; k++)
		{
		  dumc = aa[imax][k];
		  aa[imax][k] = aa[j][k];
		  aa[j][k] = dumc;
		}
	      // Change the parity of d.
	      d = -d;
	      // Interchange the scale factor.
	      vv[imax] = vv[j];
	    }  // End if(j != imax)
	  indx[j]=imax;
	  if(j != dimension - 1)
	    {
	      dumc = 1.0/aa[j][j];
	      for(i = j+1; i < dimension; i++) aa[i][j] = aa[i][j]*dumc;
	    }
	} // End for(j = 0; j < dimension; j++)
    } // END if(flag == 1)
  // Calculate the determinant using the decomposed matrix.
  if(flag == 0) determinant = 0.0;
  else
    {    // Multiply the diagonal elements.
      for(int diagonal = 0; diagonal < dimension; diagonal++)
	{
	  determinant = determinant * aa[diagonal][diagonal];
	}
      determinant=(my_real)d*determinant;
    } // End if(flag == 0) ... else
  return determinant;
}
// End of function Determinant(M,d).
