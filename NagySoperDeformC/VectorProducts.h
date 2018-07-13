//external int dimensions;

/***************** OPERATORS *****************/
template <typename T>
vector<T> operator-(const vector<T>& q1, const vector<T>& q2){
  vector<T> q(dimensions,0.);
  for(int i=0; i<dimensions; i++) q.at(i)=q1.at(i)-q2.at(i);
  return q;
}
rvector operator+(const rvector& q1, const rvector& q2){
  rvector q(dimensions,0.);
  for(int i=0; i<dimensions; i++) q.at(i)=q1.at(i)+q2.at(i);
  return q;
}
rvector operator+=(rvector& lhs,const rvector& rhs){
  for(int i=0; i<dimensions; i++) lhs.at(i)+=rhs.at(i);
  return lhs;
}
rvector operator*(const my_real c, const rvector & q1){
  rvector q(dimensions);
  for(int i=0; i<dimensions; i++) q.at(i)=c*q1.at(i);
  return q;
}

/***************** Global Function *****************/

//Scalar Product of two 4-momenta
template <typename T1,typename T2>
T1 SP(const vector<T1>& q1, const vector<T2>& q2){
  T1 sp=q1.at(0)*q2.at(0);
  for(int i=1; i<dimensions; i++)
    sp-=q1.at(i)*q2.at(i);
  return sp;
}

//Square 4-momenta
template <typename T>
T square(const vector<T>& q){
  return SP(q,q);
}

//Euclidean Scalar Product
template <typename T1,typename T2>
T1 ESP(const vector<T1>& q1, const vector<T2>& q2){
  T1 sp=q1.at(0)*q2.at(0);
  for(int i=1; i<dimensions; i++)
    sp+=q1.at(i)*q2.at(i);
  return sp;
}

//Determinant
complex<my_real> Determinant(const vector<cvector > & matrix);

