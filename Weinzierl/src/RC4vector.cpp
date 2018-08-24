#include "DCD.h"

/*==========================================================
  =                    Real 4-vector  
  ==========================================================*/

//Constructor
DIdeform::R4vector::R4vector()
{
  std::vector<my_real>(4, 0.).swap(v);
}

DIdeform::R4vector::R4vector(const std::vector<my_real> &w)
{
  if (w.size() != 4)
  {
    printf("R4vector needs to be initialized with a 4d vector\n");
    exit(1);
  }
  v = w;
}

//Vector functions
my_real &DIdeform::R4vector::operator[](const int mu)
{
  return this->v.at(mu);
}

const my_real &DIdeform::R4vector::operator()(const int mu) const
{
  return this->v.at(mu);
}

DIdeform::R4vector DIdeform::R4vector::dual() const
{
  std::vector<my_real> dual_v(4);
  dual_v[0] = (*this)(0);
  for (int mu = 1; mu < 4; mu++)
    dual_v[mu] = -(*this)(mu);
  return dual_v;
}

my_real DIdeform::R4vector::square() const
{
  return (*this) * (*this);
}

//Operators
DIdeform::R4vector DIdeform::R4vector::operator+(const DIdeform::R4vector &rhs) const
{
  std::vector<my_real> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = (*this)(mu) + rhs(mu);
  return result;
}

DIdeform::R4vector DIdeform::R4vector::operator-()
{
  std::vector<my_real> negative(4, 0.);
  for (int mu = 0; mu < 4; mu++)
    negative[mu] = -(*this)(mu);
  return negative;
}

DIdeform::R4vector DIdeform::R4vector::operator-(const DIdeform::R4vector &rhs) const
{
  std::vector<my_real> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = (*this)(mu)-rhs(mu);
  return result;
}

DIdeform::R4vector DIdeform::operator*(const my_real a, const DIdeform::R4vector &w)
{
  std::vector<my_real> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = a * w(mu);
  return result;
}

std::ostream &DIdeform::operator<<(std::ostream &os, DIdeform::R4vector rhs) 
{
  os << "[";
  for (int mu = 0; mu < 4; mu++)
  {
    if (rhs(mu) > 0.0)
      os << "+";
    os << rhs(mu);
    if (mu < 3)
      os << ", ";
  }
  os << "]";
  return os;
}

//Scalar product with metric (+,-,-,-)
my_real DIdeform::R4vector::operator*(const DIdeform::R4vector &rhs) const
{
  my_real result = (*this)(0) * rhs(0);
  for (int mu = 1; mu < 4; mu++)
    result -= (*this)(mu)*rhs(mu);
  return result;
}
/*--------------- End Real 4-vector --------------------*/

/*==========================================================
  =                  Complex 4-vector
  ==========================================================*/

//Constructor
DIdeform::C4vector::C4vector()
{
  std::vector<my_comp>(4, 0.).swap(v);
}

DIdeform::C4vector::C4vector(const std::vector<my_comp> &w)
{
  if (w.size() != 4)
  {
    printf("R4vector needs to be initialized with a 4d vector\n");
    exit(1);
  }
  v = w;
}

//Vector functions
my_comp &DIdeform::C4vector::operator[](const int mu)
{
  return this->v.at(mu);
}

const my_comp &DIdeform::C4vector::operator()(const int mu) const
{
  return this->v.at(mu);
}

DIdeform::C4vector DIdeform::C4vector::dual() const
{
  std::vector<my_comp> dual_v(4);
  dual_v[0] = (*this)(0);
  for (int mu = 1; mu < 4; mu++)
    dual_v[mu] = -(*this)(mu);
  return dual_v;
}

my_comp DIdeform::C4vector::square() const
{
  return (*this) * (*this);
}

//Operators
DIdeform::C4vector DIdeform::C4vector::operator+(const DIdeform::C4vector &rhs) const
{
  std::vector<my_comp> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = (*this)(mu) + rhs(mu);
  return result;
}

DIdeform::C4vector DIdeform::C4vector::operator-()
{
  std::vector<my_comp> negative(4, 0.);
  my_comp zero(0., 0.);
  for (int mu = 0; mu < 4; mu++)
    negative[mu] = zero - (*this)(mu);
  return negative;
}

DIdeform::C4vector DIdeform::C4vector::operator-(const DIdeform::C4vector &rhs) const
{
  std::vector<my_comp> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = (*this)(mu)-rhs(mu);
  return result;
}

DIdeform::C4vector DIdeform::operator*(const my_comp a, const DIdeform::C4vector &w) 
{
  std::vector<my_comp> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = a * w(mu);
  return result;
}

std::ostream &DIdeform::operator<<(std::ostream &os, DIdeform::C4vector rhs)
{
  os << "[";
  for (int mu = 0; mu < 4; mu++)
  {
    os << rhs(mu);
    if (mu < 3)
      os << ", ";
  }
  os << "]";
  return os;
}

//Scalar product with metric (+,-,-,-)
my_comp DIdeform::C4vector::operator*(const DIdeform::C4vector &rhs) const
{
  my_comp result = (*this)(0) * rhs(0);
  for (int mu = 1; mu < 4; mu++)
    result -= (*this)(mu)*rhs(mu);
  return result;
}

/*--------------- End Complex 4-vector --------------------*/

/*==========================================================
  =                   Mixed Operators
  ==========================================================*/
DIdeform::C4vector DIdeform::operator*(const my_comp a, const DIdeform::R4vector &w)
{
  std::vector<my_comp> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = a * w(mu);
  return result;
}

DIdeform::C4vector DIdeform::operator+(const DIdeform::C4vector &lhs, const DIdeform::R4vector &rhs)
{
  std::vector<my_comp> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = lhs(mu) + rhs(mu);
  return result;
}

DIdeform::C4vector DIdeform::operator+(const DIdeform::R4vector &lhs, const DIdeform::C4vector &rhs)
{
  return rhs + lhs;
}

DIdeform::C4vector DIdeform::operator-(const DIdeform::C4vector &lhs, const DIdeform::R4vector &rhs)
{
  std::vector<my_comp> result(4);
  for (int mu = 0; mu < 4; mu++)
    result[mu] = lhs(mu)-rhs(mu);
  return result;
}

DIdeform::C4vector DIdeform::operator-(const DIdeform::R4vector &lhs, const DIdeform::C4vector &rhs)
{
  return rhs - lhs;
}

//Scalar product with metric (+,-,-,-)
my_comp DIdeform::operator*(const DIdeform::R4vector &lhs, const DIdeform::C4vector &rhs)
{
  my_comp result = lhs(0) * rhs(0);
  for (int mu = 1; mu < 4; mu++)
    result -= lhs(mu)*rhs(mu);
  return result;
}

my_comp DIdeform::operator*(const DIdeform::C4vector &lhs, const DIdeform::R4vector &rhs)
{
  return rhs * lhs;
}