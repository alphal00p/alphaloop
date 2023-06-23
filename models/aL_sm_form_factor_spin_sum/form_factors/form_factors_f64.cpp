#include <cmath>
#include <complex>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_result.h>
#include <math.h>
#include <mp++/complex.hpp>
#include <mp++/complex128.hpp>

using namespace std;
typedef mppp::complex128 complex128;

// supported configurations:

/*
1 top quark (n_massive = 1, n_massless = 0)
1 up quark (n_massive = 0, n_massless = 1)
*/

// this is a constant in here for now, later this should be an argument of the
// functions
double FORM_FACTOR_PHASE = 1.0;

bool debug = true;

complex<double> global_prefactor(0.0, 2.0 / (M_PI * M_PI));

complex<double> zero_f64(0.0, 0.0);
double mu2 = 1000.0;

//  int n_massive = 0;
//  int n_massless = 1;

int n_massive = 4;
int n_massless = 5;

// charm, tau, bottom, top
double pref_arr_massive[4] = {3.0 * 16.0 / 81.0, 1.0, 3.0 / 81.0,
                              3.0 * 16.0 / 81.0};

// up, down, strange, electron, muon
double pref_arr_massless[5] = {3.0 * 16.0 / 81.0, 3.0 / 81.0, 3.0 / 81.0, 1.0,
                               1.0};

// charm, tau, bottom, top
double M2_arr[4] = {1.5 * 1.5, 1.777 * 1.777, 4.75 * 4.75, 173.0 * 173.0};
// double M2_arr[1] = {173.0 * 173.0};

complex<double> I(0.0, 1.0);

struct polynomial_parameters {
  complex<double> xs;
  complex<double> xt;
  complex<double> xu;
  complex<double> a02;
  complex<double> b02s;
  complex<double> b02t;
  complex<double> b02u;
  complex<double> b020;
  complex<double> c02s;
  complex<double> c02t;
  complex<double> c02u;
  complex<double> d02su;
  complex<double> d02st;
  complex<double> d02tu;
};

class c_inf;

c_inf operator+(const c_inf &cl, const c_inf &cr);
c_inf operator+(const double &xl, const c_inf &cr);
c_inf operator+(const c_inf &cl, const double &xr);

c_inf operator-(const c_inf &cl, const c_inf &cr);
c_inf operator-(const double &xl, const c_inf &cr);
c_inf operator-(const c_inf &cl, const double &xr);

c_inf operator*(const c_inf &cl, const c_inf &cr);
c_inf operator*(const double &xl, const c_inf &cr);
c_inf operator*(const c_inf &cl, const double &xr);

c_inf operator/(const c_inf &cl, const c_inf &cr);
c_inf operator/(const double &xl, const c_inf &cr);
c_inf operator/(const c_inf &cl, const double &xr);

c_inf sqrt(const c_inf &c);
c_inf log(const c_inf &c);
c_inf dilog(const c_inf &c);

class c_inf {
  // the sign of the infinitessimal part is currently not correctly carried
  // through the logs and dilogs, but only through the elementary arithmetic
  // operations and the sqrt;
public:
  complex<double> val;
  double inf;

  c_inf(complex<double> new_val, double new_inf) {
    val = new_val;
    inf = new_inf;
  }

  int sign() const {
    if (inf > 0) {
      return 1;
    } else {
      return -1;
    }
  }

  friend ostream &operator<<(ostream &os, const c_inf &c);
};

ostream &operator<<(ostream &os, const c_inf &c) {
  cout << c.val << " sign: " << c.sign();
  return os;
}

c_inf operator+(const c_inf &cl, const c_inf &cr) {
  complex<double> add_val = cl.val + cr.val;
  double add_inf = cl.inf + cr.inf;

  return c_inf(add_val, add_inf);
}

c_inf operator+(const double &xl, const c_inf &cr) {
  complex<double> add_val = xl + cr.val;
  double add_inf = cr.inf;

  return c_inf(add_val, add_inf);
}

c_inf operator+(const c_inf &cl, const double &xr) {
  complex<double> add_val = cl.val + xr;
  double add_inf = cl.inf;

  return c_inf(add_val, add_inf);
}

c_inf operator-(const c_inf &c) { return c_inf(-c.val, -c.inf); }

c_inf operator-(const c_inf &cl, const c_inf &cr) {
  complex<double> sub_val = cl.val - cr.val;
  double sub_inf = cl.inf - cr.inf;

  return c_inf(sub_val, sub_inf);
}

c_inf operator-(const double &xl, const c_inf &cr) {
  complex<double> sub_val = xl - cr.val;
  double sub_inf = -cr.inf;

  return c_inf(sub_val, sub_inf);
}

c_inf operator-(const c_inf &cl, const double &xr) {
  complex<double> sub_val = cl.val - xr;
  double sub_inf = cl.inf;

  return c_inf(sub_val, sub_inf);
}

c_inf operator*(const c_inf &cl, const c_inf &cr) {
  complex<double> mul_val = cl.val * cr.val;
  double mul_inf = cl.val.real() * cr.inf + cr.val.real() * cl.inf;

  return c_inf(mul_val, mul_inf);
}
c_inf operator*(const double &xl, const c_inf &cr) {
  complex<double> mul_val = xl * cr.val;
  double mul_inf = xl * cr.inf;

  return c_inf(mul_val, mul_inf);
}

c_inf operator*(const c_inf &cl, double &xr) {
  complex<double> mul_val = cl.val * xr;
  double mul_inf = cl.inf * xr;

  return c_inf(mul_val, mul_inf);
}

c_inf operator/(const c_inf &cl, const c_inf &cr) {
  complex<double> div_val = cl.val / cr.val;
  double div_inf =
      ((cl.inf * cr.val - cl.val * cr.inf) / (cr.val * cr.val)).real();

  return c_inf(div_val, div_inf);
}

c_inf operator/(const double &xl, const c_inf &cr) {
  complex<double> div_val = xl / cr.val;
  double div_inf = (-xl / cr.val * cr.inf).real();

  return c_inf(div_val, div_inf);
}

c_inf operator/(const c_inf &cl, const double &xr) {
  complex<double> div_val = cl.val / xr;
  double div_inf = cl.inf / xr;

  return c_inf(div_val, div_inf);
}

c_inf log(c_inf const &c) {
  complex<double> complex_res;
  if (c.val.real() > 0 || c.val.imag() != 0) {
    complex_res = log(c.val);
    return c_inf(complex_res, 1 / complex_res.real() * c.inf);
  } else {
    complex_res = log(abs(c.val)) + I * (M_PI * c.sign());
  }
  return c_inf(complex_res, c.inf);
}

c_inf sqrt(c_inf const &c) {
  complex<double> complex_res;
  if (c.val.real() > 0 || c.val.imag() != 0) {
    complex_res = sqrt(c.val);
    return c_inf(complex_res, 1 / (2.0 * sqrt(abs(c.val))) * c.inf);
  } else {
    complex_res = sqrt(abs(c.val)) * I * double(c.sign());
    return c_inf(complex_res, c.inf);
  }
}

c_inf dilog(c_inf const &c) {
  complex<double> complex_res;
  if (c.val.real() < 1 || c.val.imag() != 0) {
    gsl_sf_result re;
    gsl_sf_result im;

    gsl_sf_complex_dilog_xy_e(c.val.real(), c.val.imag(), &re, &im);
    complex_res = complex<double>(re.val, im.val);

    return c_inf(complex_res, c.inf);
  } else {
    double d_log_res = -gsl_sf_dilog(1 / c.val.real());
    c_inf addition_part = -0.5 * log(-c) * log(-c) - M_PI * M_PI / 6.0;

    return d_log_res + addition_part;
  }
}

complex<double> a0_massive(double m2, double mu2) {
  return m2 * (log(mu2 / m2) + 1.0);
}

complex<double> b0_massive(double s, double m2, double mu2) {
  c_inf m2_inf(m2, -1);
  c_inf beta = sqrt(1.0 - 4.0 * m2_inf / s);
  c_inf lambda_plus = 0.5 * (1.0 + beta);
  c_inf lambda_minus = 0.5 * (1.0 - beta);

  complex<double> ln1 = log(-lambda_plus / lambda_minus).val;
  complex<double> ln2 = log(mu2 / m2);
  return 2.0 - beta.val * ln1 + ln2;
}

complex<double> b0_massless(double s, double mu2) {
  c_inf s_inf(s, -1);
  return 2.0 + log(mu2 / (-s_inf)).val;
}

complex<double> b0_massive_0(double m2, double mu2) { return log(mu2 / m2); }

complex<double> c0_massless(double s, double mu2) {
  if (s == 0.0) {
    return 0.0;
  } else {
    c_inf s_inf = c_inf(s, 1);
    complex<double> ln = log(-s_inf / mu2).val;
    return 1.0 / (2.0 * s) * (ln * ln);
  }
}

complex<double> c0_massive(double s, double m2) {
  c_inf m2_inf(m2, -1);
  c_inf beta = sqrt(1.0 - 4.0 * m2_inf / s);
  c_inf lambda_plus = 0.5 * (1.0 + beta);
  c_inf lambda_minus = 0.5 * (1.0 - beta);
  complex<double> ln_rat = log(-lambda_minus / lambda_plus).val;
  return 1.0 / (2.0 * s) * ln_rat * ln_rat;
}

complex<double> d0_massless(double s, double t, double mu2) {
  c_inf s_inf(s, 1);
  c_inf t_inf(t, 1);

  complex<double> ln_rat = log(s_inf / t_inf).val;
  complex<double> ln_1 = log(-s_inf / mu2).val;
  complex<double> ln_2 = log(-t_inf / mu2).val;

  return (-ln_rat * ln_rat - M_PI * M_PI + ln_1 * ln_1 + ln_2 * ln_2) / (s * t);
}

complex<double> d0_massive(double s, double t, double m2) {
  c_inf m2_inf(m2, -1);
  double u = -s - t;
  c_inf x = s / (4 * m2_inf);
  c_inf y = t / (4 * m2_inf);

  c_inf beta_x = sqrt(1.0 - 1.0 / x);
  c_inf beta_y = sqrt(1.0 - 1.0 / y);
  c_inf beta_xy = sqrt(1.0 - 1.0 / x - 1.0 / y);

  complex<double> pref = (3.0 / (4.0 * x * y * beta_xy)).val / (6.0 * m2 * m2);

  complex<double> term1 = (2.0 * log((beta_xy + beta_x) / (beta_xy + beta_y)) *
                               log((beta_xy + beta_x) / (beta_xy + beta_y)) +
                           log((beta_xy - beta_x) / (beta_xy + beta_x)) *
                               log((beta_xy - beta_y) / (beta_xy + beta_y)) -
                           M_PI * M_PI / 2)
                              .val;
  complex<double> term2 = (2.0 * dilog((beta_x - 1.0) / (beta_xy + beta_x)) -
                           2.0 * dilog(-(beta_xy - beta_x) / (beta_x + 1)) -
                           log((beta_x + 1.0) / (beta_xy + beta_x)) *
                               log((beta_x + 1.0) / (beta_xy + beta_x)))
                              .val;
  complex<double> term3 = (2.0 * dilog((beta_y - 1.0) / (beta_xy + beta_y)) -
                           2.0 * dilog(-(beta_xy - beta_y) / (beta_y + 1)) -
                           log((beta_y + 1.0) / (beta_xy + beta_y)) *
                               log((beta_y + 1.0) / (beta_xy + beta_y)))
                              .val;

  return pref * (term1 + term2 + term3);
}

void fill_ff_params_massless(double s, double t, double u,
                             polynomial_parameters *p_a

) {

  p_a->xs = s;
  p_a->xt = t;
  p_a->xu = u;

  p_a->a02 = 0.0;

  p_a->b02s = b0_massless(s, mu2);
  p_a->b02t = b0_massless(t, mu2);
  p_a->b02u = b0_massless(u, mu2);
  p_a->b020 = 0.0;

  p_a->c02s = c0_massless(s, mu2);
  p_a->c02t = c0_massless(t, mu2);
  p_a->c02u = c0_massless(u, mu2);

  p_a->d02su = d0_massless(s, u, mu2);
  p_a->d02st = d0_massless(s, t, mu2);
  p_a->d02tu = d0_massless(t, u, mu2);
}

void fill_ff_params_massive(double s, double t, double u, double m2,
                            polynomial_parameters *p_a) {

  p_a->xs = s / m2;
  p_a->xt = t / m2;
  p_a->xu = u / m2;

  double m4 = m2 * m2;

  p_a->a02 = a0_massive(m2, mu2) / m2;

  p_a->b02s = b0_massive(s, m2, mu2);
  p_a->b02t = b0_massive(t, m2, mu2);
  p_a->b02u = b0_massive(u, m2, mu2);
  p_a->b020 = b0_massive_0(m2, mu2);

  p_a->c02s = c0_massive(s, m2) * m2;
  p_a->c02t = c0_massive(t, m2) * m2;
  p_a->c02u = c0_massive(u, m2) * m2;

  p_a->d02su = d0_massive(s, u, m2) * m4;
  p_a->d02st = d0_massive(s, t, m2) * m4;
  p_a->d02tu = d0_massive(t, u, m2) * m4;
}

complex<double> astu_polynomial(polynomial_parameters p_a) { 
	complex<double> astu;
	complex<double> Z[18];
	 

Z[0]=3.0*p_a.c02u;
Z[1]=Z[0] + 3.0*p_a.b02u;
Z[2]=16.0*p_a.d02tu;
Z[3]=3.0*p_a.b02s;
Z[4]=24.0*p_a.c02s - Z[2] - 16.0*p_a.d02su + Z[3] + 2.0 - Z[1];
Z[5]=4.0*p_a.c02s;
Z[6]=9.0*p_a.c02u;
Z[7]= - Z[5] - Z[6];
Z[8]=3.0*p_a.xu;
Z[7]=Z[7]*Z[8];
Z[4]=4.0*Z[4] + Z[7];
Z[7]=p_a.xu*p_a.xu;
Z[4]=Z[4]*Z[7];
Z[9]=3.0*p_a.c02t;
Z[10]=4.0*p_a.d02st;
Z[11]= - Z[9] - Z[10];
Z[12]=p_a.c02t - p_a.c02s;
Z[13]=p_a.xu*p_a.d02tu;
Z[12]=3.0*Z[12] - 5.0*Z[13];
Z[12]=p_a.xu*Z[12];
Z[11]=8.0*Z[11] + Z[12];
Z[11]=p_a.xt*Z[11];
Z[12]=9.0*p_a.c02s;
Z[2]= - Z[12] + Z[2] - 9.0*p_a.c02t;
Z[2]=p_a.xu*Z[2];
Z[2]=4.0*p_a.b02t + Z[2];
Z[2]=p_a.xu*Z[2];
Z[2]=Z[2] + Z[11];
Z[11]=2.0*p_a.xt;
Z[2]=Z[2]*Z[11];
Z[2]=Z[4] + Z[2];
Z[2]=p_a.xt*Z[2];
Z[4]=Z[3] - p_a.b02u;
Z[6]=Z[6] + Z[4];
Z[13]=p_a.xu*p_a.xu*p_a.xu;
Z[6]=Z[6]*Z[13];
Z[2]=4.0*Z[6] + Z[2];
Z[2]=p_a.xt*Z[2];
Z[6]=3.0*p_a.c02s;
Z[14]= - Z[10] - Z[6];
Z[15]= - 2.0*p_a.d02st - p_a.c02s;
Z[16]=p_a.xu*p_a.d02su;
Z[15]=12.0*Z[15] - 11.0*Z[16];
Z[15]=p_a.xu*Z[15];
Z[16]= - p_a.xt + Z[8];
Z[16]=p_a.d02st*Z[16];
Z[17]=Z[10] + p_a.c02s;
Z[16]=2.0*Z[17] + Z[16];
Z[16]=Z[16]*Z[11];
Z[14]=Z[16] + 8.0*Z[14] + Z[15];
Z[14]=p_a.xt*Z[14];
Z[15]=p_a.d02su - p_a.c02s;
Z[15]=8.0*Z[15] + Z[0];
Z[16]=2.0*Z[7];
Z[15]=Z[15]*Z[16];
Z[14]=Z[15] + Z[14];
Z[14]=p_a.xt*Z[14];
Z[15]=2.0*p_a.xs;
Z[17]= - Z[15] - 10.0*p_a.xt;
Z[17]=p_a.d02su*Z[17];
Z[5]=Z[17] - Z[5] + Z[0];
Z[5]=p_a.xs*Z[7]*Z[5];
Z[17]=Z[13]*p_a.c02u;
Z[5]=Z[5] - Z[17] + Z[14];
Z[5]=Z[5]*Z[15];
Z[9]=Z[9] + 8.0*p_a.d02st;
Z[6]= - Z[6] - Z[9];
Z[14]=p_a.d02st - p_a.d02su - p_a.d02tu;
Z[14]=p_a.xu*Z[14];
Z[9]=Z[14] - p_a.c02s - Z[9];
Z[8]=Z[9]*Z[8];
Z[9]=p_a.c02t + Z[10];
Z[9]=p_a.xt*Z[9];
Z[6]=4.0*Z[9] + 8.0*Z[6] + Z[8];
Z[6]=Z[6]*Z[11];
Z[8]=5.0*p_a.d02su + 3.0*p_a.d02tu;
Z[8]=4.0*Z[8] - Z[12];
Z[8]=4.0*Z[8] + Z[0];
Z[8]=p_a.xu*Z[8];
Z[8]=8.0*p_a.b02s + Z[8];
Z[8]=p_a.xu*Z[8];
Z[6]=Z[8] + Z[6];
Z[6]=p_a.xt*Z[6];
Z[8]=p_a.xu*p_a.c02u;
Z[1]=12.0*p_a.c02s + 5.0*p_a.b02s - Z[1];
Z[1]=2.0*Z[1] - 13.0*Z[8];
Z[1]=Z[1]*Z[16];
Z[1]=Z[1] + Z[6];
Z[1]=p_a.xt*Z[1];
Z[1]=Z[1] + Z[5];
Z[1]=p_a.xs*Z[1];
Z[1]=Z[2] + Z[1];
Z[1]=p_a.xs*Z[1];
Z[2]=6.0*p_a.c02t + p_a.b02t - Z[3];
Z[2]=8.0*Z[2] - 5.0*Z[8];
Z[2]=Z[2]*Z[7];
Z[3]=p_a.xt*p_a.d02tu;
Z[5]=Z[3] - 4.0*p_a.d02tu + p_a.c02t;
Z[0]=Z[0] - 8.0*Z[5];
Z[5]=Z[7]*p_a.xt;
Z[0]=Z[0]*Z[5];
Z[0]=Z[2] + Z[0];
Z[0]=p_a.xt*Z[0];
Z[2]=6.0*p_a.c02u - Z[4];
Z[2]=Z[2]*Z[13];
Z[0]=8.0*Z[2] + Z[0];
Z[0]=Z[0]*p_a.xt*p_a.xt;
Z[0]=Z[0] + Z[1];
Z[0]=p_a.xs*Z[0];
Z[1]= - 2.0*p_a.c02t - Z[3];
Z[1]=Z[1]*Z[5];
Z[1]= - 2.0*Z[17] + Z[1];
Z[1]=Z[1]*p_a.xt*p_a.xt*p_a.xt*p_a.xt;
Z[0]=4.0*Z[1] + Z[0];

astu = Z[0];
 
return astu; }

complex<double> atsu_polynomial(polynomial_parameters p_a) { 
	complex<double> atsu;
	complex<double> Z[18];
	 

   Z[0]=3.0*p_a.c02u;
   Z[1]=Z[0] + 3.0*p_a.b02u;
   Z[2]=16.0*p_a.d02su;
   Z[3]=3.0*p_a.b02t;
   Z[4]=24.0*p_a.c02t - Z[2] - 16.0*p_a.d02tu + Z[3] + 2.0 - Z[1];
   Z[5]=4.0*p_a.c02t;
   Z[6]=9.0*p_a.c02u;
   Z[7]= - Z[5] - Z[6];
   Z[8]=3.0*p_a.xu;
   Z[7]=Z[7]*Z[8];
   Z[4]=4.0*Z[4] + Z[7];
   Z[7]=p_a.xu*p_a.xu;
   Z[4]=Z[4]*Z[7];
   Z[9]=3.0*p_a.c02s;
   Z[10]=4.0*p_a.d02st;
   Z[11]= - Z[9] - Z[10];
   Z[12]=p_a.c02s - p_a.c02t;
   Z[13]=p_a.xu*p_a.d02su;
   Z[12]=3.0*Z[12] - 5.0*Z[13];
   Z[12]=p_a.xu*Z[12];
   Z[11]=8.0*Z[11] + Z[12];
   Z[11]=p_a.xs*Z[11];
   Z[12]=9.0*p_a.c02t;
   Z[2]= - Z[12] + Z[2] - 9.0*p_a.c02s;
   Z[2]=p_a.xu*Z[2];
   Z[2]=4.0*p_a.b02s + Z[2];
   Z[2]=p_a.xu*Z[2];
   Z[2]=Z[2] + Z[11];
   Z[11]=2.0*p_a.xs;
   Z[2]=Z[2]*Z[11];
   Z[2]=Z[4] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[4]=Z[3] - p_a.b02u;
   Z[6]=Z[6] + Z[4];
   Z[13]=p_a.xu*p_a.xu*p_a.xu;
   Z[6]=Z[6]*Z[13];
   Z[2]=4.0*Z[6] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[6]=3.0*p_a.c02t;
   Z[14]= - Z[10] - Z[6];
   Z[15]= - 2.0*p_a.d02st - p_a.c02t;
   Z[16]=p_a.xu*p_a.d02tu;
   Z[15]=12.0*Z[15] - 11.0*Z[16];
   Z[15]=p_a.xu*Z[15];
   Z[16]= - p_a.xs + Z[8];
   Z[16]=p_a.d02st*Z[16];
   Z[17]=Z[10] + p_a.c02t;
   Z[16]=2.0*Z[17] + Z[16];
   Z[16]=Z[16]*Z[11];
   Z[14]=Z[16] + 8.0*Z[14] + Z[15];
   Z[14]=p_a.xs*Z[14];
   Z[15]=p_a.d02tu - p_a.c02t;
   Z[15]=8.0*Z[15] + Z[0];
   Z[16]=2.0*Z[7];
   Z[15]=Z[15]*Z[16];
   Z[14]=Z[15] + Z[14];
   Z[14]=p_a.xs*Z[14];
   Z[15]=2.0*p_a.xt;
   Z[17]= - Z[15] - 10.0*p_a.xs;
   Z[17]=p_a.d02tu*Z[17];
   Z[5]=Z[17] - Z[5] + Z[0];
   Z[5]=p_a.xt*Z[7]*Z[5];
   Z[17]=Z[13]*p_a.c02u;
   Z[5]=Z[5] - Z[17] + Z[14];
   Z[5]=Z[5]*Z[15];
   Z[9]=Z[9] + 8.0*p_a.d02st;
   Z[6]= - Z[6] - Z[9];
   Z[14]=p_a.d02st - p_a.d02tu - p_a.d02su;
   Z[14]=p_a.xu*Z[14];
   Z[9]=Z[14] - p_a.c02t - Z[9];
   Z[8]=Z[9]*Z[8];
   Z[9]=p_a.c02s + Z[10];
   Z[9]=p_a.xs*Z[9];
   Z[6]=4.0*Z[9] + 8.0*Z[6] + Z[8];
   Z[6]=Z[6]*Z[11];
   Z[8]=5.0*p_a.d02tu + 3.0*p_a.d02su;
   Z[8]=4.0*Z[8] - Z[12];
   Z[8]=4.0*Z[8] + Z[0];
   Z[8]=p_a.xu*Z[8];
   Z[8]=8.0*p_a.b02t + Z[8];
   Z[8]=p_a.xu*Z[8];
   Z[6]=Z[8] + Z[6];
   Z[6]=p_a.xs*Z[6];
   Z[8]=p_a.xu*p_a.c02u;
   Z[1]=12.0*p_a.c02t + 5.0*p_a.b02t - Z[1];
   Z[1]=2.0*Z[1] - 13.0*Z[8];
   Z[1]=Z[1]*Z[16];
   Z[1]=Z[1] + Z[6];
   Z[1]=p_a.xs*Z[1];
   Z[1]=Z[1] + Z[5];
   Z[1]=p_a.xt*Z[1];
   Z[1]=Z[2] + Z[1];
   Z[1]=p_a.xt*Z[1];
   Z[2]=6.0*p_a.c02s + p_a.b02s - Z[3];
   Z[2]=8.0*Z[2] - 5.0*Z[8];
   Z[2]=Z[2]*Z[7];
   Z[3]=p_a.xs*p_a.d02su;
   Z[5]=Z[3] - 4.0*p_a.d02su + p_a.c02s;
   Z[0]=Z[0] - 8.0*Z[5];
   Z[5]=Z[7]*p_a.xs;
   Z[0]=Z[0]*Z[5];
   Z[0]=Z[2] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[2]=6.0*p_a.c02u - Z[4];
   Z[2]=Z[2]*Z[13];
   Z[0]=8.0*Z[2] + Z[0];
   Z[0]=Z[0]*p_a.xs*p_a.xs;
   Z[0]=Z[0] + Z[1];
   Z[0]=p_a.xt*Z[0];
   Z[1]= - 2.0*p_a.c02s - Z[3];
   Z[1]=Z[1]*Z[5];
   Z[1]= - 2.0*Z[17] + Z[1];
   Z[1]=Z[1]*p_a.xs*p_a.xs*p_a.xs*p_a.xs;
   Z[0]=4.0*Z[1] + Z[0];

atsu = Z[0];
 
return atsu; }

complex<double> aust_polynomial(polynomial_parameters p_a) { 
	complex<double> aust;
	complex<double> Z[17];
	 

   Z[0]=3.0*p_a.c02t;
   Z[1]=Z[0] + 3.0*p_a.b02t;
   Z[2]=16.0*p_a.d02st;
   Z[3]=3.0*p_a.b02u;
   Z[4]=24.0*p_a.c02u - Z[2] - 16.0*p_a.d02tu + Z[3] + 2.0 - Z[1];
   Z[5]=9.0*p_a.c02t;
   Z[6]=4.0*p_a.c02u;
   Z[7]= - Z[5] - Z[6];
   Z[8]=3.0*p_a.xt;
   Z[7]=Z[7]*Z[8];
   Z[4]=4.0*Z[4] + Z[7];
   Z[7]=p_a.xt*p_a.xt;
   Z[4]=Z[4]*Z[7];
   Z[9]=3.0*p_a.c02s;
   Z[10]=4.0*p_a.d02su;
   Z[11]= - Z[9] - Z[10];
   Z[12]=p_a.c02s - p_a.c02u;
   Z[13]=p_a.xt*p_a.d02st;
   Z[12]=3.0*Z[12] - 5.0*Z[13];
   Z[12]=p_a.xt*Z[12];
   Z[11]=8.0*Z[11] + Z[12];
   Z[11]=p_a.xs*Z[11];
   Z[12]=p_a.c02u + p_a.c02s;
   Z[2]=Z[2] - 9.0*Z[12];
   Z[2]=p_a.xt*Z[2];
   Z[2]=4.0*p_a.b02s + Z[2];
   Z[2]=p_a.xt*Z[2];
   Z[2]=Z[2] + Z[11];
   Z[11]=2.0*p_a.xs;
   Z[2]=Z[2]*Z[11];
   Z[2]=Z[4] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[4]=Z[3] - p_a.b02t;
   Z[5]=Z[5] + Z[4];
   Z[12]=p_a.xt*p_a.xt*p_a.xt;
   Z[5]=Z[5]*Z[12];
   Z[2]=4.0*Z[5] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[5]=3.0*p_a.c02u;
   Z[13]= - Z[10] - Z[5];
   Z[14]= - 2.0*p_a.d02su - p_a.c02u;
   Z[15]=p_a.xt*p_a.d02tu;
   Z[14]=12.0*Z[14] - 11.0*Z[15];
   Z[14]=p_a.xt*Z[14];
   Z[15]= - p_a.xs + Z[8];
   Z[15]=p_a.d02su*Z[15];
   Z[16]=Z[10] + p_a.c02u;
   Z[15]=2.0*Z[16] + Z[15];
   Z[15]=Z[15]*Z[11];
   Z[13]=Z[15] + 8.0*Z[13] + Z[14];
   Z[13]=p_a.xs*Z[13];
   Z[14]=p_a.c02u - p_a.d02tu;
   Z[14]=Z[0] - 8.0*Z[14];
   Z[15]=2.0*Z[7];
   Z[14]=Z[14]*Z[15];
   Z[13]=Z[14] + Z[13];
   Z[13]=p_a.xs*Z[13];
   Z[14]=2.0*p_a.xu;
   Z[16]= - Z[14] - 10.0*p_a.xs;
   Z[16]=p_a.d02tu*Z[16];
   Z[6]=Z[16] + Z[0] - Z[6];
   Z[6]=p_a.xu*Z[7]*Z[6];
   Z[16]=Z[12]*p_a.c02t;
   Z[6]=Z[6] - Z[16] + Z[13];
   Z[6]=Z[6]*Z[14];
   Z[9]=Z[9] + 8.0*p_a.d02su;
   Z[5]= - Z[5] - Z[9];
   Z[13]=p_a.d02su - p_a.d02tu - p_a.d02st;
   Z[13]=p_a.xt*Z[13];
   Z[9]=Z[13] - p_a.c02u - Z[9];
   Z[8]=Z[9]*Z[8];
   Z[9]=p_a.c02s + Z[10];
   Z[9]=p_a.xs*Z[9];
   Z[5]=4.0*Z[9] + 8.0*Z[5] + Z[8];
   Z[5]=Z[5]*Z[11];
   Z[8]=5.0*p_a.d02tu + 3.0*p_a.d02st;
   Z[8]= - 36.0*p_a.c02u + 16.0*Z[8] + Z[0];
   Z[8]=p_a.xt*Z[8];
   Z[8]=8.0*p_a.b02u + Z[8];
   Z[8]=p_a.xt*Z[8];
   Z[5]=Z[8] + Z[5];
   Z[5]=p_a.xs*Z[5];
   Z[8]=p_a.xt*p_a.c02t;
   Z[1]=12.0*p_a.c02u + 5.0*p_a.b02u - Z[1];
   Z[1]=2.0*Z[1] - 13.0*Z[8];
   Z[1]=Z[1]*Z[15];
   Z[1]=Z[1] + Z[5];
   Z[1]=p_a.xs*Z[1];
   Z[1]=Z[1] + Z[6];
   Z[1]=p_a.xu*Z[1];
   Z[1]=Z[2] + Z[1];
   Z[1]=p_a.xu*Z[1];
   Z[2]=6.0*p_a.c02s + p_a.b02s - Z[3];
   Z[2]=8.0*Z[2] - 5.0*Z[8];
   Z[2]=Z[2]*Z[7];
   Z[3]=p_a.xs*p_a.d02st;
   Z[5]=Z[3] - 4.0*p_a.d02st + p_a.c02s;
   Z[0]=Z[0] - 8.0*Z[5];
   Z[5]=Z[7]*p_a.xs;
   Z[0]=Z[0]*Z[5];
   Z[0]=Z[2] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[2]=6.0*p_a.c02t - Z[4];
   Z[2]=Z[2]*Z[12];
   Z[0]=8.0*Z[2] + Z[0];
   Z[0]=Z[0]*p_a.xs*p_a.xs;
   Z[0]=Z[0] + Z[1];
   Z[0]=p_a.xu*Z[0];
   Z[1]= - 2.0*p_a.c02s - Z[3];
   Z[1]=Z[1]*Z[5];
   Z[1]= - 2.0*Z[16] + Z[1];
   Z[1]=Z[1]*p_a.xs*p_a.xs*p_a.xs*p_a.xs;
   Z[0]=4.0*Z[1] + Z[0];

aust = Z[0];
 
return aust; }

complex<double> bstu_polynomial(polynomial_parameters p_a) { 
	complex<double> bstu;
	complex<double> Z[29];
	 

   Z[0]=2.0*p_a.d02st;
   Z[1]=3.0*p_a.c02t;
   Z[2]=Z[0] + Z[1];
   Z[3]=3.0*p_a.xu;
   Z[4]=p_a.c02s + p_a.d02st;
   Z[5]=p_a.c02t - Z[4];
   Z[5]=Z[5]*Z[3];
   Z[2]=4.0*Z[2] + Z[5];
   Z[5]=p_a.d02st*Z[3];
   Z[6]= - p_a.xt*Z[0];
   Z[4]=Z[6] + 4.0*Z[4] + Z[5];
   Z[4]=p_a.xt*Z[4];
   Z[2]=2.0*Z[2] + Z[4];
   Z[4]=2.0*p_a.xt;
   Z[2]=Z[2]*Z[4];
   Z[5]=p_a.xu*p_a.d02su;
   Z[6]=37.0*Z[5];
   Z[7]=26.0*p_a.d02su - 33.0*p_a.c02s;
   Z[7]=2.0*Z[7] + Z[6];
   Z[8]=p_a.xu*p_a.xu;
   Z[7]=Z[7]*Z[8];
   Z[2]=Z[7] + Z[2];
   Z[2]=p_a.xt*Z[2];
   Z[7]=p_a.b02s + p_a.b02u;
   Z[9]=3.0*p_a.c02s;
   Z[10]=Z[9] + Z[7];
   Z[11]=2.0*p_a.d02su;
   Z[12]=Z[11] - Z[9];
   Z[13]=7.0*p_a.c02u;
   Z[12]=4.0*Z[12] - Z[13];
   Z[12]=5.0*Z[12] + 7.0*Z[5];
   Z[12]=p_a.xu*Z[12];
   Z[10]=8.0*Z[10] + Z[12];
   Z[12]=2.0*Z[8];
   Z[10]=Z[10]*Z[12];
   Z[2]=Z[10] + Z[2];
   Z[2]=p_a.xt*Z[2];
   Z[10]=3.0*p_a.c02u;
   Z[14]= - Z[5] - 4.0*p_a.c02s - Z[10];
   Z[15]=p_a.xu*p_a.xu*p_a.xu;
   Z[16]=4.0*Z[15];
   Z[14]=Z[14]*Z[16];
   Z[17]=10.0*p_a.d02su;
   Z[18]= - 4.0*Z[5] + Z[17] - 19.0*p_a.c02s;
   Z[12]=Z[18]*Z[12];
   Z[18]=Z[8]*p_a.xt;
   Z[19]=Z[18]*p_a.d02su;
   Z[12]=Z[12] - 7.0*Z[19];
   Z[12]=p_a.xt*Z[12];
   Z[20]= - p_a.c02s - Z[5];
   Z[21]=p_a.xs*p_a.d02su;
   Z[20]= - 4.0*Z[21] + 8.0*Z[20];
   Z[20]=Z[8]*Z[20];
   Z[19]= - 11.0*Z[19] + Z[20];
   Z[19]=p_a.xs*Z[19];
   Z[12]=Z[19] + Z[14] + Z[12];
   Z[12]=p_a.xs*Z[12];
   Z[14]=2.0*p_a.c02u;
   Z[19]=Z[14] + p_a.c02s;
   Z[20]=p_a.xu*p_a.xu*p_a.xu*p_a.xu;
   Z[21]=Z[19]*Z[20];
   Z[2]=2.0*Z[12] - 16.0*Z[21] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[12]=9.0*p_a.c02t;
   Z[21]= - Z[9] + 4.0*p_a.d02st + Z[12];
   Z[22]=p_a.xu + 1.0;
   Z[22]=p_a.d02st*Z[22];
   Z[22]=p_a.c02t + Z[22];
   Z[22]=Z[22]*Z[3];
   Z[23]=Z[9] + 7.0*p_a.d02st + p_a.c02t;
   Z[24]=p_a.xt*p_a.d02st;
   Z[23]=2.0*Z[23] - 3.0*Z[24];
   Z[23]=p_a.xt*Z[23];
   Z[21]=Z[23] + 4.0*Z[21] + Z[22];
   Z[21]=p_a.xt*Z[21];
   Z[22]=Z[1] - Z[9];
   Z[23]=p_a.b02t - p_a.b02s;
   Z[24]=Z[22] - Z[23];
   Z[25]= - Z[0] + p_a.c02t;
   Z[25]=2.0*Z[25] - 5.0*p_a.c02s;
   Z[25]=6.0*Z[25] + Z[6];
   Z[25]=p_a.xu*Z[25];
   Z[24]=16.0*Z[24] + Z[25];
   Z[24]=p_a.xu*Z[24];
   Z[21]=Z[24] + 4.0*Z[21];
   Z[21]=p_a.xt*Z[21];
   Z[24]=5.0*p_a.b02s;
   Z[25]=Z[24] + 3.0*p_a.b02t;
   Z[26]=8.0*p_a.d02su;
   Z[27]=24.0*p_a.c02s - Z[26] + 9.0*p_a.b02u + Z[25];
   Z[28]=21.0*p_a.c02u;
   Z[6]=Z[6] - 46.0*p_a.c02s - Z[28];
   Z[6]=p_a.xu*Z[6];
   Z[6]=4.0*Z[27] + Z[6];
   Z[6]=p_a.xu*Z[6];
   Z[6]= - 32.0*p_a.b02t + Z[6];
   Z[6]=p_a.xu*Z[6];
   Z[6]=Z[6] + Z[21];
   Z[6]=p_a.xt*Z[6];
   Z[21]=9.0*p_a.c02u;
   Z[27]=6.0*p_a.c02s + Z[7];
   Z[27]=2.0*Z[27] + Z[21];
   Z[17]= - Z[28] + Z[17] - 7.0*p_a.c02s;
   Z[17]=p_a.xu*Z[17];
   Z[17]=2.0*Z[27] + Z[17];
   Z[16]=Z[17]*Z[16];
   Z[6]=Z[16] + Z[6];
   Z[6]=p_a.xt*Z[6];
   Z[16]=p_a.xu*p_a.xu*p_a.xu*p_a.xu*p_a.xu;
   Z[17]=p_a.c02u*Z[16];
   Z[2]=Z[2] - 8.0*Z[17] + Z[6];
   Z[2]=p_a.xs*Z[2];
   Z[6]= - 2.0*Z[23] + Z[22];
   Z[17]=Z[0] - p_a.d02tu;
   Z[11]=Z[11] + Z[17];
   Z[11]=p_a.xu*Z[11];
   Z[17]=Z[17] + p_a.c02s;
   Z[17]= - p_a.c02t - 2.0*Z[17];
   Z[11]=2.0*Z[17] + Z[11];
   Z[11]=Z[11]*Z[3];
   Z[6]=16.0*Z[6] + Z[11];
   Z[6]=p_a.xu*Z[6];
   Z[11]=p_a.d02st + Z[1];
   Z[11]=2.0*Z[11] - Z[9];
   Z[0]=Z[0] + p_a.c02t;
   Z[17]=p_a.d02st - p_a.d02tu;
   Z[17]=p_a.xu*Z[17];
   Z[17]=2.0*Z[0] + Z[17];
   Z[17]=Z[17]*Z[3];
   Z[0]=p_a.xt*Z[0];
   Z[0]=12.0*Z[0] + 8.0*Z[11] + Z[17];
   Z[0]=Z[0]*Z[4];
   Z[0]=Z[6] + Z[0];
   Z[0]=p_a.xt*Z[0];
   Z[4]= - Z[9] + 3.0*p_a.d02tu - 23.0*p_a.d02su;
   Z[4]=3.0*Z[5] + 2.0*Z[4] + 5.0*p_a.c02u;
   Z[4]=p_a.xu*Z[4];
   Z[5]=Z[24] + 1.0 + 6.0*p_a.b02u;
   Z[5]=12.0*p_a.c02s + 2.0*Z[5] - 5.0*p_a.b02t;
   Z[4]=2.0*Z[5] + Z[4];
   Z[4]=p_a.xu*Z[4];
   Z[5]= - p_a.b02s - p_a.b02t;
   Z[4]=16.0*Z[5] + Z[4];
   Z[4]=p_a.xu*Z[4];
   Z[0]=2.0*Z[4] + Z[0];
   Z[0]=p_a.xt*Z[0];
   Z[4]=2.0*p_a.b02t;
   Z[5]= - Z[4] - Z[7];
   Z[6]=18.0*p_a.c02u + 36.0*p_a.c02s - 16.0*p_a.d02su + p_a.b02u + Z[25];
   Z[9]= - 92.0*p_a.d02su - p_a.c02u;
   Z[9]=p_a.xu*Z[9];
   Z[6]=4.0*Z[6] + Z[9];
   Z[6]=p_a.xu*Z[6];
   Z[5]=32.0*Z[5] + Z[6];
   Z[5]=Z[5]*Z[8];
   Z[0]=Z[5] + Z[0];
   Z[0]=p_a.xt*Z[0];
   Z[5]=p_a.xu*p_a.c02u;
   Z[6]=24.0*Z[19] - 7.0*Z[5];
   Z[6]=Z[6]*Z[20];
   Z[0]=2.0*Z[6] + Z[0];
   Z[0]=p_a.xt*Z[0];
   Z[0]=Z[0] + Z[2];
   Z[0]=p_a.xs*Z[0];
   Z[2]=8.0*p_a.d02tu;
   Z[6]=24.0*p_a.c02u;
   Z[9]=3.0*p_a.b02u + p_a.b02s;
   Z[9]=Z[6] - Z[26] - Z[2] + 3.0*Z[9] - 4.0*p_a.b02t;
   Z[9]=4.0*Z[9] - 27.0*Z[5];
   Z[9]=p_a.xu*Z[9];
   Z[11]= - p_a.b02t - p_a.b02u - 2.0*p_a.b02s;
   Z[9]=32.0*Z[11] + Z[9];
   Z[9]=Z[9]*Z[8];
   Z[2]=12.0*p_a.c02t - Z[2] + Z[4] + 7.0*p_a.b02u + 3.0*p_a.b02s;
   Z[4]=4.0*p_a.d02tu;
   Z[11]= - Z[4] - Z[12];
   Z[3]= - p_a.d02tu*Z[3];
   Z[3]=Z[3] + 2.0*Z[11] - 39.0*p_a.c02u;
   Z[3]=p_a.xu*Z[3];
   Z[2]=4.0*Z[2] + Z[3];
   Z[2]=p_a.xu*Z[2];
   Z[2]= - 32.0*p_a.b02s + Z[2];
   Z[2]=p_a.xu*Z[2];
   Z[3]=Z[4] - Z[1];
   Z[11]=p_a.xu*p_a.d02tu;
   Z[3]=2.0*Z[3] - Z[11];
   Z[3]=Z[3]*Z[8];
   Z[12]=Z[18]*p_a.d02tu;
   Z[3]=5.0*Z[3] - 14.0*Z[12];
   Z[3]=p_a.xt*Z[3];
   Z[2]=Z[2] + Z[3];
   Z[2]=p_a.xt*Z[2];
   Z[2]=Z[9] + Z[2];
   Z[2]=p_a.xt*Z[2];
   Z[3]= - p_a.b02t - Z[7];
   Z[7]= - 4.0*p_a.d02su + Z[21];
   Z[7]=p_a.xu*Z[7];
   Z[3]=4.0*Z[3] + Z[7];
   Z[3]=Z[3]*Z[15];
   Z[2]=8.0*Z[3] + Z[2];
   Z[2]=p_a.xt*Z[2];
   Z[3]=Z[16]*Z[6];
   Z[2]=Z[3] + Z[2];
   Z[2]=p_a.xt*Z[2];
   Z[0]=Z[2] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[2]=Z[1] + p_a.b02t;
   Z[3]=Z[10] + Z[2];
   Z[1]= - Z[1] - Z[13];
   Z[1]=p_a.xu*Z[1];
   Z[1]=4.0*Z[3] + Z[1];
   Z[1]=Z[1]*Z[15];
   Z[3]=p_a.d02tu - p_a.c02t;
   Z[3]=10.0*Z[3] - Z[11];
   Z[3]=Z[3]*Z[8];
   Z[4]= - Z[18]*Z[4];
   Z[3]=Z[3] + Z[4];
   Z[3]=p_a.xt*Z[3];
   Z[4]=p_a.c02u - p_a.d02tu;
   Z[4]= - 7.0*p_a.c02t - 10.0*Z[4];
   Z[4]=p_a.xu*Z[4];
   Z[2]=4.0*Z[2] + Z[4];
   Z[2]=Z[2]*Z[8];
   Z[2]=Z[2] + Z[3];
   Z[2]=p_a.xt*Z[2];
   Z[1]=Z[1] + Z[2];
   Z[1]=p_a.xt*Z[1];
   Z[2]=4.0*p_a.c02u - Z[5];
   Z[2]=Z[2]*Z[20];
   Z[1]=3.0*Z[2] + Z[1];
   Z[1]=Z[1]*p_a.xt*p_a.xt*p_a.xt;
   Z[0]=4.0*Z[1] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[1]= - p_a.c02t - p_a.c02u;
   Z[1]=Z[1]*Z[15];
   Z[2]= - 2.0*p_a.c02t - Z[11];
   Z[2]=Z[2]*Z[8];
   Z[2]=Z[2] - Z[12];
   Z[2]=p_a.xt*Z[2];
   Z[1]=2.0*Z[1] + Z[2];
   Z[1]=p_a.xt*Z[1];
   Z[2]= - Z[20]*Z[14];
   Z[1]=Z[2] + Z[1];
   Z[1]=Z[1]*p_a.xt*p_a.xt*p_a.xt*p_a.xt*p_a.xt;
   Z[0]=8.0*Z[1] + Z[0];

bstu = Z[0];
 
return bstu; }

complex<double> btsu_polynomial(polynomial_parameters p_a) { 
	complex<double> btsu;
	complex<double> Z[28];
	 

   Z[0]=2.0*p_a.d02st;
   Z[1]=3.0*p_a.c02s;
   Z[2]=Z[0] + Z[1];
   Z[3]=3.0*p_a.xu;
   Z[4]=p_a.c02t + p_a.d02st;
   Z[5]=p_a.c02s - Z[4];
   Z[5]=Z[5]*Z[3];
   Z[2]=4.0*Z[2] + Z[5];
   Z[5]=p_a.d02st*Z[3];
   Z[6]= - p_a.xs*Z[0];
   Z[4]=Z[6] + 4.0*Z[4] + Z[5];
   Z[4]=p_a.xs*Z[4];
   Z[2]=2.0*Z[2] + Z[4];
   Z[4]=2.0*p_a.xs;
   Z[2]=Z[2]*Z[4];
   Z[5]=p_a.xu*p_a.d02tu;
   Z[6]=37.0*Z[5];
   Z[7]=26.0*p_a.d02tu - 33.0*p_a.c02t;
   Z[7]=2.0*Z[7] + Z[6];
   Z[8]=p_a.xu*p_a.xu;
   Z[7]=Z[7]*Z[8];
   Z[2]=Z[7] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[7]=p_a.b02t + p_a.b02u;
   Z[9]=3.0*p_a.c02t;
   Z[10]=Z[9] + Z[7];
   Z[11]=2.0*p_a.d02tu;
   Z[12]=Z[11] - Z[9];
   Z[13]=7.0*p_a.c02u;
   Z[12]=4.0*Z[12] - Z[13];
   Z[12]=5.0*Z[12] + 7.0*Z[5];
   Z[12]=p_a.xu*Z[12];
   Z[10]=8.0*Z[10] + Z[12];
   Z[12]=2.0*Z[8];
   Z[10]=Z[10]*Z[12];
   Z[2]=Z[10] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[10]=3.0*p_a.c02u;
   Z[14]= - Z[5] - 4.0*p_a.c02t - Z[10];
   Z[15]=p_a.xu*p_a.xu*p_a.xu;
   Z[16]=4.0*Z[15];
   Z[14]=Z[14]*Z[16];
   Z[17]=10.0*p_a.d02tu;
   Z[18]= - 4.0*Z[5] + Z[17] - 19.0*p_a.c02t;
   Z[12]=Z[18]*Z[12];
   Z[18]=Z[8]*p_a.xs;
   Z[19]=Z[18]*p_a.d02tu;
   Z[12]=Z[12] - 7.0*Z[19];
   Z[12]=p_a.xs*Z[12];
   Z[20]= - p_a.c02t - Z[5];
   Z[21]=p_a.xt*p_a.d02tu;
   Z[20]= - 4.0*Z[21] + 8.0*Z[20];
   Z[20]=Z[8]*Z[20];
   Z[19]= - 11.0*Z[19] + Z[20];
   Z[19]=p_a.xt*Z[19];
   Z[12]=Z[19] + Z[14] + Z[12];
   Z[12]=p_a.xt*Z[12];
   Z[14]=2.0*p_a.c02u;
   Z[19]=Z[14] + p_a.c02t;
   Z[20]=p_a.xu*p_a.xu*p_a.xu*p_a.xu;
   Z[21]=Z[19]*Z[20];
   Z[2]=2.0*Z[12] - 16.0*Z[21] + Z[2];
   Z[2]=p_a.xt*Z[2];
   Z[12]=9.0*p_a.c02s;
   Z[21]= - Z[9] + 4.0*p_a.d02st + Z[12];
   Z[22]=p_a.xu + 1.0;
   Z[22]=p_a.d02st*Z[22];
   Z[22]=p_a.c02s + Z[22];
   Z[22]=Z[22]*Z[3];
   Z[23]=Z[9] + 7.0*p_a.d02st + p_a.c02s;
   Z[24]=p_a.xs*p_a.d02st;
   Z[23]=2.0*Z[23] - 3.0*Z[24];
   Z[23]=p_a.xs*Z[23];
   Z[21]=Z[23] + 4.0*Z[21] + Z[22];
   Z[21]=p_a.xs*Z[21];
   Z[22]=Z[1] - Z[9];
   Z[23]=p_a.b02t - p_a.b02s;
   Z[24]=Z[22] + Z[23];
   Z[25]= - Z[0] + p_a.c02s;
   Z[25]=2.0*Z[25] - 5.0*p_a.c02t;
   Z[25]=6.0*Z[25] + Z[6];
   Z[25]=p_a.xu*Z[25];
   Z[24]=16.0*Z[24] + Z[25];
   Z[24]=p_a.xu*Z[24];
   Z[21]=Z[24] + 4.0*Z[21];
   Z[21]=p_a.xs*Z[21];
   Z[24]=5.0*p_a.b02t;
   Z[25]=8.0*p_a.d02tu;
   Z[26]=3.0*p_a.b02u + p_a.b02s;
   Z[26]=24.0*p_a.c02t - Z[25] + 3.0*Z[26] + Z[24];
   Z[27]=21.0*p_a.c02u;
   Z[6]=Z[6] - 46.0*p_a.c02t - Z[27];
   Z[6]=p_a.xu*Z[6];
   Z[6]=4.0*Z[26] + Z[6];
   Z[6]=p_a.xu*Z[6];
   Z[6]= - 32.0*p_a.b02s + Z[6];
   Z[6]=p_a.xu*Z[6];
   Z[6]=Z[6] + Z[21];
   Z[6]=p_a.xs*Z[6];
   Z[21]=9.0*p_a.c02u;
   Z[26]=6.0*p_a.c02t + Z[7];
   Z[26]=2.0*Z[26] + Z[21];
   Z[17]= - Z[27] + Z[17] - 7.0*p_a.c02t;
   Z[17]=p_a.xu*Z[17];
   Z[17]=2.0*Z[26] + Z[17];
   Z[16]=Z[17]*Z[16];
   Z[6]=Z[16] + Z[6];
   Z[6]=p_a.xs*Z[6];
   Z[16]=p_a.c02u*p_a.xu*p_a.xu*p_a.xu*p_a.xu*p_a.xu;
   Z[2]=Z[2] - 8.0*Z[16] + Z[6];
   Z[2]=p_a.xt*Z[2];
   Z[6]=2.0*Z[23] + Z[22];
   Z[17]=Z[0] - p_a.d02su;
   Z[11]=Z[11] + Z[17];
   Z[11]=p_a.xu*Z[11];
   Z[17]=Z[17] + p_a.c02t;
   Z[17]= - p_a.c02s - 2.0*Z[17];
   Z[11]=2.0*Z[17] + Z[11];
   Z[11]=Z[11]*Z[3];
   Z[6]=16.0*Z[6] + Z[11];
   Z[6]=p_a.xu*Z[6];
   Z[11]=p_a.d02st + Z[1];
   Z[11]=2.0*Z[11] - Z[9];
   Z[0]=Z[0] + p_a.c02s;
   Z[17]= - p_a.d02su + p_a.d02st;
   Z[17]=p_a.xu*Z[17];
   Z[17]=2.0*Z[0] + Z[17];
   Z[17]=Z[17]*Z[3];
   Z[0]=p_a.xs*Z[0];
   Z[0]=12.0*Z[0] + 8.0*Z[11] + Z[17];
   Z[0]=Z[0]*Z[4];
   Z[0]=Z[6] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[4]= - Z[9] + 3.0*p_a.d02su - 23.0*p_a.d02tu;
   Z[4]=3.0*Z[5] + 2.0*Z[4] + 5.0*p_a.c02u;
   Z[4]=p_a.xu*Z[4];
   Z[5]=1.0 + 6.0*p_a.b02u;
   Z[5]=12.0*p_a.c02t + 10.0*p_a.b02t + 2.0*Z[5] - 5.0*p_a.b02s;
   Z[4]=2.0*Z[5] + Z[4];
   Z[4]=p_a.xu*Z[4];
   Z[5]= - p_a.b02s - p_a.b02t;
   Z[4]=16.0*Z[5] + Z[4];
   Z[4]=p_a.xu*Z[4];
   Z[0]=2.0*Z[4] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[4]=2.0*p_a.b02s;
   Z[5]= - Z[4] - Z[7];
   Z[6]=18.0*p_a.c02u + 36.0*p_a.c02t - 16.0*p_a.d02tu + Z[24] + p_a.b02u + 3.0*p_a.b02s;
   Z[9]= - 92.0*p_a.d02tu - p_a.c02u;
   Z[9]=p_a.xu*Z[9];
   Z[6]=4.0*Z[6] + Z[9];
   Z[6]=p_a.xu*Z[6];
   Z[5]=32.0*Z[5] + Z[6];
   Z[5]=Z[5]*Z[8];
   Z[0]=Z[5] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[5]=p_a.xu*p_a.c02u;
   Z[6]=24.0*Z[19] - 7.0*Z[5];
   Z[6]=Z[6]*Z[20];
   Z[0]=2.0*Z[6] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[0]=Z[0] + Z[2];
   Z[0]=p_a.xt*Z[0];
   Z[2]= - 8.0*p_a.d02su + 3.0*p_a.b02t;
   Z[4]=12.0*p_a.c02s + 7.0*p_a.b02u + Z[4] + Z[2];
   Z[6]=4.0*p_a.d02su;
   Z[9]= - Z[6] - Z[12];
   Z[3]= - p_a.d02su*Z[3];
   Z[3]=Z[3] + 2.0*Z[9] - 39.0*p_a.c02u;
   Z[3]=p_a.xu*Z[3];
   Z[3]=4.0*Z[4] + Z[3];
   Z[3]=p_a.xu*Z[3];
   Z[3]= - 32.0*p_a.b02t + Z[3];
   Z[3]=p_a.xu*Z[3];
   Z[4]=Z[6] - Z[1];
   Z[9]=p_a.xu*p_a.d02su;
   Z[4]=2.0*Z[4] - Z[9];
   Z[4]=Z[4]*Z[8];
   Z[11]=Z[18]*p_a.d02su;
   Z[4]=5.0*Z[4] - 14.0*Z[11];
   Z[4]=p_a.xs*Z[4];
   Z[3]=Z[3] + Z[4];
   Z[3]=p_a.xs*Z[3];
   Z[2]= - Z[25] + 9.0*p_a.b02u - 4.0*p_a.b02s + Z[2] + 24.0*p_a.c02u;
   Z[2]=4.0*Z[2] - 27.0*Z[5];
   Z[2]=p_a.xu*Z[2];
   Z[4]= - 2.0*p_a.b02t - p_a.b02u - p_a.b02s;
   Z[2]=32.0*Z[4] + Z[2];
   Z[2]=Z[2]*Z[8];
   Z[2]=Z[2] + Z[3];
   Z[2]=p_a.xs*Z[2];
   Z[3]= - p_a.b02s - Z[7];
   Z[4]= - 4.0*p_a.d02tu + Z[21];
   Z[4]=p_a.xu*Z[4];
   Z[3]=4.0*Z[3] + Z[4];
   Z[3]=Z[3]*Z[15];
   Z[2]=8.0*Z[3] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[2]=24.0*Z[16] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[0]=Z[2] + Z[0];
   Z[0]=p_a.xt*Z[0];
   Z[2]=Z[1] + p_a.b02s;
   Z[3]=Z[10] + Z[2];
   Z[1]= - Z[1] - Z[13];
   Z[1]=p_a.xu*Z[1];
   Z[1]=4.0*Z[3] + Z[1];
   Z[1]=Z[1]*Z[15];
   Z[3]=p_a.d02su - p_a.c02s;
   Z[3]=10.0*Z[3] - Z[9];
   Z[3]=Z[3]*Z[8];
   Z[4]= - Z[18]*Z[6];
   Z[3]=Z[3] + Z[4];
   Z[3]=p_a.xs*Z[3];
   Z[4]=p_a.c02u - p_a.d02su;
   Z[4]= - 7.0*p_a.c02s - 10.0*Z[4];
   Z[4]=p_a.xu*Z[4];
   Z[2]=4.0*Z[2] + Z[4];
   Z[2]=Z[2]*Z[8];
   Z[2]=Z[2] + Z[3];
   Z[2]=p_a.xs*Z[2];
   Z[1]=Z[1] + Z[2];
   Z[1]=p_a.xs*Z[1];
   Z[2]=4.0*p_a.c02u - Z[5];
   Z[2]=Z[2]*Z[20];
   Z[1]=3.0*Z[2] + Z[1];
   Z[1]=Z[1]*p_a.xs*p_a.xs*p_a.xs;
   Z[0]=4.0*Z[1] + Z[0];
   Z[0]=p_a.xt*Z[0];
   Z[1]= - p_a.c02s - p_a.c02u;
   Z[1]=Z[1]*Z[15];
   Z[2]= - 2.0*p_a.c02s - Z[9];
   Z[2]=Z[2]*Z[8];
   Z[2]=Z[2] - Z[11];
   Z[2]=p_a.xs*Z[2];
   Z[1]=2.0*Z[1] + Z[2];
   Z[1]=p_a.xs*Z[1];
   Z[2]= - Z[20]*Z[14];
   Z[1]=Z[2] + Z[1];
   Z[1]=Z[1]*p_a.xs*p_a.xs*p_a.xs*p_a.xs*p_a.xs;
   Z[0]=8.0*Z[1] + Z[0];

btsu = Z[0];
 
return btsu; }

complex<double> bust_polynomial(polynomial_parameters p_a) { 
	complex<double> bust;
	complex<double> Z[28];
	 

   Z[0]=2.0*p_a.d02su;
   Z[1]=3.0*p_a.c02s;
   Z[2]=Z[0] + Z[1];
   Z[3]=3.0*p_a.xt;
   Z[4]=p_a.c02u + p_a.d02su;
   Z[5]=p_a.c02s - Z[4];
   Z[5]=Z[5]*Z[3];
   Z[2]=4.0*Z[2] + Z[5];
   Z[5]=p_a.d02su*Z[3];
   Z[6]= - p_a.xs*Z[0];
   Z[4]=Z[6] + 4.0*Z[4] + Z[5];
   Z[4]=p_a.xs*Z[4];
   Z[2]=2.0*Z[2] + Z[4];
   Z[4]=2.0*p_a.xs;
   Z[2]=Z[2]*Z[4];
   Z[5]=p_a.xt*p_a.d02tu;
   Z[6]=37.0*Z[5];
   Z[7]=26.0*p_a.d02tu - 33.0*p_a.c02u;
   Z[7]=2.0*Z[7] + Z[6];
   Z[8]=p_a.xt*p_a.xt;
   Z[7]=Z[7]*Z[8];
   Z[2]=Z[7] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[7]=8.0*p_a.d02tu;
   Z[9]=7.0*p_a.c02t;
   Z[10]=12.0*p_a.c02u;
   Z[11]= - Z[10] + Z[7] - Z[9];
   Z[11]=5.0*Z[11] + 7.0*Z[5];
   Z[11]=p_a.xt*Z[11];
   Z[12]=p_a.b02u + p_a.b02t;
   Z[13]=3.0*p_a.c02u;
   Z[14]=Z[13] + Z[12];
   Z[11]=8.0*Z[14] + Z[11];
   Z[14]=2.0*Z[8];
   Z[11]=Z[11]*Z[14];
   Z[2]=Z[11] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[11]=3.0*p_a.c02t;
   Z[15]= - Z[5] - Z[11] - 4.0*p_a.c02u;
   Z[16]=p_a.xt*p_a.xt*p_a.xt;
   Z[17]=4.0*Z[16];
   Z[15]=Z[15]*Z[17];
   Z[18]=10.0*p_a.d02tu;
   Z[19]= - 4.0*Z[5] + Z[18] - 19.0*p_a.c02u;
   Z[14]=Z[19]*Z[14];
   Z[19]=Z[8]*p_a.xs;
   Z[20]=Z[19]*p_a.d02tu;
   Z[14]=Z[14] - 7.0*Z[20];
   Z[14]=p_a.xs*Z[14];
   Z[21]= - p_a.c02u - Z[5];
   Z[22]=p_a.xu*p_a.d02tu;
   Z[21]= - 4.0*Z[22] + 8.0*Z[21];
   Z[21]=Z[8]*Z[21];
   Z[20]= - 11.0*Z[20] + Z[21];
   Z[20]=p_a.xu*Z[20];
   Z[14]=Z[20] + Z[15] + Z[14];
   Z[14]=p_a.xu*Z[14];
   Z[15]=2.0*p_a.c02t;
   Z[20]=Z[15] + p_a.c02u;
   Z[21]=p_a.xt*p_a.xt*p_a.xt*p_a.xt;
   Z[22]=Z[20]*Z[21];
   Z[2]=2.0*Z[14] - 16.0*Z[22] + Z[2];
   Z[2]=p_a.xu*Z[2];
   Z[14]=9.0*p_a.c02s;
   Z[22]= - Z[13] + 4.0*p_a.d02su + Z[14];
   Z[23]=p_a.xt + 1.0;
   Z[23]=p_a.d02su*Z[23];
   Z[23]=p_a.c02s + Z[23];
   Z[23]=Z[23]*Z[3];
   Z[24]=Z[13] + 7.0*p_a.d02su + p_a.c02s;
   Z[25]=p_a.xs*p_a.d02su;
   Z[24]=2.0*Z[24] - 3.0*Z[25];
   Z[24]=p_a.xs*Z[24];
   Z[22]=Z[24] + 4.0*Z[22] + Z[23];
   Z[22]=p_a.xs*Z[22];
   Z[23]=Z[13] - Z[1];
   Z[24]=p_a.b02u - p_a.b02s;
   Z[25]= - Z[23] + Z[24];
   Z[26]= - Z[0] + p_a.c02s;
   Z[26]=2.0*Z[26] - 5.0*p_a.c02u;
   Z[26]=6.0*Z[26] + Z[6];
   Z[26]=p_a.xt*Z[26];
   Z[25]=16.0*Z[25] + Z[26];
   Z[25]=p_a.xt*Z[25];
   Z[22]=Z[25] + 4.0*Z[22];
   Z[22]=p_a.xs*Z[22];
   Z[25]=5.0*p_a.b02u;
   Z[26]=3.0*p_a.b02t + p_a.b02s;
   Z[26]=24.0*p_a.c02u - Z[7] + 3.0*Z[26] + Z[25];
   Z[27]=21.0*p_a.c02t;
   Z[6]=Z[6] - Z[27] - 46.0*p_a.c02u;
   Z[6]=p_a.xt*Z[6];
   Z[6]=4.0*Z[26] + Z[6];
   Z[6]=p_a.xt*Z[6];
   Z[6]= - 32.0*p_a.b02s + Z[6];
   Z[6]=p_a.xt*Z[6];
   Z[6]=Z[6] + Z[22];
   Z[6]=p_a.xs*Z[6];
   Z[22]=9.0*p_a.c02t;
   Z[26]=Z[10] + 2.0*Z[12] + Z[22];
   Z[18]= - 7.0*p_a.c02u + Z[18] - Z[27];
   Z[18]=p_a.xt*Z[18];
   Z[18]=2.0*Z[26] + Z[18];
   Z[17]=Z[18]*Z[17];
   Z[6]=Z[17] + Z[6];
   Z[6]=p_a.xs*Z[6];
   Z[17]=p_a.c02t*p_a.xt*p_a.xt*p_a.xt*p_a.xt*p_a.xt;
   Z[2]=Z[2] - 8.0*Z[17] + Z[6];
   Z[2]=p_a.xu*Z[2];
   Z[6]=2.0*Z[24] - Z[23];
   Z[18]=Z[0] - p_a.d02st;
   Z[23]=Z[18] + p_a.c02u;
   Z[23]= - p_a.c02s - 2.0*Z[23];
   Z[18]=2.0*p_a.d02tu + Z[18];
   Z[18]=p_a.xt*Z[18];
   Z[18]=2.0*Z[23] + Z[18];
   Z[18]=Z[18]*Z[3];
   Z[6]=16.0*Z[6] + Z[18];
   Z[6]=p_a.xt*Z[6];
   Z[18]=p_a.d02su + Z[1];
   Z[13]=2.0*Z[18] - Z[13];
   Z[0]=Z[0] + p_a.c02s;
   Z[18]=p_a.d02su - p_a.d02st;
   Z[18]=p_a.xt*Z[18];
   Z[18]=2.0*Z[0] + Z[18];
   Z[18]=Z[18]*Z[3];
   Z[0]=p_a.xs*Z[0];
   Z[0]=12.0*Z[0] + 8.0*Z[13] + Z[18];
   Z[0]=Z[0]*Z[4];
   Z[0]=Z[6] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[4]=1.0 + 6.0*p_a.b02t;
   Z[4]=Z[10] + 10.0*p_a.b02u + 2.0*Z[4] - 5.0*p_a.b02s;
   Z[6]=3.0*p_a.d02st - 23.0*p_a.d02tu;
   Z[5]=3.0*Z[5] - 6.0*p_a.c02u + 2.0*Z[6] + 5.0*p_a.c02t;
   Z[5]=p_a.xt*Z[5];
   Z[4]=2.0*Z[4] + Z[5];
   Z[4]=p_a.xt*Z[4];
   Z[5]= - p_a.b02s - p_a.b02u;
   Z[4]=16.0*Z[5] + Z[4];
   Z[4]=p_a.xt*Z[4];
   Z[0]=2.0*Z[4] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[4]=2.0*p_a.b02s;
   Z[5]= - Z[4] - Z[12];
   Z[6]=36.0*p_a.c02u + 18.0*p_a.c02t - 16.0*p_a.d02tu + Z[25] + p_a.b02t + 3.0*p_a.b02s;
   Z[10]= - 92.0*p_a.d02tu - p_a.c02t;
   Z[10]=p_a.xt*Z[10];
   Z[6]=4.0*Z[6] + Z[10];
   Z[6]=p_a.xt*Z[6];
   Z[5]=32.0*Z[5] + Z[6];
   Z[5]=Z[5]*Z[8];
   Z[0]=Z[5] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[5]=p_a.xt*p_a.c02t;
   Z[6]=24.0*Z[20] - 7.0*Z[5];
   Z[6]=Z[6]*Z[21];
   Z[0]=2.0*Z[6] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[0]=Z[0] + Z[2];
   Z[0]=p_a.xu*Z[0];
   Z[2]= - 8.0*p_a.d02st + 3.0*p_a.b02u;
   Z[4]=12.0*p_a.c02s + 7.0*p_a.b02t + Z[4] + Z[2];
   Z[6]=4.0*p_a.d02st;
   Z[10]= - Z[6] - Z[14];
   Z[3]= - p_a.d02st*Z[3];
   Z[3]=Z[3] + 2.0*Z[10] - 39.0*p_a.c02t;
   Z[3]=p_a.xt*Z[3];
   Z[3]=4.0*Z[4] + Z[3];
   Z[3]=p_a.xt*Z[3];
   Z[3]= - 32.0*p_a.b02u + Z[3];
   Z[3]=p_a.xt*Z[3];
   Z[4]=Z[6] - Z[1];
   Z[6]=p_a.xt*p_a.d02st;
   Z[4]=2.0*Z[4] - Z[6];
   Z[4]=Z[4]*Z[8];
   Z[10]=Z[19]*p_a.d02st;
   Z[4]=5.0*Z[4] - 14.0*Z[10];
   Z[4]=p_a.xs*Z[4];
   Z[3]=Z[3] + Z[4];
   Z[3]=p_a.xs*Z[3];
   Z[2]= - Z[7] + 9.0*p_a.b02t - 4.0*p_a.b02s + Z[2] + 24.0*p_a.c02t;
   Z[2]=4.0*Z[2] - 27.0*Z[5];
   Z[2]=p_a.xt*Z[2];
   Z[4]= - 2.0*p_a.b02u - p_a.b02t - p_a.b02s;
   Z[2]=32.0*Z[4] + Z[2];
   Z[2]=Z[2]*Z[8];
   Z[2]=Z[2] + Z[3];
   Z[2]=p_a.xs*Z[2];
   Z[3]= - p_a.b02s - Z[12];
   Z[4]= - 4.0*p_a.d02tu + Z[22];
   Z[4]=p_a.xt*Z[4];
   Z[3]=4.0*Z[3] + Z[4];
   Z[3]=Z[3]*Z[16];
   Z[2]=8.0*Z[3] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[2]=24.0*Z[17] + Z[2];
   Z[2]=p_a.xs*Z[2];
   Z[0]=Z[2] + Z[0];
   Z[0]=p_a.xu*Z[0];
   Z[2]=Z[1] + p_a.b02s;
   Z[3]=Z[11] + Z[2];
   Z[1]= - Z[1] - Z[9];
   Z[1]=p_a.xt*Z[1];
   Z[1]=4.0*Z[3] + Z[1];
   Z[1]=Z[1]*Z[16];
   Z[3]=p_a.c02t - p_a.d02st;
   Z[3]= - 7.0*p_a.c02s - 10.0*Z[3];
   Z[3]=p_a.xt*Z[3];
   Z[2]=4.0*Z[2] + Z[3];
   Z[2]=Z[2]*Z[8];
   Z[3]=p_a.d02st - p_a.c02s;
   Z[3]=10.0*Z[3] - Z[6];
   Z[3]=Z[3]*Z[8];
   Z[3]=Z[3] - 4.0*Z[10];
   Z[3]=p_a.xs*Z[3];
   Z[2]=Z[2] + Z[3];
   Z[2]=p_a.xs*Z[2];
   Z[1]=Z[1] + Z[2];
   Z[1]=p_a.xs*Z[1];
   Z[2]=4.0*p_a.c02t - Z[5];
   Z[2]=Z[2]*Z[21];
   Z[1]=3.0*Z[2] + Z[1];
   Z[1]=Z[1]*p_a.xs*p_a.xs*p_a.xs;
   Z[0]=4.0*Z[1] + Z[0];
   Z[0]=p_a.xu*Z[0];
   Z[1]= - p_a.c02s - p_a.c02t;
   Z[1]=Z[1]*Z[16];
   Z[2]= - 2.0*p_a.c02s - Z[6];
   Z[2]=Z[2]*Z[8];
   Z[2]=Z[2] - Z[10];
   Z[2]=p_a.xs*Z[2];
   Z[1]=2.0*Z[1] + Z[2];
   Z[1]=p_a.xs*Z[1];
   Z[2]= - Z[21]*Z[15];
   Z[1]=Z[2] + Z[1];
   Z[1]=Z[1]*p_a.xs*p_a.xs*p_a.xs*p_a.xs*p_a.xs;
   Z[0]=8.0*Z[1] + Z[0];

bust = Z[0];
 
return bust; }

complex<double> cstu_polynomial(polynomial_parameters p_a) { 
	complex<double> cstu;
	complex<double> Z[35];
	 

   Z[0]=3.0*p_a.xu;
   Z[1]= - p_a.d02su - p_a.c02u;
   Z[1]=Z[1]*Z[0];
   Z[2]=3.0*p_a.c02u;
   Z[3]=Z[2] + p_a.b02t;
   Z[4]=3.0*p_a.d02su;
   Z[5]=p_a.a02 - 1.0;
   Z[6]= - p_a.c02t - Z[5];
   Z[1]=Z[1] - Z[4] + 3.0*Z[6] + 7.0*p_a.b02u - Z[3];
   Z[1]=p_a.xu*Z[1];
   Z[6]=3.0*p_a.b02u;
   Z[7]=Z[6] - p_a.a02;
   Z[1]=Z[1] - 2.0 - Z[7];
   Z[8]=p_a.xu*p_a.xu;
   Z[9]=2.0*Z[8];
   Z[1]=Z[1]*Z[9];
   Z[10]=3.0*p_a.b02t;
   Z[11]=4.0*p_a.c02s;
   Z[12]=2.0*p_a.c02t;
   Z[13]= - Z[11] - Z[12] + Z[10];
   Z[13]=p_a.xu*Z[13];
   Z[14]=p_a.xu + p_a.xt;
   Z[14]=p_a.d02st*Z[14];
   Z[15]=2.0*p_a.c02s;
   Z[16]=p_a.d02st - p_a.c02t;
   Z[16]=3.0*Z[16] - Z[15] + Z[14];
   Z[17]=2.0*p_a.xt;
   Z[16]=Z[16]*Z[17];
   Z[13]=Z[13] + Z[16];
   Z[13]=p_a.xt*Z[13];
   Z[16]=p_a.xu*p_a.c02u;
   Z[18]=p_a.b02u - Z[5] + p_a.b02t;
   Z[18]=2.0*Z[18] + Z[16];
   Z[18]=p_a.xu*Z[18];
   Z[19]=8.0*p_a.b02t;
   Z[18]=Z[19] + Z[18];
   Z[18]=p_a.xu*Z[18];
   Z[13]=Z[18] + Z[13];
   Z[18]=3.0*p_a.xt;
   Z[13]=Z[13]*Z[18];
   Z[1]=Z[1] + Z[13];
   Z[1]=p_a.xt*Z[1];
   Z[13]=2.0*p_a.b020;
   Z[19]= - Z[19] + Z[13] - Z[6];
   Z[20]=4.0*p_a.c02u;
   Z[21]= - Z[4] - Z[20];
   Z[21]=Z[21]*Z[0];
   Z[21]=Z[21] + 30.0*p_a.c02u - 18.0*p_a.d02su + 13.0*p_a.b02u + 12.0*p_a.c02s;
   Z[21]=p_a.xu*Z[21];
   Z[19]=2.0*Z[19] + Z[21];
   Z[21]=p_a.xu*p_a.xu*p_a.xu;
   Z[19]=Z[19]*Z[21];
   Z[1]=Z[19] + Z[1];
   Z[1]=p_a.xt*Z[1];
   Z[19]=Z[2] + p_a.c02s;
   Z[22]=4.0*Z[19] - Z[16];
   Z[23]=p_a.xu*p_a.xu*p_a.xu*p_a.xu*p_a.xu;
   Z[22]=Z[22]*Z[23];
   Z[1]=3.0*Z[22] + Z[1];
   Z[1]=p_a.xt*Z[1];
   Z[22]=p_a.xu*p_a.d02su;
   Z[24]=Z[22] + Z[2] - Z[11] + 7.0*p_a.d02su;
   Z[24]=p_a.xu*Z[24];
   Z[25]=Z[5] - Z[6];
   Z[26]=Z[25] - p_a.b02t;
   Z[24]=Z[24] - Z[26];
   Z[24]=p_a.xu*Z[24];
   Z[24]=10.0*p_a.b02t + Z[24];
   Z[24]=p_a.xu*Z[24];
   Z[27]= - Z[15] - 4.0*p_a.c02t + Z[10];
   Z[27]=p_a.xu*Z[27];
   Z[28]= - p_a.c02s + 3.0*p_a.d02st - 8.0*p_a.c02t;
   Z[28]=2.0*Z[28] + Z[14];
   Z[28]=p_a.xt*Z[28];
   Z[27]=Z[27] + Z[28];
   Z[27]=p_a.xt*Z[27];
   Z[24]=Z[24] + Z[27];
   Z[24]=Z[24]*Z[18];
   Z[27]=3.0*p_a.c02s;
   Z[28]= - Z[27] + 4.0*p_a.d02su;
   Z[28]=Z[22] + 2.0*Z[28] - 9.0*p_a.c02u;
   Z[28]=Z[28]*Z[0];
   Z[29]=6.0*p_a.c02s;
   Z[30]=5.0*p_a.b02u + Z[29];
   Z[28]=Z[28] - Z[2] + 5.0*Z[30] - 12.0*p_a.d02su;
   Z[28]=p_a.xu*Z[28];
   Z[28]=Z[13] + Z[28];
   Z[28]=Z[28]*Z[8];
   Z[24]=Z[28] + Z[24];
   Z[24]=p_a.xt*Z[24];
   Z[28]=5.0*p_a.c02u;
   Z[30]=Z[28] + Z[29];
   Z[31]=2.0*p_a.d02su;
   Z[28]= - Z[28] - p_a.c02s + Z[31];
   Z[28]=p_a.xu*Z[28];
   Z[28]=Z[28] + p_a.b02u + Z[30];
   Z[32]=p_a.xu*p_a.xu*p_a.xu*p_a.xu;
   Z[28]=Z[28]*Z[32];
   Z[24]=6.0*Z[28] + Z[24];
   Z[24]=p_a.xt*Z[24];
   Z[27]= - Z[27] + Z[31];
   Z[27]= - 5.0*Z[22] + 4.0*Z[27] + p_a.c02u;
   Z[27]=Z[27]*Z[21];
   Z[4]= - Z[21]*Z[4];
   Z[28]= - p_a.xt*Z[12];
   Z[4]=Z[4] + Z[28];
   Z[4]=Z[4]*Z[17];
   Z[4]=Z[27] + Z[4];
   Z[4]=p_a.xt*Z[4];
   Z[27]= - p_a.xs*Z[31];
   Z[27]=Z[27] - Z[11] - 3.0*Z[22];
   Z[27]=Z[21]*Z[27];
   Z[28]=6.0*p_a.d02su;
   Z[33]=Z[21]*p_a.xt;
   Z[34]= - Z[33]*Z[28];
   Z[27]=Z[34] + Z[27];
   Z[27]=p_a.xs*Z[27];
   Z[30]= - Z[22] - Z[30];
   Z[30]=Z[30]*Z[32];
   Z[4]=Z[27] + Z[30] + Z[4];
   Z[4]=p_a.xs*Z[4];
   Z[27]= - Z[8]*Z[31];
   Z[27]=Z[27] - Z[12] + p_a.b02t;
   Z[27]=p_a.xu*Z[27];
   Z[30]=p_a.d02st - 7.0*p_a.c02t;
   Z[17]=Z[30]*Z[17];
   Z[17]=Z[27] + Z[17];
   Z[17]=p_a.xt*Z[17];
   Z[27]=p_a.c02u - Z[11] + 5.0*p_a.d02su;
   Z[27]=3.0*Z[27] - Z[22];
   Z[27]=p_a.xu*Z[27];
   Z[30]=2.0*p_a.b02u;
   Z[27]=Z[30] + Z[27];
   Z[27]=p_a.xu*Z[27];
   Z[27]=4.0*p_a.b02t + Z[27];
   Z[27]=p_a.xu*Z[27];
   Z[17]=Z[27] + Z[17];
   Z[17]=p_a.xt*Z[17];
   Z[27]=p_a.c02u - p_a.d02su;
   Z[34]= - p_a.c02s - Z[27];
   Z[34]=p_a.xu*Z[34];
   Z[11]=6.0*Z[34] + p_a.b02u + Z[11];
   Z[34]=2.0*Z[21];
   Z[11]=Z[11]*Z[34];
   Z[11]=Z[11] + Z[17];
   Z[11]=p_a.xt*Z[11];
   Z[17]=Z[19]*Z[23];
   Z[4]=Z[4] - 2.0*Z[17] + Z[11];
   Z[4]=p_a.xs*Z[4];
   Z[11]=p_a.xu*p_a.xu*p_a.xu*p_a.xu*p_a.xu*p_a.xu;
   Z[17]=Z[11]*Z[2];
   Z[4]=3.0*Z[4] - Z[17] + Z[24];
   Z[4]=p_a.xs*Z[4];
   Z[1]=Z[1] + Z[4];
   Z[1]=p_a.xs*Z[1];
   Z[4]=6.0*Z[22];
   Z[19]=2.0*p_a.b02t;
   Z[5]=Z[4] + Z[28] - 18.0*p_a.c02s - Z[19] - 3.0*Z[5] - Z[30];
   Z[5]=p_a.xu*Z[5];
   Z[6]=Z[19] + Z[6];
   Z[22]=3.0*p_a.b020 - p_a.a02 - Z[6];
   Z[5]=2.0*Z[22] + Z[5];
   Z[5]=Z[5]*Z[21];
   Z[22]=3.0 - 2.0*p_a.a02;
   Z[3]=Z[4] - Z[29] + 3.0*Z[22] + p_a.b02u + Z[3];
   Z[3]=p_a.xu*Z[3];
   Z[4]= - 3.0 + p_a.b020 - Z[7];
   Z[3]=2.0*Z[4] + Z[3];
   Z[3]=Z[3]*Z[8];
   Z[4]= - p_a.xu*Z[26];
   Z[4]=Z[19] + Z[4];
   Z[4]=p_a.xu*Z[4];
   Z[19]=p_a.d02st - p_a.c02s;
   Z[14]=2.0*Z[19] + Z[14];
   Z[14]=p_a.xt*Z[14];
   Z[15]=p_a.b02t - Z[15];
   Z[15]=p_a.xu*Z[15];
   Z[14]=Z[15] + Z[14];
   Z[14]=p_a.xt*Z[14];
   Z[4]=Z[4] + Z[14];
   Z[4]=Z[4]*Z[18];
   Z[3]=Z[3] + Z[4];
   Z[3]=p_a.xt*Z[3];
   Z[3]=Z[5] + Z[3];
   Z[3]=p_a.xt*Z[3];
   Z[4]= - p_a.c02s + Z[27];
   Z[4]=Z[4]*Z[0];
   Z[4]=Z[4] + p_a.b020 - Z[6];
   Z[5]=2.0*Z[32];
   Z[4]=Z[4]*Z[5];
   Z[3]=Z[4] + Z[3];
   Z[3]=p_a.xt*Z[3];
   Z[4]=p_a.c02u*Z[11];
   Z[3]=6.0*Z[4] + Z[3];
   Z[3]=p_a.xt*Z[3];
   Z[1]=Z[3] + Z[1];
   Z[1]=p_a.xs*Z[1];
   Z[3]=2.0*p_a.d02tu;
   Z[4]= - Z[20] + Z[31] + Z[3] + p_a.b02t;
   Z[4]=Z[4]*Z[0];
   Z[6]=Z[10] + 4.0*p_a.b020 + Z[7];
   Z[4]=2.0*Z[6] + Z[4];
   Z[4]=Z[4]*Z[21];
   Z[6]=p_a.c02t - p_a.d02tu;
   Z[7]=p_a.b02u - Z[6];
   Z[11]=2.0*p_a.c02u;
   Z[14]= - p_a.d02tu + Z[11];
   Z[14]=p_a.xu*Z[14];
   Z[7]=Z[14] + 2.0*Z[7] - p_a.b02t;
   Z[7]=Z[7]*Z[0];
   Z[14]= - 1.0 + Z[13];
   Z[7]=2.0*Z[14] + Z[7];
   Z[7]=Z[7]*Z[8];
   Z[8]=Z[3] - p_a.c02t;
   Z[14]= - p_a.xu*Z[8];
   Z[14]=p_a.b02u + Z[14];
   Z[9]=Z[14]*Z[9];
   Z[14]=Z[33]*p_a.d02tu;
   Z[9]=Z[9] + Z[14];
   Z[9]=Z[9]*Z[18];
   Z[7]=Z[7] + Z[9];
   Z[7]=p_a.xt*Z[7];
   Z[4]=Z[4] + Z[7];
   Z[4]=p_a.xt*Z[4];
   Z[7]=Z[10] + Z[13] - Z[25];
   Z[2]=Z[31] - Z[2];
   Z[0]=Z[2]*Z[0];
   Z[0]=2.0*Z[7] + Z[0];
   Z[0]=Z[0]*Z[32];
   Z[0]=Z[0] + Z[4];
   Z[0]=p_a.xt*Z[0];
   Z[0]= - Z[17] + Z[0];
   Z[0]=Z[0]*p_a.xt*p_a.xt;
   Z[0]=Z[0] + Z[1];
   Z[0]=p_a.xs*Z[0];
   Z[1]=Z[12] + p_a.b02t;
   Z[2]=Z[16] - Z[11] - Z[1];
   Z[2]=Z[2]*Z[5];
   Z[4]=p_a.xu*p_a.d02tu;
   Z[6]=4.0*Z[6] + Z[4];
   Z[6]=Z[6]*Z[21];
   Z[3]=Z[33]*Z[3];
   Z[3]=Z[6] + Z[3];
   Z[3]=p_a.xt*Z[3];
   Z[6]=Z[11] - Z[8];
   Z[6]=p_a.xu*Z[6];
   Z[1]=Z[6] - Z[1];
   Z[1]=Z[1]*Z[34];
   Z[1]=Z[1] + Z[3];
   Z[1]=p_a.xt*Z[1];
   Z[1]=Z[2] + Z[1];
   Z[1]=p_a.xt*Z[1];
   Z[2]= - Z[23]*Z[20];
   Z[1]=Z[2] + Z[1];
   Z[1]=Z[1]*p_a.xt*p_a.xt*p_a.xt*p_a.xt;
   Z[0]=3.0*Z[1] + Z[0];
   Z[0]=p_a.xs*Z[0];
   Z[1]=p_a.c02t + p_a.c02u;
   Z[1]=Z[1]*Z[5];
   Z[2]=Z[12] + Z[4];
   Z[2]=Z[2]*Z[21];
   Z[2]=Z[2] + Z[14];
   Z[2]=p_a.xt*Z[2];
   Z[1]=Z[1] + Z[2];
   Z[1]=p_a.xt*Z[1];
   Z[2]=Z[23]*Z[11];
   Z[1]=Z[2] + Z[1];
   Z[1]=Z[1]*p_a.xt*p_a.xt*p_a.xt*p_a.xt*p_a.xt*p_a.xt;

cstu = Z[0] + 3.0*Z[1];
 
return cstu; }

complex<double> astu_polynomial_massless(polynomial_parameters p_a) {

  complex<double> shat = p_a.xs;
  complex<double> that = p_a.xt;
  complex<double> uhat = p_a.xu;

  complex<double> shat2 = p_a.xs * p_a.xs;
  complex<double> that2 = p_a.xt * p_a.xt;
  complex<double> uhat2 = p_a.xu * p_a.xu;

  complex<double> shat3 = shat2 * shat;
  complex<double> that3 = that2 * that;
  complex<double> uhat3 = uhat2 * uhat;

  complex<double> A1stu_FormFactor_Massless =
      8. - (8. * (shat2 - 3. * that * uhat) * p_a.b02s) / (that * uhat) -
      (8. * that2 * p_a.b02t) / (shat * uhat) -
      (4. * (3. * shat2 - 2. * that * uhat + shat * (3. * that + uhat)) *
       p_a.b02u) /
          (shat * that) -
      (2. *
       (3. * shat * that2 * uhat * (that + 6. * uhat) +
        3. * that2 * uhat * (that2 + 3. * that * uhat + 2. * uhat2) -
        4. * shat2 * (that3 - 3. * that2 * uhat - 3. * that * uhat2 + uhat3)) *
       p_a.c02s) /
          (that2 * uhat2) +
      (2. * that *
       (shat3 * (4. * that - 9. * uhat) +
        3. * shat2 * (that - 3. * uhat) * uhat + 4. * that * uhat3) *
       p_a.c02t) /
          (shat2 * uhat2) -
      ((-8. * that3 * uhat2 + 3. * shat2 * that2 * (that + 9. * uhat) +
        shat3 * (3. * that2 + 24. * that * uhat - 8. * uhat2)) *
       p_a.c02u) /
          (shat2 * that2) -
      (2. * shat * that * (2. * shat * (that - 3. * uhat) - 3. * uhat2) *
       p_a.d02st) /
          uhat2 -
      (2. * shat *
       (3. * that3 + shat * (3. * that2 - 6. * that * uhat + 2. * uhat2)) *
       p_a.d02su) /
          that2 +
      (2. * that *
       (3. * shat3 - 2. * that * uhat2 + 3. * shat2 * (that + 2. * uhat)) *
       p_a.d02tu) /
          shat2;

  return A1stu_FormFactor_Massless;
}

complex<double> atsu_polynomial_massless(polynomial_parameters p_a) {
  complex<double> shat = p_a.xs;
  complex<double> that = p_a.xt;
  complex<double> uhat = p_a.xu;

  complex<double> shat2 = p_a.xs * p_a.xs;
  complex<double> that2 = p_a.xt * p_a.xt;
  complex<double> uhat2 = p_a.xu * p_a.xu;

  complex<double> shat3 = shat2 * shat;
  complex<double> that3 = that2 * that;
  complex<double> uhat3 = uhat2 * uhat;

  complex<double> A1tsu_FormFactor_Massless =
      8. - (8. * shat2 * p_a.b02s) / (that * uhat) -
      (8. * (that2 - 3. * shat * uhat) * p_a.b02t) / (shat * uhat) -
      (4. * (shat * (3. * that - 2. * uhat) + that * (3. * that + uhat)) *
       p_a.b02u) /
          (shat * that) +
      (8. * shat2 * (that3 + 3. * that2 * uhat + uhat3) * p_a.c02s) /
          (that2 * uhat2) -
      (2. *
       (3. * shat2 * shat2 * uhat + 12. * shat * that2 * uhat2 -
        4. * that2 * uhat3 +
        6. * shat2 * uhat * (2. * that2 + 3. * that * uhat + uhat2) +
        shat3 * (-4. * that2 + 3. * that * uhat + 9. * uhat2)) *
       p_a.c02t) /
          (shat2 * uhat2) -
      ((24. * shat * that3 * uhat - 8. * that3 * uhat2 +
        3. * shat2 * that2 * (that + 9. * uhat) +
        shat3 * (3. * that2 - 8. * uhat2)) *
       p_a.c02u) /
          (shat2 * that2) -
      (2. * shat * that * (2. * shat * that - 3. * uhat * (2. * that + uhat)) *
       p_a.d02st) /
          uhat2 +
      (2. * shat *
       (3. * that2 * (that + 2. * uhat) + shat * (3. * that2 - 2. * uhat2)) *
       p_a.d02su) /
          that2 -
      (2. * that *
       (3. * shat3 + 3. * shat2 * that - 6. * shat * that * uhat +
        2. * that * uhat2) *
       p_a.d02tu) /
          shat2;
  return A1tsu_FormFactor_Massless;
}

complex<double> aust_polynomial_massless(polynomial_parameters p_a) {
  complex<double> shat = p_a.xs;
  complex<double> that = p_a.xt;
  complex<double> uhat = p_a.xu;

  complex<double> shat2 = p_a.xs * p_a.xs;
  complex<double> that2 = p_a.xt * p_a.xt;
  complex<double> uhat2 = p_a.xu * p_a.xu;

  complex<double> shat3 = shat2 * shat;
  complex<double> that3 = that2 * that;
  complex<double> uhat3 = uhat2 * uhat;

  complex<double> A1ust_FormFactor_Massless =
      8. - (8. * shat2 * p_a.b02s) / (that * uhat) -
      (4. * (shat * (-2. * that + 3. * uhat) + uhat * (that + 3. * uhat)) *
       p_a.b02t) /
          (shat * uhat) +
      (8. * (3. * shat * that - uhat2) * p_a.b02u) / (shat * that) +
      (8. * shat2 * (that3 + 3. * that * uhat2 + uhat3) * p_a.c02s) /
          (that2 * uhat2) +
      ((-24. * shat * that * uhat3 + 8. * that2 * uhat3 -
        3. * shat2 * uhat2 * (9. * that + uhat) +
        shat3 * (8. * that2 - 3. * uhat2)) *
       p_a.c02t) /
          (shat2 * uhat2) -
      (2. *
       (3. * shat2 * shat2 * that + 12. * shat * that2 * uhat2 -
        4. * that3 * uhat2 +
        shat3 * (9. * that2 + 3. * that * uhat - 4. * uhat2) +
        6. * shat2 * that * (that2 + 3. * that * uhat + 2. * uhat2)) *
       p_a.c02u) /
          (shat2 * that2) +
      (2. * shat *
       (3. * uhat2 * (2. * that + uhat) + shat * (-2. * that2 + 3. * uhat2)) *
       p_a.d02st) /
          uhat2 -
      (2. * shat * uhat * (-3. * that2 + 2. * shat * uhat - 6. * that * uhat) *
       p_a.d02su) /
          that2 -
      (2. * uhat *
       (3. * shat3 + 3. * shat2 * uhat - 6. * shat * that * uhat +
        2. * that2 * uhat) *
       p_a.d02tu) /
          shat2;
  return A1ust_FormFactor_Massless;
}

complex<double> bstu_polynomial_massless(polynomial_parameters p_a) {
  complex<double> shat = p_a.xs;
  complex<double> that = p_a.xt;
  complex<double> uhat = p_a.xu;

  complex<double> shat2 = p_a.xs * p_a.xs;
  complex<double> that2 = p_a.xt * p_a.xt;
  complex<double> uhat2 = p_a.xu * p_a.xu;

  complex<double> shat3 = shat2 * shat;
  complex<double> that3 = that2 * that;
  complex<double> uhat3 = uhat2 * uhat;

  complex<double> shat4 = shat3 * shat;
  complex<double> that4 = that3 * that;
  complex<double> uhat4 = uhat3 * uhat;

  complex<double> DB111stu_FormFactor_Massless =
      8. / uhat +
      (4. *
       (8. * that3 + 9. * that2 * uhat + 11. * that * uhat2 + 4. * uhat3 +
        shat * (4. * that2 + 3. * that * uhat + 4. * uhat2)) *
       p_a.b02s) /
          (that2 * uhat2) -
      (4. *
       (shat2 * (4. * that - 3. * uhat) + 2. * that * uhat2 +
        shat * (8. * that2 + 3. * that * uhat - 3. * uhat2)) *
       p_a.b02t) /
          (shat * that * uhat2) +
      (4. *
       (shat2 * (3. * that - 4. * uhat) + that2 * (3. * that + 5. * uhat) +
        shat * (6. * that2 - 5. * that * uhat - 4. * uhat2)) *
       p_a.b02u) /
          (shat * that2 * uhat) -
      (2. * shat *
       (-12. * that4 + 6. * shat2 * uhat2 + 9. * that2 * uhat2 -
        17. * that * uhat3 - 26. * uhat4 +
        shat * (-4. * that3 + 6. * that2 * uhat + 9. * that * uhat2 -
                12. * uhat3)) *
       p_a.c02s) /
          (that2 * uhat3) +
      (2. *
       (6. * shat4 * uhat - 3. * shat * that * uhat2 * (that + 3. * uhat) -
        2. * that * uhat3 * (that + 3. * uhat) +
        3. * shat2 * that * (4. * that2 + 2. * that * uhat - uhat2) +
        shat3 * (4. * that2 + 6. * that * uhat + 6. * uhat2)) *
       p_a.c02t) /
          (shat2 * uhat3) -
      ((2. * shat4 * (9. * that - 8. * uhat) +
        4. * that3 * uhat * (that + 3. * uhat) +
        3. * shat * that3 * (5. * that + 9. * uhat) +
        shat3 * (21. * that2 - 24. * that * uhat - 32. * uhat2) +
        shat2 * (18. * that3 - 27. * that2 * uhat - 58. * that * uhat2 -
                 16. * uhat3)) *
       p_a.c02u) /
          (shat2 * that3) -
      (2. * shat * that *
       (2. * shat * that + 6. * that2 - 3. * shat * uhat - 3. * uhat2) *
       p_a.d02st) /
          uhat3 +
      (shat *
       (3. * that3 + shat2 * (6. * that - 8. * uhat) - 9. * that2 * uhat -
        26. * that * uhat2 - 8. * uhat3 +
        shat * (9. * that2 - 12. * that * uhat - 16. * uhat2)) *
       p_a.d02su) /
          that3 +
      (that *
       (3. * shat2 + 3. * shat * (that + 3. * uhat) +
        2. * uhat * (that + 3. * uhat)) *
       p_a.d02tu) /
          shat2;

  return DB111stu_FormFactor_Massless;
}

complex<double> btsu_polynomial_massless(polynomial_parameters p_a) {
  complex<double> shat = p_a.xs;
  complex<double> that = p_a.xt;
  complex<double> uhat = p_a.xu;

  complex<double> shat2 = p_a.xs * p_a.xs;
  complex<double> that2 = p_a.xt * p_a.xt;
  complex<double> uhat2 = p_a.xu * p_a.xu;

  complex<double> shat3 = shat2 * shat;
  complex<double> that3 = that2 * that;
  complex<double> uhat3 = uhat2 * uhat;

  complex<double> shat4 = shat3 * shat;
  complex<double> that4 = that3 * that;
  complex<double> uhat4 = uhat3 * uhat;

  complex<double> DB111tsu_FormFactor_Massless =
      8. / uhat -
      (8. * (4. * shat * that + 2. * that2 + 3. * that * uhat + uhat2) *
       p_a.b02s) /
          (that * uhat2) +
      (4. *
       (8. * shat2 + 4. * shat * that + 9. * shat * uhat + 3. * that * uhat +
        7. * uhat2) *
       p_a.b02t) /
          (shat * uhat2) +
      (4. *
       (3. * shat2 + 6. * shat * that + 3. * that2 + 5. * shat * uhat -
        that * uhat) *
       p_a.b02u) /
          (shat * that * uhat) +
      (2. * shat *
       (-3. * that * uhat2 * (that2 + 3. * that * uhat + 2. * uhat2) +
        4. * shat2 * (3. * that3 + 2. * uhat3) +
        shat * (4. * that4 + 6. * that3 * uhat - 3. * that2 * uhat2 +
                6. * that * uhat3 + 8. * uhat4)) *
       p_a.c02s) /
          (that3 * uhat3) +
      (2. * that *
       (12. * shat4 + 4. * shat3 * that -
        3. * shat2 * uhat * (2. * that + 3. * uhat) +
        shat * uhat2 * (-9. * that + 17. * uhat) +
        2. * uhat2 * (-3. * that2 + 6. * that * uhat + 13. * uhat2)) *
       p_a.c02t) /
          (shat2 * uhat3) +
      ((16. * shat4 * uhat + shat * that3 * (-21. * that + 43. * uhat) -
        3. * shat2 * that * (6. * that2 + 9. * that * uhat + 4. * uhat2) +
        shat3 * (-15. * that2 + 12. * that * uhat + 16. * uhat2) +
        2. * that3 * (-9. * that2 + 12. * that * uhat + 29. * uhat2)) *
       p_a.c02u) /
          (shat2 * that3) -
      (2. * shat2 * that * (6. * shat + 2. * that + 3. * uhat) * p_a.d02st) /
          uhat3 -
      (shat *
       (8. * shat2 * uhat -
        3. * that * (that2 + 3. * that * uhat + 2. * uhat2) +
        shat * (-3. * that2 + 6. * that * uhat + 8. * uhat2)) *
       p_a.d02su) /
          that3 +
      (that *
       (3. * shat2 + 9. * shat * that + 6. * that2 - 17. * shat * uhat -
        12. * that * uhat - 26. * uhat2) *
       p_a.d02tu) /
          shat2;
  return DB111tsu_FormFactor_Massless;
}

complex<double> bust_polynomial_massless(polynomial_parameters p_a) {
  complex<double> shat = p_a.xs;
  complex<double> that = p_a.xt;
  complex<double> uhat = p_a.xu;

  complex<double> shat2 = p_a.xs * p_a.xs;
  complex<double> that2 = p_a.xt * p_a.xt;
  complex<double> uhat2 = p_a.xu * p_a.xu;

  complex<double> shat3 = shat2 * shat;
  complex<double> that3 = that2 * that;
  complex<double> uhat3 = uhat2 * uhat;

  complex<double> shat4 = shat3 * shat;
  complex<double> that4 = that3 * that;
  complex<double> uhat4 = uhat3 * uhat;
  complex<double> DB111ust_FormFactor_Massless =
      8. / that -
      (8. * (that2 + 3. * that * uhat + 2. * uhat * (2. * shat + uhat)) *
       p_a.b02s) /
          (that2 * uhat) +
      (4. *
       (3. * shat2 + 5. * shat * that + 6. * shat * uhat - that * uhat +
        3. * uhat2) *
       p_a.b02t) /
          (shat * that * uhat) +
      (4. *
       (8. * shat2 + 9. * shat * that + 7. * that2 + 4. * shat * uhat +
        3. * that * uhat) *
       p_a.b02u) /
          (shat * that2) +
      (2. * shat *
       (-3. * that2 * uhat * (2. * that2 + 3. * that * uhat + uhat2) +
        4. * shat2 * (2. * that3 + 3. * uhat3) +
        shat * (8. * that4 + 6. * that3 * uhat - 3. * that2 * uhat2 +
                6. * that * uhat3 + 4. * uhat4)) *
       p_a.c02s) /
          (that3 * uhat3) +
      ((16. * shat4 * that + shat * (43. * that - 21. * uhat) * uhat3 +
        shat3 * (16. * that2 + 12. * that * uhat - 15. * uhat2) +
        2. * uhat3 * (29. * that2 + 12. * that * uhat - 9. * uhat2) -
        3. * shat2 * uhat * (4. * that2 + 9. * that * uhat + 6. * uhat2)) *
       p_a.c02t) /
          (shat2 * uhat3) +
      (2. * uhat *
       (12. * shat4 + shat * that2 * (17. * that - 9. * uhat) +
        4. * shat3 * uhat - 3. * shat2 * that * (3. * that + 2. * uhat) +
        2. * that2 * (13. * that2 + 6. * that * uhat - 3. * uhat2)) *
       p_a.c02u) /
          (shat2 * that3) -
      (shat *
       (8. * shat2 * that +
        shat * (8. * that2 + 6. * that * uhat - 3. * uhat2) -
        3. * uhat * (2. * that2 + 3. * that * uhat + uhat2)) *
       p_a.d02st) /
          uhat3 -
      (2. * shat2 * uhat * (6. * shat + 3. * that + 2. * uhat) * p_a.d02su) /
          that3 +
      (uhat *
       (3. * shat2 - 17. * shat * that - 26. * that2 + 9. * shat * uhat -
        12. * that * uhat + 6. * uhat2) *
       p_a.d02tu) /
          shat2;
  return DB111ust_FormFactor_Massless;
}

complex<double> cstu_polynomial_massless(polynomial_parameters p_a) {
  return -(4.0 * p_a.xs * p_a.xs * p_a.xt + 8. * p_a.xs * p_a.xt * p_a.xt +
           4.0 * p_a.xt * p_a.xt * p_a.xt + 13. * p_a.xs * p_a.xt * p_a.xu +
           10. * p_a.xt * p_a.xt * p_a.xu + 3.0 * p_a.xs * p_a.xu * p_a.xu +
           12. * p_a.xt * p_a.xu * p_a.xu + 3. * p_a.xu * p_a.xu * p_a.xu);
}

void APHOAMPFFSTU_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E42;
  complex<double> dynamic_prefactor;
  polynomial_parameters p_a;

  complex<double> tot_res = 0.0;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E42 = E4 * E4;

  // massless loops

  double pref_sum = 0.0;

  for (int n = 0; n < n_massless; n++) {
    pref_sum += pref_arr_massless[n];
  }

  fill_ff_params_massless(s.real(), t.real(), u.real(), &p_a);

  complex<double> expr_massless;

  expr_massless = astu_polynomial_massless(p_a);

  dynamic_prefactor = 1. / (96.0 * E1 * E2 * E3 * E42 * s * s);

  tot_res += pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {

    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = astu_polynomial(p_a);

    dynamic_prefactor = 1. / (96.0 * p_a.xs * p_a.xs * p_a.xt * p_a.xt *
                              p_a.xu * p_a.xu * E1 * E2 * E3 * E42 * s * s);

    tot_res += pref_arr_massive[n] * dynamic_prefactor * expr_massive;
  }

  *out = global_prefactor * tot_res;
}

void APHOAMPFFTSU_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E42;
  complex<double> dynamic_prefactor;
  polynomial_parameters p_a;

  complex<double> tot_res = 0.0;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E42 = E4 * E4;

  // massless loops
  double pref_sum = 0.0;

  for (int n = 0; n < n_massless; n++) {
    pref_sum += pref_arr_massless[n];
  }

  fill_ff_params_massless(s.real(), t.real(), u.real(), &p_a);

  complex<double> expr_massless = atsu_polynomial_massless(p_a);

  dynamic_prefactor = 1. / (96. * E1 * E2 * E3 * E42 * t * t);

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {
    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = atsu_polynomial(p_a);

    dynamic_prefactor = 1. / (96.0 * p_a.xs * p_a.xs * p_a.xt * p_a.xt *
                              p_a.xu * p_a.xu * E1 * E2 * E3 * E42 * t * t);

    tot_res += pref_arr_massive[n] * dynamic_prefactor * expr_massive;
  }

  *out = global_prefactor * tot_res;
}

void APHOAMPFFUST_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E42;
  complex<double> dynamic_prefactor;
  polynomial_parameters p_a;

  complex<double> tot_res = 0.0;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E42 = E4 * E4;

  // massless loops

  double pref_sum = 0.0;

  for (int n = 0; n < n_massless; n++) {
    pref_sum += pref_arr_massless[n];
  }

  fill_ff_params_massless(s.real(), t.real(), u.real(), &p_a);

  complex<double> expr_massless = aust_polynomial_massless(p_a);

  // aust_massless

  dynamic_prefactor = 1. / (96.0 * E1 * E2 * E3 * E42 * u * u);

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {
    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = aust_polynomial(p_a);

    dynamic_prefactor = 1. / (96.0 * p_a.xs * p_a.xs * p_a.xt * p_a.xt *
                              p_a.xu * p_a.xu * E1 * E2 * E3 * E42 * u * u);

    tot_res += pref_arr_massive[n] * dynamic_prefactor * expr_massive;
  }

  *out = global_prefactor * tot_res;
}

void BPHOAMPFFSTU_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E32, E42;
  complex<double> dynamic_prefactor;
  polynomial_parameters p_a;

  complex<double> tot_res = 0.0;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E32 = E3 * E3;
  E42 = E4 * E4;

  // massless loops

  double pref_sum = 0.0;

  for (int n = 0; n < n_massless; n++) {
    pref_sum += pref_arr_massless[n];
  }

  fill_ff_params_massless(s.real(), t.real(), u.real(), &p_a);

  complex<double> expr_massless = bstu_polynomial_massless(p_a);

  dynamic_prefactor = 1. / (48.0 * (E1 * E2 * E32 * E42 * s * s * t * t * u));

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {
    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = bstu_polynomial(p_a);

    dynamic_prefactor =
        1. / (48.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt * p_a.xt *
              p_a.xu * p_a.xu * p_a.xu *
              (E1 * E2 * E32 * E42 * s * s * t * t * u) * M2_arr[n]);

    tot_res += pref_arr_massive[n] * dynamic_prefactor * expr_massive;
  }

  *out = global_prefactor * tot_res;
}

void BPHOAMPFFTSU_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E12, E42;
  complex<double> dynamic_prefactor;
  polynomial_parameters p_a;

  complex<double> tot_res = 0.0;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E12 = E1 * E1;
  E42 = E4 * E4;

  // massless loops

  double pref_sum = 0.0;

  for (int n = 0; n < n_massless; n++) {
    pref_sum += pref_arr_massless[n];
  }

  fill_ff_params_massless(s.real(), t.real(), u.real(), &p_a);

  complex<double> expr_massless = btsu_polynomial_massless(p_a);

  dynamic_prefactor = 1. / (48.0 * (E12 * E2 * E3 * E42 * s * s * t * t * u));

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {

    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = btsu_polynomial(p_a);

    dynamic_prefactor =
        1. / (48.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt * p_a.xt *
              p_a.xu * p_a.xu * p_a.xu *
              (E12 * E2 * E3 * E42 * s * s * t * t * u) * M2_arr[n]);

    tot_res += pref_arr_massive[n] * dynamic_prefactor * expr_massive;
  }

  *out = global_prefactor * tot_res;
}

void BPHOAMPFFUST_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E22, E42;
  complex<double> dynamic_prefactor;
  polynomial_parameters p_a;

  complex<double> tot_res = 0.0;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E22 = E2 * E2;
  E42 = E4 * E4;

  // massless loops

  double pref_sum = 0.0;

  for (int n = 0; n < n_massless; n++) {
    pref_sum += pref_arr_massless[n];
  }

  fill_ff_params_massless(s.real(), t.real(), u.real(), &p_a);

  complex<double> expr_massless = bust_polynomial_massless(p_a);

  dynamic_prefactor = 1. / (48.0 * (E1 * E22 * E3 * E42 * s * s * t * u * u));

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {
    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = bust_polynomial(p_a);

    dynamic_prefactor =
        1. / (48.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt * p_a.xt *
              p_a.xu * p_a.xu * p_a.xu *
              (E1 * E22 * E3 * E42 * s * s * t * u * u * M2_arr[n]));

    tot_res += pref_arr_massive[n] * dynamic_prefactor * expr_massive;
  }

  *out = global_prefactor * tot_res;
}

void CPHOAMPFFSTU_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {

  if (debug) {
    cout << "start evaluation CPHOAMPFFSTU ---------------- " << endl;
  }
  complex<double> s, t, u, E4;
  complex<double> dynamic_prefactor;
  polynomial_parameters p_a;

  complex<double> tot_res = 0.0;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;

  if (debug) {
    cout << "s: " << s << endl;
    cout << "t: " << t << endl;
    cout << "u: " << u << endl;
    cout << "E4: " << E4 << endl;
  }

  // massless loops

  double pref_sum = 0.0;

  for (int n = 0; n < n_massless; n++) {
    pref_sum += pref_arr_massless[n];
  }

  if (debug) {
    cout << "massless prefactor: " << pref_sum << endl;
  }

  fill_ff_params_massless(s.real(), t.real(), u.real(), &p_a);

  if (debug) {
    cout << "massless parameters: " << endl;
    cout << "xs: " << p_a.xs << endl;
    cout << "xt: " << p_a.xt << endl;
    cout << "xu: " << p_a.xu << endl;
    cout << "a02: " << p_a.a02 << endl;
    cout << "b02s: " << p_a.b02s << endl;
    cout << "b02t: " << p_a.b02t << endl;
    cout << "b02u: " << p_a.b02u << endl;
    cout << "b020: " << p_a.b020 << endl;
    cout << "c02s: " << p_a.c02s << endl;
    cout << "c02t: " << p_a.c02t << endl;
    cout << "c02u: " << p_a.c02u << endl;
    cout << "d02su: " << p_a.d02su << endl;
    cout << "d02st: " << p_a.d02st << endl;
    cout << "d02tu: " << p_a.d02tu << endl;
  }

  complex<double> expr_massless = cstu_polynomial_massless(p_a);

  if (debug) {
    cout << "expr_massless: " << expr_massless << endl;
  }

  dynamic_prefactor = 1. / (3.0 * p_a.xs * p_a.xt * p_a.xt * p_a.xu * p_a.xu *
                            (E1 * E2 * E3 * E4 * s * t * t * u));

  if (debug) {
    cout << "dynamic_prefacor: " << dynamic_prefactor << endl;
  }

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  if (debug) {
    cout << "total massless contribution " << tot_res << endl;
  }

  // massive loops

  for (int n = 0; n < n_massive; n++) {
    if (debug) {
      cout << "massless loop #: " << n + 1 << "/" << n_massive << endl;
      cout << "M2: " << M2_arr[n] << endl;
      cout << "pref: " << pref_arr_massive[n] << endl;
    }

    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    if (debug) {
      cout << "massless parameters: " << endl;
      cout << "xs: " << p_a.xs << endl;
      cout << "xt: " << p_a.xt << endl;
      cout << "xu: " << p_a.xu << endl;
      cout << "a02: " << p_a.a02 << endl;
      cout << "b02s: " << p_a.b02s << endl;
      cout << "b02t: " << p_a.b02t << endl;
      cout << "b02u: " << p_a.b02u << endl;
      cout << "b020: " << p_a.b020 << endl;
      cout << "c02s: " << p_a.c02s << endl;
      cout << "c02t: " << p_a.c02t << endl;
      cout << "c02u: " << p_a.c02u << endl;
      cout << "d02su: " << p_a.d02su << endl;
      cout << "d02st: " << p_a.d02st << endl;
      cout << "d02tu: " << p_a.d02tu << endl;
    }

    complex<double> expr_massive = cstu_polynomial(p_a);

    if (debug) {
      cout << "expr_massive: " << expr_massive << endl;
    }

    double M4 = M2_arr[n] * M2_arr[n];

    dynamic_prefactor =
        1. / (3.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt *
              p_a.xt * p_a.xt * p_a.xu * p_a.xu * p_a.xu * p_a.xu *
              (E1 * E2 * E3 * E4 * s * t * t * u) * M4);

    if (debug) {
      cout << "dynamic_prefactor: " << dynamic_prefactor << endl;
    }

    complex<double> massive_contr =
        pref_arr_massive[n] * dynamic_prefactor * expr_massive;

    if (debug) {
      cout << "massive contribution: " << massive_contr << endl;
    }

    tot_res += massive_contr;
  }

  *out = global_prefactor * tot_res;
  if (debug) {
    cout << "out: " << *out << endl;
    cout << "end evaluation CPHOAMPFFSTU --------------------------" << endl;
  }
}