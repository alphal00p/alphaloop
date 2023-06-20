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

bool debug = false;

double global_prefactor = 2.0 / (M_PI * M_PI);

complex<double> zero_f64(0.0, 0.0);
double mu2 = 100.0;
double scale = 1000.0;

int n_massive = 1;
int n_massless = 0;

double pref_arr_massless[0];
double pref_arr_massive[1] = {3.0 * 16.0 / 81.0};

double M2_arr[1] = {173.0 * 173.0};

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

complex<double> astu_polynomial(polynomial_parameters p_a) { int astu; }

complex<double> atsu_polynomial(polynomial_parameters p_a) { int atsu; }

complex<double> aust_polynomial(polynomial_parameters p_a) { int aust; }

complex<double> bstu_polynomial(polynomial_parameters p_a) { int bstu; }

complex<double> btsu_polynomial(polynomial_parameters p_a) { int btsu; }

complex<double> bust_polynomial(polynomial_parameters p_a) { int bust; }

complex<double> cstu_polynomial(polynomial_parameters p_a) { int cstu; }

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

  expr_massless = astu_polynomial(p_a);

  dynamic_prefactor = 1. / (96.0 * p_a.xs * p_a.xs * p_a.xt * p_a.xt * p_a.xu *
                            p_a.xu * E1 * E2 * E3 * E42 * s * s);

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

  complex<double> expr_massless = atsu_polynomial(p_a);

  dynamic_prefactor = 1. / (96.0 * p_a.xs * p_a.xs * p_a.xt * p_a.xt * p_a.xu *
                            p_a.xu * E1 * E2 * E3 * E42 * t * t);

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {
    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = atsu_polynomial(p_a);

    dynamic_prefactor = 1. / (96.0 * p_a.xs * p_a.xs * p_a.xt * p_a.xt *
                              p_a.xu * p_a.xu * E1 * E2 * E3 * E42 * t * t);

    tot_res = pref_arr_massive[n] * dynamic_prefactor * expr_massive;
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

  complex<double> expr_massless = aust_polynomial(p_a);

  // aust_massless

  dynamic_prefactor = 1. / (96.0 * p_a.xs * p_a.xs * p_a.xt * p_a.xt * p_a.xu *
                            p_a.xu * E1 * E2 * E3 * E42 * u * u);

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {
    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = aust_polynomial(p_a);

    dynamic_prefactor = 1. / (96.0 * p_a.xs * p_a.xs * p_a.xt * p_a.xt *
                              p_a.xu * p_a.xu * E1 * E2 * E3 * E42 * u * u);

    tot_res = pref_arr_massive[n] * dynamic_prefactor * expr_massive;
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

  complex<double> expr_massless = bstu_polynomial(p_a);

  dynamic_prefactor = 1. / (48.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt *
                            p_a.xt * p_a.xu * p_a.xu * p_a.xu *
                            (E1 * E2 * E32 * E42 * s * s * t * t * u));

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {
    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = bstu_polynomial(p_a);

    dynamic_prefactor =
        1. / (48.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt * p_a.xt *
              p_a.xu * p_a.xu * p_a.xu *
              (E1 * E2 * E32 * E42 * s * s * t * t * u) * M2_arr[n]);

    tot_res = pref_arr_massive[n] * dynamic_prefactor * expr_massive;
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

  complex<double> expr_massless = btsu_polynomial(p_a);

  dynamic_prefactor = 1. / (48.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt *
                            p_a.xt * p_a.xu * p_a.xu * p_a.xu *
                            (E12 * E2 * E3 * E42 * s * s * t * t * u));

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {

    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = btsu_polynomial(p_a);

    dynamic_prefactor =
        1. / (48.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt * p_a.xt *
              p_a.xu * p_a.xu * p_a.xu *
              (E12 * E2 * E3 * E42 * s * s * t * t * u) * M2_arr[n]);

    tot_res = pref_arr_massive[n] * dynamic_prefactor * expr_massive;
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

  complex<double> expr_massless = bust_polynomial(p_a);

  dynamic_prefactor = 1. / (48.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt *
                            p_a.xt * p_a.xu * p_a.xu * p_a.xu *
                            (E1 * E22 * E3 * E42 * s * s * t * u * u));

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {
    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = bust_polynomial(p_a);

    dynamic_prefactor =
        1. / (48.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt * p_a.xt *
              p_a.xu * p_a.xu * p_a.xu *
              (E1 * E22 * E3 * E42 * s * s * t * u * u * M2_arr[n]));

    tot_res = pref_arr_massive[n] * dynamic_prefactor * expr_massive;
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

  complex<double> expr_massless = cstu_polynomial(p_a);

  dynamic_prefactor =
      1. / (3.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt * p_a.xt *
            p_a.xt * p_a.xu * p_a.xu * p_a.xu * p_a.xu *
            (E1 * E2 * E3 * E4 * s * t * t * u));

  tot_res = pref_sum * dynamic_prefactor * expr_massless;

  // massive loops

  for (int n = 0; n < n_massive; n++) {
    fill_ff_params_massive(s.real(), t.real(), u.real(), M2_arr[n], &p_a);

    complex<double> expr_massive = cstu_polynomial(p_a);

    double M4 = M2_arr[n] * M2_arr[n];

    dynamic_prefactor =
        1. / (3.0 * p_a.xs * p_a.xs * p_a.xs * p_a.xs * p_a.xt * p_a.xt *
              p_a.xt * p_a.xt * p_a.xu * p_a.xu * p_a.xu * p_a.xu *
              (E1 * E2 * E3 * E4 * s * t * t * u) * M4);

    tot_res = pref_arr_massive[n] * dynamic_prefactor * expr_massive;
  }

  *out = global_prefactor * tot_res;
  if (debug) {
    cout << "out: " << out << endl;
    cout << "end evaluation CPHOAMPFFSTU --------------------------" << endl;
  }
}