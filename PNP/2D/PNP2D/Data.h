/***********************************
 **@file BH2D.cpp
 **@author Peimeng Yin
 **date 2017 4 3
 **
 **brief   the exact solution u and the function f.
 **
 **********************************/


double _alpha_coef_();
double _beta_coef_();
double _gamma_coef_();
double _epsilon_coef_();
double _g_coef_();
double _B_coef_();
double _L_coef_();
double _L1_coef_();
double _mu_(const double u);
double _phi_(const double u);
double _H_(const double u);

double _U_(const double t, const double * p);
double _U_(const double u);
double _U_0_(const double * p);
double _phi_(const double t);

double _u_(const double t, const double * p);

std::vector<double> _u_grad_(const double t, const double * p);
std::vector<double> _c1_grad_(const double t, const double * p);
std::vector<double> _c2_grad_(const double t, const double * p);

double _u_laplace_(const double t, const double * p);

double _q_(const double t, const double * p);

std::vector<double> _q_grad_(const double t, const double * p);

double _u_0_(const double * p);

double _q_0_(const double * p);

std::vector<double> _u_grad_0_(const double * p);

std::vector<double> _q_grad_0_(const double * p);

double _rho_(const double t, const double * p);
double _f1_(const double t, const double * p);
double _f2_(const double t, const double * p);

double _c1_(const double t, const double * p);
double _c1_0_(const double * p);

double _c2_(const double t, const double * p);
double _c2_0_(const double * p);

double _g_(const double * p);