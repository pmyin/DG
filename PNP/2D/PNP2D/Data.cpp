/***********************************
 **@file BH2D.cpp
 **@author Peimeng Yin
 **date 2017 4 3
 **
 **brief   This is an example in 2d.
 **
 **********************************/

#include "BH2d.h"
#include "Data.h"


#define PI (4.0*atan(1.0))

// Example 2D f=0 1D data Example modified Zhongming 2D positivity  nonzero Neumann, has not finished
/*

#define a_coef 1.0
#define t_coef (1.0)
#define c1_coef 1.0    // (4.8e-4)   // C1,C2 coefficients
#define c2_coef 1.0    // (1.0e-4)   // C1,C2 coefficients
#define u_coef (1.0/420.0)   // phi coefficient 


double _L_coef_()       // used in readData for Pnt[j][0]=Pnt[j][0]*_L_coef_()+_L1_coef_();
{
	double val= 1.0;
	return val;
}

double _L1_coef_()
{
	double val= 0.0;
	return val;
}


double _rho_(const double t, const double * p)
{
	//double val = -u_coef*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	//      val += -u_coef*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	double      val = 0.0+0.0*(_c2_(t,p)-_c1_(t,p));
	return val;
}

double _u_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = -u_coef*exp(-t_coef*t)*( 10.0*pow(p[0],7)-28.0*pow(p[0],6)+21.0*pow(p[0],5) );
	      val += -u_coef*exp(-t_coef*t)*( 10.0*pow(p[1],7)-28.0*pow(p[1],6)+21.0*pow(p[1],5) );
		  val = 0.0*val;
	return val;
}

double _v_laplace_(const double t, const double * p)
{
	double val = 0.0*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	      val += 0.0*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	return val;
}

std::vector<double> _u_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*u_coef*exp(-t_coef*t)*( 0.0 );
	if (p[1] == 1.0)
		val[1]=1.0;
	if (p[1] == 1.0)
		
    return val;
}

double _c1_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c1_coef*( PI*sin(PI*p[0])+PI*sin(PI*p[1]) );
	val = 0.1/2*val;
	return val;
}

double _c2_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c2_coef*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) + p[1]*p[1]*(1.0-p[1])*(1.0-p[1]) );
	val = 15*0.2*val;
	return val;
}

std::vector<double> _c1_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

std::vector<double> _c2_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

double _u_0_(const double * p)
{
	double val = _u_(0.0,p);
	return val;
}

double _c1_0_(const double * p)
{
	double val = _c1_(0.0,p);
	return val;
}

double _c2_0_(const double * p)
{
	double val = _c2_(0.0,p);
	return val;
}

double _f1_(const double t, const double * p)
{
	double val0 = ( 50*pow(p[0],9)-198*pow(p[0],8)+292*pow(p[0],7)-189*pow(p[0],6)+45*pow(p[0],5)  )/(30*exp(2*t));
	double val1 = (-pow(p[0],4)+2*pow(p[0],3)-13*pow(p[0],2)+12*p[0]-2 )/(exp(t));
	double val2 = ( 50*pow(p[1],9)-198*pow(p[1],8)+292*pow(p[1],7)-189*pow(p[1],6)+45*pow(p[1],5)  )/(30*exp(2*t));
	double val3 = (-pow(p[1],4)+2*pow(p[1],3)-13*pow(p[1],2)+12*p[1]-2 )/(exp(t));
	double val = val0+val1+val2+val3; 
	val = 0.0*val;
	return val;
}

double _f2_(const double t, const double * p)
{
	double val0 = ( (p[0]-1)*( 110*pow(p[0],9)-430*pow(p[0],8)+623*pow(p[0],7)-393*pow(p[0],6)+90*pow(p[0],5) ) )/(60*exp(2*t));
	double val1 = (p[0]-1)*(pow(p[0],4)-2*pow(p[0],3)+21*pow(p[0],2)-16*p[0]+2 )/(exp(t));
	double val2 = ( (p[1]-1)*( 110*pow(p[1],9)-430*pow(p[1],8)+623*pow(p[1],7)-393*pow(p[1],6)+90*pow(p[1],5) ) )/(60*exp(2*t));
	double val3 = (p[1]-1)*(pow(p[1],4)-2*pow(p[1],3)+21*pow(p[1],2)-16*p[1]+2 )/(exp(t));
	double val = val0+val1+val2+val3; 
	val = 0.0*val;
	return val;
}

double _g_(const double * p)
{
	double val = exp(p[0]*p[1]);
	val = 1.0*val;
	return val;
}

// */

// Example 2D f=0 1D data Example modified Zhongming 2D positivity, zero Neumann BC
/*

#define a_coef 1.0
#define t_coef (1.0)
#define c1_coef 1.0    // (4.8e-4)   // C1,C2 coefficients
#define c2_coef 1.0    // (1.0e-4)   // C1,C2 coefficients
#define u_coef (1.0/420.0)   // phi coefficient 


double _L_coef_()       // used in readData for Pnt[j][0]=Pnt[j][0]*_L_coef_()+_L1_coef_();
{
	double val= 1.0;
	return val;
}

double _L1_coef_()
{
	double val= 0.0;
	return val;
}


double _rho_(const double t, const double * p)
{
	//double val = -u_coef*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	//      val += -u_coef*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	double      val = 0.0+0.0*(_c2_(t,p)-_c1_(t,p));
	return val;
}

double _u_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = -u_coef*exp(-t_coef*t)*( 10.0*pow(p[0],7)-28.0*pow(p[0],6)+21.0*pow(p[0],5) );
	      val += -u_coef*exp(-t_coef*t)*( 10.0*pow(p[1],7)-28.0*pow(p[1],6)+21.0*pow(p[1],5) );
		  val = 0.0*val;
	return val;
}

double _v_laplace_(const double t, const double * p)
{
	double val = 0.0*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	      val += 0.0*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	return val;
}

std::vector<double> _u_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*u_coef*exp(-t_coef*t)*( 0.0 );
    val[1]=0.0*u_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

double _c1_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c1_coef*( PI*sin(PI*p[0])+PI*sin(PI*p[1]) );
	val = 0.1/2*val;
	return val;
}

double _c2_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c2_coef*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) + p[1]*p[1]*(1.0-p[1])*(1.0-p[1]) );
	val = 15*0.2*val;
	return val;
}

std::vector<double> _c1_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

std::vector<double> _c2_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

double _u_0_(const double * p)
{
	double val = _u_(0.0,p);
	return val;
}

double _c1_0_(const double * p)
{
	double val = _c1_(0.0,p);
	return val;
}

double _c2_0_(const double * p)
{
	double val = _c2_(0.0,p);
	return val;
}

double _f1_(const double t, const double * p)
{
	double val0 = ( 50*pow(p[0],9)-198*pow(p[0],8)+292*pow(p[0],7)-189*pow(p[0],6)+45*pow(p[0],5)  )/(30*exp(2*t));
	double val1 = (-pow(p[0],4)+2*pow(p[0],3)-13*pow(p[0],2)+12*p[0]-2 )/(exp(t));
	double val2 = ( 50*pow(p[1],9)-198*pow(p[1],8)+292*pow(p[1],7)-189*pow(p[1],6)+45*pow(p[1],5)  )/(30*exp(2*t));
	double val3 = (-pow(p[1],4)+2*pow(p[1],3)-13*pow(p[1],2)+12*p[1]-2 )/(exp(t));
	double val = val0+val1+val2+val3; 
	val = 0.0*val;
	return val;
}

double _f2_(const double t, const double * p)
{
	double val0 = ( (p[0]-1)*( 110*pow(p[0],9)-430*pow(p[0],8)+623*pow(p[0],7)-393*pow(p[0],6)+90*pow(p[0],5) ) )/(60*exp(2*t));
	double val1 = (p[0]-1)*(pow(p[0],4)-2*pow(p[0],3)+21*pow(p[0],2)-16*p[0]+2 )/(exp(t));
	double val2 = ( (p[1]-1)*( 110*pow(p[1],9)-430*pow(p[1],8)+623*pow(p[1],7)-393*pow(p[1],6)+90*pow(p[1],5) ) )/(60*exp(2*t));
	double val3 = (p[1]-1)*(pow(p[1],4)-2*pow(p[1],3)+21*pow(p[1],2)-16*p[1]+2 )/(exp(t));
	double val = val0+val1+val2+val3; 
	val = 0.0*val;
	return val;
}

double _g_(const double * p)
{
	double val = exp(p[0]*p[1]);
	val = 1.0*val;
	return val;
}

// */

// Example 2D f=0 1D data Example modified Zhongming 2D mass, energy 
/*

#define a_coef 1.0
#define t_coef (1.0)
#define c1_coef 1.0    // (4.8e-4)   // C1,C2 coefficients
#define c2_coef 1.0    // (1.0e-4)   // C1,C2 coefficients
#define u_coef (1.0/420.0)   // phi coefficient 


double _L_coef_()       // used in readData for Pnt[j][0]=Pnt[j][0]*_L_coef_()+_L1_coef_();
{
	double val= 1.0;
	return val;
}

double _L1_coef_()
{
	double val= 0.0;
	return val;
}


double _rho_(const double t, const double * p)
{
	//double val = -u_coef*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	//      val += -u_coef*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	double      val = 0.0*(_c2_(t,p)-_c1_(t,p));
	return val;
}

double _u_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = -u_coef*exp(-t_coef*t)*( 10.0*pow(p[0],7)-28.0*pow(p[0],6)+21.0*pow(p[0],5) );
	      val += -u_coef*exp(-t_coef*t)*( 10.0*pow(p[1],7)-28.0*pow(p[1],6)+21.0*pow(p[1],5) );
		  val = 0.0*val;
	return val;
}

double _v_laplace_(const double t, const double * p)
{
	double val = 0.0*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	      val += 0.0*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	return val;
}

std::vector<double> _u_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*u_coef*exp(-t_coef*t)*( 0.0 );
    val[1]=0.0*u_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

double _c1_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c1_coef*( 1.0+0.5*PI*sin(PI*p[0])+0.5*PI*sin(PI*p[1]) );
	val = 0.1*val;
	return val;
}

double _c2_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c2_coef*( 4.0-p[0]-p[1] );
	val = 0.1*val;
	return val;
}

std::vector<double> _c1_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

std::vector<double> _c2_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

double _u_0_(const double * p)
{
	double val = _u_(0.0,p);
	return val;
}

double _c1_0_(const double * p)
{
	double val = _c1_(0.0,p);
	return val;
}

double _c2_0_(const double * p)
{
	double val = _c2_(0.0,p);
	return val;
}

double _f1_(const double t, const double * p)
{
	double val0 = ( 50*pow(p[0],9)-198*pow(p[0],8)+292*pow(p[0],7)-189*pow(p[0],6)+45*pow(p[0],5)  )/(30*exp(2*t));
	double val1 = (-pow(p[0],4)+2*pow(p[0],3)-13*pow(p[0],2)+12*p[0]-2 )/(exp(t));
	double val2 = ( 50*pow(p[1],9)-198*pow(p[1],8)+292*pow(p[1],7)-189*pow(p[1],6)+45*pow(p[1],5)  )/(30*exp(2*t));
	double val3 = (-pow(p[1],4)+2*pow(p[1],3)-13*pow(p[1],2)+12*p[1]-2 )/(exp(t));
	double val = val0+val1+val2+val3; 
	val = 0.0*val;
	return val;
}

double _f2_(const double t, const double * p)
{
	double val0 = ( (p[0]-1)*( 110*pow(p[0],9)-430*pow(p[0],8)+623*pow(p[0],7)-393*pow(p[0],6)+90*pow(p[0],5) ) )/(60*exp(2*t));
	double val1 = (p[0]-1)*(pow(p[0],4)-2*pow(p[0],3)+21*pow(p[0],2)-16*p[0]+2 )/(exp(t));
	double val2 = ( (p[1]-1)*( 110*pow(p[1],9)-430*pow(p[1],8)+623*pow(p[1],7)-393*pow(p[1],6)+90*pow(p[1],5) ) )/(60*exp(2*t));
	double val3 = (p[1]-1)*(pow(p[1],4)-2*pow(p[1],3)+21*pow(p[1],2)-16*p[1]+2 )/(exp(t));
	double val = val0+val1+val2+val3; 
	val = 0.0*val;
	return val;
}

double _g_(const double * p)
{
	double val = exp(p[0]*p[1]);
	val = 1.0*val;
	return val;
}

// */

// Example 2D f=0 1D data Example modified Zhongming 2D, added rho have not finished 
/*

#define a_coef 1.0
#define t_coef (1.0)
#define c1_coef 1.0    // (4.8e-4)   // C1,C2 coefficients
#define c2_coef 1.0    // (1.0e-4)   // C1,C2 coefficients
#define u_coef (1.0/420.0)   // phi coefficient 


double _L_coef_()       // used in readData for Pnt[j][0]=Pnt[j][0]*_L_coef_()+_L1_coef_();
{
	double val= 1.0;
	return val;
}

double _L1_coef_()
{
	double val= 0.0;
	return val;
}


double _rho_(const double t, const double * p)
{
	double val = -u_coef*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	      val += -u_coef*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
//	double      val = 1.0*(_c2_(t,p)-_c1_(t,p));
	return val;
}

double _u_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = -u_coef*exp(-t_coef*t)*( 10.0*pow(p[0],7)-28.0*pow(p[0],6)+21.0*pow(p[0],5) );
	      val += -u_coef*exp(-t_coef*t)*( 10.0*pow(p[1],7)-28.0*pow(p[1],6)+21.0*pow(p[1],5) );
		  val = 0.5*val;
	return val;
}

double _v_laplace_(const double t, const double * p)
{
	double val = 0.0*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	      val += 0.0*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	return val;
}

std::vector<double> _u_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*u_coef*exp(-t_coef*t)*( 0.0 );
    val[1]=0.0*u_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

double _c1_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c1_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) );
	      val += c1_coef*exp(-t_coef*t)*( p[1]*p[1]*(1.0-p[1])*(1.0-p[1]) );
	val = 0.5*val;
	return val;
}

double _c2_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c2_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0])*(1.0-p[0]) );
	      val += c2_coef*exp(-t_coef*t)*( p[1]*p[1]*(1.0-p[1])*(1.0-p[1])*(1.0-p[1]) );
	val = 0.5*val;
	return val;
}

std::vector<double> _c1_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

std::vector<double> _c2_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

double _u_0_(const double * p)
{
	double val = _u_(0.0,p);
	return val;
}

double _c1_0_(const double * p)
{
	double val = _c1_(0.0,p);
	return val;
}

double _c2_0_(const double * p)
{
	double val = _c2_(0.0,p);
	return val;
}

double _f1_(const double t, const double * p)
{
	double val0 = ( 50*pow(p[0],9)-198*pow(p[0],8)+292*pow(p[0],7)-189*pow(p[0],6)+45*pow(p[0],5)  )/(30*exp(2*t));
	double val1 = (-pow(p[0],4)+2*pow(p[0],3)-13*pow(p[0],2)+12*p[0]-2 )/(exp(t));
	double val2 = ( 50*pow(p[1],9)-198*pow(p[1],8)+292*pow(p[1],7)-189*pow(p[1],6)+45*pow(p[1],5)  )/(30*exp(2*t));
	double val3 = (-pow(p[1],4)+2*pow(p[1],3)-13*pow(p[1],2)+12*p[1]-2 )/(exp(t));
	double val = val0+val1+val2+val3; 
	val = 0.5*val;
	return val;
}

double _f2_(const double t, const double * p)
{
	double val0 = ( (p[0]-1)*( 110*pow(p[0],9)-430*pow(p[0],8)+623*pow(p[0],7)-393*pow(p[0],6)+90*pow(p[0],5) ) )/(60*exp(2*t));
	double val1 = (p[0]-1)*(pow(p[0],4)-2*pow(p[0],3)+21*pow(p[0],2)-16*p[0]+2 )/(exp(t));
	double val2 = ( (p[1]-1)*( 110*pow(p[1],9)-430*pow(p[1],8)+623*pow(p[1],7)-393*pow(p[1],6)+90*pow(p[1],5) ) )/(60*exp(2*t));
	double val3 = (p[1]-1)*(pow(p[1],4)-2*pow(p[1],3)+21*pow(p[1],2)-16*p[1]+2 )/(exp(t));
	double val = val0+val1+val2+val3; 
	val = 0.5*val;
	return val;
}

double _g_(const double * p)
{
	double val = exp(p[0]*p[1]);
	val = 1.0*val;
	return val;
}

// */

// Example 2D f=0 1D data Example modified Zhongming 2D /*

#define a_coef 1.0
#define t_coef (1.0e-1)
#define c1_coef (1.0e-1)    // (4.8e-4)   // C1,C2 coefficients
#define c2_coef (1.0e-1)    // (1.0e-4)   // C1,C2 coefficients
#define c1_up (1e-2*c1_coef)   // move c up
#define c2_up (1e-2*c2_coef)   // move c up
#define u_coef (1.0e-1)   // phi coefficient 


double _L_coef_()       // used in readData for Pnt[j][0]=Pnt[j][0]*_L_coef_()+_L1_coef_();
{
	double val= 1.0;
	return val;
}

double _L1_coef_()
{
	double val= 0.0;
	return val;
}


double _rho_(const double t, const double * p)
{
	double val = -u_coef*exp(-t_coef*t)*(  20.0*p[0]*p[0]*p[0]-30*p[0]*p[0]+12*p[0]-1.0  );
	      val += -u_coef*exp(-t_coef*t)*(  20.0*p[1]*p[1]*p[1]-30*p[1]*p[1]+12*p[1]-1.0  );
	      val += (_c2_(t,p)-_c1_(t,p));
	return val;
}

double _u_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = u_coef*exp(-t_coef*t)*( (2*p[0]-1.0)*p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) );
	      val += u_coef*exp(-t_coef*t)*( (2*p[1]-1.0)*p[1]*p[1]*(1.0-p[1])*(1.0-p[1]) );
		  val = 0.5*val;
	return val;
}

double _v_laplace_(const double t, const double * p)
{
	double val = 0.0*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	      val += 0.0*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	return val;
}

std::vector<double> _u_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*u_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0])+2*p[0]*(p[0]-0.5)*(1.0-p[0])*(1.0-p[0])-2*(p[0]-0.5)*p[0]*p[0]*(1.0-p[0]) );
    val[1]=0.0*u_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0])+2*p[0]*(p[0]-0.5)*(1.0-p[0])*(1.0-p[0])-2*(p[0]-0.5)*p[0]*p[0]*(1.0-p[0]) );
    return val;
}

double _c1_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c1_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) );
	      val += c1_coef*exp(-t_coef*t)*( p[1]*p[1]*(1.0-p[1])*(1.0-p[1]) );
	val = 1.0*val+c1_up;
	return val;
}

double _c2_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c2_coef*exp(-t_coef*t)*( p[0]*p[0]*p[0]*(1.0-p[0])*(1.0-p[0])*(1.0-p[0]) );
	      val += c2_coef*exp(-t_coef*t)*( p[1]*p[1]*p[1]*(1.0-p[1])*(1.0-p[1])*(1.0-p[1]) );
	val = 1.0*val+c2_up;
	return val;
}

std::vector<double> _c1_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

std::vector<double> _c2_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

double _u_0_(const double * p)
{
	double val = _u_(0.0,p);
	return val;
}

double _c1_0_(const double * p)
{
	double val = _c1_(0.0,p);
	return val;
}

double _c2_0_(const double * p)
{
	double val = _c2_(0.0,p);
	return val;
}

double _f1_(const double t, const double * p)
{
	double val0 = -t_coef*(_c1_(t,p)-c1_up);
	double val1 = -(12*p[0]*p[0]-12*p[0]+2)-(12*p[1]*p[1]-12*p[1]+2);
	double val2 = _c1_(t,p)*_rho_(t,p);
	double val3 = - 2*p[0]*p[0]*(1.0-p[0])*(1.0-p[0])*(2*p[0]-1)*(5*p[0]*p[0]-5*p[0]+1);
	double val4 = - 2*p[1]*p[1]*(1.0-p[1])*(1.0-p[1])*(2*p[1]-1)*(5*p[1]*p[1]-5*p[1]+1);
	double val = val0+val1+val2+val3+val4; 
	return val;
}

double _f2_(const double t, const double * p)
{
	double val0 = -t_coef*(_c2_(t,p)-c2_up);
	double val1 = 6*(p[0]-1)*p[0]*(5*p[0]*p[0]-5*p[0]+1) + 6*(p[1]-1)*p[1]*(5*p[1]*p[1]-5*p[1]+1);
	double val2 = -_c2_(t,p)*_rho_(t,p);
	double val3 = - 3.0*pow(p[0],3)*pow((p[0]-1),3)*(2*p[0]-1)*(5*p[0]*p[0]-5*p[0]+1);
	double val4 = - 3.0*pow(p[1],3)*pow((p[1]-1),3)*(2*p[1]-1)*(5*p[1]*p[1]-5*p[1]+1);
	double val = val0+val1+val2+val3+val4; 
	return val;
}

double _g_(const double * p)
{
	double val = exp(p[0]*p[1]);
	val = 1.0*val;
	return val;
}

// */

// Example 2D f=0 1D data
/*

#define a_coef 1.0
#define t_coef (1.0)
#define c1_coef 1.0    // (4.8e-4)   // C1,C2 coefficients
#define c2_coef 1.0    // (1.0e-4)   // C1,C2 coefficients
#define u_coef (1.0/420.0)   // phi coefficient 


double _L_coef_()       // used in readData for Pnt[j][0]=Pnt[j][0]*_L_coef_()+_L1_coef_();
{
	double val= 1.0;
	return val;
}

double _L1_coef_()
{
	double val= 0.0;
	return val;
}


double _rho_(const double t, const double * p)
{
	//double val = -u_coef*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	//      val += -u_coef*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	double      val = 0.0*(_c2_(t,p)-_c1_(t,p));
	return val;
}

double _u_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = -u_coef*exp(-t_coef*t)*( 10.0*pow(p[0],7)-28.0*pow(p[0],6)+21.0*pow(p[0],5) );
	return val;
}

double _v_laplace_(const double t, const double * p)
{
	double val = 0.0*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	      val += 0.0*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	return val;
}

std::vector<double> _u_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*u_coef*exp(-t_coef*t)*( 0.0 );
    val[1]=0.0*u_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

double _c1_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c1_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) );
	val = 1.0*val;
	return val;
}

double _c2_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c2_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0])*(1.0-p[0]) );
	val = 1.0*val;
	return val;
}

std::vector<double> _c1_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c1_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

std::vector<double> _c2_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=0.0*c2_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}



double _u_0_(const double * p)
{
	double val = _u_(0.0,p);
	return val;
}

double _c1_0_(const double * p)
{
	double val = _c1_(0.0,p);
	return val;
}

double _c2_0_(const double * p)
{
	double val = _c2_(0.0,p);
	return val;
}

double _f1_(const double t, const double * p)
{
	double val0 = ( 50*pow(p[0],9)-198*pow(p[0],8)+292*pow(p[0],7)-189*pow(p[0],6)+45*pow(p[0],5)  )/(30*exp(2*t));
	double val1 = (-pow(p[0],4)+2*pow(p[0],3)-13*pow(p[0],2)+12*p[0]-2 )/(exp(t));
	double val = val0+val1; 
	return val;
}

double _f2_(const double t, const double * p)
{
	double val0 = ( (p[0]-1)*( 110*pow(p[0],9)-430*pow(p[0],8)+623*pow(p[0],7)-393*pow(p[0],6)+90*pow(p[0],5) ) )/(60*exp(2*t));
	double val1 = (p[0]-1)*(pow(p[0],4)-2*pow(p[0],3)+21*pow(p[0],2)-16*p[0]+2 )/(exp(t));
	double val = val0+val1;
	return val;
}

double _g_(const double * p)
{
	double val = exp(p[0]*p[1]);
	val = 1.0*val;
	return val;
}

// */

// Example 2D f=0    
/*

#define a_coef 1.0
#define t_coef (1.0)
#define c1_coef 1.0    // (4.8e-4)   // C1,C2 coefficients
#define c2_coef 1.0    // (1.0e-4)   // C1,C2 coefficients
#define u_coef 1.0e-2   // phi coefficient 


double _L_coef_()       // used in readData for Pnt[j][0]=Pnt[j][0]*_L_coef_()+_L1_coef_();
{
	double val= 1.0;
	return val;
}

double _L1_coef_()
{
	double val= 0.0;
	return val;
}


double _rho_(const double t, const double * p)
{
	double val = -u_coef*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	      val += -u_coef*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	      val += _c2_(t,p)-_c1_(t,p);
	return val;
}

double _u_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = u_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) + p[1]*p[1]*(1.0-p[1])*(1.0-p[1]) );
	return val;
}

double _v_laplace_(const double t, const double * p)
{
	double val = exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	      val += exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	return val;
}

std::vector<double> _u_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=u_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=u_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

double _c1_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c1_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) + p[1]*p[1]*(1.0-p[1])*(1.0-p[1]) );
	val = 1.0*val;
	return val;
}

double _c2_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c2_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) + p[1]*p[1]*(1.0-p[1])*(1.0-p[1]) );
	val = 1.0*val;
	return val;
}

std::vector<double> _c1_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=c1_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=c1_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

std::vector<double> _c2_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=c2_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=c2_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}



double _u_0_(const double * p)
{
	double val = _u_(0.0,p);
	return val;
}

double _c1_0_(const double * p)
{
	double val = _c1_(0.0,p);
	return val;
}

double _c2_0_(const double * p)
{
	double val = _c2_(0.0,p);
	return val;
}

double _f1_(const double t, const double * p)
{
	double val0 = -t_coef*c1_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) + p[1]*p[1]*(1.0-p[1])*(1.0-p[1]) );
	double val1 = c1_coef*_v_laplace_(t,p);
	double val2 = _c1_(t,p)*_rho_(t,p);
	double val3 = _c1_grad_(t,p)[0]*_u_grad_(t,p)[0]+_c1_grad_(t,p)[1]*_u_grad_(t,p)[1];
	double val = val0-val1-val2-val3; 
	return val;
}

double _f2_(const double t, const double * p)
{
	double val0 = -t_coef*c2_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) + p[1]*p[1]*(1.0-p[1])*(1.0-p[1]) );
	double val1 = c2_coef*_v_laplace_(t,p);
	double val2 = _c2_(t,p)*_rho_(t,p);
	double val3 = _c2_grad_(t,p)[0]*_u_grad_(t,p)[0]+_c2_grad_(t,p)[1]*_u_grad_(t,p)[1];
	double val = val0-val1+val2+val3; 
	return val;
}

double _g_(const double * p)
{
	double val = exp(p[0]*p[1]);
	val = 1.0*val;
	return val;
}

// */

// Example 2D f=0 Example 3 old editing 
/*

#define a_coef 1.0
#define t_coef (1.0e-4)
#define c1_coef (1.5e-4)   // C1,C2 coefficients
#define c2_coef (1.0e-4)   // C1,C2 coefficients
#define c1_up (c1_coef)   // move c up
#define c2_up (c2_coef)   // move c up
#define q_coef 1.0    // to test q_i
#define u_coef (2.0e-4)   // phi coefficient 


double _L_coef_()       // used in readData for Pnt[j][0]=Pnt[j][0]*_L_coef_()+_L1_coef_();
{
	double val= 1.0*PI;
	return val;
}

double _L1_coef_()
{
	double val= 0.0;
	return val;
}


double _rho_(const double t, const double * p)
{
	double val = -u_coef*exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	      val += -u_coef*exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	      val += _c2_(t,p)-_c1_(t,p);
	return val;
}

double _u_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = u_coef*exp(-t_coef*t)*( p[0]*p[0]*(1.0-p[0])*(1.0-p[0]) + p[1]*p[1]*(1.0-p[1])*(1.0-p[1]) );
	return val;
}

double _v_laplace_(const double t, const double * p)
{
	double val = exp(-t_coef*t)*(  2.0*p[0]*p[0]+2.0*(1.0-p[0])*(1.0-p[0]) - 8.0*p[0]*(1.0-p[0])  );
	      val += exp(-t_coef*t)*(  2.0*p[1]*p[1]+2.0*(1.0-p[1])*(1.0-p[1]) - 8.0*p[1]*(1.0-p[1])  );
	return val;
}

std::vector<double> _u_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=u_coef*exp(-t_coef*t)*( 2.0*p[0]*(1.0-p[0])*(1.0-2.0*p[0]) );
    val[1]=u_coef*exp(-t_coef*t)*( 2.0*p[1]*(1.0-p[1])*(1.0-2.0*p[1]) );
    return val;
}

double _c1_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c1_coef*a_coef*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1])+c1_up;
	val = 1.0*val;
	return val;
}

double _c2_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c2_coef*a_coef*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1])+c2_up;
	val = 1.0*val;
	return val;
}

std::vector<double> _c1_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=-c1_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
    val[1]=-c1_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
    return val;
}

std::vector<double> _c2_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=-c2_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
    val[1]=-c2_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
    return val;
}



double _u_0_(const double * p)
{
	double val = _u_(0.0,p);
	return val;
}

double _c1_0_(const double * p)
{
	double val = _c1_(0.0,p);
	return val;
}

double _c2_0_(const double * p)
{
	double val = _c2_(0.0,p);
	return val;
}

double _f1_(const double t, const double * p)
{
	double val0 = exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1]);
	double val1 = exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
	double val2 = exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
	double val = (2.0*pow(a_coef,4)-t_coef*a_coef*a_coef)*c1_coef*val0+q_coef*_c1_(t,p)*_rho_(t,p) - q_coef*u_coef*c1_coef*pow(a_coef,4)*(val1*val1+val2*val2); 
	return val;
}

double _f2_(const double t, const double * p)
{
	double val0 = exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1]);
	double val1 = exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
	double val2 = exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
	double val = (2.0*pow(a_coef,4)-t_coef*a_coef*a_coef)*c2_coef*val0-q_coef*_c2_(t,p)*_rho_(t,p) + q_coef*u_coef*c2_coef*pow(a_coef,4)*(val1*val1+val2*val2); 
	return val;
}

double _g_(const double * p)
{
	double val = exp(p[0]*p[1]);
	val = 1.0*val;
	return val;
}

// */

// Example 2D f=0 Example 3 not good
/*

#define a_coef 1.0
#define t_coef (1.0e-4)
#define c1_coef (2.0e-4)   // C1,C2 coefficients
#define c2_coef (1.0e-4)   // C1,C2 coefficients
#define c1_up (0.0*c1_coef)   // move c up
#define c2_up (0.0*c2_coef)   // move c up
#define q_coef 1.0    // to test q_i
#define u_coef (1.0e-4)   // phi coefficient 


double _L_coef_()       // used in readData for Pnt[j][0]=Pnt[j][0]*_L_coef_()+_L1_coef_();
{
	double val= 1.0*PI;
	return val;
}

double _L1_coef_()
{
	double val= 0.0;
	return val;
}


double _rho_(const double t, const double * p)
{
	double val= u_coef*8.0*(a_coef*a_coef)*exp(-t_coef*t)*cos(2.0*a_coef*p[0])*cos(2.0*a_coef*p[1])+_c2_(t,p)-_c1_(t,p);
	return val;
}

double _u_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = u_coef*exp(-t_coef*t)*cos(2.0*a_coef*p[0])*cos(2.0*a_coef*p[1]);
	return val;
}

std::vector<double> _u_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=-u_coef*2.0*a_coef*exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
    val[1]=-u_coef*2.0*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
    return val;
}

double _c1_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = exp(-t_coef*t)*c1_coef*a_coef*a_coef*(cos(a_coef*p[0])*cos(a_coef*p[1]))+c1_coef*a_coef*a_coef;
	val = 1.0*val;
	return val;
}

double _c2_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = exp(-t_coef*t)*c2_coef*a_coef*a_coef*(cos(a_coef*p[0])*cos(a_coef*p[1]))+c2_coef*a_coef*a_coef;
	val = 1.0*val;
	return val;
}

std::vector<double> _c1_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=-c1_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
    val[1]=-c1_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
    return val;
}

std::vector<double> _c2_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=-c2_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
    val[1]=-c2_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
    return val;
}



double _u_0_(const double * p)
{
	double val = _u_(0.0,p);
	return val;
}

double _c1_0_(const double * p)
{
	double val = _c1_(0.0,p);
	return val;
}

double _c2_0_(const double * p)
{
	double val = _c2_(0.0,p);
	return val;
}

double _f1_(const double t, const double * p)
{
	double val0 = exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1]);
	double val1 = exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
	double val2 = exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
//	double val = (2.0*pow(a_coef,4))*c1_coef*val0-t_coef*_c1_(t,p)  +q_coef*_c1_(t,p)*_rho_(t,p) - q_coef*u_coef*c1_coef*pow(a_coef,4)*(val1*val1+val2*val2); 
	double val = (2.0*pow(a_coef,4)-t_coef*a_coef*a_coef)*c1_coef*val0+q_coef*_c1_(t,p)*_rho_(t,p) - q_coef*u_coef*c1_coef*pow(a_coef,4)*(val1*val1+val2*val2); 
	return val;
}

double _f2_(const double t, const double * p)
{
	double val0 = exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1]);
	double val1 = exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
	double val2 = exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
	double val = (2.0*pow(a_coef,4)-t_coef*a_coef*a_coef)*c2_coef*val0-q_coef*_c2_(t,p)*_rho_(t,p) + q_coef*u_coef*c2_coef*pow(a_coef,4)*(val1*val1+val2*val2); 
	return val;
}

double _g_(const double * p)
{
	double val = exp(p[0]*p[1]);
	val = 1.0*val;
	return val;
}

// */

// Example 2D f=0 Example 3 old 
/*

#define a_coef 1.0
#define t_coef (1.0e-3)
#define c1_coef (1.0e-3)   // C1,C2 coefficients
#define c2_coef (1.0e-3)   // C1,C2 coefficients
#define c1_up (c1_coef)   // move c up
#define c2_up (c2_coef)   // move c up
#define q_coef 1.0    // to test q_i
#define u_coef (1.0e-3)   // phi coefficient 


double _L_coef_()       // used in readData for Pnt[j][0]=Pnt[j][0]*_L_coef_()+_L1_coef_();
{
	double val= 1.0*PI;
	return val;
}

double _L1_coef_()
{
	double val= 0.0;
	return val;
}


double _rho_(const double t, const double * p)
{
	double val= u_coef*2.0*(a_coef*a_coef)*exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1])+_c2_(t,p)-_c1_(t,p);
	return val;
}

double _u_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = u_coef*exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1]);
	return val;
}

std::vector<double> _u_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=-u_coef*a_coef*exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
    val[1]=-u_coef*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
    return val;
}

double _c1_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c1_coef*a_coef*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1])+c1_up;
	val = 1.0*val;
	return val;
}

double _c2_(const double t, const double * p)
{
//	std::cout<<"c_coef= "<<c_coef<<std::endl;
	double val = c2_coef*a_coef*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1])+c2_up;
	val = 1.0*val;
	return val;
}

std::vector<double> _c1_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=-c1_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
    val[1]=-c1_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
    return val;
}

std::vector<double> _c2_grad_(const double t, const double * p)
{
    std::vector<double> val(2);
    val[0]=-c2_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
    val[1]=-c2_coef*a_coef*a_coef*a_coef*exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
    return val;
}



double _u_0_(const double * p)
{
	double val = _u_(0.0,p);
	return val;
}

double _c1_0_(const double * p)
{
	double val = _c1_(0.0,p);
	return val;
}

double _c2_0_(const double * p)
{
	double val = _c2_(0.0,p);
	return val;
}

double _f1_(const double t, const double * p)
{
	double val0 = exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1]);
	double val1 = exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
	double val2 = exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
	double val = (2.0*pow(a_coef,4)-t_coef*a_coef*a_coef)*c1_coef*val0+q_coef*_c1_(t,p)*_rho_(t,p) - q_coef*u_coef*c1_coef*pow(a_coef,4)*(val1*val1+val2*val2); 
	return val;
}

double _f2_(const double t, const double * p)
{
	double val0 = exp(-t_coef*t)*cos(a_coef*p[0])*cos(a_coef*p[1]);
	double val1 = exp(-t_coef*t)*sin(a_coef*p[0])*cos(a_coef*p[1]);
	double val2 = exp(-t_coef*t)*cos(a_coef*p[0])*sin(a_coef*p[1]);
	double val = (2.0*pow(a_coef,4)-t_coef*a_coef*a_coef)*c2_coef*val0-q_coef*_c2_(t,p)*_rho_(t,p) + q_coef*u_coef*c2_coef*pow(a_coef,4)*(val1*val1+val2*val2); 
	return val;
}

double _g_(const double * p)
{
	double val = exp(p[0]*p[1]);
	val = 1.0*val;
	return val;
}

// */

