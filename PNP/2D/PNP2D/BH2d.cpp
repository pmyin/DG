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

///////////////////////////////////////////////////////////////
MyPoint::MyPoint()
{
  for (int i = 0;i < 2;i ++)
    x[i] = 0;
}

MyPoint::MyPoint(const double * data)
{
  for (int i = 0;i < 2;i ++)
    x[i] = data[i];
}

MyPoint::MyPoint(const MyPoint& p)
{
  int i;
  for (i = 0;i < 2;i ++)
    x[i] = p.x[i];
}

MyPoint::MyPoint(double d, ...)
{
  x[0] = d;
  va_list ap;
  va_start(ap, d);
  for (int i = 1;i < 2;i ++) {
    x[i] = va_arg(ap, double);
  }
  va_end(ap);
}

MyPoint::~MyPoint()
{}

MyPoint::operator const double *() const
{
  return x;
}

MyPoint::operator double *()
{
  return x;
}

MyPoint& MyPoint::operator=(const MyPoint& p)
{
  int i;
  for (i = 0;i < 2;i ++)
    x[i] = p.x[i];
  return *this;
}

const double& MyPoint::operator[](int i) const
{
    return x[i];
}

double& MyPoint::operator[](int i)
{
  return x[i];
}

double MyPoint::length() const
{
  double v = 0.;
  for (int i = 0;i < 2;i ++)
    v += x[i]*x[i];
  return sqrt(v);
}

MyPoint& MyPoint::operator+=(const MyPoint& p)
{
  for (int i = 0;i < 2;i ++)
    x[i] += p.x[i];
  return *this;
}

MyPoint& MyPoint::operator-=(const MyPoint& p)
{
  for (int i = 0;i < 2;i ++)
    x[i] -= p.x[i];
  return *this;
}

MyPoint& MyPoint::operator*=(const double& s)
{
  for (int i = 0;i < 2;i ++)
    x[i] *= s;
  return *this;
}

MyPoint& MyPoint::operator/=(const double& s)
{
  for (int i = 0;i < 2;i ++)
    x[i] /= s;
  return *this;
}

MyPoint midpoint(const MyPoint& p1, const MyPoint& p2)
{
  double p[2];
  for (int i = 0;i < 2;i ++)
    p[i] = (p1.x[i] + p2.x[i])/2.;
  return p;
}

double distance(const MyPoint& p1, const MyPoint& p2)
{
  double d = 0.0;
  for (int i = 0;i < 2;i ++)
    d += (p1.x[i] - p2.x[i])*(p1.x[i] - p2.x[i]);
  return sqrt(d);
}

MyPoint barycenter(const std::vector<MyPoint>& p, const double * w)
{
  int i, j, k;
  double bc[2];
  k = p.size();
  if (w == NULL) {
    for (i = 0;i < 2;i ++) {
      bc[i] = 0;
      for (j = 0;j < k;j ++)
	bc[i] += p[j][i];
      bc[i] /= k;
    }
  }
  else {
    double sw = 0;
    for (j = 0;j < k;j ++) sw += w[j];
    for (i = 0;i < 2;k ++) {
      bc[i] = 0;
      for (j = 0;j < k;j ++)
	bc[i] += w[j]*p[j][i];
      bc[i] /= sw;
    }
  }
  return bc;
}

MyPoint operator+(const MyPoint& p1, const MyPoint& p2)
{
  double p[2];
  for (int i = 0;i < 2;i ++)
    p[i] = p1.x[i] + p2.x[i];
  return p;
}

MyPoint operator-(const MyPoint& p1, const MyPoint& p2)
{
  double p[2];
  for (int i = 0;i < 2;i ++)
    p[i] = p1.x[i] - p2.x[i];
  return p;
}

std::istream& operator>>(std::istream& is, MyPoint& p)
{
  for (int i = 0;i < 2;i ++)
    is >> p.x[i];
  return is;
}

std::ostream& operator<<(std::ostream& os, const MyPoint& p)
{
  for (int i = 0;i < 2;i ++)
    os << p.x[i] << "\t";
  return os;
}

//////////////////////////////////////////单元类，即二维情况
void TmpEle::buildTE(int i)/*i与Gauss精度有关*/
{
    Pnt.resize(4);/*anticlockwise*/
    Pnt[0][0]=-1.0;
    Pnt[0][1]=-1.0;
    Pnt[1][0]=1.0;
    Pnt[1][1]=-1.0;
    Pnt[2][0]=1.0;
    Pnt[2][1]=1.0;
    Pnt[3][0]=-1.0;
    Pnt[3][1]=1.0;
    
    /*n点Gauss求积公式具有2*n-1次代数精度*/
    
    switch(i)
    {
        case 1:/*单方向上1个Gauss点，代数精度为1*/
            GaussPnt.resize(1);
            GaussWeight.resize(1);
            GaussPnt[0][0]=0.0;
            GaussPnt[0][1]=0.0;
            
            GaussWeight[0]=4.0;
            break;
        case 2:/*单方向上2个Gauss点，代数精度为3*/
            GaussPnt.resize(4);
            GaussWeight.resize(4);
            GaussPnt[0][0]=-0.5773502691896257;
            GaussPnt[1][0]=-0.5773502691896257;
            GaussPnt[2][0]=0.5773502691896257;
            GaussPnt[3][0]=0.5773502691896257;
            GaussPnt[0][1]=-0.5773502691896257;
            GaussPnt[1][1]=0.5773502691896257;
            GaussPnt[2][1]=-0.5773502691896257;
            GaussPnt[3][1]=0.5773502691896257;
            
            GaussWeight[0]=1.0;
            GaussWeight[1]=1.0;
            GaussWeight[2]=1.0;
            GaussWeight[3]=1.0;
            break;
        case 3:/*单方向上3个Gauss点，代数精度为5*/
            GaussPnt.resize(9);
            GaussWeight.resize(9);
            GaussPnt[0][0]=-0.7745966692414834;
            GaussPnt[1][0]=-0.7745966692414834;
            GaussPnt[2][0]=-0.7745966692414834;
            GaussPnt[3][0]=0.0000000000000000;
            GaussPnt[4][0]=0.0000000000000000;
            GaussPnt[5][0]=0.0000000000000000;
            GaussPnt[6][0]=0.7745966692414834;
            GaussPnt[7][0]=0.7745966692414834;
            GaussPnt[8][0]=0.7745966692414834;
            
            GaussPnt[0][1]=-0.7745966692414834;
            GaussPnt[1][1]=0.0000000000000000;
            GaussPnt[2][1]=0.7745966692414834;
            GaussPnt[3][1]=-0.7745966692414834;
            GaussPnt[4][1]=0.0000000000000000;
            GaussPnt[5][1]=0.7745966692414834;
            GaussPnt[6][1]=-0.7745966692414834;
            GaussPnt[7][1]=0.0000000000000000;
            GaussPnt[8][1]=0.7745966692414834;
            
            GaussWeight[0]=0.5555555555555556*0.5555555555555556;
            GaussWeight[1]=0.5555555555555556*0.8888888888888888;
            GaussWeight[2]=0.5555555555555556*0.5555555555555556;
            GaussWeight[3]=0.8888888888888888*0.5555555555555556;
            GaussWeight[4]=0.8888888888888888*0.8888888888888888;
            GaussWeight[5]=0.8888888888888888*0.5555555555555556; 
            GaussWeight[6]=0.5555555555555556*0.5555555555555556;  
            GaussWeight[7]=0.5555555555555556*0.8888888888888888; 
            GaussWeight[8]=0.5555555555555556*0.5555555555555556;
            break;
        case 4:/*单方向上4个Gauss点，代数精度为7*/
            GaussPnt.resize(16);
            GaussWeight.resize(16);
            GaussPnt[0][0]=-0.3399810435848563;
            GaussPnt[1][0]=-0.3399810435848563;
            GaussPnt[2][0]=-0.3399810435848563;
            GaussPnt[3][0]=-0.3399810435848563;
            
            GaussPnt[4][0]=-0.8611363115940526;
            GaussPnt[5][0]=-0.8611363115940526;
            GaussPnt[6][0]=-0.8611363115940526;
            GaussPnt[7][0]=-0.8611363115940526;
            
            GaussPnt[8][0]=0.3399810435848563;
            GaussPnt[9][0]=0.3399810435848563;
            GaussPnt[10][0]=0.3399810435848563;
            GaussPnt[11][0]=0.3399810435848563;
            
            GaussPnt[12][0]=0.8611363115940526;
            GaussPnt[13][0]=0.8611363115940526;
            GaussPnt[14][0]=0.8611363115940526;
            GaussPnt[15][0]=0.8611363115940526;
            
            GaussPnt[0][1]=-0.3399810435848563;
            GaussPnt[1][1]=-0.8611363115940526;
            GaussPnt[2][1]=0.3399810435848563;
            GaussPnt[3][1]=0.8611363115940526;
            
            GaussPnt[4][1]=-0.3399810435848563;
            GaussPnt[5][1]=-0.8611363115940526;
            GaussPnt[6][1]=0.3399810435848563;
            GaussPnt[7][1]=0.8611363115940526;
            
            GaussPnt[8][1]=-0.3399810435848563;
            GaussPnt[9][1]=-0.8611363115940526;
            GaussPnt[10][1]=0.3399810435848563;
            GaussPnt[11][1]=0.8611363115940526;
            
            GaussPnt[12][1]=-0.3399810435848563;
            GaussPnt[13][1]=-0.8611363115940526;
            GaussPnt[14][1]=0.3399810435848563;
            GaussPnt[15][1]=0.8611363115940526;
            
            GaussWeight[0]=0.6521451548625461*0.6521451548625461;
            GaussWeight[1]=0.6521451548625461*0.3478548451374538;
            GaussWeight[2]=0.6521451548625461*0.6521451548625461;
            GaussWeight[3]=0.6521451548625461*0.3478548451374538;
            
            GaussWeight[4]=0.3478548451374538*0.6521451548625461;
            GaussWeight[5]=0.3478548451374538*0.3478548451374538; 
            GaussWeight[6]=0.3478548451374538*0.6521451548625461;  
            GaussWeight[7]=0.3478548451374538*0.3478548451374538; 
            
            GaussWeight[8]=0.6521451548625461*0.6521451548625461; 
            GaussWeight[9]=0.6521451548625461*0.3478548451374538;
            GaussWeight[10]=0.6521451548625461*0.6521451548625461;
            GaussWeight[11]=0.6521451548625461*0.3478548451374538;
            
            GaussWeight[12]=0.3478548451374538*0.6521451548625461; 
            GaussWeight[13]=0.3478548451374538*0.3478548451374538;  
            GaussWeight[14]=0.3478548451374538*0.6521451548625461; 
            GaussWeight[15]=0.3478548451374538*0.3478548451374538; 
            break;
        case 5:/*单方向上5个Gauss点，代数精度为9*/
            GaussPnt.resize(25);
            GaussWeight.resize(25);
            GaussPnt[0][0]=-0.5384693101056831;
            GaussPnt[1][0]=-0.5384693101056831;
            GaussPnt[2][0]=-0.5384693101056831;
            GaussPnt[3][0]=-0.5384693101056831;
            GaussPnt[4][0]=-0.5384693101056831;
            
            GaussPnt[5][0]=-0.90617984593866401;
            GaussPnt[6][0]=-0.9061798459386640;
            GaussPnt[7][0]=-0.9061798459386640;
            GaussPnt[8][0]=-0.9061798459386640;
            GaussPnt[9][0]=-0.9061798459386640;
            
            GaussPnt[10][0]=0.0000000000000000;
            GaussPnt[11][0]=0.0000000000000000;
            GaussPnt[12][0]=0.0000000000000000;
            GaussPnt[13][0]=0.0000000000000000;
            GaussPnt[14][0]=0.0000000000000000;
            
            GaussPnt[15][0]=0.5384693101056831;
            GaussPnt[16][0]=0.5384693101056831;
            GaussPnt[17][0]=0.5384693101056831;
            GaussPnt[18][0]=0.5384693101056831;
            GaussPnt[19][0]=0.5384693101056831;
            
            GaussPnt[20][0]=0.9061798459386640;
            GaussPnt[21][0]=0.9061798459386640;
            GaussPnt[22][0]=0.9061798459386640;
            GaussPnt[23][0]=0.9061798459386640;
            GaussPnt[24][0]=0.9061798459386640;
            
            GaussPnt[0][1]=-0.5384693101056831;
            GaussPnt[1][1]=-0.9061798459386640;
            GaussPnt[2][1]=0.0000000000000000;
            GaussPnt[3][1]=0.5384693101056831;
            GaussPnt[4][1]=0.9061798459386640;
            
            GaussPnt[5][1]=-0.5384693101056831;
            GaussPnt[6][1]=-0.9061798459386640;
            GaussPnt[7][1]=0.0000000000000000;
            GaussPnt[8][1]=0.5384693101056831;
            GaussPnt[9][1]=0.9061798459386640;
            
            GaussPnt[10][1]=-0.5384693101056831;
            GaussPnt[11][1]=-0.9061798459386640;
            GaussPnt[12][1]=0.0000000000000000;
            GaussPnt[13][1]=0.5384693101056831;
            GaussPnt[14][1]=0.9061798459386640;
            
            GaussPnt[15][1]=-0.5384693101056831;
            GaussPnt[16][1]=-0.9061798459386640;
            GaussPnt[17][1]=0.0000000000000000;
            GaussPnt[18][1]=0.5384693101056831;
            GaussPnt[19][1]=0.9061798459386640;
            
            GaussPnt[20][1]=-0.5384693101056831;
            GaussPnt[21][1]=-0.9061798459386640;
            GaussPnt[22][1]=0.0000000000000000;
            GaussPnt[23][1]=0.5384693101056831;
            GaussPnt[24][1]=0.9061798459386640;
            
            GaussWeight[0]=0.4786286704993665*0.4786286704993665;
            GaussWeight[1]=0.4786286704993665*0.2369268850561891;
            GaussWeight[2]=0.4786286704993665*0.5688888888888889;
            GaussWeight[3]=0.4786286704993665*0.4786286704993665;
            GaussWeight[4]=0.4786286704993665*0.2369268850561891;
            
            GaussWeight[5]=0.2369268850561891*0.4786286704993665;
            GaussWeight[6]=0.2369268850561891*0.2369268850561891;
            GaussWeight[7]=0.2369268850561891*0.5688888888888889; 
            GaussWeight[8]=0.2369268850561891*0.4786286704993665;
            GaussWeight[9]=0.2369268850561891*0.2369268850561891;  
            
            GaussWeight[10]=0.5688888888888889*0.4786286704993665;
            GaussWeight[11]=0.5688888888888889*0.2369268850561891;
            GaussWeight[12]=0.5688888888888889*0.5688888888888889;
            GaussWeight[13]=0.5688888888888889*0.4786286704993665;
            GaussWeight[14]=0.5688888888888889*0.2369268850561891;
            
            GaussWeight[15]=0.4786286704993665*0.4786286704993665; 
            GaussWeight[16]=0.4786286704993665*0.2369268850561891;  
            GaussWeight[17]=0.4786286704993665*0.5688888888888889; 
            GaussWeight[18]=0.4786286704993665*0.4786286704993665;
            GaussWeight[19]=0.4786286704993665*0.2369268850561891; 
            
            GaussWeight[20]=0.2369268850561891*0.4786286704993665;
            GaussWeight[21]=0.2369268850561891*0.2369268850561891;
            GaussWeight[22]=0.2369268850561891*0.5688888888888889;
            GaussWeight[23]=0.2369268850561891*0.4786286704993665;
            GaussWeight[24]=0.2369268850561891*0.2369268850561891;  
            break;

    case  10:
    GaussPnt.resize(100);
    GaussWeight.resize(100);
    GaussPnt[0][0] = -0.1488743389816312;
    GaussPnt[1][0] = -0.1488743389816312;
    GaussPnt[2][0] = -0.1488743389816312;
    GaussPnt[3][0] = -0.1488743389816312;
    GaussPnt[4][0] = -0.1488743389816312;
    GaussPnt[5][0] = -0.1488743389816312;
    GaussPnt[6][0] = -0.1488743389816312;
    GaussPnt[7][0] = -0.1488743389816312;
    GaussPnt[8][0] = -0.1488743389816312;
    GaussPnt[9][0] = -0.1488743389816312;

    GaussPnt[10][0] = 0.1488743389816312;
    GaussPnt[11][0] = 0.1488743389816312;
    GaussPnt[12][0] = 0.1488743389816312;
    GaussPnt[13][0] = 0.1488743389816312;
    GaussPnt[14][0] = 0.1488743389816312;
    GaussPnt[15][0] = 0.1488743389816312;
    GaussPnt[16][0] = 0.1488743389816312;
    GaussPnt[17][0] = 0.1488743389816312;
    GaussPnt[18][0] = 0.1488743389816312;
    GaussPnt[19][0] = 0.1488743389816312;

    GaussPnt[20][0] = -0.4333953941292472;
    GaussPnt[21][0] = -0.4333953941292472;
    GaussPnt[22][0] = -0.4333953941292472;
    GaussPnt[23][0] = -0.4333953941292472;
    GaussPnt[24][0] = -0.4333953941292472;
    GaussPnt[25][0] = -0.4333953941292472;
    GaussPnt[26][0] = -0.4333953941292472;
    GaussPnt[27][0] = -0.4333953941292472;
    GaussPnt[28][0] = -0.4333953941292472;
    GaussPnt[29][0] = -0.4333953941292472;

    GaussPnt[30][0] = 0.4333953941292472;
    GaussPnt[31][0] = 0.4333953941292472;
    GaussPnt[32][0] = 0.4333953941292472;
    GaussPnt[33][0] = 0.4333953941292472;
    GaussPnt[34][0] = 0.4333953941292472;
    GaussPnt[35][0] = 0.4333953941292472;
    GaussPnt[36][0] = 0.4333953941292472;
    GaussPnt[37][0] = 0.4333953941292472;
    GaussPnt[38][0] = 0.4333953941292472;
    GaussPnt[39][0] = 0.4333953941292472;

    GaussPnt[40][0] = -0.6794095682990244;
    GaussPnt[41][0] = -0.6794095682990244;
    GaussPnt[42][0] = -0.6794095682990244;
    GaussPnt[43][0] = -0.6794095682990244;
    GaussPnt[44][0] = -0.6794095682990244;
    GaussPnt[45][0] = -0.6794095682990244;
    GaussPnt[46][0] = -0.6794095682990244;
    GaussPnt[47][0] = -0.6794095682990244;
    GaussPnt[48][0] = -0.6794095682990244;
    GaussPnt[49][0] = -0.6794095682990244;

    GaussPnt[50][0] = 0.6794095682990244;
    GaussPnt[51][0] = 0.6794095682990244;
    GaussPnt[52][0] = 0.6794095682990244;
    GaussPnt[53][0] = 0.6794095682990244;
    GaussPnt[54][0] = 0.6794095682990244;
    GaussPnt[55][0] = 0.6794095682990244;
    GaussPnt[56][0] = 0.6794095682990244;
    GaussPnt[57][0] = 0.6794095682990244;
    GaussPnt[58][0] = 0.6794095682990244;
    GaussPnt[59][0] = 0.6794095682990244;

    GaussPnt[60][0] = -0.8650633666889845;
    GaussPnt[61][0] = -0.8650633666889845;
    GaussPnt[62][0] = -0.8650633666889845;
    GaussPnt[63][0] = -0.8650633666889845;
    GaussPnt[64][0] = -0.8650633666889845;
    GaussPnt[65][0] = -0.8650633666889845;
    GaussPnt[66][0] = -0.8650633666889845;
    GaussPnt[67][0] = -0.8650633666889845;
    GaussPnt[68][0] = -0.8650633666889845;
    GaussPnt[69][0] = -0.8650633666889845;

    GaussPnt[70][0] = 0.8650633666889845;
    GaussPnt[71][0] = 0.8650633666889845;
    GaussPnt[72][0] = 0.8650633666889845;
    GaussPnt[73][0] = 0.8650633666889845;
    GaussPnt[74][0] = 0.8650633666889845;
    GaussPnt[75][0] = 0.8650633666889845;
    GaussPnt[76][0] = 0.8650633666889845;
    GaussPnt[77][0] = 0.8650633666889845;
    GaussPnt[78][0] = 0.8650633666889845;
    GaussPnt[79][0] = 0.8650633666889845;

    GaussPnt[80][0] = -0.9739065285171717;
    GaussPnt[81][0] = -0.9739065285171717;
    GaussPnt[82][0] = -0.9739065285171717;
    GaussPnt[83][0] = -0.9739065285171717;
    GaussPnt[84][0] = -0.9739065285171717;
    GaussPnt[85][0] = -0.9739065285171717;
    GaussPnt[86][0] = -0.9739065285171717;
    GaussPnt[87][0] = -0.9739065285171717;
    GaussPnt[88][0] = -0.9739065285171717;
    GaussPnt[89][0] = -0.9739065285171717;

    GaussPnt[90][0] = 0.9739065285171717;
    GaussPnt[91][0] = 0.9739065285171717;
    GaussPnt[92][0] = 0.9739065285171717;
    GaussPnt[93][0] = 0.9739065285171717;
    GaussPnt[94][0] = 0.9739065285171717;
    GaussPnt[95][0] = 0.9739065285171717;
    GaussPnt[96][0] = 0.9739065285171717;
    GaussPnt[97][0] = 0.9739065285171717;
    GaussPnt[98][0] = 0.9739065285171717;
    GaussPnt[99][0] = 0.9739065285171717;

    GaussPnt[0][1] = -0.1488743389816312;
    GaussPnt[1][1] = 0.1488743389816312;
    GaussPnt[2][1] = -0.4333953941292472;
    GaussPnt[3][1] = 0.4333953941292472;
    GaussPnt[4][1] = -0.6794095682990244;
    GaussPnt[5][1] = 0.6794095682990244;
    GaussPnt[6][1] = -0.8650633666889845;
    GaussPnt[7][1] = 0.8650633666889845;
    GaussPnt[8][1] = -0.9739065285171717;
    GaussPnt[9][1] = 0.9739065285171717;

    GaussPnt[10][1] = -0.1488743389816312;
    GaussPnt[11][1] = 0.1488743389816312;
    GaussPnt[12][1] = -0.4333953941292472;
    GaussPnt[13][1] = 0.4333953941292472;
    GaussPnt[14][1] = -0.6794095682990244;
    GaussPnt[15][1] = 0.6794095682990244;
    GaussPnt[16][1] = -0.8650633666889845;
    GaussPnt[17][1] = 0.8650633666889845;
    GaussPnt[18][1] = -0.9739065285171717;
    GaussPnt[19][1] = 0.9739065285171717;

    GaussPnt[20][1] = -0.1488743389816312;
    GaussPnt[21][1] = 0.1488743389816312;
    GaussPnt[22][1] = -0.4333953941292472;
    GaussPnt[23][1] = 0.4333953941292472;
    GaussPnt[24][1] = -0.6794095682990244;
    GaussPnt[25][1] = 0.6794095682990244;
    GaussPnt[26][1] = -0.8650633666889845;
    GaussPnt[27][1] = 0.8650633666889845;
    GaussPnt[28][1] = -0.9739065285171717;
    GaussPnt[29][1] = 0.9739065285171717;

    GaussPnt[30][1] = -0.1488743389816312;
    GaussPnt[31][1] = 0.1488743389816312;
    GaussPnt[32][1] = -0.4333953941292472;
    GaussPnt[33][1] = 0.4333953941292472;
    GaussPnt[34][1] = -0.6794095682990244;
    GaussPnt[35][1] = 0.6794095682990244;
    GaussPnt[36][1] = -0.8650633666889845;
    GaussPnt[37][1] = 0.8650633666889845;
    GaussPnt[38][1] = -0.9739065285171717;
    GaussPnt[39][1] = 0.9739065285171717;

    GaussPnt[40][1] = -0.1488743389816312;
    GaussPnt[41][1] = 0.1488743389816312;
    GaussPnt[42][1] = -0.4333953941292472;
    GaussPnt[43][1] = 0.4333953941292472;
    GaussPnt[44][1] = -0.6794095682990244;
    GaussPnt[45][1] = 0.6794095682990244;
    GaussPnt[46][1] = -0.8650633666889845;
    GaussPnt[47][1] = 0.8650633666889845;
    GaussPnt[48][1] = -0.9739065285171717;
    GaussPnt[49][1] = 0.9739065285171717;

    GaussPnt[50][1] = -0.1488743389816312;
    GaussPnt[51][1] = 0.1488743389816312;
    GaussPnt[52][1] = -0.4333953941292472;
    GaussPnt[53][1] = 0.4333953941292472;
    GaussPnt[54][1] = -0.6794095682990244;
    GaussPnt[55][1] = 0.6794095682990244;
    GaussPnt[56][1] = -0.8650633666889845;
    GaussPnt[57][1] = 0.8650633666889845;
    GaussPnt[58][1] = -0.9739065285171717;
    GaussPnt[59][1] = 0.9739065285171717;

    GaussPnt[60][1] = -0.1488743389816312;
    GaussPnt[61][1] = 0.1488743389816312;
    GaussPnt[62][1] = -0.4333953941292472;
    GaussPnt[63][1] = 0.4333953941292472;
    GaussPnt[64][1] = -0.6794095682990244;
    GaussPnt[65][1] = 0.6794095682990244;
    GaussPnt[66][1] = -0.8650633666889845;
    GaussPnt[67][1] = 0.8650633666889845;
    GaussPnt[68][1] = -0.9739065285171717;
    GaussPnt[69][1] = 0.9739065285171717;

    GaussPnt[70][1] = -0.1488743389816312;
    GaussPnt[71][1] = 0.1488743389816312;
    GaussPnt[72][1] = -0.4333953941292472;
    GaussPnt[73][1] = 0.4333953941292472;
    GaussPnt[74][1] = -0.6794095682990244;
    GaussPnt[75][1] = 0.6794095682990244;
    GaussPnt[76][1] = -0.8650633666889845;
    GaussPnt[77][1] = 0.8650633666889845;
    GaussPnt[78][1] = -0.9739065285171717;
    GaussPnt[79][1] = 0.9739065285171717;

    GaussPnt[80][1] = -0.1488743389816312;
    GaussPnt[81][1] = 0.1488743389816312;
    GaussPnt[82][1] = -0.4333953941292472;
    GaussPnt[83][1] = 0.4333953941292472;
    GaussPnt[84][1] = -0.6794095682990244;
    GaussPnt[85][1] = 0.6794095682990244;
    GaussPnt[86][1] = -0.8650633666889845;
    GaussPnt[87][1] = 0.8650633666889845;
    GaussPnt[88][1] = -0.9739065285171717;
    GaussPnt[89][1] = 0.9739065285171717;

    GaussPnt[90][1] = -0.1488743389816312;
    GaussPnt[91][1] = 0.1488743389816312;
    GaussPnt[92][1] = -0.4333953941292472;
    GaussPnt[93][1] = 0.4333953941292472;
    GaussPnt[94][1] = -0.6794095682990244;
    GaussPnt[95][1] = 0.6794095682990244;
    GaussPnt[96][1] = -0.8650633666889845;
    GaussPnt[97][1] = 0.8650633666889845;
    GaussPnt[98][1] = -0.9739065285171717;
    GaussPnt[99][1] = 0.9739065285171717;

    GaussWeight[0] = 0.08733456739325577;
    GaussWeight[1] = 0.08733456739325577;
    GaussWeight[2] = 0.07957483846557164;
    GaussWeight[3] = 0.07957483846557164;
    GaussWeight[4] = 0.06474532742811088;
    GaussWeight[5] = 0.06474532742811088;
    GaussWeight[6] = 0.04416649409029918;
    GaussWeight[7] = 0.04416649409029918;
    GaussWeight[8] = 0.0197029973375154;
    GaussWeight[9] = 0.0197029973375154;
    GaussWeight[10] = 0.08733456739325577;
    GaussWeight[11] = 0.08733456739325577;
    GaussWeight[12] = 0.07957483846557164;
    GaussWeight[13] = 0.07957483846557164;
    GaussWeight[14] = 0.06474532742811088;
    GaussWeight[15] = 0.06474532742811088;
    GaussWeight[16] = 0.04416649409029918;
    GaussWeight[17] = 0.04416649409029918;
    GaussWeight[18] = 0.0197029973375154;
    GaussWeight[19] = 0.0197029973375154;
    GaussWeight[20] = 0.07957483846557164;
    GaussWeight[21] = 0.07957483846557164;
    GaussWeight[22] = 0.07250456612796833;
    GaussWeight[23] = 0.07250456612796833;
    GaussWeight[24] = 0.05899266608023902;
    GaussWeight[25] = 0.05899266608023902;
    GaussWeight[26] = 0.04024227448222964;
    GaussWeight[27] = 0.04024227448222964;
    GaussWeight[28] = 0.01795237415398764;
    GaussWeight[29] = 0.01795237415398764;
    GaussWeight[30] = 0.07957483846557164;
    GaussWeight[31] = 0.07957483846557164;
    GaussWeight[32] = 0.07250456612796833;
    GaussWeight[33] = 0.07250456612796833;
    GaussWeight[34] = 0.05899266608023902;
    GaussWeight[35] = 0.05899266608023902;
    GaussWeight[36] = 0.04024227448222964;
    GaussWeight[37] = 0.04024227448222964;
    GaussWeight[38] = 0.01795237415398764;
    GaussWeight[39] = 0.01795237415398764;
    GaussWeight[40] = 0.06474532742811088;
    GaussWeight[41] = 0.06474532742811088;
    GaussWeight[42] = 0.05899266608023902;
    GaussWeight[43] = 0.05899266608023902;
    GaussWeight[44] = 0.04799883424048428;
    GaussWeight[45] = 0.04799883424048428;
    GaussWeight[46] = 0.0327427524585067;
    GaussWeight[47] = 0.0327427524585067;
    GaussWeight[48] = 0.01460678230864109;
    GaussWeight[49] = 0.01460678230864109;
    GaussWeight[50] = 0.06474532742811088;
    GaussWeight[51] = 0.06474532742811088;
    GaussWeight[52] = 0.05899266608023902;
    GaussWeight[53] = 0.05899266608023902;
    GaussWeight[54] = 0.04799883424048428;
    GaussWeight[55] = 0.04799883424048428;
    GaussWeight[56] = 0.0327427524585067;
    GaussWeight[57] = 0.0327427524585067;
    GaussWeight[58] = 0.01460678230864109;
    GaussWeight[59] = 0.01460678230864109;
    GaussWeight[60] = 0.04416649409029918;
    GaussWeight[61] = 0.04416649409029918;
    GaussWeight[62] = 0.04024227448222964;
    GaussWeight[63] = 0.04024227448222964;
    GaussWeight[64] = 0.0327427524585067;
    GaussWeight[65] = 0.0327427524585067;
    GaussWeight[66] = 0.02233570576292875;
    GaussWeight[67] = 0.02233570576292875;
    GaussWeight[68] = 0.00996412235661632;
    GaussWeight[69] = 0.00996412235661632;
    GaussWeight[70] = 0.04416649409029918;
    GaussWeight[71] = 0.04416649409029918;
    GaussWeight[72] = 0.04024227448222964;
    GaussWeight[73] = 0.04024227448222964;
    GaussWeight[74] = 0.0327427524585067;
    GaussWeight[75] = 0.0327427524585067;
    GaussWeight[76] = 0.02233570576292875;
    GaussWeight[77] = 0.02233570576292875;
    GaussWeight[78] = 0.00996412235661632;
    GaussWeight[79] = 0.00996412235661632;
    GaussWeight[80] = 0.0197029973375154;
    GaussWeight[81] = 0.0197029973375154;
    GaussWeight[82] = 0.01795237415398764;
    GaussWeight[83] = 0.01795237415398764;
    GaussWeight[84] = 0.01460678230864109;
    GaussWeight[85] = 0.01460678230864109;
    GaussWeight[86] = 0.00996412235661632;
    GaussWeight[87] = 0.00996412235661632;
    GaussWeight[88] = 0.004445068151927637;
    GaussWeight[89] = 0.004445068151927637;
    GaussWeight[90] = 0.0197029973375154;
    GaussWeight[91] = 0.0197029973375154;
    GaussWeight[92] = 0.01795237415398764;
    GaussWeight[93] = 0.01795237415398764;
    GaussWeight[94] = 0.01460678230864109;
    GaussWeight[95] = 0.01460678230864109;
    GaussWeight[96] = 0.00996412235661632;
    GaussWeight[97] = 0.00996412235661632;
    GaussWeight[98] = 0.004445068151927637;
    GaussWeight[99] = 0.004445068151927637;
    break;					
			
			
        default:
            std::cout<<"There are no Gauss quadrature points with accuracy of"<<i<<"-th order"<<std::endl;

    }
    
	int xPntNO=5;
	int PntNO=xPntNO*xPntNO;
	MeshPnt.resize(PntNO);
	double dx = 2.0/(xPntNO-1);
	for(int i=0; i<xPntNO; i++)
		for(int j=0; j<xPntNO; j++)
		{
			MeshPnt[xPntNO*i+j][0]=-1.0+j*dx;
			MeshPnt[xPntNO*i+j][1]=-1.0+i*dx;
		}    
    
    
}
 
MyPoint TmpEle::getPnt(int i)
{
  return Pnt[i];
} 
 
double TmpEle::getVolume()
{
  return 4.0;
} 

std::vector<MyPoint> TmpEle::getGaussPnt()
{
   return GaussPnt;
}

std::vector<double> TmpEle::getGaussWeight()
{ 
   return GaussWeight;
}

std::vector<MyPoint> TmpEle::getMeshPnt()
{
   return MeshPnt;
}

MyPoint TmpEle::Local_to_Global(const MyPoint & lp, const std::vector<MyPoint> & gv) const
{
  MyPoint gp;
  gp[0] = (gv[0][0] + gv[2][0])/2.0 + lp[0]*(gv[2][0]-gv[0][0])/2.0;
  gp[1] = (gv[0][1] + gv[2][1])/2.0 + lp[1]*(gv[2][1]-gv[0][1])/2.0;
	return gp;
}

MyPoint TmpEle::Global_to_Local(const MyPoint & gp, const std::vector<MyPoint> & gv) const
{
  MyPoint lp;
  lp[0] = (2.0*gp[0]-gv[0][0]-gv[2][0])/(gv[2][0]-gv[0][0]);
  lp[1] = (2.0*gp[1]-gv[0][1]-gv[2][1])/(gv[2][1]-gv[0][1]);
	return lp;
}

double TmpEle::Local_to_Global_jacobian(const MyPoint & lp, const std::vector<MyPoint> & gv) const
{
  double larea = (Pnt[2][0]-Pnt[0][0])*(Pnt[2][1]-Pnt[0][1]);
  double garea = (gv[2][0]-gv[0][0])*(gv[2][1]-gv[0][1]);
  return garea/larea;
}

double TmpEle::Global_to_Local_jacobian(const MyPoint& gp, const std::vector<MyPoint> & gv) const
{
  double larea = (Pnt[2][0]-Pnt[0][0])*(Pnt[2][1]-Pnt[0][1]);
  double garea = (gv[2][0]-gv[0][0])*(gv[2][1]-gv[0][1]);
  return larea/garea;
}

std::vector<MyPoint> TmpEle::Local_to_Global(const std::vector<MyPoint>& lp, const std::vector<MyPoint> & gv) const
{
        std::vector<MyPoint> gp(lp.size());
	for(int i=0; i<gp.size(); i++){
 	 	gp[i][0] = (gv[0][0] + gv[2][0])/2.0 + lp[i][0]*(gv[2][0]-gv[0][0])/2.0;
  		gp[i][1] = (gv[0][1] + gv[2][1])/2.0 + lp[i][1]*(gv[2][1]-gv[0][1])/2.0;
	}
	return gp;
}

std::vector<MyPoint> TmpEle::Global_to_Local(const std::vector<MyPoint>& gp, const std::vector<MyPoint> & gv) const
{
        std::vector<MyPoint> lp(gp.size());
	for(int i=0; i<gp.size(); i++){
  		lp[i][0] = (2.0*gp[i][0]-gv[0][0]-gv[2][0])/(gv[2][0]-gv[0][0]);
  		lp[i][1] = (2.0*gp[i][1]-gv[0][1]-gv[2][1])/(gv[2][1]-gv[0][1]);
	}
	return lp;
}

std::vector<double> TmpEle::Local_to_Global_jacobian(const std::vector<MyPoint>& lp, const std::vector<MyPoint> & gv) const
{
        std::vector<double> gj(lp.size());     
  	double larea = (Pnt[2][0]-Pnt[0][0])*(Pnt[2][1]-Pnt[0][1]);
  	double garea = (gv[2][0]-gv[0][0])*(gv[2][1]-gv[0][1]);
	for(int i=0; i<gj.size(); i++){
  		gj[i]=garea/larea;
	}
	return gj;
}

std::vector<double> TmpEle::Global_to_Local_jacobian(const std::vector<MyPoint>& gp, const std::vector<MyPoint> & gv) const
{
        std::vector<double> lj(gp.size());
  	double larea = (Pnt[2][0]-Pnt[0][0])*(Pnt[2][1]-Pnt[0][1]);
  	double garea = (gv[2][0]-gv[0][0])*(gv[2][1]-gv[0][1]);
        for(int i=0; i<lj.size(); i++){
	   lj[i]=larea/garea;
	}
	return lj;	
}

//////////////////////////////////////////////////////////////

void TmpEdge::buildTE(int i)
{
  Pnt.resize(2);
  Pnt[0]=-1.0;
  Pnt[1]=1.0; 
  switch(i){
     case 1:
        GaussPnt.resize(2);
	GaussWeight.resize(2);
	GaussPnt[0]=0.500000000;
	GaussPnt[1]=-0.500000000;
	GaussWeight[0]=1.000000000;
	GaussWeight[1]=1.000000000;
        break;

    case  2:
    GaussPnt.resize(2);
    GaussWeight.resize(2);
    GaussPnt[0] = -0.5773502691896257;
    GaussPnt[1] = 0.5773502691896257;
    GaussWeight[0] = 1;
    GaussWeight[1] = 1;
    break;
    case  3:
    GaussPnt.resize(3);
    GaussWeight.resize(3);
    GaussPnt[0] = 0;
    GaussPnt[1] = -0.7745966692414834;
    GaussPnt[2] = 0.7745966692414834;
    GaussWeight[0] = 0.8888888888888888;
    GaussWeight[1] = 0.5555555555555556;
    GaussWeight[2] = 0.5555555555555556;
    break;
    case  4:
    GaussPnt.resize(4);
    GaussWeight.resize(4);
    GaussPnt[0] = -0.3399810435848563;
    GaussPnt[1] = 0.3399810435848563;
    GaussPnt[2] = -0.8611363115940526;
    GaussPnt[3] = 0.8611363115940526;
    GaussWeight[0] = 0.6521451548625461;
    GaussWeight[1] = 0.6521451548625461;
    GaussWeight[2] = 0.3478548451374538;
    GaussWeight[3] = 0.3478548451374538;
    break;
    case  5:
    GaussPnt.resize(5);
    GaussWeight.resize(5);
    GaussPnt[0] = 0;
    GaussPnt[1] = -0.5384693101056831;
    GaussPnt[2] = 0.5384693101056831;
    GaussPnt[3] = -0.906179845938664;
    GaussPnt[4] = 0.906179845938664;
    GaussWeight[0] = 0.5688888888888889;
    GaussWeight[1] = 0.4786286704993665;
    GaussWeight[2] = 0.4786286704993665;
    GaussWeight[3] = 0.2369268850561891;
    GaussWeight[4] = 0.2369268850561891;
    break;
    case  6:
    GaussPnt.resize(6);
    GaussWeight.resize(6);
    GaussPnt[0] = 0.6612093864662645;
    GaussPnt[1] = -0.6612093864662645;
    GaussPnt[2] = -0.2386191860831969;
    GaussPnt[3] = 0.2386191860831969;
    GaussPnt[4] = -0.9324695142031521;
    GaussPnt[5] = 0.9324695142031521;
    GaussWeight[0] = 0.3607615730481386;
    GaussWeight[1] = 0.3607615730481386;
    GaussWeight[2] = 0.467913934572691;
    GaussWeight[3] = 0.467913934572691;
    GaussWeight[4] = 0.1713244923791704;
    GaussWeight[5] = 0.1713244923791704;
    break;
	
    case  10:
    GaussPnt.resize(10);
    GaussWeight.resize(10);
    GaussPnt[0] = -0.1488743389816312;
    GaussPnt[1] = 0.1488743389816312;
    GaussPnt[2] = -0.4333953941292472;
    GaussPnt[3] = 0.4333953941292472;
    GaussPnt[4] = -0.6794095682990244;
    GaussPnt[5] = 0.6794095682990244;
    GaussPnt[6] = -0.8650633666889845;
    GaussPnt[7] = 0.8650633666889845;
    GaussPnt[8] = -0.9739065285171717;
    GaussPnt[9] = 0.9739065285171717;
    GaussWeight[0] = 0.2955242247147529;
    GaussWeight[1] = 0.2955242247147529;
    GaussWeight[2] = 0.2692667193099963;
    GaussWeight[3] = 0.2692667193099963;
    GaussWeight[4] = 0.219086362515982;
    GaussWeight[5] = 0.219086362515982;
    GaussWeight[6] = 0.1494513491505806;
    GaussWeight[7] = 0.1494513491505806;
    GaussWeight[8] = 0.0666713443086881;
    GaussWeight[9] = 0.0666713443086881;
    break;	
	
	
     case 12:
        GaussPnt.resize(12);
	GaussWeight.resize(12);
	GaussPnt[0]=0.033765243;
	GaussPnt[1]=0.169395307;
	GaussPnt[2]=0.380690407;
	GaussPnt[3]=0.619309593;
	GaussPnt[4]=0.830604693;
	GaussPnt[5]=0.966234757;
	GaussPnt[6]=-0.033765243;
	GaussPnt[7]=-0.169395307;
	GaussPnt[8]=-0.380690407;
	GaussPnt[9]=-0.619309593;
	GaussPnt[10]=-0.830604693;
	GaussPnt[11]=-0.966234757;
	GaussWeight[0]=0.042831123*2.0;
	GaussWeight[1]=0.090190393*2.0;
	GaussWeight[2]=0.116978484*2.0;
	GaussWeight[3]=0.116978484*2.0;
	GaussWeight[4]=0.090190393*2.0;
	GaussWeight[5]=0.042831123*2.0;
	GaussWeight[6]=0.042831123*2.0;
	GaussWeight[7]=0.090190393*2.0; 
	GaussWeight[8]=0.116978484*2.0;
	GaussWeight[9]=0.116978484*2.0;
	GaussWeight[10]=0.090190393*2.0;
	GaussWeight[11]=0.042831123*2.0;
        break;
     case 14:
        GaussPnt.resize(14);
	GaussWeight.resize(14);
	GaussPnt[0]=0.025446044;
	GaussPnt[1]=0.129234407;
	GaussPnt[2]=0.297077424;
	GaussPnt[3]=0.500000000;
	GaussPnt[4]=0.702922576;
	GaussPnt[5]=0.870765593;
	GaussPnt[6]=0.974553956;
	GaussPnt[7]=-0.025446044;
	GaussPnt[8]=-0.129234407;
	GaussPnt[9]=-0.297077424;
	GaussPnt[10]=-0.500000000;
	GaussPnt[11]=-0.702922576;
	GaussPnt[12]=-0.870765593;
	GaussPnt[13]=-0.974553956;
	GaussWeight[0]=0.032371242*2.0;
	GaussWeight[1]=0.069926348*2.0;
	GaussWeight[2]=0.095457513*2.0;
	GaussWeight[3]=0.104489796*2.0;
	GaussWeight[4]=0.095457513*2.0;
	GaussWeight[5]=0.069926348*2.0;
	GaussWeight[6]=0.032371242*2.0;
	GaussWeight[7]=0.032371242*2.0;  
	GaussWeight[8]=0.069926348*2.0;
	GaussWeight[9]=0.095457513*2.0; 
	GaussWeight[10]=0.104489796*2.0;
	GaussWeight[11]=0.095457513*2.0; 
	GaussWeight[12]=0.069926348*2.0;
	GaussWeight[13]=0.032371242*2.0;       
        break;
     case 16:
        GaussPnt.resize(16);
	GaussWeight.resize(16);
	GaussPnt[0]=0.019855072;
	GaussPnt[1]=0.101666761;
	GaussPnt[2]=0.237233795;
	GaussPnt[3]=0.408282679;
	GaussPnt[4]=0.591717321;
	GaussPnt[5]=0.762766205;
	GaussPnt[6]=0.898333239;
	GaussPnt[7]=0.980144928;
	GaussPnt[8]=-0.019855072;
	GaussPnt[9]=-0.101666761;
	GaussPnt[10]=-0.237233795;
	GaussPnt[11]=-0.408282679;
	GaussPnt[12]=-0.591717321;
	GaussPnt[13]=-0.762766205;
	GaussPnt[14]=-0.898333239;
	GaussPnt[15]=-0.980144928;
	GaussWeight[0]=0.025307134*2.0;
	GaussWeight[1]=0.055595259*2.0;
	GaussWeight[2]=0.078426661*2.0;
	GaussWeight[3]=0.090670946*2.0;
	GaussWeight[4]=0.090670946*2.0;
	GaussWeight[5]=0.078426661*2.0;
	GaussWeight[6]=0.055595259*2.0;
	GaussWeight[7]=0.025307134*2.0;  
	GaussWeight[8]=0.025307134*2.0;
	GaussWeight[9]=0.055595259*2.0; 
	GaussWeight[10]=0.078426661*2.0;
	GaussWeight[11]=0.090670946*2.0;  
	GaussWeight[12]=0.090670946*2.0;
	GaussWeight[13]=0.078426661*2.0;    
	GaussWeight[14]=0.055595259*2.0;
	GaussWeight[15]=0.025307134*2.0;        
        break;
     case 18:
        GaussPnt.resize(18);
	GaussWeight.resize(18);
	GaussPnt[0]=0.015919880;
	GaussPnt[1]=0.081984446;
	GaussPnt[2]=0.193314284;
	GaussPnt[3]=0.337873288;
	GaussPnt[4]=0.500000000;
	GaussPnt[5]=0.662126712;
	GaussPnt[6]=0.806685716;
	GaussPnt[7]=0.918015554;
	GaussPnt[8]=0.984080120;
	GaussPnt[9]=-0.015919880;
	GaussPnt[10]=-0.081984446;
	GaussPnt[11]=-0.193314284;
	GaussPnt[12]=-0.337873288;
	GaussPnt[13]=-0.500000000;
	GaussPnt[14]=-0.662126712;
	GaussPnt[15]=-0.806685716;
	GaussPnt[16]=-0.918015554;
	GaussPnt[17]=-0.984080120;
	GaussWeight[0]=0.020318597*2.0;
	GaussWeight[1]=0.045162040*2.0;
	GaussWeight[2]=0.065152674*2.0;
	GaussWeight[3]=0.078086769*2.0;
	GaussWeight[4]=0.082559839*2.0;
	GaussWeight[5]=0.078086769*2.0;
	GaussWeight[6]=0.065152674*2.0;
	GaussWeight[7]=0.045162040*2.0;  
	GaussWeight[8]=0.020318597*2.0;
	GaussWeight[9]=0.020318597*2.0; 
	GaussWeight[10]=0.045162040*2.0;
	GaussWeight[11]=0.065152674*2.0;  
	GaussWeight[12]=0.078086769*2.0;
	GaussWeight[13]=0.082559839*2.0;    
	GaussWeight[14]=0.078086769*2.0;
	GaussWeight[15]=0.065152674*2.0;   
	GaussWeight[16]=0.045162040*2.0;
	GaussWeight[17]=0.020318597*2.0; 
        break;
     case 20:
        GaussPnt.resize(20);
	GaussWeight.resize(20);
	GaussPnt[0]=0.013046736;
	GaussPnt[1]=0.067468317;
	GaussPnt[2]=0.160295216;
	GaussPnt[3]=0.283302303;
	GaussPnt[4]=0.425562830;
	GaussPnt[5]=0.574437170;
	GaussPnt[6]=0.716697697;
	GaussPnt[7]=0.839704784;
	GaussPnt[8]=0.932531683;
	GaussPnt[9]=0.986953264;
	GaussPnt[10]=-0.013046736;
	GaussPnt[11]=-0.067468317;
	GaussPnt[12]=-0.160295216;
	GaussPnt[13]=-0.283302303;
	GaussPnt[14]=-0.425562830;
	GaussPnt[15]=-0.574437170;
	GaussPnt[16]=-0.716697697;
	GaussPnt[17]=-0.839704784;
	GaussPnt[18]=-0.932531683;
	GaussPnt[19]=-0.986953264;
	GaussWeight[0]=0.016667836*2.0;
	GaussWeight[1]=0.037362837*2.0;
	GaussWeight[2]=0.054771591*2.0;
	GaussWeight[3]=0.067316680*2.0;
	GaussWeight[4]=0.073881056*2.0;
	GaussWeight[5]=0.073881056*2.0;
	GaussWeight[6]=0.067316680*2.0;
	GaussWeight[7]=0.054771591*2.0;  
	GaussWeight[8]=0.037362837*2.0;
	GaussWeight[9]=0.016667836*2.0; 
	GaussWeight[10]=0.016667836*2.0;
	GaussWeight[11]=0.037362837*2.0;  
	GaussWeight[12]=0.054771591*2.0;
	GaussWeight[13]=0.067316680*2.0;    
	GaussWeight[14]=0.073881056*2.0;
	GaussWeight[15]=0.073881056*2.0;   
	GaussWeight[16]=0.067316680*2.0;
	GaussWeight[17]=0.054771591*2.0; 
	GaussWeight[18]=0.037362837*2.0;
	GaussWeight[19]=0.016667836*2.0;
        break;

        case  30:
    GaussPnt.resize(30);
    GaussWeight.resize(30);
    GaussPnt[0] = -0.0514718425553177;
    GaussPnt[1] = 0.0514718425553177;
    GaussPnt[2] = -0.1538699136085835;
    GaussPnt[3] = 0.1538699136085835;
    GaussPnt[4] = -0.2546369261678899;
    GaussPnt[5] = 0.2546369261678899;
    GaussPnt[6] = -0.3527047255308781;
    GaussPnt[7] = 0.3527047255308781;
    GaussPnt[8] = -0.4470337695380892;
    GaussPnt[9] = 0.4470337695380892;
    GaussPnt[10] = -0.5366241481420199;
    GaussPnt[11] = 0.5366241481420199;
    GaussPnt[12] = -0.6205261829892429;
    GaussPnt[13] = 0.6205261829892429;
    GaussPnt[14] = -0.6978504947933158;
    GaussPnt[15] = 0.6978504947933158;
    GaussPnt[16] = -0.7677774321048262;
    GaussPnt[17] = 0.7677774321048262;
    GaussPnt[18] = -0.8295657623827684;
    GaussPnt[19] = 0.8295657623827684;
    GaussPnt[20] = -0.8825605357920527;
    GaussPnt[21] = 0.8825605357920527;
    GaussPnt[22] = -0.9262000474292743;
    GaussPnt[23] = 0.9262000474292743;
    GaussPnt[24] = -0.9600218649683075;
    GaussPnt[25] = 0.9600218649683075;
    GaussPnt[26] = -0.9836681232797472;
    GaussPnt[27] = 0.9836681232797472;
    GaussPnt[28] = -0.9968934840746495;
    GaussPnt[29] = 0.9968934840746495;
    GaussWeight[0] = 0.1028526528935588;
    GaussWeight[1] = 0.1028526528935588;
    GaussWeight[2] = 0.1017623897484055;
    GaussWeight[3] = 0.1017623897484055;
    GaussWeight[4] = 0.0995934205867953;
    GaussWeight[5] = 0.0995934205867953;
    GaussWeight[6] = 0.0963687371746443;
    GaussWeight[7] = 0.0963687371746443;
    GaussWeight[8] = 0.0921225222377861;
    GaussWeight[9] = 0.0921225222377861;
    GaussWeight[10] = 0.086899787201083;
    GaussWeight[11] = 0.086899787201083;
    GaussWeight[12] = 0.0807558952294202;
    GaussWeight[13] = 0.0807558952294202;
    GaussWeight[14] = 0.0737559747377052;
    GaussWeight[15] = 0.0737559747377052;
    GaussWeight[16] = 0.0659742298821805;
    GaussWeight[17] = 0.0659742298821805;
    GaussWeight[18] = 0.0574931562176191;
    GaussWeight[19] = 0.0574931562176191;
    GaussWeight[20] = 0.0484026728305941;
    GaussWeight[21] = 0.0484026728305941;
    GaussWeight[22] = 0.0387991925696271;
    GaussWeight[23] = 0.0387991925696271;
    GaussWeight[24] = 0.0287847078833234;
    GaussWeight[25] = 0.0287847078833234;
    GaussWeight[26] = 0.018466468311091;
    GaussWeight[27] = 0.018466468311091;
    GaussWeight[28] = 0.0079681924961666;
    GaussWeight[29] = 0.0079681924961666;
    break;
    
     default:
        std::cout<<"There are no Gauss quadrature points with accuracy of "<<i<<"-th order"<<std::endl;
  }
  
}
 
double TmpEdge::getPnt(int i)
{
  return Pnt[i];
} 
 
double TmpEdge::getVolume()
{
  return Pnt[1]-Pnt[0];
} 

std::vector<double> TmpEdge::getGaussPnt()
{
   return GaussPnt;
}

std::vector<double> TmpEdge::getGaussWeight()
{ 
   return GaussWeight;
}
 
double TmpEdge::Local_to_Global(const double & lp, const std::vector<double> & gv) const
{
        double gp;
	double lambda[2];
	lambda[0] = (1.0 - lp)/2.0;
	lambda[1] = (lp + 1.0)/2.0;
	gp = lambda[0]*gv[0] + lambda[1]*gv[1];
	return gp;
}

double TmpEdge::Global_to_Local(const double & gp, const std::vector<double> & gv) const
{
        double lp;
	double lambda[2];
	lambda[0] = (gv[1] - gp)/(gv[1] - gv[0]);
	lambda[1] = (gp - gv[0])/(gv[1] - gv[0]);
	lp = -lambda[0] + lambda[1];
	return lp;
}

double TmpEdge::Local_to_Global_jacobian(const double & lp, const std::vector<double> & gv) const
{
        return (gv[1]-gv[0])/2.0;
}

double TmpEdge::Global_to_Local_jacobian(const double& gp, const std::vector<double> & gv) const
{
        return 2.0/(gv[1]-gv[0]);
}

std::vector<double> TmpEdge::Local_to_Global(const std::vector<double>& lp, const std::vector<double> & gv) const
{
        std::vector<double> gp(lp.size());
	double lambda[2];
	for(int i=0; i<gp.size(); i++){
	    lambda[0] = (1.0 - lp[i])/2.0;
	    lambda[1] = (lp[i] + 1.0)/2.0;
            gp[i] = lambda[0]*gv[0] + lambda[1]*gv[1]; 
	}
	return gp;
}

std::vector<double> TmpEdge::Global_to_Local(const std::vector<double>& gp, const std::vector<double> & gv) const
{
        std::vector<double> lp(gp.size());
	double lambda[2];
	for(int i=0; i<gp.size(); i++){
	    lambda[0] = (gv[1] - gp[i])/(gv[1] - gv[0]);
	    lambda[1] = (gp[i] - gv[0])/(gv[1] - gv[0]);
	    lp[i] = -lambda[0] + lambda[1]; 
	}
	return lp;
}

std::vector<double> TmpEdge::Local_to_Global_jacobian(const std::vector<double>& lp, const std::vector<double> & gv) const
{
        std::vector<double> gj(lp.size());
	for(int i=0; i<gj.size(); i++){
	   gj[i]=(gv[1]-gv[0])/2.0;
	}
	return gj;
}

std::vector<double> TmpEdge::Global_to_Local_jacobian(const std::vector<double>& gp, const std::vector<double> & gv) const
{
        std::vector<double> lj(gp.size());
        for(int i=0; i<lj.size(); i++){
	   lj[i]=2.0/(gv[1]-gv[0]);
	}
	return lj;	
}

///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////

void Mesh::readData(const std::string& f)
{
  int i;
  std::ifstream is(f.c_str());
  is >> i;
  Pnt.resize(i);
  for(int j=0; j<i; j++){
    is>>Pnt[j][0]>>Pnt[j][1];
// 
/*
    Pnt[j][0]=Pnt[j][0]*40;
    Pnt[j][1]=Pnt[j][1]*40;
// */
// 
/*
    Pnt[j][0]=Pnt[j][0]*8.0*PI-4.0*PI;
    Pnt[j][1]=Pnt[j][1]*8.0*PI-4.0*PI;
// */
// /*
    Pnt[j][0]=Pnt[j][0]*_L_coef_()+_L1_coef_();
    Pnt[j][1]=Pnt[j][1]*_L_coef_()+_L1_coef_();
// */
  }

  int en;
  is>>en;
  Edg.resize(en);
  for(int j=0; j<en; j++){
    Edg[j].resize(3);
    is>>Edg[j][0]>>Edg[j][1]>>Edg[j][2];
  }

  int n;
  is>>n;
  Ele.resize(n);
  for(int j=0; j<n; j++){
    Ele[j].resize(9);
    is>>Ele[j][0]>>Ele[j][1]>>Ele[j][2]>>Ele[j][3]>>Ele[j][4]>>Ele[j][5]>>Ele[j][6]>>Ele[j][7]>>Ele[j][8];
  }
  is.close();   
}

int Mesh::getEdgVtx(int i, int j)
{
   return Edg[i][j];
}

int Mesh::getEleVtx(int i, int j)
{
   if(j>3){
      std::cout<<"Error of index j in element\t"<<j<<std::endl;
      exit(0);
   }
   return Ele[i][j];
}

int Mesh::getEleEdg(int i, int j)
{
   if(j>3){
      std::cout<<"Error of index j in edge\t"<<j<<std::endl;
      exit(0);
   }
   return Ele[i][j+4];
}

std::vector<int> Mesh::getEdgVtx(int i)
{
   return Edg[i];
}

std::vector<int> Mesh::getEleVtx(int i)
{
   std::vector<int> EleVtx(4);
   for(int j=0; j<4; j++){
      EleVtx[j]=Ele[i][j];
   }
   return EleVtx;
}

std::vector<int> Mesh::getEleEdg(int i)
{
   std::vector<int> EleEdg(4);
   for(int j=0; j<4; j++){
      EleEdg[j]=Ele[i][j+4];
   }
   return EleEdg;
}

MyPoint Mesh::getPnt(int i)
{
   return Pnt[i];
}

std::vector<MyPoint> Mesh::getPnt(std::vector<int> i)
{
   std::vector<MyPoint> pt;
   pt.resize(i.size());
   for(int j=0; j<pt.size(); j++){
      pt[j]=Pnt[i[j]];
   }
   return pt;
}

int Mesh::n_edge()
{
   return Edg.size();
}

int Mesh::n_element()
{
   return Ele.size();
}

int Mesh::n_point()
{
   return Pnt.size();
}

int Mesh::getEleBndInfo(int i)
{
   return Ele[i][8];
}

int Mesh::getEdgBndInfo(int i)
{
   return Edg[i][2];
}
////////////////////////////////////////////////////////////////////////////////////

BH2D::BH2D(const std::string& file)
{
  mesh.readData(file);
}

void BH2D::init()
{
   tmpEle.buildTE(5);
   tmpEdge.buildTE(5);
   build_Edge_Patch();
   build_Element_Patch();
//   buildSparsityPattern();

}

// 
/*
void BH2D::basisValue(MyPoint &p, std::vector<MyPoint> &v, std::vector<double> & val)
{
  val.resize(9);
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
  val[0] = 1.0;
  val[1] = p[0]-Pc[0];
  val[2] = p[1]-Pc[1];
  val[3] = pow(p[0]-Pc[0], 2);
  val[4] = (p[0]-Pc[0])*(p[1]-Pc[1]);
  val[5] = pow(p[1]-Pc[1], 2);
  val[6] = pow(p[0]-Pc[0], 2)*pow(p[1]-Pc[1], 1);
  val[7] = pow(p[0]-Pc[0], 1)*pow(p[1]-Pc[1], 2);
  val[8] = pow(p[0]-Pc[0], 2)*pow(p[1]-Pc[1], 2);
}

void BH2D::basisGrad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
   val.resize(9);
   for(int i=0; i<val.size(); i++){
       val[i].resize(2);
   }

  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

  val[0][0] = 0.0;
  val[0][1] = 0.0;

  val[1][0] = 1.0;
  val[1][1] = 0.0;

  val[2][0] = 0.0;
  val[2][1] = 1.0;

  val[3][0] = 2.0*(p[0]-Pc[0]);
  val[3][1] = 0.0;

  val[4][0] = p[1]-Pc[1];
  val[4][1] = p[0]-Pc[0];

  val[5][0] = 0.0;
  val[5][1] = 2.0*(p[1]-Pc[1]);

  val[6][0] = 2.0*(p[0]-Pc[0])*(p[1]-Pc[1]);
  val[6][1] = pow(p[0]-Pc[0], 2);

  val[7][0] = pow(p[1]-Pc[1], 2);
  val[7][1] = 2.0*(p[0]-Pc[0])*(p[1]-Pc[1]);

  val[8][0] = 2.0*pow(p[0]-Pc[0], 1)*pow(p[1]-Pc[1], 2);
  val[8][1] = 2.0*pow(p[0]-Pc[0], 2)*pow(p[1]-Pc[1], 1);
}

void BH2D::basisValue(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
//  std::cout<<" Pc[0]= "<<Pc[0]<<" Pc[1]= "<<Pc[1]<<std::endl;

   int n=p.size();
   val.resize(9);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<n; i++){
      val[0][i] = 1.0;
      val[1][i] = p[i][0]-Pc[0];
      val[2][i] = p[i][1]-Pc[1];
      val[3][i] = pow(p[i][0]-Pc[0], 2);
      val[4][i] = (p[i][0]-Pc[0])*(p[i][1]-Pc[1]);
      val[5][i] = pow(p[i][1]-Pc[1], 2);
      val[6][i] = pow(p[i][0]-Pc[0], 2)*pow(p[i][1]-Pc[1], 1);
      val[7][i] = pow(p[i][0]-Pc[0], 1)*pow(p[i][1]-Pc[1], 2);
      val[8][i] = pow(p[i][0]-Pc[0], 2)*pow(p[i][1]-Pc[1], 2);
      
//      if(v[0][0]==p[0][0]&&v[0][1]==v[0][1])
//     std::cout<<"val[i]= "<<val[8][i]<<std::endl;
   }
     
}

void BH2D::basisGrad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

   int n=p.size();
   val.resize(9);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(2);
      }
   }
   for(int i=0; i<n; i++){

  	val[0][i][0] = 0.0;
  	val[0][i][1] = 0.0;

  	val[1][i][0] = 1.0;
 	val[1][i][1] = 0.0;

  	val[2][i][0] = 0.0;
  	val[2][i][1] = 1.0;

  	val[3][i][0] = 2.0*(p[i][0]-Pc[0]);
  	val[3][i][1] = 0.0;

  	val[4][i][0] = p[i][1]-Pc[1];
  	val[4][i][1] = p[i][0]-Pc[0];

  	val[5][i][0] = 0.0;
  	val[5][i][1] = 2.0*(p[i][1]-Pc[1]);

  	val[6][i][0] = 2.0*(p[i][0]-Pc[0])*(p[i][1]-Pc[1]);
  	val[6][i][1] = pow(p[i][0]-Pc[0], 2);

  	val[7][i][0] = pow(p[i][1]-Pc[1], 2);
  	val[7][i][1] = 2.0*(p[i][0]-Pc[0])*(p[i][1]-Pc[1]);

  	val[8][i][0] = 2.0*pow(p[i][0]-Pc[0], 1)*pow(p[i][1]-Pc[1], 2);
  	val[8][i][1] = 2.0*pow(p[i][0]-Pc[0], 2)*pow(p[i][1]-Pc[1], 1);
   }
}

// */

// 
/*
void BH2D::basisValue(MyPoint &p, std::vector<MyPoint> &v, std::vector<double> & val)
{
  val.resize(6);
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
  val[0] = 1.0;
  val[1] = p[0]-Pc[0];
  val[2] = p[1]-Pc[1];
  val[3] = power(p[0]-Pc[0], 2);
  val[4] = (p[0]-Pc[0])*(p[1]-Pc[1]);
  val[5] = power(p[1]-Pc[1], 2);
}

void BH2D::basisGrad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
   val.resize(6);
   for(int i=0; i<val.size(); i++){
       val[i].resize(2);
   }

  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

  val[0][0] = 0.0;
  val[0][1] = 0.0;

  val[1][0] = 1.0;
  val[1][1] = 0.0;

  val[2][0] = 0.0;
  val[2][1] = 1.0;

  val[3][0] = 2.0*(p[0]-Pc[0]);
  val[3][1] = 0.0;

  val[4][0] = p[1]-Pc[1];
  val[4][1] = p[0]-Pc[0];

  val[5][0] = 0.0;
  val[5][1] = 2.0*(p[1]-Pc[1]);
}

void BH2D::basisValue(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
//  std::cout<<" Pc[0]= "<<Pc[0]<<" Pc[1]= "<<Pc[1]<<std::endl;

   int n=p.size();
   val.resize(6);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<n; i++){
      val[0][i] = 1.0;
      val[1][i] = p[i][0]-Pc[0];
      val[2][i] = p[i][1]-Pc[1];
      val[3][i] = power(p[i][0]-Pc[0], 2);
      val[4][i] = (p[i][0]-Pc[0])*(p[i][1]-Pc[1]);
      val[5][i] = power(p[i][1]-Pc[1], 2);
      
//      if(v[0][0]==p[0][0]&&v[0][1]==v[0][1])
//     std::cout<<"val[i]= "<<val[8][i]<<std::endl;
   }
     
}

void BH2D::basisGrad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

   int n=p.size();
   val.resize(6);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(2);
      }
   }
   for(int i=0; i<n; i++){

  	val[0][i][0] = 0.0;
  	val[0][i][1] = 0.0;

  	val[1][i][0] = 1.0;
 	val[1][i][1] = 0.0;

  	val[2][i][0] = 0.0;
  	val[2][i][1] = 1.0;

  	val[3][i][0] = 2.0*(p[i][0]-Pc[0]);
  	val[3][i][1] = 0.0;

  	val[4][i][0] = p[i][1]-Pc[1];
  	val[4][i][1] = p[i][0]-Pc[0];

  	val[5][i][0] = 0.0;
  	val[5][i][1] = 2.0*(p[i][1]-Pc[1]);

   }
}

// */


// 
/*
void BH2D::basisValue(MyPoint &p, std::vector<MyPoint> &v, std::vector<double> & val)
{
  val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

  for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++)
        {
           val[i*Msize1+j]= power(p[0]-Pc[0], i)*power(p[1]-Pc[1], j);
        }
}

void BH2D::basisGrad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
	val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
       val[i].resize(2);
   }

  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
  
    for(int i=0; i<Msize1; i++)
        for(int j=0; j<Msize1; j++)
		{
            val[i*Msize1+j][0] = i*power(p[0]-Pc[0], i-1)*power(p[1]-Pc[1], j);
            val[i*Msize1+j][1] = j*power(p[0]-Pc[0], i)*power(p[1]-Pc[1], j-1);
		}

}

void BH2D::basis2Grad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
       val[i].resize(3);
   }


  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

  for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++){
        val[i*Msize1+j][0] =i*(i-1)*power(p[0]-Pc[0], i-2)*power(p[1]-Pc[1], j);
        val[i*Msize1+j][1] =i*j*power(p[0]-Pc[0], i-1)*power(p[1]-Pc[1], j-1);
        val[i*Msize1+j][2] =j*(j-1)*power(p[0]-Pc[0], i)*power(p[1]-Pc[1], j-2);
	 }
}

void BH2D::basisValue(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
   int n=p.size();
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }

   for(int l=0; l<n; l++){
   for(int i=0; i<Msize1; i++){
     for(int j=0; j<Msize1; j++)
        {
           val[i*Msize1+j][l]= power(p[l][0]-Pc[0], i)*power(p[l][1]-Pc[1], j);
        }
   }
   }
}

void BH2D::basisGrad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

   int n=p.size();
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(2);
      }
   }
   for(int l=0; l<n; l++){
   for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++)
	 {
        val[i*Msize1+j][l][0] = i*power(p[l][0]-Pc[0], i-1)*power(p[l][1]-Pc[1], j);
        val[i*Msize1+j][l][1] = j*power(p[l][0]-Pc[0], i)*power(p[l][1]-Pc[1], j-1);
	 }
//	std::cout<<val[7][l][0]<<"   "<<2.0*(p[l][0]-Pc[0])*(p[l][1]-Pc[1])<<std::endl;
   }

}

void BH2D::basis2Grad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

   int n=p.size();
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(3);
      }
   }

   for(int l=0; l<n; l++){
   for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++)
	 {
        val[i*Msize1+j][l][0] =i*(i-1)*power(p[l][0]-Pc[0], i-2)*power(p[l][1]-Pc[1], j);
        val[i*Msize1+j][l][1] =i*j*power(p[l][0]-Pc[0], i-1)*power(p[l][1]-Pc[1], j-1);
        val[i*Msize1+j][l][2] =j*(j-1)*power(p[l][0]-Pc[0], i)*power(p[l][1]-Pc[1], j-2);
	 }
   }
}
// */

//  Polynomials 1, x,y,xy,.... 
/*

void BH2D::basisValue(MyPoint &p, std::vector<MyPoint> &v, std::vector<double> & val)
{
  val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

  for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++)
        {
           val[i*Msize1+j]= power(p[0]-Pc[0], i)*power(p[1]-Pc[1], j);
        }

}

void BH2D::basisGrad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
	val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
       val[i].resize(2);
   }

  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
  
    for(int i=0; i<Msize1; i++)
        for(int j=0; j<Msize1; j++)
		{
            val[i*Msize1+j][0] = i*power(p[0]-Pc[0], i-1)*power(p[1]-Pc[1], j);
            val[i*Msize1+j][1] = j*power(p[0]-Pc[0], i)*power(p[1]-Pc[1], j-1);
		}

}

void BH2D::basis2Grad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
       val[i].resize(3);
   }


  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

  for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++){
        val[i*Msize1+j][0] =i*(i-1)*power(p[0]-Pc[0], i-2)*power(p[1]-Pc[1], j);
        val[i*Msize1+j][1] =i*j*power(p[0]-Pc[0], i-1)*power(p[1]-Pc[1], j-1);
        val[i*Msize1+j][2] =j*(j-1)*power(p[0]-Pc[0], i)*power(p[1]-Pc[1], j-2);
	 }
              
   
}

void BH2D::basisValue(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
   int n=p.size();
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }

   for(int l=0; l<n; l++){
   for(int i=0; i<Msize1; i++){
     for(int j=0; j<Msize1; j++)
        {
           val[i*Msize1+j][l]= power(p[l][0]-Pc[0], i)*power(p[l][1]-Pc[1], j);
        }
   }
   }

}

void BH2D::basisGrad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

   int n=p.size();
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(2);
      }
   }
   for(int l=0; l<n; l++){
   for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++)
	 {
        val[i*Msize1+j][l][0] = i*power(p[l][0]-Pc[0], i-1)*power(p[l][1]-Pc[1], j);
        val[i*Msize1+j][l][1] = j*power(p[l][0]-Pc[0], i)*power(p[l][1]-Pc[1], j-1);
	 }
//	std::cout<<val[7][l][0]<<"   "<<2.0*(p[l][0]-Pc[0])*(p[l][1]-Pc[1])<<std::endl;
   }

}

void BH2D::basis2Grad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

   int n=p.size();
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(3);
      }
   }

   for(int l=0; l<n; l++){
   for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++)
	 {
        val[i*Msize1+j][l][0] =i*(i-1)*power(p[l][0]-Pc[0], i-2)*power(p[l][1]-Pc[1], j);
        val[i*Msize1+j][l][1] =i*j*power(p[l][0]-Pc[0], i-1)*power(p[l][1]-Pc[1], j-1);
        val[i*Msize1+j][l][2] =j*(j-1)*power(p[l][0]-Pc[0], i)*power(p[l][1]-Pc[1], j-2);
	 }
   }
}

// End of the polynomails 1,x,y,xy, ... */



// Legendre polynomial in 2D 
/*

void BH2D::basisValue(MyPoint &p, std::vector<MyPoint> &v, std::vector<double> & val)
{
    val.resize(ele_dof);
//   int Msize1 = sqrt(ele_dof);
    MyPoint Pc;
    Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
    Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
	double hx=v[2][0]-v[0][0];
	double hy=v[2][1]-v[0][1];
	
	


    std::vector<double> stdval;
    stdval.resize(16);

    stdval[0] = 1.0/2.0;
    stdval[1] = sqrt(3.0)/hx*(p[0]-Pc[0]);
    stdval[2] = sqrt(3.0)/hy*(p[1]-Pc[1]);
    stdval[3] = 6.0/hx/hy*(p[0]-Pc[0])*(p[1]-Pc[1]);
    stdval[4] = sqrt(5.0)/hx/hx*(3.0*power(p[0]-Pc[0], 2)-hx*hx/4.0);
    stdval[5] = sqrt(5.0)/hy/hy*(3.0*power(p[1]-Pc[1], 2)-hy*hy/4.0);
    stdval[6] = sqrt(15.0)/4.0*(3.0*power(p[0]-Pc[0], 2)-1.0)*(p[1]-Pc[1]);
    stdval[7] = sqrt(15.0)/4.0*(3.0*power(p[1]-Pc[1], 2)-1.0)*(p[0]-Pc[0]);
    stdval[8] = sqrt(7.0)/4.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]));
    stdval[9] = sqrt(7.0)/4.0*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1])); 
    stdval[10] = 5.0/8.0*(3.0*power(p[0]-Pc[0], 2)-1.0)*(3.0*power(p[1]-Pc[1], 2)-1.0);	
    stdval[11] = sqrt(21.0)/4.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]))*(p[1]-Pc[1]);
    stdval[12] = sqrt(21.0)/4.0*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]))*(p[0]-Pc[0]);
    stdval[13] = sqrt(35.0)/8.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]))*(3.0*power(p[1]-Pc[1], 2)-1.0);
    stdval[14] = sqrt(35.0)/8.0*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]))*(3.0*power(p[0]-Pc[0], 2)-1.0);
    stdval[15] = 7.0/8.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]))*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]));

    for(int j=0; j<ele_dof; j++)
        val[j]=stdval[j];

}

void BH2D::basisGrad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
	val.resize(ele_dof);
//    int Msize1 = sqrt(ele_dof);
    for(int i=0; i<val.size(); i++){
        val[i].resize(2);
    }

    MyPoint Pc;
    Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
    Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
	double hx=v[2][0]-v[0][0];
	double hy=v[2][1]-v[0][1];

    std::vector<std::vector<double> > stdval;
    stdval.resize(16);
    for(int i=0; i<stdval.size(); i++){
        stdval[i].resize(2);
    }

    stdval[0][0] = 0.0;
    stdval[0][1] = 0.0;

    stdval[1][0] = sqrt(3.0)/hx;
    stdval[1][1] = 0.0;

    stdval[2][0] = 0.0;
    stdval[2][1] = sqrt(3.0)/hy;

    stdval[3][0] = 6.0/hx/hy*(p[1]-Pc[1]);
    stdval[3][1] = 6.0/hx/hy*(p[0]-Pc[0]);

    stdval[4][0] = sqrt(5.0)/hx/hx*(6.0*(p[0]-Pc[0]));
    stdval[4][1] = 0.0;

    stdval[5][0] = 0.0;
    stdval[5][1] = sqrt(5.0)/hy/hy*(6.0*(p[1]-Pc[1]));

    stdval[6][0] = sqrt(15.0)/4.0*(6.0*(p[0]-Pc[0]))*(p[1]-Pc[1]);
    stdval[6][1] = sqrt(15.0)/4.0*(3.0*power(p[0]-Pc[0], 2)-1.0);

    stdval[7][0] = sqrt(15.0)/4.0*(3.0*power(p[1]-Pc[1], 2)-1.0);
    stdval[7][1] = sqrt(15.0)/4.0*(6.0*(p[1]-Pc[1]))*(p[0]-Pc[0]);

    stdval[8][0] = sqrt(7.0)/4.0*(15.0*power(p[0]-Pc[0], 2)-3.0);
    stdval[8][1] = 0.0;

    stdval[9][0] = 0.0;
    stdval[9][1] = sqrt(7.0)/4.0*(15.0*power(p[1]-Pc[1], 2)-3.0);
	
    stdval[10][0] = 5.0/8.0*(6.0*(p[0]-Pc[0]))*(3.0*power(p[1]-Pc[1], 2)-1.0);
    stdval[10][1] = 5.0/8.0*(6.0*(p[1]-Pc[1]))*(3.0*power(p[0]-Pc[0], 2)-1.0);	

    stdval[11][0] = sqrt(21.0)/4.0*(15.0*power(p[0]-Pc[0], 2)-3.0)*(p[1]-Pc[1]);
    stdval[11][1] = sqrt(21.0)/4.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]));

    stdval[12][0] = sqrt(21.0)/4.0*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]));
    stdval[12][1] = sqrt(21.0)/4.0*(15.0*power(p[1]-Pc[1], 2)-3.0)*(p[0]-Pc[0]);

    stdval[13][0] = sqrt(35.0)/8.0*(15.0*power(p[0]-Pc[0], 2)-3.0)*(3.0*power(p[1]-Pc[1], 2)-1.0);
    stdval[13][1] = sqrt(35.0)/8.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]))*(6.0*(p[1]-Pc[1]));

    stdval[14][0] = sqrt(35.0)/8.0*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]))*(6.0*(p[0]-Pc[0]));
    stdval[14][1] = sqrt(35.0)/8.0*(15.0*power(p[1]-Pc[1], 2)-3.0)*(3.0*power(p[0]-Pc[0], 2)-1.0);

    stdval[15][0] = 7.0/8.0*(15.0*power(p[0]-Pc[0], 2)-3.0)*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]));
    stdval[15][1] = 7.0/8.0*(15.0*power(p[1]-Pc[1], 2)-3.0)*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]));
  
    for(int j=0; j<ele_dof; j++)
        for(int k=0; k<2; k++)
            val[j][k] = stdval[j][k];
 
}

void BH2D::basisValue(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
    MyPoint Pc;
    Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
    Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
	double hx=v[2][0]-v[0][0];
	double hy=v[2][1]-v[0][1];
	
    int n=p.size();
    val.resize(ele_dof);
//    int Msize1 = sqrt(ele_dof);
    for(int i=0; i<val.size(); i++){
       val[i].resize(n);
    }

    std::vector<double> stdval;
    stdval.resize(16);

    for(int i=0; i<n; i++){

    stdval[0] = 1.0/2.0;
    stdval[1] = sqrt(3.0)/hx*(p[i][0]-Pc[0]);
    stdval[2] = sqrt(3.0)/hy*(p[i][1]-Pc[1]);
    stdval[3] = 6.0/hx/hy*(p[i][0]-Pc[0])*(p[i][1]-Pc[1]);
    stdval[4] = sqrt(5.0)/hx/hx*(3.0*power(p[i][0]-Pc[0], 2)-hx*hx/4.0);
    stdval[5] = sqrt(5.0)/hy/hy*(3.0*power(p[i][1]-Pc[1], 2)-hy*hy/4.0);
    stdval[6] = sqrt(15.0)/4.0*(3.0*power(p[i][0]-Pc[0], 2)-1.0)*(p[i][1]-Pc[1]);
    stdval[7] = sqrt(15.0)/4.0*(3.0*power(p[i][1]-Pc[1], 2)-1.0)*(p[i][0]-Pc[0]);
    stdval[8] = sqrt(7.0)/4.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]));
    stdval[9] = sqrt(7.0)/4.0*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]));
    stdval[10] = 5.0/8.0*(3.0*power(p[i][0]-Pc[0], 2)-1.0)*(3.0*power(p[i][1]-Pc[1], 2)-1.0);	
    stdval[11] = sqrt(21.0)/4.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]))*(p[i][1]-Pc[1]);
    stdval[12] = sqrt(21.0)/4.0*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]))*(p[i][0]-Pc[0]);
    stdval[13] = sqrt(35.0)/8.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]))*(3.0*power(p[i][1]-Pc[1], 2)-1.0);
    stdval[14] = sqrt(35.0)/8.0*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]))*(3.0*power(p[i][0]-Pc[0], 2)-1.0);
    stdval[15] = 7.0/8.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]))*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]));

    for(int j=0; j<ele_dof; j++)
        val[j][i]=stdval[j];
        
    }

}

void BH2D::basisGrad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
  double hx=v[2][0]-v[0][0];
  double hy=v[2][1]-v[0][1];

   int n=p.size();
   val.resize(ele_dof);

   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(2);
      }
   }
   
    std::vector<std::vector<double> > stdval;
    stdval.resize(16);
    for(int i=0; i<stdval.size(); i++){
        stdval[i].resize(2);
    }
    
    for(int i=0; i<n; i++){

    stdval[0][0] = 0.0;
    stdval[0][1] = 0.0;

    stdval[1][0] = sqrt(3.0)/hx;
    stdval[1][1] = 0.0;

    stdval[2][0] = 0.0;
    stdval[2][1] = sqrt(3.0)/hy;

    stdval[3][0] = 6.0/hx/hy*(p[i][1]-Pc[1]);
    stdval[3][1] = 6.0/hx/hy*(p[i][0]-Pc[0]);

    stdval[4][0] = sqrt(5.0)/hx/hx*(6.0*(p[i][0]-Pc[0]));
    stdval[4][1] = 0.0;

    stdval[5][0] = 0.0;
    stdval[5][1] = sqrt(5.0)/hy/hy*(6.0*(p[i][1]-Pc[1]));

    stdval[6][0] = sqrt(15.0)/4.0*(6.0*(p[i][0]-Pc[0]))*(p[i][1]-Pc[1]);
    stdval[6][1] = sqrt(15.0)/4.0*(3.0*power(p[i][0]-Pc[0], 2)-1.0);

    stdval[7][0] = sqrt(15.0)/4.0*(3.0*power(p[i][1]-Pc[1], 2)-1.0);
    stdval[7][1] = sqrt(15.0)/4.0*(6.0*(p[i][1]-Pc[1]))*(p[i][0]-Pc[0]);

    stdval[8][0] = sqrt(7.0)/4.0*(15.0*power(p[i][0]-Pc[0], 2)-3.0);
    stdval[8][1] = 0.0;

    stdval[9][0] = 0.0;
    stdval[9][1] = sqrt(7.0)/4.0*(15.0*power(p[i][1]-Pc[1], 2)-3.0);
	
    stdval[10][0] = 5.0/8.0*(6.0*(p[i][0]-Pc[0]))*(3.0*power(p[i][1]-Pc[1], 2)-1.0);
    stdval[10][1] = 5.0/8.0*(6.0*(p[i][1]-Pc[1]))*(3.0*power(p[i][0]-Pc[0], 2)-1.0);	

    stdval[11][0] = sqrt(21.0)/4.0*(15.0*power(p[i][0]-Pc[0], 2)-3.0)*(p[i][1]-Pc[1]);
    stdval[11][1] = sqrt(21.0)/4.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]));

    stdval[12][0] = sqrt(21.0)/4.0*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]));
    stdval[12][1] = sqrt(21.0)/4.0*(15.0*power(p[i][1]-Pc[1], 2)-3.0)*(p[i][0]-Pc[0]);

    stdval[13][0] = sqrt(35.0)/8.0*(15.0*power(p[i][0]-Pc[0], 2)-3.0)*(3.0*power(p[i][1]-Pc[1], 2)-1.0);
    stdval[13][1] = sqrt(35.0)/8.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]))*(6.0*(p[i][1]-Pc[1]));

    stdval[14][0] = sqrt(35.0)/8.0*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]))*(6.0*(p[i][0]-Pc[0]));
    stdval[14][1] = sqrt(35.0)/8.0*(15.0*power(p[i][1]-Pc[1], 2)-3.0)*(3.0*power(p[i][0]-Pc[0], 2)-1.0);

    stdval[15][0] = 7.0/8.0*(15.0*power(p[i][0]-Pc[0], 2)-3.0)*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]));
    stdval[15][1] = 7.0/8.0*(15.0*power(p[i][1]-Pc[1], 2)-3.0)*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]));

    for(int j=0; j<ele_dof; j++)
        for(int k=0; k<2; k++)
            val[j][i][k] = stdval[j][k];
  	
  	
   }

}

// End of Legendre polynomial in 2D  */

double BH2D::power(double &f, int &n)
{
	double val;
	if(n==0)
		val=1.0;
	else if(n<0)
		val=0.0;
	else
	{
		val=1.0;
		for(int i=0; i<n; i++)
		    val=val*f;
	}
	return val;
}

double BH2D::power(double f, int n)
{
	double val;
	if(n==0)
		val=1.0;
	else if(n<0)
		val=0.0;
	else
	{
		val=1.0;
		for(int i=0; i<n; i++)
		    val=val*f;
	}
	return val;
}

//  Polynomials 1, x,y,xy,.... 
/*

void BH2D::basisValue(MyPoint &p, std::vector<MyPoint> &v, std::vector<double> & val)
{
  val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

  for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++)
        {
           val[i*Msize1+j]= power(p[0]-Pc[0], i)*power(p[1]-Pc[1], j);
        }

}

void BH2D::basisGrad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
	val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
       val[i].resize(2);
   }

  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
  
    for(int i=0; i<Msize1; i++)
        for(int j=0; j<Msize1; j++)
		{
            val[i*Msize1+j][0] = i*power(p[0]-Pc[0], i-1)*power(p[1]-Pc[1], j);
            val[i*Msize1+j][1] = j*power(p[0]-Pc[0], i)*power(p[1]-Pc[1], j-1);
		}

}

void BH2D::basis2Grad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
       val[i].resize(3);
   }


  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

  for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++){
        val[i*Msize1+j][0] =i*(i-1)*power(p[0]-Pc[0], i-2)*power(p[1]-Pc[1], j);
        val[i*Msize1+j][1] =i*j*power(p[0]-Pc[0], i-1)*power(p[1]-Pc[1], j-1);
        val[i*Msize1+j][2] =j*(j-1)*power(p[0]-Pc[0], i)*power(p[1]-Pc[1], j-2);
	 }
              
   
}

void BH2D::basisValue(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
   int n=p.size();
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }

   for(int l=0; l<n; l++){
   for(int i=0; i<Msize1; i++){
     for(int j=0; j<Msize1; j++)
        {
           val[i*Msize1+j][l]= power(p[l][0]-Pc[0], i)*power(p[l][1]-Pc[1], j);
        }
   }
   }

}

void BH2D::basisGrad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

   int n=p.size();
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(2);
      }
   }
   for(int l=0; l<n; l++){
   for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++)
	 {
        val[i*Msize1+j][l][0] = i*power(p[l][0]-Pc[0], i-1)*power(p[l][1]-Pc[1], j);
        val[i*Msize1+j][l][1] = j*power(p[l][0]-Pc[0], i)*power(p[l][1]-Pc[1], j-1);
	 }
//	std::cout<<val[7][l][0]<<"   "<<2.0*(p[l][0]-Pc[0])*(p[l][1]-Pc[1])<<std::endl;
   }

}

void BH2D::basis2Grad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

   int n=p.size();
   val.resize(ele_dof);
  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(3);
      }
   }

   for(int l=0; l<n; l++){
   for(int i=0; i<Msize1; i++)
     for(int j=0; j<Msize1; j++)
	 {
        val[i*Msize1+j][l][0] =i*(i-1)*power(p[l][0]-Pc[0], i-2)*power(p[l][1]-Pc[1], j);
        val[i*Msize1+j][l][1] =i*j*power(p[l][0]-Pc[0], i-1)*power(p[l][1]-Pc[1], j-1);
        val[i*Msize1+j][l][2] =j*(j-1)*power(p[l][0]-Pc[0], i)*power(p[l][1]-Pc[1], j-2);
	 }
   }
}

// End of the polynomails 1,x,y,xy, ... */



// Legendre polynomial in 2D 
/*

void BH2D::basisValue(MyPoint &p, std::vector<MyPoint> &v, std::vector<double> & val)
{
    val.resize(ele_dof);
//   int Msize1 = sqrt(ele_dof);
    MyPoint Pc;
    Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
    Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;


    std::vector<double> stdval;
    stdval.resize(16);

    stdval[0] = 1.0/2.0;
    stdval[1] = sqrt(3.0)/2.0*(p[0]-Pc[0]);
    stdval[2] = sqrt(3.0)/2.0*(p[1]-Pc[1]);
    stdval[3] = 3.0/2.0*(p[0]-Pc[0])*(p[1]-Pc[1]);
    stdval[4] = sqrt(5.0)/4.0*(3.0*power(p[0]-Pc[0], 2)-1.0);
    stdval[5] = sqrt(5.0)/4.0*(3.0*power(p[1]-Pc[1], 2)-1.0);
    stdval[6] = sqrt(15.0)/4.0*(3.0*power(p[0]-Pc[0], 2)-1.0)*(p[1]-Pc[1]);
    stdval[7] = sqrt(15.0)/4.0*(3.0*power(p[1]-Pc[1], 2)-1.0)*(p[0]-Pc[0]);
    stdval[8] = 5.0/8.0*(3.0*power(p[0]-Pc[0], 2)-1.0)*(3.0*power(p[1]-Pc[1], 2)-1.0);
    stdval[9] = sqrt(7.0)/4.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]));
    stdval[10] = sqrt(7.0)/4.0*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]));
    stdval[11] = sqrt(21.0)/4.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]))*(p[1]-Pc[1]);
    stdval[12] = sqrt(21.0)/4.0*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]))*(p[0]-Pc[0]);
    stdval[13] = sqrt(35.0)/8.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]))*(3.0*power(p[1]-Pc[1], 2)-1.0);
    stdval[14] = sqrt(35.0)/8.0*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]))*(3.0*power(p[0]-Pc[0], 2)-1.0);
    stdval[15] = 7.0/8.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]))*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]));

    for(int j=0; j<ele_dof; j++)
        val[j]=stdval[j];

}

void BH2D::basisGrad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
	val.resize(ele_dof);
//    int Msize1 = sqrt(ele_dof);
    for(int i=0; i<val.size(); i++){
        val[i].resize(2);
    }

    MyPoint Pc;
    Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
    Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

    std::vector<std::vector<double> > stdval;
    stdval.resize(16);
    for(int i=0; i<stdval.size(); i++){
        stdval[i].resize(2);
    }

    stdval[0][0] = 0.0;
    stdval[0][1] = 0.0;

    stdval[1][0] = sqrt(3.0)/2.0;
    stdval[1][1] = 0.0;

    stdval[2][0] = 0.0;
    stdval[2][1] = sqrt(3.0)/2.0;

    stdval[3][0] = 3.0/2.0*(p[1]-Pc[1]);
    stdval[3][1] = 3.0/2.0*(p[0]-Pc[0]);

    stdval[4][0] = sqrt(5.0)/4.0*(6.0*(p[0]-Pc[0]));
    stdval[4][1] = 0.0;

    stdval[5][0] = 0.0;
    stdval[5][1] = sqrt(5.0)/4.0*(6.0*(p[1]-Pc[1]));

    stdval[6][0] = sqrt(15.0)/4.0*(6.0*(p[0]-Pc[0]))*(p[1]-Pc[1]);
    stdval[6][1] = sqrt(15.0)/4.0*(3.0*power(p[0]-Pc[0], 2)-1.0);

    stdval[7][0] = sqrt(15.0)/4.0*(3.0*power(p[1]-Pc[1], 2)-1.0);
    stdval[7][1] = sqrt(15.0)/4.0*(6.0*(p[1]-Pc[1]))*(p[0]-Pc[0]);

    stdval[8][0] = 5.0/8.0*(6.0*(p[0]-Pc[0]))*(3.0*power(p[1]-Pc[1], 2)-1.0);
    stdval[8][1] = 5.0/8.0*(6.0*(p[1]-Pc[1]))*(3.0*power(p[0]-Pc[0], 2)-1.0);

    stdval[9][0] = sqrt(7.0)/4.0*(15.0*power(p[0]-Pc[0], 2)-3.0);
    stdval[9][1] = 0.0;

    stdval[10][0] = 0.0;
    stdval[10][1] = sqrt(7.0)/4.0*(15.0*power(p[1]-Pc[1], 2)-3.0);

    stdval[11][0] = sqrt(21.0)/4.0*(15.0*power(p[0]-Pc[0], 2)-3.0)*(p[1]-Pc[1]);
    stdval[11][1] = sqrt(21.0)/4.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]));

    stdval[12][0] = sqrt(21.0)/4.0*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]));
    stdval[12][1] = sqrt(21.0)/4.0*(15.0*power(p[1]-Pc[1], 2)-3.0)*(p[0]-Pc[0]);

    stdval[13][0] = sqrt(35.0)/8.0*(15.0*power(p[0]-Pc[0], 2)-3.0)*(3.0*power(p[1]-Pc[1], 2)-1.0);
    stdval[13][1] = sqrt(35.0)/8.0*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]))*(6.0*(p[1]-Pc[1]));

    stdval[14][0] = sqrt(35.0)/8.0*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]))*(6.0*(p[0]-Pc[0]));
    stdval[14][1] = sqrt(35.0)/8.0*(15.0*power(p[1]-Pc[1], 2)-3.0)*(3.0*power(p[0]-Pc[0], 2)-1.0);

    stdval[15][0] = 7.0/8.0*(15.0*power(p[0]-Pc[0], 2)-3.0)*(5.0*power(p[1]-Pc[1], 3)-3.0*(p[1]-Pc[1]));
    stdval[15][1] = 7.0/8.0*(15.0*power(p[1]-Pc[1], 2)-3.0)*(5.0*power(p[0]-Pc[0], 3)-3.0*(p[0]-Pc[0]));
  
    for(int j=0; j<ele_dof; j++)
        for(int k=0; k<2; k++)
            val[j][k] = stdval[j][k];
 
}

void BH2D::basisValue(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
    MyPoint Pc;
    Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
    Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
    int n=p.size();
    val.resize(ele_dof);
//    int Msize1 = sqrt(ele_dof);
    for(int i=0; i<val.size(); i++){
       val[i].resize(n);
    }

    std::vector<double> stdval;
    stdval.resize(16);

    for(int i=0; i<n; i++){

    stdval[0] = 1.0/2.0;
    stdval[1] = sqrt(3.0)/2.0*(p[i][0]-Pc[0]);
    stdval[2] = sqrt(3.0)/2.0*(p[i][1]-Pc[1]);
    stdval[3] = 3.0/2.0*(p[i][0]-Pc[0])*(p[i][1]-Pc[1]);
    stdval[4] = sqrt(5.0)/4.0*(3.0*power(p[i][0]-Pc[0], 2)-1.0);
    stdval[5] = sqrt(5.0)/4.0*(3.0*power(p[i][1]-Pc[1], 2)-1.0);
    stdval[6] = sqrt(15.0)/4.0*(3.0*power(p[i][0]-Pc[0], 2)-1.0)*(p[i][1]-Pc[1]);
    stdval[7] = sqrt(15.0)/4.0*(3.0*power(p[i][1]-Pc[1], 2)-1.0)*(p[i][0]-Pc[0]);
    stdval[8] = 5.0/8.0*(3.0*power(p[i][0]-Pc[0], 2)-1.0)*(3.0*power(p[i][1]-Pc[1], 2)-1.0);
    stdval[9] = sqrt(7.0)/4.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]));
    stdval[10] = sqrt(7.0)/4.0*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]));
    stdval[11] = sqrt(21.0)/4.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]))*(p[i][1]-Pc[1]);
    stdval[12] = sqrt(21.0)/4.0*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]))*(p[i][0]-Pc[0]);
    stdval[13] = sqrt(35.0)/8.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]))*(3.0*power(p[i][1]-Pc[1], 2)-1.0);
    stdval[14] = sqrt(35.0)/8.0*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]))*(3.0*power(p[i][0]-Pc[0], 2)-1.0);
    stdval[15] = 7.0/8.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]))*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]));

    for(int j=0; j<ele_dof; j++)
        val[j][i]=stdval[j];
        
    }

}

void BH2D::basisGrad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

   int n=p.size();
   val.resize(ele_dof);

   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(2);
      }
   }
   
    std::vector<std::vector<double> > stdval;
    stdval.resize(16);
    for(int i=0; i<stdval.size(); i++){
        stdval[i].resize(2);
    }
    
    for(int i=0; i<n; i++){

    stdval[0][0] = 0.0;
    stdval[0][1] = 0.0;

    stdval[1][0] = sqrt(3.0)/2.0;
    stdval[1][1] = 0.0;

    stdval[2][0] = 0.0;
    stdval[2][1] = sqrt(3.0)/2.0;

    stdval[3][0] = 3.0/2.0*(p[i][1]-Pc[1]);
    stdval[3][1] = 3.0/2.0*(p[i][0]-Pc[0]);

    stdval[4][0] = sqrt(5.0)/4.0*(6.0*(p[i][0]-Pc[0]));
    stdval[4][1] = 0.0;

    stdval[5][0] = 0.0;
    stdval[5][1] = sqrt(5.0)/4.0*(6.0*(p[i][1]-Pc[1]));

    stdval[6][0] = sqrt(15.0)/4.0*(6.0*(p[i][0]-Pc[0]))*(p[i][1]-Pc[1]);
    stdval[6][1] = sqrt(15.0)/4.0*(3.0*power(p[i][0]-Pc[0], 2)-1.0);

    stdval[7][0] = sqrt(15.0)/4.0*(3.0*power(p[i][1]-Pc[1], 2)-1.0);
    stdval[7][1] = sqrt(15.0)/4.0*(6.0*(p[i][1]-Pc[1]))*(p[i][0]-Pc[0]);

    stdval[8][0] = 5.0/8.0*(6.0*(p[i][0]-Pc[0]))*(3.0*power(p[i][1]-Pc[1], 2)-1.0);
    stdval[8][1] = 5.0/8.0*(6.0*(p[i][1]-Pc[1]))*(3.0*power(p[i][0]-Pc[0], 2)-1.0);

    stdval[9][0] = sqrt(7.0)/4.0*(15.0*power(p[i][0]-Pc[0], 2)-3.0);
    stdval[9][1] = 0.0;

    stdval[10][0] = 0.0;
    stdval[10][1] = sqrt(7.0)/4.0*(15.0*power(p[i][1]-Pc[1], 2)-3.0);

    stdval[11][0] = sqrt(21.0)/4.0*(15.0*power(p[i][0]-Pc[0], 2)-3.0)*(p[i][1]-Pc[1]);
    stdval[11][1] = sqrt(21.0)/4.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]));

    stdval[12][0] = sqrt(21.0)/4.0*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]));
    stdval[12][1] = sqrt(21.0)/4.0*(15.0*power(p[i][1]-Pc[1], 2)-3.0)*(p[i][0]-Pc[0]);

    stdval[13][0] = sqrt(35.0)/8.0*(15.0*power(p[i][0]-Pc[0], 2)-3.0)*(3.0*power(p[i][1]-Pc[1], 2)-1.0);
    stdval[13][1] = sqrt(35.0)/8.0*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]))*(6.0*(p[i][1]-Pc[1]));

    stdval[14][0] = sqrt(35.0)/8.0*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]))*(6.0*(p[i][0]-Pc[0]));
    stdval[14][1] = sqrt(35.0)/8.0*(15.0*power(p[i][1]-Pc[1], 2)-3.0)*(3.0*power(p[i][0]-Pc[0], 2)-1.0);

    stdval[15][0] = 7.0/8.0*(15.0*power(p[i][0]-Pc[0], 2)-3.0)*(5.0*power(p[i][1]-Pc[1], 3)-3.0*(p[i][1]-Pc[1]));
    stdval[15][1] = 7.0/8.0*(15.0*power(p[i][1]-Pc[1], 2)-3.0)*(5.0*power(p[i][0]-Pc[0], 3)-3.0*(p[i][0]-Pc[0]));

    for(int j=0; j<ele_dof; j++)
        for(int k=0; k<2; k++)
            val[j][i][k] = stdval[j][k];
  	
  	
   }

}


// End of Legendre polynomial in 2D  */

// P^k polynomial in 2D /*

void BH2D::basisValue(MyPoint &p, std::vector<MyPoint> &v, std::vector<double> & val)
{
    val.resize(ele_dof);
//   int Msize1 = sqrt(ele_dof);
    MyPoint Pc;
    Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
    Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;


    std::vector<double> stdval;
    stdval.resize(10);

    stdval[0] = 1.0;
    stdval[1] = p[0]-Pc[0];
    stdval[2] = p[1]-Pc[1];
    stdval[3] = (p[0]-Pc[0])*(p[1]-Pc[1]);
    stdval[4] = pow(p[0]-Pc[0], 2);
    stdval[5] = pow(p[1]-Pc[1], 2);
    stdval[6] = pow(p[0]-Pc[0], 2)*(p[1]-Pc[1]);
    stdval[7] = pow(p[1]-Pc[1], 2)*(p[0]-Pc[0]);
    stdval[8] = pow(p[0]-Pc[0], 3);
    stdval[9] = pow(p[1]-Pc[1], 3);

    for(int j=0; j<ele_dof; j++)
        val[j]=stdval[j];

}

void BH2D::basisGrad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
	val.resize(ele_dof);
//    int Msize1 = sqrt(ele_dof);
    for(int i=0; i<val.size(); i++){
        val[i].resize(2);
    }

    MyPoint Pc;
    Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
    Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

    std::vector<std::vector<double> > stdval;
    stdval.resize(10);
    for(int i=0; i<stdval.size(); i++){
        stdval[i].resize(2);
    }

    stdval[0][0] = 0.0;
    stdval[0][1] = 0.0;

    stdval[1][0] = 1.0;
    stdval[1][1] = 0.0;

    stdval[2][0] = 0.0;
    stdval[2][1] = 1.0;

    stdval[3][0] = (p[1]-Pc[1]);
    stdval[3][1] = (p[0]-Pc[0]);

    stdval[4][0] = (2.0*(p[0]-Pc[0]));
    stdval[4][1] = 0.0;

    stdval[5][0] = 0.0;
    stdval[5][1] = (2.0*(p[1]-Pc[1]));

    stdval[6][0] = 2.0*(p[0]-Pc[0])*(p[1]-Pc[1]);
    stdval[6][1] = pow(p[0]-Pc[0], 2);

    stdval[7][0] = pow(p[1]-Pc[1], 2);
    stdval[7][1] = 2.0*(p[1]-Pc[1])*(p[0]-Pc[0]);

    stdval[8][0] = 3.0*pow(p[0]-Pc[0], 2);
    stdval[8][1] = 0.0;

    stdval[9][0] = 0.0;
    stdval[9][1] = 3.0*pow(p[1]-Pc[1], 2);

    for(int j=0; j<ele_dof; j++)
        for(int k=0; k<2; k++)
            val[j][k] = stdval[j][k];
 
}

void BH2D::basis2Grad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
   val.resize(ele_dof);
//  int Msize1 = sqrt(ele_dof);
   for(int i=0; i<val.size(); i++){
       val[i].resize(3);
   }


  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
  
  std::vector<std::vector<double> > stdval;
  stdval.resize(10);
  for(int i=0; i<stdval.size(); i++){
      stdval[i].resize(2);
  }
  

  stdval[0][0] = 0.0;
  stdval[0][1] = 0.0;
  stdval[0][2] = 0.0;

  stdval[1][0] = 0.0;
  stdval[1][1] = 0.0;
  stdval[1][2] = 0.0;

  stdval[2][0] = 0.0;
  stdval[2][1] = 0.0;
  stdval[2][2] = 0.0;

  stdval[3][0] = 0.0;
  stdval[3][1] = 1.0;
  stdval[3][2] = 0.0;

  stdval[4][0] = 2.0;
  stdval[4][1] = 0.0;
  stdval[4][2] = 0.0;

  stdval[5][0] = 0.0;
  stdval[5][1] = 0.0;
  stdval[5][2] = 2.0;

  stdval[6][0] = 2.0*(p[1]-Pc[1]);
  stdval[6][1] = 2.0*(p[0]-Pc[0]);
  stdval[6][2] = 0.0;

  stdval[7][0] = 0.0;
  stdval[7][1] = 2.0*(p[1]-Pc[1]);
  stdval[7][2] = 2.0*(p[0]-Pc[0]);

  stdval[8][0] = 6.0*(p[0]-Pc[0]);
  stdval[8][1] = 0.0;
  stdval[8][2] = 0.0;

  stdval[9][0] = 0.0;
  stdval[9][1] = 0.0;
  stdval[9][2] = 6.0*(p[1]-Pc[1]);              

    for(int j=0; j<ele_dof; j++)
        for(int k=0; k<3; k++)
            val[j][k] = stdval[j][k];

}

void BH2D::basisValue(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val)
{
    MyPoint Pc;
    Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
    Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;
    int n=p.size();
    val.resize(ele_dof);
//    int Msize1 = sqrt(ele_dof);
    for(int i=0; i<val.size(); i++){
       val[i].resize(n);
    }

    std::vector<double> stdval;
    stdval.resize(10);

    for(int i=0; i<n; i++){

    stdval[0] = 1.0;
    stdval[1] = p[i][0]-Pc[0];
    stdval[2] = p[i][1]-Pc[1];
    stdval[3] = (p[i][0]-Pc[0])*(p[i][1]-Pc[1]);
    stdval[4] = pow(p[i][0]-Pc[0], 2);
    stdval[5] = pow(p[i][1]-Pc[1], 2);
    stdval[6] = pow(p[i][0]-Pc[0], 2)*(p[i][1]-Pc[1]);
    stdval[7] = pow(p[i][1]-Pc[1], 2)*(p[i][0]-Pc[0]);
    stdval[8] = pow(p[i][0]-Pc[0], 3);
    stdval[9] = pow(p[i][1]-Pc[1], 3);

    for(int j=0; j<ele_dof; j++)
        val[j][i]=stdval[j];
        
    }

}

void BH2D::basisGrad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

   int n=p.size();
   val.resize(ele_dof);

   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(2);
      }
   }
   
    std::vector<std::vector<double> > stdval;
    stdval.resize(10);
    for(int i=0; i<stdval.size(); i++){
        stdval[i].resize(2);
    }
    
    for(int i=0; i<n; i++){

    stdval[0][0] = 0.0;
    stdval[0][1] = 0.0;

    stdval[1][0] = 1.0;
    stdval[1][1] = 0.0;

    stdval[2][0] = 0.0;
    stdval[2][1] = 1.0;

    stdval[3][0] = p[i][1]-Pc[1];
    stdval[3][1] = p[i][0]-Pc[0];

    stdval[4][0] = 2.0*(p[i][0]-Pc[0]);
    stdval[4][1] = 0.0;

    stdval[5][0] = 0.0;
    stdval[5][1] = 2.0*(p[i][1]-Pc[1]);

    stdval[6][0] = 2.0*(p[i][0]-Pc[0])*(p[i][1]-Pc[1]);
    stdval[6][1] = pow(p[i][0]-Pc[0], 2);

    stdval[7][0] = pow(p[i][1]-Pc[1], 2);
    stdval[7][1] = 2.0*(p[i][1]-Pc[1])*(p[i][0]-Pc[0]);

    stdval[8][0] = 3.0*pow(p[i][0]-Pc[0], 2);
    stdval[8][1] = 0.0;

    stdval[9][0] = 0.0;
    stdval[9][1] = 3.0*pow(p[i][1]-Pc[1], 2);

    for(int j=0; j<ele_dof; j++)
        for(int k=0; k<2; k++)
            val[j][i][k] = stdval[j][k];
  	
  	
   }

}

void BH2D::basis2Grad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val)
{
  MyPoint Pc;
  Pc[0]=(v[0][0]+v[1][0]+v[2][0]+v[3][0])/4.0;
  Pc[1]=(v[0][1]+v[1][1]+v[2][1]+v[3][1])/4.0;

   int n=p.size();
   val.resize(ele_dof);

   for(int i=0; i<val.size(); i++){
      val[i].resize(n);
   }
   for(int i=0; i<val.size(); i++){
      for(int j=0; j<n; j++){
         val[i][j].resize(3);
      }
   }
   
    std::vector<std::vector<double> > stdval;
    stdval.resize(10);
    for(int i=0; i<stdval.size(); i++){
        stdval[i].resize(3);
    }
   
  for(int i=0; i<n; i++){
  stdval[0][0] = 0.0;
  stdval[0][1] = 0.0;
  stdval[0][2] = 0.0;

  stdval[1][0] = 0.0;
  stdval[1][1] = 0.0;
  stdval[1][2] = 0.0;

  stdval[2][0] = 0.0;
  stdval[2][1] = 0.0;
  stdval[2][2] = 0.0;

  stdval[3][0] = 0.0;
  stdval[3][1] = 1.0;
  stdval[3][2] = 0.0;

  stdval[4][0] = 2.0;
  stdval[4][1] = 0.0;
  stdval[4][2] = 0.0;

  stdval[5][0] = 0.0;
  stdval[5][1] = 0.0;
  stdval[5][2] = 2.0;

  stdval[6][0] = 2.0*(p[i][1]-Pc[1]);
  stdval[6][1] = 2.0*(p[i][0]-Pc[0]);
  stdval[6][2] = 0.0;

  stdval[7][0] = 0.0;
  stdval[7][1] = 2.0*(p[i][1]-Pc[1]);
  stdval[7][2] = 2.0*(p[i][0]-Pc[0]);

  stdval[8][0] = 6.0*(p[i][0]-Pc[0]);
  stdval[8][1] = 0.0;
  stdval[8][2] = 0.0;

  stdval[9][0] = 0.0;
  stdval[9][1] = 0.0;
  stdval[9][2] = 6.0*(p[i][1]-Pc[1]);
  
    for(int j=0; j<ele_dof; j++)
        for(int k=0; k<3; k++)
            val[j][i][k] = stdval[j][k];
  
  }
   
}


// */

void BH2D::build_Edge_Patch()
{
   int n_edge=mesh.n_edge();
   edge_patch.resize(n_edge);
   NPedge_patch.resize(n_edge);
   int n_ele=mesh.n_element();
   std::vector<int> bnd_ele_set;
   for(int i=0; i<n_ele; i++){
      int eleBndInfo=mesh.getEleBndInfo(i);
      if(eleBndInfo==1){
		 bnd_ele_set.push_back(i);
	  }
   }
   for(int i=0; i<n_ele; i++){
      std::vector<int> eleEdg=mesh.getEleEdg(i);
      for(int j=0; j<eleEdg.size(); j++){
         edge_patch[eleEdg[j]].push_back(i);
         NPedge_patch[eleEdg[j]].push_back(i);
         int edgBndInfo=mesh.getEdgBndInfo(eleEdg[j]);
         if(edgBndInfo==1){//deal the boundary edge for the periodic boundary condition
			std::vector<int> edgjVtx=mesh.getEdgVtx(eleEdg[j]);
			std::vector<MyPoint> edgjPnt=mesh.getPnt(edgjVtx);
			MyPoint edgjmid=midpoint(edgjPnt[0], edgjPnt[1]);
			double dfx=fabs(edgjPnt[0][0]-edgjPnt[1][0]);
			double dfy=fabs(edgjPnt[0][1]-edgjPnt[1][1]);
           if(dfy<1.0e-08){//边平行与x轴
		    for(int k=0; k<bnd_ele_set.size(); k++){
				if(edge_patch[eleEdg[j]].size()==2) break;
				int bk=bnd_ele_set[k];
                if(bk!=i){//如果第k个边界单元不是i, 则找出该单元的边界边, 并比较边界边与第i个单元的边界边j
					std::vector<int> elebkEdg=mesh.getEleEdg(bk);
					for(int m=0; m<elebkEdg.size(); m++){
         				int edgBndInfobk=mesh.getEdgBndInfo(elebkEdg[m]);
						if(edgBndInfobk==1){
							std::vector<int> edgbkVtx=mesh.getEdgVtx(elebkEdg[m]);
							std::vector<MyPoint> edgbkPnt=mesh.getPnt(edgbkVtx);
							MyPoint edgbkmid=midpoint(edgbkPnt[0], edgbkPnt[1]);
							double dxjbk=fabs(edgjmid[0]-edgbkmid[0]);
							if(dxjbk<1.0e-08){
//								edge_patch[eleEdg[j]].push_back(bk);
//								edge_patch[elebkEdg[m]].push_back(i);	
								break;
							}
						}
					}
				}
			}
		   }
           else{//边平行与y轴
		    for(int k=0; k<bnd_ele_set.size(); k++){
				if(edge_patch[eleEdg[j]].size()==2) break;
				int bk=bnd_ele_set[k];
                if(bk!=i){//如果第k个边界单元不是i, 则找出该单元的边界边, 并比较边界边与第i个单元的边界边j
					std::vector<int> elebkEdg=mesh.getEleEdg(bk);
					for(int m=0; m<elebkEdg.size(); m++){
         				int edgBndInfobk=mesh.getEdgBndInfo(elebkEdg[m]);
						if(edgBndInfobk==1){
							std::vector<int> edgbkVtx=mesh.getEdgVtx(elebkEdg[m]);
							std::vector<MyPoint> edgbkPnt=mesh.getPnt(edgbkVtx);
							MyPoint edgbkmid=midpoint(edgbkPnt[0], edgbkPnt[1]);
							double dyjbk=fabs(edgjmid[1]-edgbkmid[1]);
							if(dyjbk<1.0e-08){
//								edge_patch[eleEdg[j]].push_back(bk);	
//								edge_patch[elebkEdg[m]].push_back(i);
								break;
							}
						}
					}
				}
			}
		   }
	     }
      }
   }
// 
/*
   for(int i=0; i<mesh.n_edge(); i++){
		std::cout<<i<<" : ";
      for(int j=0; j<edge_patch[i].size(); j++){
	      std::cout<<"\t"<<edge_patch[i][j];
	  }

		std::cout<<"\n";
	}
   getchar();
// */
}

void BH2D::build_Element_Patch()
{
   int n_ele=mesh.n_element();
   element_patch.resize(n_ele);
   element_patch_edge.resize(n_ele);
   for(int i=0; i<n_ele; i++){
      std::vector<int> eleEdg=mesh.getEleEdg(i);
      for(int j=0; j<eleEdg.size(); j++){
         for(int k=0; k<edge_patch[eleEdg[j]].size(); k++){
            if(edge_patch[eleEdg[j]][k]!=i){
                element_patch[i].push_back(edge_patch[eleEdg[j]][k]);
                element_patch_edge[i].push_back(j);
            }
         }
      }
   }
// 
/*
   for(int i=0; i<n_ele; i++){
		std::cout<<i<<" : ";
      for(int j=0; j<element_patch[i].size(); j++){
	      std::cout<<"\t"<<element_patch[i][j];
	  }

		std::cout<<"\n";
	}
   getchar();
// */
}

void BH2D::buildSparsityPattern()
{
   int n_ele=mesh.n_element();
//   int ele_dof=9;
//   int v_size=ele_dof*n_ele;
   int n_dof=n_ele*ele_dof;
   int n_max_couple_dof=ele_dof*5;
   sparsity_pattern.reinit(n_dof, n_dof, n_max_couple_dof);

	//For Matrix A
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
			    sparsity_pattern.add(ele_dof*i+j, ele_dof*i+k);
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
					    sparsity_pattern.add(ele_dof*i+j, ele_dof*im+k);
				}
			} 
		}
	}
// 
/*
	//For Matrix B
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
			    sparsity_pattern.add(ele_dof*i+j, v_size+ele_dof*i+k);
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
				    sparsity_pattern.add(ele_dof*i+j, v_size+ele_dof*im+k);
				}
			} 
		}
	}
	//For Matrix C
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
			    sparsity_pattern.add(v_size+ele_dof*i+j, ele_dof*i+k);
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
				    sparsity_pattern.add(v_size+ele_dof*i+j, ele_dof*im+k);
				}
			} 
		}
	}
	//For Matrix D
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
			    sparsity_pattern.add(v_size+ele_dof*i+j, v_size+ele_dof*i+k);
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
//				    sparsity_pattern.add(v_size+ele_dof*i+j, v_size+ele_dof*im+k);
				}
			} 
		}
	}
// */	
	sparsity_pattern.compress();
	

// /*
	std::ofstream os("stiff_matrixUW.html");
	sparsity_pattern.print_svg(os);
// */
}

void BH2D::buildSparseMatrix()
{
    int n_ele=mesh.n_element();
//    int ele_dof=9;
//	int v_size=ele_dof*n_ele;
	StiffMatrix.reinit(sparsity_pattern);
	for(int i=0; i<n_ele; i++){
		std::vector< std::vector< std::vector<double> > > Ele_matrix;
		
	    getElementMatrix(i, Ele_matrix);
// 
	/*
		if(i==0)
for(int ii=0; ii<Ele_matrix[0].size();ii++){
    for(int jj=0; jj<Ele_matrix[0].size();jj++)
		std::cout<<Ele_matrix[0][ii][ele_dof*0+jj]<<"  ";
	std::cout<<std::endl;
}
// */

// 	/*
		for(int j=0; j<ele_dof; j++){
			for(int m=0; m<ele_dof; m++){
				StiffMatrix.add(ele_dof*i+j, ele_dof*i+m, Ele_matrix[0][j][m]);
//				StiffMatrix.add(ele_dof*i+j, ele_dof*i+m, sqrt(-_alpha_coef_())*theta*Ele_matrix[1][j][m]);
//				StiffMatrix.add(v_size+ele_dof*i+j, ele_dof*i+m, sqrt(-_alpha_coef_())*theta*Ele_matrix[2][j][m]);           
//				StiffMatrix.add(v_size+ele_dof*i+j, v_size+ele_dof*i+m, theta*Ele_matrix[3][j][m]);    
			}
		}
        std::vector<int> eleEdg=mesh.getEleEdg(i);
// /*        
		for(int k=0; k<eleEdg.size(); k++)
		{
			if(edge_patch[eleEdg[k]].size()>1)
			{
			    int ik=edge_patch[eleEdg[k]][0];
			    if(ik==i)
				ik=edge_patch[eleEdg[k]][1];
			    for(int j=0; j<ele_dof; j++)
			    {
				    for(int m=0; m<ele_dof; m++)
				    {
					    StiffMatrix.add(ele_dof*i+j, ele_dof*ik+m, Ele_matrix[0][j][ele_dof*(k+1)+m]);
//					    StiffMatrix.add(v_size+ele_dof*i+j, ele_dof*ik+m, sqrt(-_alpha_coef_())*theta*Ele_matrix[2][j][ele_dof*(k+1)+m]);     
				    }
			    } 
		    }
		}
		
// */ 
	}      
//	std::ofstream os("stiff_matrixUW.gnuplot");
//	StiffMatrix.get_sparsity_pattern().print_svg(os);
}

// /*
void BH2D::getElementMatrix(int n, std::vector<std::vector< std::vector<double> > > &Ele_matrix)
{
    int MatrixNO = 4;
    Ele_matrix.resize(MatrixNO);
//    int np=element_patch[n].size();
    int np=4;
    for(int i=0; i<MatrixNO; i++){
        Ele_matrix[i].resize(ele_dof);
        for(int j=0; j<ele_dof; j++)
            Ele_matrix[i][j].resize(ele_dof*(np+1), 0.0);
    }
//	int ele_dof=9;
   	std::vector<int> EV=mesh.getEleVtx(n);
   	std::vector<MyPoint> EP=mesh.getPnt(EV);
  	std::vector<MyPoint> GP=tmpEle.getGaussPnt();

	std::vector<double> GW=tmpEle.getGaussWeight();
	std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	std::vector<std::vector<double> > basis_value;
    basisValue(q_point, EP, basis_value);
	std::vector<std::vector<std::vector<double> > > basis_grad;
    basisGrad(q_point, EP, basis_grad);
    
    double h=EP[2][0]-EP[0][0];
    double h3=power(EP[2][0]-EP[0][0],1);
    
	for(int i = 0; i < q_point.size(); i++)
	{
		double Jxw = GW[i]*jacobian[i];
		for(int j=0; j<basis_grad.size(); j++){
			for(int k=0; k<basis_grad.size(); k++){
               	Ele_matrix[0][j][k]+=Jxw*(basis_grad[k][i][0]*basis_grad[j][i][0]+basis_grad[k][i][1]*basis_grad[j][i][1]); //B Matrix
			}
		}
	}
   	std::vector<double> ec(2);
   	ec[0]=(EP[0][0]+EP[1][0]+EP[2][0]+EP[3][0])/4.0;
   	ec[1]=(EP[0][1]+EP[1][1]+EP[2][1]+EP[3][1])/4.0;
    //     /*
    std::vector<int> eleEdg=mesh.getEleEdg(n);
    for(int i=0; i<eleEdg.size(); i++)
	{
		std::vector<int> edgVtx=mesh.getEdgVtx(eleEdg[i]);
// /*
		std::vector<MyPoint> edgi;
   		edgi=mesh.getPnt(edgVtx);
//		std::cout<<edgi[1][0]<<"   "<<edgi[1][1]<<std::endl;
// */
//  
/*
        std::vector<MyPoint> edgi(2);
   		edgi[0]=mesh.getPnt(edgVtx[0]);
   		edgi[1]=mesh.getPnt(edgVtx[1]);
//		std::cout<<edgi[1][0]<<"   "<<edgi[1][1]<<std::endl;
// */
        MyPoint midEdg=midpoint(edgi[0],edgi[1]);
		
//		std::cout<<ec[0]<<"   "<<ec[1]<<std::endl;

		int EdgLocation=1;
        if(fabs(midEdg[0]-ec[0])>1.0e-08){//该边平行y轴
			if(midEdg[0]-ec[0]>1.0e-08)
				EdgLocation=1;//该边是单元的右边
			else
				EdgLocation=3;//该边是单元的左边
		}
		else{//该边平行x轴
			if(midEdg[1]-ec[1]>1.0e-08)
				EdgLocation=2;//该边是单元的上边
			else
				EdgLocation=4;//该边是单元的下边Matrix
		}
//		int ik=element_patch[n][i];
     
	    if(edge_patch[eleEdg[i]].size()>1)     // if the edge has two neighbor.
		{
	 
		int ik=edge_patch[eleEdg[i]][0];
		if(ik==n)
			ik=edge_patch[eleEdg[i]][1];
// 
/*
        if(n==0 && i==0){
		    std::cout<<midEdg[0]<<"   "<<midEdg[1]<<std::endl;
			std::cout<<ec[0]<<"   "<<ec[1]<<std::endl;
			std::cout<< EdgLocation <<std::endl;
//			std::cout<<edgi[0][0]<<"   "<<edgi[0][1]<<"   "<<edgi[1][0]<<"   "<<edgi[1][1]<<std::endl;
		}
// */

	  	std::vector<int> EVk=mesh.getEleVtx(ik);
// /*
        std::vector<MyPoint> EPk;
        EPk=mesh.getPnt(EVk);
//		std::cout<<EPk[1][0]<<"   "<<EPk[1][1]<<std::endl;
// */
// 
/*
        std::vector<MyPoint> EPk(4);
        EPk[0]=mesh.getPnt(EVk[0]);
        EPk[1]=mesh.getPnt(EVk[1]);
        EPk[2]=mesh.getPnt(EVk[2]);
        EPk[3]=mesh.getPnt(EVk[3]);
		std::cout<<EPk[1][0]<<"   "<<EPk[1][1]<<std::endl;
// */		
		
		switch(EdgLocation){
			case 1: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();
	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointEE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
             		q_pointEE[ie][1]=GPy[ie];
          		}
///////////////////////////////////////////////////////////////////
// test close 
/*
          		for(int ie=0; ie<GPE.size(); ie++){
             		if(fabs(q_pointE[ie][0]-EP[1][0])<1.0e-07 || fabs(q_pointE[ie][0]-EP[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"Error1!!!"<<std::endl;
						getchar();
					}
             		if(fabs(q_pointEE[ie][0]-EPk[1][0])<1.0e-07 || fabs(q_pointEE[ie][0]-EPk[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"error1!!!"<<std::endl;
						getchar();
					}
          		}
// */ 
///////////////////////////////////////////////////////////////////
	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);
				std::vector<std::vector<std::vector<double> > > basis_grad2I;
    	  		basis2Grad(q_pointE, EP, basis_grad2I);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
				std::vector<std::vector<std::vector<double> > > basis_grad2E;
    	  		basis2Grad(q_pointEE, EPk, basis_grad2E);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
							Ele_matrix[0][m][j]+=Jxw*beta_0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_matrix[0][m][ele_dof*(i+1)+j]+=-Jxw*beta_0/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);																			
							
                  			Ele_matrix[0][m][j]+=-Jxw/2*(basis_gradI[j][ie][0]*basis_valueI[m][ie]
											+ basis_valueI[j][ie]*basis_gradI[m][ie][0]);

                  		    Ele_matrix[0][m][ele_dof*(i+1)+j]+=Jxw/2*(-basis_gradE[j][ie][0]*basis_valueI[m][ie]
											+ basis_valueE[j][ie]*basis_gradI[m][ie][0]);
											
							Ele_matrix[0][m][j]+=Jxw*beta1*h*(basis_grad2I[j][ie][0]*basis_valueI[m][ie]);

                  		    Ele_matrix[0][m][ele_dof*(i+1)+j]+=-Jxw*beta1*h*(basis_grad2E[j][ie][0]*basis_valueI[m][ie]);

//							Ele_matrix[2][m][j]-=Jxw/2*(basis_gradI[j][ie][0]*basis_valueI[m][ie]
//											+ basis_valueI[j][ie]*basis_gradI[m][ie][0]);

//                  		    Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw/2*(-basis_gradE[j][ie][0]*basis_valueI[m][ie]
//											+ basis_valueE[j][ie]*basis_gradI[m][ie][0]);
											
//							Ele_matrix[1][m][j]+=Jxw*beta0x/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[1][m][ele_dof*(i+1)+j]-=Jxw*beta0x/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);

//							Ele_matrix[2][m][j]-=Jxw*beta0x/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw*beta0x/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);								
						}
					}
          		}
				break;	
			}
			case 2: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointEE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
             		q_pointEE[ie][0]=GPx[ie];
          		}
////////////////////////////////////////////////////////////////////
// test close 
/*
          		for(int ie=0; ie<GPE.size(); ie++){
             		if(fabs(q_pointE[ie][1]-EP[1][1])<1.0e-07 || fabs(q_pointE[ie][1]-EP[3][1])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"Error2!!!"<<std::endl;
						getchar();
					}
             		if(fabs(q_pointEE[ie][1]-EPk[1][1])<1.0e-07 || fabs(q_pointEE[ie][1]-EPk[3][1])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"error2!!!"<<std::endl;
						getchar();
					}
          		}
// */ 
////////////////////////////////////////////////////////////////////
	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);
				std::vector<std::vector<std::vector<double> > > basis_grad2I;
    	  		basis2Grad(q_pointE, EP, basis_grad2I);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
				std::vector<std::vector<std::vector<double> > > basis_grad2E;
    	  		basis2Grad(q_pointEE, EPk, basis_grad2E);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
							
							Ele_matrix[0][m][j]+=Jxw*beta_0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_matrix[0][m][ele_dof*(i+1)+j]-=Jxw*beta_0/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);								
							
                  			Ele_matrix[0][m][j]-=Jxw/2*(basis_gradI[j][ie][1]*basis_valueI[m][ie]
											+basis_valueI[j][ie]*basis_gradI[m][ie][1]);

                  		    Ele_matrix[0][m][ele_dof*(i+1)+j]+=Jxw/2*(-basis_gradE[j][ie][1]*basis_valueI[m][ie]
											+basis_valueE[j][ie]*basis_gradI[m][ie][1]);
											
							Ele_matrix[0][m][j]+=Jxw*beta1*h*(basis_grad2I[j][ie][2]*basis_valueI[m][ie]);
											
							Ele_matrix[0][m][ele_dof*(i+1)+j]+=-Jxw*beta1*h*(basis_grad2E[j][ie][2]*basis_valueI[m][ie]);	
											
//							Ele_matrix[2][m][j]-=Jxw/2*(basis_gradI[j][ie][1]*basis_valueI[m][ie]
//											+basis_valueI[j][ie]*basis_gradI[m][ie][1]);

 //                 		    Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw/2*(-basis_gradE[j][ie][1]*basis_valueI[m][ie]
//											+basis_valueE[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_matrix[1][m][j]+=Jxw*beta0y/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[1][m][ele_dof*(i+1)+j]-=Jxw*beta0y/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);

//							Ele_matrix[2][m][j]-=Jxw*beta0y/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw*beta0y/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);												
						}
					}
          		}
				break;	
			}
			case 3: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointEE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
             		q_pointEE[ie][1]=GPy[ie];
          		}
////////////////////////////////////////////////////////////////////
// test close 
/*
          		for(int ie=0; ie<GPE.size(); ie++){
             		if(fabs(q_pointE[ie][0]-EP[1][0])<1.0e-07 || fabs(q_pointE[ie][0]-EP[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"Error3!!!"<<std::endl;
						getchar();
					}
             		if(fabs(q_pointEE[ie][0]-EPk[1][0])<1.0e-07 || fabs(q_pointEE[ie][0]-EPk[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"error3!!!"<<std::endl;
						getchar();
					}
          		}
// */
////////////////////////////////////////////////////////////////////
	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);
				std::vector<std::vector<std::vector<double> > > basis_grad2I;
    	  		basis2Grad(q_pointE, EP, basis_grad2I);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
				std::vector<std::vector<std::vector<double> > > basis_grad2E;
    	  		basis2Grad(q_pointEE, EPk, basis_grad2E);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
							Ele_matrix[0][m][j]+=Jxw*beta_0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_matrix[0][m][ele_dof*(i+1)+j]+=-Jxw*beta_0/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);															
							
                  			Ele_matrix[0][m][j]+=Jxw/2*(basis_gradI[j][ie][0]*basis_valueI[m][ie]
											+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
											
                  		    Ele_matrix[0][m][ele_dof*(i+1)+j]+=Jxw/2*(basis_gradE[j][ie][0]*basis_valueI[m][ie]
											-basis_valueE[j][ie]*basis_gradI[m][ie][0]);
							
							Ele_matrix[0][m][j]+=Jxw*beta1*h*(basis_grad2I[j][ie][0]*basis_valueI[m][ie]);

                  		    Ele_matrix[0][m][ele_dof*(i+1)+j]+=-Jxw*beta1*h*(basis_grad2E[j][ie][0]*basis_valueI[m][ie]);

							
//							Ele_matrix[2][m][j]+=Jxw/2*(basis_gradI[j][ie][0]*basis_valueI[m][ie]
//											+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
											
//                  		    Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw/2*(basis_gradE[j][ie][0]*basis_valueI[m][ie]
//											-basis_valueE[j][ie]*basis_gradI[m][ie][0]);
											
//							Ele_matrix[1][m][j]+=Jxw*beta0x/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[1][m][ele_dof*(i+1)+j]-=Jxw*beta0x/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);

//							Ele_matrix[2][m][j]-=Jxw*beta0x/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw*beta0x/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);																							
						}
					}
          		}
				break;	
			}
			case 4: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointEE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
             		q_pointEE[ie][0]=GPx[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);
				std::vector<std::vector<std::vector<double> > > basis_grad2I;
    	  		basis2Grad(q_pointE, EP, basis_grad2I);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
				std::vector<std::vector<std::vector<double> > > basis_grad2E;
    	  		basis2Grad(q_pointEE, EPk, basis_grad2E);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
							Ele_matrix[0][m][j]+=Jxw*beta_0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_matrix[0][m][ele_dof*(i+1)+j]+=-Jxw*beta_0/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);		
							
                  			Ele_matrix[0][m][j]+=Jxw/2*(basis_gradI[j][ie][1]*basis_valueI[m][ie]
											+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
                  		    Ele_matrix[0][m][ele_dof*(i+1)+j]+=Jxw/2*(basis_gradE[j][ie][1]*basis_valueI[m][ie]
											-basis_valueE[j][ie]*basis_gradI[m][ie][1]);
							
							Ele_matrix[0][m][j]+=Jxw*beta1*h*(basis_grad2I[j][ie][2]*basis_valueI[m][ie]);
											
							Ele_matrix[0][m][ele_dof*(i+1)+j]+=-Jxw*beta1*h*(basis_grad2E[j][ie][2]*basis_valueI[m][ie]);	
							
//							Ele_matrix[2][m][j]+=Jxw/2*(basis_gradI[j][ie][1]*basis_valueI[m][ie]
//											+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
//                  		    Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw/2*(basis_gradE[j][ie][1]*basis_valueI[m][ie]
//											-basis_valueE[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_matrix[1][m][j]+=Jxw*beta0y/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[1][m][ele_dof*(i+1)+j]-=Jxw*beta0y/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);

//							Ele_matrix[2][m][j]-=Jxw*beta0y/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw*beta0y/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);												
						}
					}
          		}
				break;
			}
			default:
			{
				std::cout<<"Some problem in location the edg, please chech carefully!"<<std::endl;
				break;	
			}
		}
		
		
		}
		else     // if the edge has only one neighbor.
		{
		
		
		switch(EdgLocation){
			case 1: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();
	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);
				
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_matrix[0][m][j]+=theta*Jxw*beta_0/h3*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[0][m][ele_dof*(i+1)+j]-=Jxw*beta_0/h3*(basis_valueE[j][ie]*basis_valueI[m][ie]);																			
							
                  			Ele_matrix[0][m][j]+=-Jxw*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);

//							Ele_matrix[2][m][j]-=Jxw*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
							
							Ele_matrix[0][m][j]+=Jxw*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_matrix[2][m][j]+=Jxw*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);

											
//							Ele_matrix[1][m][j]+=Jxw*beta0x/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[1][m][ele_dof*(i+1)+j]-=Jxw*beta0x/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);

//							Ele_matrix[2][m][j]-=Jxw*beta0x/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw*beta0x/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);								
						}
					}
          		}
				break;	
			}
			case 2: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
							
//							Ele_matrix[0][m][j]+=theta*Jxw*beta_0/h3*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[0][m][ele_dof*(i+1)+j]-=Jxw*beta_0/h3*(basis_valueE[j][ie]*basis_valueI[m][ie]);								
							
//                  			Ele_matrix[0][m][j]+=-Jxw*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_matrix[2][m][j]-=Jxw*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
							
//							Ele_matrix[0][m][j]+=Jxw*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_matrix[2][m][j]+=Jxw*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);							
											
//							Ele_matrix[1][m][j]+=Jxw*beta0y/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[1][m][ele_dof*(i+1)+j]-=Jxw*beta0y/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);

//							Ele_matrix[2][m][j]-=Jxw*beta0y/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw*beta0y/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);												
						}
					}
          		}
				break;	
			}
			case 3: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_matrix[0][m][j]+=theta*Jxw*beta_0/h3*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[0][m][ele_dof*(i+1)+j]-=Jxw*beta_0/h3*(basis_valueE[j][ie]*basis_valueI[m][ie]);															
							
                  			Ele_matrix[0][m][j]+=Jxw*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
											
//							Ele_matrix[2][m][j]+=Jxw*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
							
							Ele_matrix[0][m][j]+=Jxw*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_matrix[2][m][j]+=Jxw*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);											
											
//							Ele_matrix[1][m][j]+=Jxw*beta0x/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[1][m][ele_dof*(i+1)+j]-=Jxw*beta0x/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);

//							Ele_matrix[2][m][j]-=Jxw*beta0x/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw*beta0x/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);																							
						}
					}
          		}
				break;	
			}
			case 4: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 

          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_matrix[0][m][j]+=theta*Jxw*beta_0/h3*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[0][m][ele_dof*(i+1)+j]-=Jxw*beta_0/h3*(basis_valueE[j][ie]*basis_valueI[m][ie]);		
							
//                  			Ele_matrix[0][m][j]+=Jxw*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_matrix[2][m][j]+=Jxw*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
							
//							Ele_matrix[0][m][j]+=Jxw*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_matrix[2][m][j]+=Jxw*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);	
											
//							Ele_matrix[1][m][j]+=Jxw*beta0y/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[1][m][ele_dof*(i+1)+j]-=Jxw*beta0y/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);

//							Ele_matrix[2][m][j]-=Jxw*beta0y/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_matrix[2][m][ele_dof*(i+1)+j]+=Jxw*beta0y/h*(basis_valueE[j][ie]*basis_valueI[m][ie]);												
						}
					}
          		}
				break;
			}
			default:
			{
				std::cout<<"Some problem in location the edg, please chech carefully!"<<std::endl;
				break;	
			}
		}		
		
		
		
		}
	}
// */
}
// */

void BH2D::buildRHS()
{
    int n_ele=mesh.n_element();
//    int ele_dof=9;
    RHS.reinit(ele_dof*n_ele);
	int v_size=ele_dof*n_ele;
    for(int i=0; i<n_ele; i++){
		std::vector<std::vector<double> > Ele_Rhs;
		getElementRHS(i, Ele_Rhs);
		for(int j=0; j<ele_dof; j++){
			RHS(ele_dof*i+j)=Ele_Rhs[0][j];
//			RHS(v_size+ele_dof*i+j)=Ele_Rhs[1][j]; 
        }
    }
}

void BH2D::getElementRHS(int n, std::vector<std::vector<double> > &Ele_Rhs)
{
//	int ele_dof=9;
    int RHSNO = 1;
    Ele_Rhs.resize(RHSNO);
    for(int l=0; l<RHSNO; l++){
        Ele_Rhs[l].resize(ele_dof, 0.0);
    }
 //   Ele_Rhs.resize(ele_dof, 0.0);
    std::vector<int> EV=mesh.getEleVtx(n);
    std::vector<MyPoint> EP=mesh.getPnt(EV);
    std::vector<MyPoint> GP=tmpEle.getGaussPnt();

    std::vector<double> GW=tmpEle.getGaussWeight();
    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
    std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
    std::vector<std::vector<double> > basis_value;
    basisValue(q_point, EP, basis_value);
	std::vector<std::vector<std::vector<double> > > basis_grad;
    basisGrad(q_point, EP, basis_grad);
		
//    std::vector<double> u_coef;
//    Read_u_Coef(n, u_coef);
//    Read_q_Coef(n, q_coef);
    
    double h=EP[2][0]-EP[0][0];
    double h3=power(EP[2][0]-EP[0][0],1);
	
	for(int j=0; j<q_point.size(); j++){
		double Jxw = GW[j]*jacobian[j];
		
		for(int k = 0; k < basis_value.size(); k++){
			for(int l = 0; l < basis_value.size(); l++){
                Ele_Rhs[0][k] += Jxw*(q_i[0]*c1_h(n*ele_dof+l)+q_i[1]*c2_h(n*ele_dof+l))*basis_value[l][j]*basis_value[k][j];
			}
			Ele_Rhs[0][k] += Jxw*( _rho_(t,q_point[j]) )*basis_value[k][j];
			
	  	}
	}
	
   	std::vector<double> ec(2);
   	ec[0]=(EP[0][0]+EP[1][0]+EP[2][0]+EP[3][0])/4.0;
   	ec[1]=(EP[0][1]+EP[1][1]+EP[2][1]+EP[3][1])/4.0;
    //     /*
    std::vector<int> eleEdg=mesh.getEleEdg(n);
    for(int i=0; i<eleEdg.size(); i++)
	{
		std::vector<int> edgVtx=mesh.getEdgVtx(eleEdg[i]);
// /*
		std::vector<MyPoint> edgi;
   		edgi=mesh.getPnt(edgVtx);
//		std::cout<<edgi[1][0]<<"   "<<edgi[1][1]<<std::endl;
// */
//  
/*
        std::vector<MyPoint> edgi(2);
   		edgi[0]=mesh.getPnt(edgVtx[0]);
   		edgi[1]=mesh.getPnt(edgVtx[1]);
//		std::cout<<edgi[1][0]<<"   "<<edgi[1][1]<<std::endl;
// */
        MyPoint midEdg=midpoint(edgi[0],edgi[1]);
		
//		std::cout<<ec[0]<<"   "<<ec[1]<<std::endl;

		int EdgLocation=1;
        if(fabs(midEdg[0]-ec[0])>1.0e-08){//该边平行y轴
			if(midEdg[0]-ec[0]>1.0e-08)
				EdgLocation=1;//该边是单元的右边
			else
				EdgLocation=3;//该边是单元的左边
		}
		else{//该边平行x轴
			if(midEdg[1]-ec[1]>1.0e-08)
				EdgLocation=2;//该边是单元的上边
			else
				EdgLocation=4;//该边是单元的下边
		}
//		int ik=element_patch[n][i];


        if(edge_patch[eleEdg[i]].size()>1)        // if the edge has two neighbor.
		{

		int ik=edge_patch[eleEdg[i]][0];
		if(ik==n)
			ik=edge_patch[eleEdg[i]][1];
		
//		std::vector<double> u_coef0, q_coef0;
 //       Read_u_Coef(ik, u_coef0);
//        Read_q_Coef(ik, q_coef0);
// 
/*
        if(n==0 && i==0){
		    std::cout<<midEdg[0]<<"   "<<midEdg[1]<<std::endl;
			std::cout<<ec[0]<<"   "<<ec[1]<<std::endl;
			std::cout<< EdgLocation <<std::endl;
//			std::cout<<edgi[0][0]<<"   "<<edgi[0][1]<<"   "<<edgi[1][0]<<"   "<<edgi[1][1]<<std::endl;
		}
// */

	  	std::vector<int> EVk=mesh.getEleVtx(ik);
// /*
        std::vector<MyPoint> EPk;
        EPk=mesh.getPnt(EVk);
//		std::cout<<EPk[1][0]<<"   "<<EPk[1][1]<<std::endl;
// */
// 
/*
        std::vector<MyPoint> EPk(4);
        EPk[0]=mesh.getPnt(EVk[0]);
        EPk[1]=mesh.getPnt(EVk[1]);
        EPk[2]=mesh.getPnt(EVk[2]);
        EPk[3]=mesh.getPnt(EVk[3]);
		std::cout<<EPk[1][0]<<"   "<<EPk[1][1]<<std::endl;
// */		
		
		switch(EdgLocation){
			case 1: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();
	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointEE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
             		q_pointEE[ie][1]=GPy[ie];
          		}
///////////////////////////////////////////////////////////////////
// test close 
/*
          		for(int ie=0; ie<GPE.size(); ie++){
             		if(fabs(q_pointE[ie][0]-EP[1][0])<1.0e-07 || fabs(q_pointE[ie][0]-EP[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"Error1!!!"<<std::endl;
						getchar();
					}
             		if(fabs(q_pointEE[ie][0]-EPk[1][0])<1.0e-07 || fabs(q_pointEE[ie][0]-EPk[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"error1!!!"<<std::endl;
						getchar();
					}
          		}
// */ 
///////////////////////////////////////////////////////////////////
	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*beta_0/h3*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);								
							
 //                 			Ele_Rhs[0][m]+=Jxw/2*antitheta*sqrt(-_alpha_coef_())*(q_coef[j]*basis_gradI[j][ie][0]*basis_valueI[m][ie]
//											+q_coef[j]*basis_valueI[j][ie]*basis_gradI[m][ie][0]);

//                  		    Ele_Rhs[0][m]-=Jxw/2*antitheta*sqrt(-_alpha_coef_())*(-q_coef0[j]*basis_gradE[j][ie][0]*basis_valueI[m][ie]
//											+q_coef0[j]*basis_valueE[j][ie]*basis_gradI[m][ie][0]);

//							Ele_Rhs[1][m]+=0.0*Jxw/2*antitheta*sqrt(-_alpha_coef_())*(u_coef[j]*basis_gradI[j][ie][0]*basis_valueI[m][ie]
//											+u_coef[j]*basis_valueI[j][ie]*basis_gradI[m][ie][0]);

//                  		    Ele_Rhs[1][m]-=0.0*Jxw/2*antitheta*sqrt(-_alpha_coef_())*(-u_coef0[j]*basis_gradE[j][ie][0]*basis_valueI[m][ie]
//											+ u_coef0[j]*basis_valueE[j][ie]*basis_gradI[m][ie][0]);

//							Ele_Rhs[0][m]-=Jxw*beta0x/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0x/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                            Ele_Rhs[1][m]+=Jxw*beta0x/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0x/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);							
						}
					}
          		}
				break;	
			}
			case 2: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointEE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
             		q_pointEE[ie][0]=GPx[ie];
          		}
////////////////////////////////////////////////////////////////////
// test close 
/*
          		for(int ie=0; ie<GPE.size(); ie++){
             		if(fabs(q_pointE[ie][1]-EP[1][1])<1.0e-07 || fabs(q_pointE[ie][1]-EP[3][1])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"Error2!!!"<<std::endl;
						getchar();
					}
             		if(fabs(q_pointEE[ie][1]-EPk[1][1])<1.0e-07 || fabs(q_pointEE[ie][1]-EPk[3][1])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"error2!!!"<<std::endl;
						getchar();
					}
          		}
// */ 
////////////////////////////////////////////////////////////////////
	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*beta_0/h3*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
//                  			Ele_Rhs[0][m]+=Jxw/2*antitheta*sqrt(-_alpha_coef_())*(q_coef[j]*basis_gradI[j][ie][1]*basis_valueI[m][ie]
//											+q_coef[j]*basis_valueI[j][ie]*basis_gradI[m][ie][1]);
//
//                  		    Ele_Rhs[0][m]-=Jxw/2*antitheta*sqrt(-_alpha_coef_())*(-q_coef0[j]*basis_gradE[j][ie][1]*basis_valueI[m][ie]
//											+q_coef0[j]*basis_valueE[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_Rhs[1][m]+=0.0*Jxw/2*antitheta*sqrt(-_alpha_coef_())*(u_coef[j]*basis_gradI[j][ie][1]*basis_valueI[m][ie]
//											+u_coef[j]*basis_valueI[j][ie]*basis_gradI[m][ie][1]);

//                  		    Ele_Rhs[1][m]-=0.0*Jxw/2*antitheta*sqrt(-_alpha_coef_())*(-u_coef0[j]*basis_gradE[j][ie][1]*basis_valueI[m][ie]
//											+u_coef0[j]*basis_valueE[j][ie]*basis_gradI[m][ie][1]);	

//							Ele_Rhs[0][m]-=Jxw*beta0y/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0y/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                           Ele_Rhs[1][m]+=Jxw*beta0y/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0y/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
					}
          		}
				break;	
			}
			case 3: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointEE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
             		q_pointEE[ie][1]=GPy[ie];
          		}
////////////////////////////////////////////////////////////////////
// test close 
/*
          		for(int ie=0; ie<GPE.size(); ie++){
             		if(fabs(q_pointE[ie][0]-EP[1][0])<1.0e-07 || fabs(q_pointE[ie][0]-EP[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"Error3!!!"<<std::endl;
						getchar();
					}
             		if(fabs(q_pointEE[ie][0]-EPk[1][0])<1.0e-07 || fabs(q_pointEE[ie][0]-EPk[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"error3!!!"<<std::endl;
						getchar();
					}
          		}
// */
////////////////////////////////////////////////////////////////////
	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*beta_0/h3*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
//                  			Ele_Rhs[0][m]-=Jxw/2*antitheta*sqrt(-_alpha_coef_())*(q_coef[j]*basis_gradI[j][ie][0]*basis_valueI[m][ie]
//											+q_coef[j]*basis_valueI[j][ie]*basis_gradI[m][ie][0]);
											
 //                 		    Ele_Rhs[0][m]-=Jxw/2*antitheta*sqrt(-_alpha_coef_())*(q_coef0[j]*basis_gradE[j][ie][0]*basis_valueI[m][ie]
//											-q_coef0[j]*basis_valueE[j][ie]*basis_gradI[m][ie][0]);
											
//							Ele_Rhs[1][m]-=0.0*Jxw/2*antitheta*sqrt(-_alpha_coef_())*(u_coef[j]*basis_gradI[j][ie][0]*basis_valueI[m][ie]
//											+u_coef[j]*basis_valueI[j][ie]*basis_gradI[m][ie][0]);
											
//                  		    Ele_Rhs[1][m]-=0.0*Jxw/2*antitheta*sqrt(-_alpha_coef_())*(u_coef0[j]*basis_gradE[j][ie][0]*basis_valueI[m][ie]
//											-u_coef0[j]*basis_valueE[j][ie]*basis_gradI[m][ie][0]);
											
//							Ele_Rhs[0][m]-=Jxw*beta0x/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0x/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

 //                           Ele_Rhs[1][m]+=Jxw*beta0x/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0x/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
					}
          		}
				break;	
			}
			case 4: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointEE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
             		q_pointEE[ie][0]=GPx[ie];
          		}
////////////////////////////////////////////////////////////////////
// test close 
/*
          		for(int ie=0; ie<GPE.size(); ie++){
             		if(fabs(q_pointE[ie][1]-EP[1][1])<1.0e-07 || fabs(q_pointE[ie][1]-EP[3][1])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"Error4!!!"<<std::endl;
						getchar();
					}
             		if(fabs(q_pointEE[ie][1]-EPk[1][1])<1.0e-07 || fabs(q_pointEE[ie][1]-EPk[3][1])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"error4!!!"<<std::endl;
						getchar();
					}
          		}
//  */
////////////////////////////////////////////////////////////////////
	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*beta_0/h3*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
 //                 			Ele_Rhs[0][m]-=Jxw/2*antitheta*sqrt(-_alpha_coef_())*(q_coef[j]*basis_gradI[j][ie][1]*basis_valueI[m][ie]
//											+q_coef[j]*basis_valueI[j][ie]*basis_gradI[m][ie][1]);
 //                 		    Ele_Rhs[0][m]-=Jxw/2*antitheta*sqrt(-_alpha_coef_())*(q_coef0[j]*basis_gradE[j][ie][1]*basis_valueI[m][ie]
//											-q_coef0[j]*basis_valueE[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_Rhs[1][m]-=0.0*Jxw/2*antitheta*sqrt(-_alpha_coef_())*(u_coef[j]*basis_gradI[j][ie][1]*basis_valueI[m][ie]
//											+u_coef[j]*basis_valueI[j][ie]*basis_gradI[m][ie][1]);
//                  		    Ele_Rhs[1][m]-=0.0*Jxw/2*antitheta*sqrt(-_alpha_coef_())*(u_coef0[j]*basis_gradE[j][ie][1]*basis_valueI[m][ie]
//											-u_coef0[j]*basis_valueE[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_Rhs[0][m]-=Jxw*beta0y/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0y/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                            Ele_Rhs[1][m]+=Jxw*beta0y/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0y/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
					}
          		}
				break;
			}
			default:
			{
				std::cout<<"Some problem in location the edg, please chech carefully!"<<std::endl;
				break;	
			}
		}
		
		
	    }// end of the two neighbor case.
		else
		{
			

			
		switch(EdgLocation){
			case 1: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();
	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);								
							
//                  			Ele_Rhs[0][m]+=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);

//							Ele_Rhs[1][m]+=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);

//							Ele_Rhs[0][m]-=Jxw*beta0x/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0x/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                            Ele_Rhs[1][m]+=Jxw*beta0x/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0x/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);							
						}
//						std::cout<<"_u_grad_(t,q_pointE[ie])[0]= "<<_u_grad_(t,q_pointE[ie])[0]<<", _u_grad_(t,q_pointE[ie])[1]= "<<_u_grad_(t,q_pointE[ie])[1]<<std::endl;
						Ele_Rhs[0][m]+=Jxw*beta_0/h*(_u_(t,q_pointE[ie])*basis_valueI[m][ie]);
						Ele_Rhs[0][m]+=-Jxw*(_u_(t,q_pointE[ie])*basis_gradI[m][ie][0]);

//						Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[0][m]-=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][0]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
						
//						Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]-=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][0]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
									
//						Ele_Rhs[1][m]+=Jxw*theta*(_u_grad_(t,q_pointE[ie])[0]*basis_valueI[m][ie]- _u_(t,q_pointE[ie])*basis_gradI[m][ie][0])
//						               +Jxw*antitheta*(_u_grad_(t-dt,q_pointE[ie])[0]*basis_valueI[m][ie]-_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
			
					}
          		}
				break;	
			}
			case 2: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
//                  			Ele_Rhs[0][m]+=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_Rhs[1][m]+=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);

//							Ele_Rhs[0][m]-=Jxw*beta0y/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0y/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                           Ele_Rhs[1][m]+=Jxw*beta0y/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0y/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
					
//					    Ele_Rhs[0][m]+=Jxw*beta_0/h*(_u_(t,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[0][m]+=-Jxw*(_u_(t,q_pointE[ie])*basis_gradI[m][ie][1]);
						Ele_Rhs[0][m]+=Jxw*(_u_grad_(t,q_pointE[ie])[1]*basis_valueI[m][ie]);
					    
//					    Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//					    Ele_Rhs[0][m]-=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][1]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
					    
//					    Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//					    Ele_Rhs[1][m]-=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][1]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);


//						Ele_Rhs[1][m]+=Jxw*theta*(_u_grad_(t,q_pointE[ie])[1]*basis_valueI[m][ie]-_u_(t,q_pointE[ie])*basis_gradI[m][ie][1])
//						               +Jxw*antitheta*(_u_grad_(t-dt,q_pointE[ie])[1]*basis_valueI[m][ie]-_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
					
					}
          		}
				break;	
			}
			case 3: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
//                  			Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
											
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
//							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]-=Jxw*beta0x/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0x/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

 //                           Ele_Rhs[1][m]+=Jxw*beta0x/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0x/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
						
						Ele_Rhs[0][m]+=Jxw*beta_0/h*(_u_(t,q_pointE[ie])*basis_valueI[m][ie]);
						Ele_Rhs[0][m]+=Jxw*(_u_(t,q_pointE[ie])*basis_gradI[m][ie][0]);
						
//						Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[0][m]+=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][0]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
						
//						Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]+=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][0]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);


//						Ele_Rhs[1][m]+=Jxw*theta*(-_u_grad_(t,q_pointE[ie])[0]*basis_valueI[m][ie]+_u_(t,q_pointE[ie])*basis_gradI[m][ie][0])
//						               +Jxw*antitheta*(-_u_grad_(t-dt,q_pointE[ie])[0]*basis_valueI[m][ie]+_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
					}
          		}
				break;	
			}
			case 4: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
//                  			Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]-=Jxw*beta0y/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0y/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                            Ele_Rhs[1][m]+=Jxw*beta0y/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0y/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
//						std::cout<<"t= "<<t<<", dt= "<<dt<<", q_pointE[ie]= "<<q_pointE[ie][0]<<" , q_pointE[ie][1]= "<<q_pointE[ie][1]<<",  _u_(t,q_pointE(ie))= "<<_u_(t,q_pointE[ie])<<std::endl;

//						Ele_Rhs[0][m]+=Jxw*beta_0/h*(_u_(t,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[0][m]+=Jxw*(_u_(t,q_pointE[ie])*basis_gradI[m][ie][1]);
						Ele_Rhs[0][m]+=-Jxw*(_u_grad_(t,q_pointE[ie])[1]*basis_valueI[m][ie]);
						
//						Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[0][m]+=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][1]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
//						
//						Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]+=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][1]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
						
//						Ele_Rhs[1][m]+=Jxw*theta*(-_u_grad_(t,q_pointE[ie])[1]*basis_valueI[m][ie]+_u_(t,q_pointE[ie])*basis_gradI[m][ie][1])
//						               +Jxw*antitheta*(-_u_grad_(t-dt,q_pointE[ie])[1]*basis_valueI[m][ie]+_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
						               
					}
          		}
				break;
			}
			default:
			{
				std::cout<<"Some problem in location the edg, please chech carefully!"<<std::endl;
				break;	
			}
		}
				
		} // end of the one neighbor case.
		
			
		
		
	}	
	
}

void BH2D::buildRHS_C()
{
    int n_ele=mesh.n_element();
//    int ele_dof=9;
    RHS0.reinit(2*ele_dof*n_ele);
    for(int i=0; i<n_ele; i++){
		std::vector<std::vector<double> > Ele_Rhs, Ele_RhsC;
		getElementRHS_g(i, Ele_Rhs);
		getElementRHS_C(i, Ele_RhsC);
		for(int j=0; j<ele_dof; j++){
			RHS0(ele_dof*i+j)=dt*Ele_Rhs[0][j]+Ele_RhsC[0][j];
			RHS0(v_size+ele_dof*i+j)=dt*Ele_Rhs[1][j]+Ele_RhsC[1][j]; 
        }
    }
}

void BH2D::getElementRHS_C(int n, std::vector<std::vector<double> > &Ele_Rhs)
{
//	int ele_dof=9;
    int RHSNO = 2;
    Ele_Rhs.resize(RHSNO);
    for(int l=0; l<RHSNO; l++){
        Ele_Rhs[l].resize(ele_dof, 0.0);
    }
 //   Ele_Rhs.resize(ele_dof, 0.0);
    std::vector<int> EV=mesh.getEleVtx(n);
    std::vector<MyPoint> EP=mesh.getPnt(EV);
    std::vector<MyPoint> GP=tmpEle.getGaussPnt();

    std::vector<double> GW=tmpEle.getGaussWeight();
    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
    std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
    std::vector<std::vector<double> > basis_value;
    basisValue(q_point, EP, basis_value);
	std::vector<std::vector<std::vector<double> > > basis_grad;
    basisGrad(q_point, EP, basis_grad);
		
//    std::vector<double> u_coef;
//    Read_u_Coef(n, u_coef);
//    Read_q_Coef(n, q_coef);
    
    double h=EP[2][0]-EP[0][0];
    double h3=power(EP[2][0]-EP[0][0],1);
	
	for(int j=0; j<q_point.size(); j++){
		double Jxw = GW[j]*jacobian[j];
		
		double phi_val=0.0;
		for(int k=0; k<basis_value.size(); k++)
		{
			phi_val += u_h(n*ele_dof+k)*basis_value[k][j];
		}
		
		for(int k = 0; k < basis_value.size(); k++){
			for(int l = 0; l < basis_value.size(); l++){
                Ele_Rhs[0][k] += Jxw*c1_h(n*ele_dof+l)*basis_value[l][j]*basis_value[k][j];
				Ele_Rhs[1][k] += Jxw*c2_h(n*ele_dof+l)*basis_value[l][j]*basis_value[k][j];
			}
			
	  	}
	}
	
}

void BH2D::getElementRHS_g(int n, std::vector<std::vector<double> > &Ele_Rhs)
{
//	int ele_dof=9;
    int RHSNO = 2;
    Ele_Rhs.resize(RHSNO);
    for(int l=0; l<RHSNO; l++){
        Ele_Rhs[l].resize(ele_dof, 0.0);
    }
 //   Ele_Rhs.resize(ele_dof, 0.0);
    std::vector<int> EV=mesh.getEleVtx(n);
    std::vector<MyPoint> EP=mesh.getPnt(EV);
    std::vector<MyPoint> GP=tmpEle.getGaussPnt();

    std::vector<double> GW=tmpEle.getGaussWeight();
    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
    std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
    std::vector<std::vector<double> > basis_value;
    basisValue(q_point, EP, basis_value);
	std::vector<std::vector<std::vector<double> > > basis_grad;
    basisGrad(q_point, EP, basis_grad);
		
//    std::vector<double> u_coef;
//    Read_u_Coef(n, u_coef);
//    Read_q_Coef(n, q_coef);
    
    double h=EP[2][0]-EP[0][0];
    double h3=power(EP[2][0]-EP[0][0],1);
	
	for(int j=0; j<q_point.size(); j++){
		double Jxw = GW[j]*jacobian[j];
		
		double phi_val=0.0;
		for(int k=0; k<basis_value.size(); k++)
		{
			phi_val += u_h(n*ele_dof+k)*basis_value[k][j];
		}
		
		double M1_val=exp(-q_i[0]*phi_val);
		double M2_val=exp(-q_i[1]*phi_val);
		
		for(int k = 0; k < basis_value.size(); k++){
			for(int l = 0; l < basis_value.size(); l++){
//                Ele_Rhs[0][k] += Jxw*c1_h(n*ele_dof+l)*basis_value[l][j]*basis_value[k][j];
//				Ele_Rhs[1][k] += Jxw*c2_h(n*ele_dof+l)*basis_value[l][j]*basis_value[k][j];
				
				Ele_Rhs[0][k] -= Jxw*M1_val*g1_h(n*ele_dof+l)*(basis_grad[l][j][0]*basis_grad[k][j][0]+basis_grad[l][j][1]*basis_grad[k][j][1]);
				Ele_Rhs[1][k] -= Jxw*M2_val*g2_h(n*ele_dof+l)*(basis_grad[l][j][0]*basis_grad[k][j][0]+basis_grad[l][j][1]*basis_grad[k][j][1]);
			}
			
			Ele_Rhs[0][k] += Jxw*( _f1_(t,q_point[j]) )*basis_value[k][j];
			Ele_Rhs[1][k] += Jxw*( _f2_(t,q_point[j]) )*basis_value[k][j];
			
	  	}
	}
	
   	std::vector<double> ec(2);
   	ec[0]=(EP[0][0]+EP[1][0]+EP[2][0]+EP[3][0])/4.0;
   	ec[1]=(EP[0][1]+EP[1][1]+EP[2][1]+EP[3][1])/4.0;
    //     /*
    std::vector<int> eleEdg=mesh.getEleEdg(n);
    for(int i=0; i<eleEdg.size(); i++)
	{
		std::vector<int> edgVtx=mesh.getEdgVtx(eleEdg[i]);
// /*
		std::vector<MyPoint> edgi;
   		edgi=mesh.getPnt(edgVtx);
//		std::cout<<edgi[1][0]<<"   "<<edgi[1][1]<<std::endl;
// */
//  
/*
        std::vector<MyPoint> edgi(2);
   		edgi[0]=mesh.getPnt(edgVtx[0]);
   		edgi[1]=mesh.getPnt(edgVtx[1]);
//		std::cout<<edgi[1][0]<<"   "<<edgi[1][1]<<std::endl;
// */
        MyPoint midEdg=midpoint(edgi[0],edgi[1]);
		
//		std::cout<<ec[0]<<"   "<<ec[1]<<std::endl;

		int EdgLocation=1;
        if(fabs(midEdg[0]-ec[0])>1.0e-08){//该边平行y轴
			if(midEdg[0]-ec[0]>1.0e-08)
				EdgLocation=1;//该边是单元的右边
			else
				EdgLocation=3;//该边是单元的左边
		}
		else{//该边平行x轴
			if(midEdg[1]-ec[1]>1.0e-08)
				EdgLocation=2;//该边是单元的上边
			else
				EdgLocation=4;//该边是单元的下边
		}
//		int ik=element_patch[n][i];


        if(edge_patch[eleEdg[i]].size()>1)        // if the edge has two neighbor.
		{

		int ik=edge_patch[eleEdg[i]][0];
		if(ik==n)
			ik=edge_patch[eleEdg[i]][1];
		
//		std::vector<double> u_coef0, q_coef0;
 //       Read_u_Coef(ik, u_coef0);
//        Read_q_Coef(ik, q_coef0);
// 
/*
        if(n==0 && i==0){
		    std::cout<<midEdg[0]<<"   "<<midEdg[1]<<std::endl;
			std::cout<<ec[0]<<"   "<<ec[1]<<std::endl;
			std::cout<< EdgLocation <<std::endl;
//			std::cout<<edgi[0][0]<<"   "<<edgi[0][1]<<"   "<<edgi[1][0]<<"   "<<edgi[1][1]<<std::endl;
		}
// */

	  	std::vector<int> EVk=mesh.getEleVtx(ik);
// /*
        std::vector<MyPoint> EPk;
        EPk=mesh.getPnt(EVk);
//		std::cout<<EPk[1][0]<<"   "<<EPk[1][1]<<std::endl;
// */
// 
/*
        std::vector<MyPoint> EPk(4);
        EPk[0]=mesh.getPnt(EVk[0]);
        EPk[1]=mesh.getPnt(EVk[1]);
        EPk[2]=mesh.getPnt(EVk[2]);
        EPk[3]=mesh.getPnt(EVk[3]);
		std::cout<<EPk[1][0]<<"   "<<EPk[1][1]<<std::endl;
// */		
		
		switch(EdgLocation){
			case 1: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();
	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointEE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
             		q_pointEE[ie][1]=GPy[ie];
          		}
///////////////////////////////////////////////////////////////////
// test close 
/*
          		for(int ie=0; ie<GPE.size(); ie++){
             		if(fabs(q_pointE[ie][0]-EP[1][0])<1.0e-07 || fabs(q_pointE[ie][0]-EP[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"Error1!!!"<<std::endl;
						getchar();
					}
             		if(fabs(q_pointEE[ie][0]-EPk[1][0])<1.0e-07 || fabs(q_pointEE[ie][0]-EPk[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"error1!!!"<<std::endl;
						getchar();
					}
          		}
// */ 
///////////////////////////////////////////////////////////////////
	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);
				std::vector<std::vector<std::vector<double> > > basis_grad2I;
    	  		basis2Grad(q_pointE, EP, basis_grad2I);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
				std::vector<std::vector<std::vector<double> > > basis_grad2E;
    	  		basis2Grad(q_pointEE, EPk, basis_grad2E);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double u_n_valE= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						u_n_valE += u_h(n*ele_dof+j)*basis_valueI[j][ie];
					}
				
					double u_n_valEE= 0.0;
					for(int j=0; j<basis_valueE.size(); j++)
					{
						u_n_valEE += u_h(ik*ele_dof+j)*basis_valueE[j][ie];
					}
				
					double phi_val = (u_n_valE+u_n_valEE)/2.0;
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
							Ele_Rhs[0][m]-=Jxw*beta_0/h*M1_val*g1_h(n*ele_dof+j)*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_Rhs[0][m]-=-Jxw*beta_0/h*M1_val*g1_h(ik*ele_dof+j)*(basis_valueE[j][ie]*basis_valueI[m][ie]);																			
							
                  			Ele_Rhs[0][m]-=-Jxw/2*M1_val*g1_h(n*ele_dof+j)*(basis_gradI[j][ie][0]*basis_valueI[m][ie]
											+ basis_valueI[j][ie]*basis_gradI[m][ie][0]);

                  		    Ele_Rhs[0][m]-=Jxw/2*M1_val*g1_h(ik*ele_dof+j)*(-basis_gradE[j][ie][0]*basis_valueI[m][ie]
											+ basis_valueE[j][ie]*basis_gradI[m][ie][0]);
											
							Ele_Rhs[0][m]-=Jxw*beta1*h*M1_val*g1_h(n*ele_dof+j)*(basis_grad2I[j][ie][0]*basis_valueI[m][ie]);

                  		    Ele_Rhs[0][m]-=-Jxw*beta1*h*M1_val*g1_h(ik*ele_dof+j)*(basis_grad2E[j][ie][0]*basis_valueI[m][ie]);	

							Ele_Rhs[1][m]-=Jxw*beta_0/h*M2_val*g2_h(n*ele_dof+j)*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_Rhs[1][m]-=-Jxw*beta_0/h*M2_val*g2_h(ik*ele_dof+j)*(basis_valueE[j][ie]*basis_valueI[m][ie]);																			
							
                  			Ele_Rhs[1][m]-=-Jxw/2*M2_val*g2_h(n*ele_dof+j)*(basis_gradI[j][ie][0]*basis_valueI[m][ie]
											+ basis_valueI[j][ie]*basis_gradI[m][ie][0]);

                  		    Ele_Rhs[1][m]-=Jxw/2*M2_val*g2_h(ik*ele_dof+j)*(-basis_gradE[j][ie][0]*basis_valueI[m][ie]
											+ basis_valueE[j][ie]*basis_gradI[m][ie][0]);
											
							Ele_Rhs[1][m]-=Jxw*beta1*h*M2_val*g2_h(n*ele_dof+j)*(basis_grad2I[j][ie][0]*basis_valueI[m][ie]);

                  		    Ele_Rhs[1][m]-=-Jxw*beta1*h*M2_val*g2_h(ik*ele_dof+j)*(basis_grad2E[j][ie][0]*basis_valueI[m][ie]);							
						}
					}
          		}
				break;	
			}
			case 2: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointEE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
             		q_pointEE[ie][0]=GPx[ie];
          		}
////////////////////////////////////////////////////////////////////
// test close 
/*
          		for(int ie=0; ie<GPE.size(); ie++){
             		if(fabs(q_pointE[ie][1]-EP[1][1])<1.0e-07 || fabs(q_pointE[ie][1]-EP[3][1])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"Error2!!!"<<std::endl;
						getchar();
					}
             		if(fabs(q_pointEE[ie][1]-EPk[1][1])<1.0e-07 || fabs(q_pointEE[ie][1]-EPk[3][1])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"error2!!!"<<std::endl;
						getchar();
					}
          		}
// */ 
////////////////////////////////////////////////////////////////////
	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);
				std::vector<std::vector<std::vector<double> > > basis_grad2I;
    	  		basis2Grad(q_pointE, EP, basis_grad2I);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
				std::vector<std::vector<std::vector<double> > > basis_grad2E;
    	  		basis2Grad(q_pointEE, EPk, basis_grad2E);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double u_n_valE= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						u_n_valE += u_h(n*ele_dof+j)*basis_valueI[j][ie];
					}
				
					double u_n_valEE= 0.0;
					for(int j=0; j<basis_valueE.size(); j++)
					{
						u_n_valEE += u_h(ik*ele_dof+j)*basis_valueE[j][ie];
					}
				
					double phi_val = (u_n_valE+u_n_valEE)/2.0;
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
							Ele_Rhs[0][m]-=Jxw*beta_0/h*M1_val*g1_h(n*ele_dof+j)*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_Rhs[0][m]+=Jxw*beta_0/h*M1_val*g1_h(ik*ele_dof+j)*(basis_valueE[j][ie]*basis_valueI[m][ie]);								
							
                  			Ele_Rhs[0][m]+=Jxw/2*M1_val*g1_h(n*ele_dof+j)*(basis_gradI[j][ie][1]*basis_valueI[m][ie]
											+basis_valueI[j][ie]*basis_gradI[m][ie][1]);

                  		    Ele_Rhs[0][m]-=Jxw/2*M1_val*g1_h(ik*ele_dof+j)*(-basis_gradE[j][ie][1]*basis_valueI[m][ie]
											+basis_valueE[j][ie]*basis_gradI[m][ie][1]);
											
							Ele_Rhs[0][m]-=Jxw*beta1*h*M1_val*g1_h(n*ele_dof+j)*(basis_grad2I[j][ie][2]*basis_valueI[m][ie]);
											
							Ele_Rhs[0][m]-=-Jxw*beta1*h*M1_val*g1_h(ik*ele_dof+j)*(basis_grad2E[j][ie][2]*basis_valueI[m][ie]);	

							Ele_Rhs[1][m]-=Jxw*beta_0/h*M2_val*g2_h(n*ele_dof+j)*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_Rhs[1][m]+=Jxw*beta_0/h*M2_val*g2_h(ik*ele_dof+j)*(basis_valueE[j][ie]*basis_valueI[m][ie]);								
							
                  			Ele_Rhs[1][m]+=Jxw/2*M2_val*g2_h(n*ele_dof+j)*(basis_gradI[j][ie][1]*basis_valueI[m][ie]
											+basis_valueI[j][ie]*basis_gradI[m][ie][1]);

                  		    Ele_Rhs[1][m]-=Jxw/2*M2_val*g2_h(ik*ele_dof+j)*(-basis_gradE[j][ie][1]*basis_valueI[m][ie]
											+basis_valueE[j][ie]*basis_gradI[m][ie][1]);
											
							Ele_Rhs[1][m]-=Jxw*beta1*h*M2_val*g2_h(n*ele_dof+j)*(basis_grad2I[j][ie][2]*basis_valueI[m][ie]);
											
							Ele_Rhs[1][m]-=-Jxw*beta1*h*M2_val*g2_h(ik*ele_dof+j)*(basis_grad2E[j][ie][2]*basis_valueI[m][ie]);
							
						}
					}
          		}
				break;	
			}
			case 3: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointEE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
             		q_pointEE[ie][1]=GPy[ie];
          		}
////////////////////////////////////////////////////////////////////
// test close 
/*
          		for(int ie=0; ie<GPE.size(); ie++){
             		if(fabs(q_pointE[ie][0]-EP[1][0])<1.0e-07 || fabs(q_pointE[ie][0]-EP[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"Error3!!!"<<std::endl;
						getchar();
					}
             		if(fabs(q_pointEE[ie][0]-EPk[1][0])<1.0e-07 || fabs(q_pointEE[ie][0]-EPk[3][0])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"error3!!!"<<std::endl;
						getchar();
					}
          		}
// */
////////////////////////////////////////////////////////////////////
	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);
				std::vector<std::vector<std::vector<double> > > basis_grad2I;
    	  		basis2Grad(q_pointE, EP, basis_grad2I);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
				std::vector<std::vector<std::vector<double> > > basis_grad2E;
    	  		basis2Grad(q_pointEE, EPk, basis_grad2E);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double u_n_valE= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						u_n_valE += u_h(n*ele_dof+j)*basis_valueI[j][ie];
					}
				
					double u_n_valEE= 0.0;
					for(int j=0; j<basis_valueE.size(); j++)
					{
						u_n_valEE += u_h(ik*ele_dof+j)*basis_valueE[j][ie];
					}
				
					double phi_val = (u_n_valE+u_n_valEE)/2.0;
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
							Ele_Rhs[0][m]-=Jxw*beta_0/h*M1_val*g1_h(n*ele_dof+j)*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_Rhs[0][m]-=-Jxw*beta_0/h*M1_val*g1_h(ik*ele_dof+j)*(basis_valueE[j][ie]*basis_valueI[m][ie]);															
							
                  			Ele_Rhs[0][m]-=Jxw/2*M1_val*g1_h(n*ele_dof+j)*(basis_gradI[j][ie][0]*basis_valueI[m][ie]
											+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
											
                  		    Ele_Rhs[0][m]-=Jxw/2*M1_val*g1_h(ik*ele_dof+j)*(basis_gradE[j][ie][0]*basis_valueI[m][ie]
											-basis_valueE[j][ie]*basis_gradI[m][ie][0]);
							
							Ele_Rhs[0][m]-=Jxw*beta1*h*M1_val*g1_h(n*ele_dof+j)*(basis_grad2I[j][ie][0]*basis_valueI[m][ie]);

                  		    Ele_Rhs[0][m]-=-Jxw*beta1*h*M1_val*g1_h(ik*ele_dof+j)*(basis_grad2E[j][ie][0]*basis_valueI[m][ie]);	

							Ele_Rhs[1][m]-=Jxw*beta_0/h*M2_val*g2_h(n*ele_dof+j)*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_Rhs[1][m]-=-Jxw*beta_0/h*M2_val*g2_h(ik*ele_dof+j)*(basis_valueE[j][ie]*basis_valueI[m][ie]);															
							
                  			Ele_Rhs[1][m]-=Jxw/2*M2_val*g2_h(n*ele_dof+j)*(basis_gradI[j][ie][0]*basis_valueI[m][ie]
											+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
											
                  		    Ele_Rhs[1][m]-=Jxw/2*M2_val*g2_h(ik*ele_dof+j)*(basis_gradE[j][ie][0]*basis_valueI[m][ie]
											-basis_valueE[j][ie]*basis_gradI[m][ie][0]);
							
							Ele_Rhs[1][m]-=Jxw*beta1*h*M2_val*g2_h(n*ele_dof+j)*(basis_grad2I[j][ie][0]*basis_valueI[m][ie]);

                  		    Ele_Rhs[1][m]-=-Jxw*beta1*h*M2_val*g2_h(ik*ele_dof+j)*(basis_grad2E[j][ie][0]*basis_valueI[m][ie]);	
							
						}
					}
          		}
				break;	
			}
			case 4: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		std::vector<MyPoint> q_pointEE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointEE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
             		q_pointEE[ie][0]=GPx[ie];
          		}
////////////////////////////////////////////////////////////////////
// test close 
/*
          		for(int ie=0; ie<GPE.size(); ie++){
             		if(fabs(q_pointE[ie][1]-EP[1][1])<1.0e-07 || fabs(q_pointE[ie][1]-EP[3][1])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"Error4!!!"<<std::endl;
						getchar();
					}
             		if(fabs(q_pointEE[ie][1]-EPk[1][1])<1.0e-07 || fabs(q_pointEE[ie][1]-EPk[3][1])<1.0e-07 ){
						int ooo=0;
					}
					else{
						std::cout<<"error4!!!"<<std::endl;
						getchar();
					}
          		}
//  */
////////////////////////////////////////////////////////////////////
	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);
				std::vector<std::vector<std::vector<double> > > basis_grad2I;
    	  		basis2Grad(q_pointE, EP, basis_grad2I);

	  			std::vector<std::vector<double> > basis_valueE;
    	  		basisValue(q_pointEE, EPk, basis_valueE);
	  			std::vector<std::vector<std::vector<double> > > basis_gradE;
    	  		basisGrad(q_pointEE, EPk, basis_gradE);
				std::vector<std::vector<std::vector<double> > > basis_grad2E;
    	  		basis2Grad(q_pointEE, EPk, basis_grad2E);
	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double u_n_valE= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						u_n_valE += u_h(n*ele_dof+j)*basis_valueI[j][ie];
					}
				
					double u_n_valEE= 0.0;
					for(int j=0; j<basis_valueE.size(); j++)
					{
						u_n_valEE += u_h(ik*ele_dof+j)*basis_valueE[j][ie];
					}
				
					double phi_val = (u_n_valE+u_n_valEE)/2.0;
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
							Ele_Rhs[0][m]-=Jxw*beta_0/h*M1_val*g1_h(n*ele_dof+j)*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_Rhs[0][m]-=-Jxw*beta_0/h*M1_val*g1_h(ik*ele_dof+j)*(basis_valueE[j][ie]*basis_valueI[m][ie]);		
							
                  			Ele_Rhs[0][m]-=Jxw/2*M1_val*g1_h(n*ele_dof+j)*(basis_gradI[j][ie][1]*basis_valueI[m][ie]
											+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
                  		    Ele_Rhs[0][m]-=Jxw/2*M1_val*g1_h(ik*ele_dof+j)*(basis_gradE[j][ie][1]*basis_valueI[m][ie]
											-basis_valueE[j][ie]*basis_gradI[m][ie][1]);
							
							Ele_Rhs[0][m]-=Jxw*beta1*h*M1_val*g1_h(n*ele_dof+j)*(basis_grad2I[j][ie][2]*basis_valueI[m][ie]);
											
							Ele_Rhs[0][m]-=-Jxw*beta1*h*M1_val*g1_h(ik*ele_dof+j)*(basis_grad2E[j][ie][2]*basis_valueI[m][ie]);		

							Ele_Rhs[1][m]-=Jxw*beta_0/h*M2_val*g2_h(n*ele_dof+j)*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
							Ele_Rhs[1][m]-=-Jxw*beta_0/h*M2_val*g2_h(ik*ele_dof+j)*(basis_valueE[j][ie]*basis_valueI[m][ie]);		
							
                  			Ele_Rhs[1][m]-=Jxw/2*M2_val*g2_h(n*ele_dof+j)*(basis_gradI[j][ie][1]*basis_valueI[m][ie]
											+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
                  		    Ele_Rhs[1][m]-=Jxw/2*M2_val*g2_h(ik*ele_dof+j)*(basis_gradE[j][ie][1]*basis_valueI[m][ie]
											-basis_valueE[j][ie]*basis_gradI[m][ie][1]);
							
							Ele_Rhs[1][m]-=Jxw*beta1*h*M2_val*g2_h(n*ele_dof+j)*(basis_grad2I[j][ie][2]*basis_valueI[m][ie]);
											
							Ele_Rhs[1][m]-=-Jxw*beta1*h*M2_val*g2_h(ik*ele_dof+j)*(basis_grad2E[j][ie][2]*basis_valueI[m][ie]);	
							
						}
					}
          		}
				break;
			}
			default:
			{
				std::cout<<"Some problem in location the edg, please chech carefully!"<<std::endl;
				break;	
			}
		}
		
		
	    }// end of the two neighbor case.
		else
		{
			

			
		switch(EdgLocation){
			case 1: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();
	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double phi_val= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						phi_val += u_h(n*ele_dof+j)*basis_valueI[j][ie];
					}
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);								
							
//                  			Ele_Rhs[0][m]+=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);

//							Ele_Rhs[1][m]+=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);

//							Ele_Rhs[0][m]-=Jxw*beta0x/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0x/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                            Ele_Rhs[1][m]+=Jxw*beta0x/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0x/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);							
						}
//						std::cout<<"_u_grad_(t,q_pointE[ie])[0]= "<<_u_grad_(t,q_pointE[ie])[0]<<", _u_grad_(t,q_pointE[ie])[1]= "<<_u_grad_(t,q_pointE[ie])[1]<<std::endl;
//						Ele_Rhs[0][m]+=Jxw*beta_0/h*(_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);

						Ele_Rhs[0][m]+=Jxw*zeroflux*M1_val*(exp(q_i[0]*_u_(t,q_pointE[ie]))*( _c1_grad_(t,q_pointE[ie])[0] + _c1_(t,q_pointE[ie])*q_i[0]*_u_grad_(t,q_pointE[ie])[0] )*basis_valueI[m][ie]);
						Ele_Rhs[1][m]+=Jxw*zeroflux*M2_val*(exp(q_i[1]*_u_(t,q_pointE[ie]))*( _c2_grad_(t,q_pointE[ie])[0] + _c2_(t,q_pointE[ie])*q_i[1]*_u_grad_(t,q_pointE[ie])[0] )*basis_valueI[m][ie]);

//						Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[0][m]-=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][0]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
						
//						Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]-=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][0]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
									
//						Ele_Rhs[1][m]+=Jxw*theta*(_u_grad_(t,q_pointE[ie])[0]*basis_valueI[m][ie]- _u_(t,q_pointE[ie])*basis_gradI[m][ie][0])
//						               +Jxw*antitheta*(_u_grad_(t-dt,q_pointE[ie])[0]*basis_valueI[m][ie]-_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
			
					}
          		}
				break;	
			}
			case 2: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double phi_val= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						phi_val += u_h(n*ele_dof+j)*basis_valueI[j][ie];
					}
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
//                  			Ele_Rhs[0][m]+=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_Rhs[1][m]+=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);

//							Ele_Rhs[0][m]-=Jxw*beta0y/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0y/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                           Ele_Rhs[1][m]+=Jxw*beta0y/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0y/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
				
						Ele_Rhs[0][m]+=Jxw*zeroflux*M1_val*(exp(q_i[0]*_u_(t,q_pointE[ie]))*( _c1_grad_(t,q_pointE[ie])[1] + _c1_(t,q_pointE[ie])*q_i[0]*_u_grad_(t,q_pointE[ie])[1] )*basis_valueI[m][ie]);
						Ele_Rhs[1][m]+=Jxw*zeroflux*M2_val*(exp(q_i[1]*_u_(t,q_pointE[ie]))*( _c2_grad_(t,q_pointE[ie])[1] + _c2_(t,q_pointE[ie])*q_i[1]*_u_grad_(t,q_pointE[ie])[1] )*basis_valueI[m][ie]);
						
					    
//					    Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//					    Ele_Rhs[0][m]-=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][1]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
					    
//					    Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//					    Ele_Rhs[1][m]-=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][1]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);


//						Ele_Rhs[1][m]+=Jxw*theta*(_u_grad_(t,q_pointE[ie])[1]*basis_valueI[m][ie]-_u_(t,q_pointE[ie])*basis_gradI[m][ie][1])
//						               +Jxw*antitheta*(_u_grad_(t-dt,q_pointE[ie])[1]*basis_valueI[m][ie]-_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
					
					}
          		}
				break;	
			}
			case 3: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double phi_val= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						phi_val += u_h(n*ele_dof+j)*basis_valueI[j][ie];
					}
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
//                  			Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
											
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
//							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]-=Jxw*beta0x/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0x/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

 //                           Ele_Rhs[1][m]+=Jxw*beta0x/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0x/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
						
						Ele_Rhs[0][m]-=Jxw*zeroflux*M1_val*(exp(q_i[0]*_u_(t,q_pointE[ie]))*( _c1_grad_(t,q_pointE[ie])[0] + _c1_(t,q_pointE[ie])*q_i[0]*_u_grad_(t,q_pointE[ie])[0] )*basis_valueI[m][ie]);
						Ele_Rhs[1][m]-=Jxw*zeroflux*M2_val*(exp(q_i[1]*_u_(t,q_pointE[ie]))*( _c2_grad_(t,q_pointE[ie])[0] + _c2_(t,q_pointE[ie])*q_i[1]*_u_grad_(t,q_pointE[ie])[0] )*basis_valueI[m][ie]);
						
//						Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[0][m]+=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][0]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
						
//						Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]+=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][0]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);


//						Ele_Rhs[1][m]+=Jxw*theta*(-_u_grad_(t,q_pointE[ie])[0]*basis_valueI[m][ie]+_u_(t,q_pointE[ie])*basis_gradI[m][ie][0])
//						               +Jxw*antitheta*(-_u_grad_(t-dt,q_pointE[ie])[0]*basis_valueI[m][ie]+_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
					}
          		}
				break;	
			}
			case 4: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double phi_val= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						phi_val += u_h(n*ele_dof+j)*basis_valueI[j][ie];
					}
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
//                  			Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]-=Jxw*beta0y/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0y/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                            Ele_Rhs[1][m]+=Jxw*beta0y/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0y/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
//						std::cout<<"t= "<<t<<", dt= "<<dt<<", q_pointE[ie]= "<<q_pointE[ie][0]<<" , q_pointE[ie][1]= "<<q_pointE[ie][1]<<",  _u_(t,q_pointE(ie))= "<<_u_(t,q_pointE[ie])<<std::endl;

						Ele_Rhs[0][m]-=Jxw*zeroflux*M1_val*(exp(q_i[0]*_u_(t,q_pointE[ie]))*( _c1_grad_(t,q_pointE[ie])[1] + _c1_(t,q_pointE[ie])*q_i[0]*_u_grad_(t,q_pointE[ie])[1] )*basis_valueI[m][ie]);
						Ele_Rhs[1][m]-=Jxw*zeroflux*M2_val*(exp(q_i[1]*_u_(t,q_pointE[ie]))*( _c2_grad_(t,q_pointE[ie])[1] + _c2_(t,q_pointE[ie])*q_i[1]*_u_grad_(t,q_pointE[ie])[1] )*basis_valueI[m][ie]);
						
//						Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[0][m]+=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][1]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
//						
//						Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]+=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][1]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
						
//						Ele_Rhs[1][m]+=Jxw*theta*(-_u_grad_(t,q_pointE[ie])[1]*basis_valueI[m][ie]+_u_(t,q_pointE[ie])*basis_gradI[m][ie][1])
//						               +Jxw*antitheta*(-_u_grad_(t-dt,q_pointE[ie])[1]*basis_valueI[m][ie]+_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
						               
					}
          		}
				break;
			}
			default:
			{
				std::cout<<"Some problem in location the edg, please chech carefully!"<<std::endl;
				break;	
			}
		}
				
		} // end of the one neighbor case.
		
			
		
		
	}	
	
}

void BH2D::Solve()
{
// 
/*
    if(t==dt){
//	    writeData_Mat();
//	    writeData_RHS();
	    fpwriteData_Mat();
	    fpwriteData_RHS();
	}
// */

// 
/*
	SparseILU<double> preconditioner;
	preconditioner.initialize(StiffMatrix, SparseILU<double>::AdditionalData());

	SolverControl solver_control(10000, 1e-10);
	PrimitiveVectorMemory<double> vector_memory;
	SolverBicgstab<Vector<double> > solver(solver_control);
	solver.solve(StiffMatrix, u_h, RHS, preconditioner);
	
// */

// /*
	SparseILU<double> preconditioner;
	preconditioner.initialize(StiffMatrix, SparseILU<double>::AdditionalData());

	SolverControl solver_control(100000, 1e-12);
	PrimitiveVectorMemory<double> vector_memory;
	SolverFGMRES<Vector<double> > solver(solver_control);      // Can change between GMRES and FGMRES
	solver.solve(StiffMatrix, u_h, RHS, preconditioner);
	
// */

// 
/*
	SparseMIC<double> preconditioner;
	preconditioner.initialize(StiffMatrix, SparseMIC<double>::AdditionalData());

	SolverControl solver_control(10000, 1e-10);
	PrimitiveVectorMemory<double> vector_memory;
	SolverBicgstab<Vector<double> > solver(solver_control);
	solver.solve(StiffMatrix, u_h, RHS, preconditioner);
	
// */

// 
/*
	SparseMIC<double> preconditioner;
	preconditioner.initialize(StiffMatrix, SparseMIC<double>::AdditionalData());

	SolverControl solver_control(10000, 1e-10);
	PrimitiveVectorMemory<double> vector_memory;
	SolverFGMRES<Vector<double> > solver(solver_control);      // Can change between GMRES and FGMRES
	solver.solve(StiffMatrix, u_h, RHS, preconditioner);
	
// */



// 
/*
	SparseMIC<double> preconditioner;
	preconditioner.initialize(StiffMatrix, SparseMIC<double>::AdditionalData());

	SolverControl solver_control(10000, 1e-10);
	PrimitiveVectorMemory<double> vector_memory;
	SolverCG<Vector<double> > solver(solver_control);
	solver.solve(StiffMatrix, u_h, RHS, preconditioner);
	
// */

// 
/*	
  SolverControl           solver_control (10000, 1e-10);
  SolverBicgstab<>        bicgstab (solver_control);
  PreconditionJacobi<>    preconditioner;
  preconditioner.initialize(StiffMatrix, 2.0);
  bicgstab.solve(StiffMatrix, u_h, RHS, preconditioner);
// */

// 
/*
  SolverControl           solver_control (10000, 1e-12);
  SolverFGMRES<>          solver (solver_control, SolverFGMRES<>::AdditionalData());
  PreconditionSSOR<>      preconditioner;
  preconditioner.initialize(StiffMatrix, 1.0);
  solver.solve (StiffMatrix, u_h, RHS, preconditioner);
// */

// 
/*
  SolverControl           solver_control (10000, 1e-10);
  SolverFGMRES<>           solver (solver_control);
  PreconditionSSOR<>      preconditioner;
  preconditioner.initialize(StiffMatrix, 1.0);
  solver.solve (StiffMatrix, u_h, RHS, preconditioner);
// */

// 
/*
  SolverControl           solver_control (10000, 1e-10);
  SolverGMRES<>           solver (solver_control);
  PreconditionSSOR<>      preconditioner;
  preconditioner.initialize(StiffMatrix, 1.2);
  solver.solve (StiffMatrix, u_h, RHS, preconditioner);
// */

// 
/*
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              solver (solver_control);
  PreconditionSSOR<>      preconditioner;
  preconditioner.initialize(StiffMatrix, 1.2);
  solver.solve(StiffMatrix, u_h, RHS, preconditioner);
// */

// 
/*
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              solver (solver_control);
  solver.solve (StiffMatrix, u_h, RHS,  PreconditionIdentity());
// */

////////////////////////////////////////////////////////////////////////////////////////////////////

// 
/*
  SolverControl           solver_control (10000, 1e-10);
  SolverBicgstab<>        bicgstab (solver_control);
  PreconditionJacobi<>    preconditioner;
  preconditioner.initialize(StiffMatrix, 1.0);
  bicgstab.solve (StiffMatrix, u_h, RHS, preconditioner);
// */

// 
/*
  SolverControl           solver_control (10000, 1e-10);
  SolverGMRES<>           solver (solver_control);
  PreconditionSSOR<>      preconditioner;
  preconditioner.initialize(StiffMatrix, 1.2);
  solver.solve (StiffMatrix, u_h, RHS, preconditioner);
// */

// 
/*
       PreconditionSelector<SparseMatrix<double>, Vector<double> >  preconditioning("jacobi", .6);
       preconditioning.use_matrix(StiffMatrix);

       SolverControl solver_control(10000, 1e-10);
       PrimitiveVectorMemory<> vector_memory;
       SolverGMRES<Vector<double> > solver(solver_control);
       solver.solve(StiffMatrix, u_h, RHS, preconditioning);
// */

// 
/*
       SparseILU<double> preconditioner;
       preconditioner.initialize(StiffMatrix, SparseILU<double>::AdditionalData());
//       SparseILU<double> ilu(sparsity_pattern);
//       ilu.decompose(StiffMatrix);

       SolverControl solver_control(10000, 1e-10);
       PrimitiveVectorMemory<double> vector_memory;
       SolverBicgstab<Vector<double> > solver(solver_control);
       solver.solve(StiffMatrix, u_h, RHS, preconditioner);
// */	


// not use
/*
//  PreconditionSelector<SparseMatrix<double>, Vector<double> >  preconditioning("jacobi", 0.6);
  PreconditionSelector<SparseMatrix<double>, Vector<double> >  preconditioning("jacobi", Jacobiintital);
  preconditioning.use_matrix(StiffMatrix);

  SolverControl solver_control(10000,TOL);
  PrimitiveVectorMemory<> vector_memory;
  SolverGMRES<Vector<double> > solver(solver_control, vector_memory,10000);
  solver.solve(StiffMatrix, u_h, RHS, preconditioning);
// */

// 
/*
     for(int i=0; i<StiffMatrix.m(); i++){
        for(int j=0; j<StiffMatrix.n(); j++){
           if(StiffMatrix.el(i,j)!=0)
              StiffMatrix.set(i,j,0.0);
        }
     }

     for(int i=0; i<RHS.size(); i++){
        RHS(i)=0;
     }
// */

/*
  SparseILU<double> ilu(sparsity_pattern);
  ilu.decompose(StiffMatrix);

  SolverControl solver_control(1000, TOL);
  PrimitiveVectorMemory<> vector_memory;
  SolverBicgstab<> solver(solver_control, vector_memory);
  solver.solve(StiffMatrix, u_h, RHS,
	       PreconditionUseMatrix<SparseILU<double>,Vector<double> >
	       (ilu, &SparseILU<double>::apply_decomposition<double>));
*/ 
// 
/*
   std::ofstream outdata;
   outdata.open("RHS.txt");
   outdata<<2*v_size<<"\n";
   for(int i=0; i<n_ele; i++){
      outdata<<RHS(9*i)<<"\t"<<RHS(9*i+1)<<"\t"<<RHS(9*i+2)<<"\t"<<RHS(9*i+3)<<"\n";
   }
   outdata.close();
// */	
	
	
}

double BH2D::ComputeL2Error(Vector<double> &f, double (*g)(const double *))
{
//	int ele_dof=9;
   double err=0.0;
   int n_ele=mesh.n_element();
   for(int j=0; j<n_ele; j++){
      std::vector<int> NV=mesh.getEleVtx(j);
      std::vector<MyPoint> EP(4);
      EP[0]=mesh.getPnt(NV[0]);
      EP[1]=mesh.getPnt(NV[1]);
      EP[2]=mesh.getPnt(NV[2]);
      EP[3]=mesh.getPnt(NV[3]);
      std::vector<MyPoint> GP=tmpEle.getGaussPnt();

	  std::vector<double> GW=tmpEle.getGaussWeight();
	  std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	  std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	  std::vector<std::vector<double> > basis_value;
      basisValue(q_point, EP, basis_value);
	  
	  for(int i = 0; i < q_point.size(); i++){
	      double Jxw = GW[i]*jacobian[i];
	      double f_value = g(q_point[i]);
              double u_h_val=0.0;
              for(int m=0; m<basis_value.size(); m++){
                  u_h_val+=f(ele_dof*j+m)*basis_value[m][i];
              }
              double df_val=f_value-u_h_val;
              err+=Jxw*df_val*df_val;
	  }
   }
   err=sqrt(err);
   return err;
}

double BH2D::ComputeL1Error(Vector<double> &f, double (*g)(const double *))
{
//	int ele_dof=9;
   double err=0.0;
   int n_ele=mesh.n_element();
   for(int j=0; j<n_ele; j++){
      std::vector<int> NV=mesh.getEleVtx(j);
      std::vector<MyPoint> EP(4);
      EP[0]=mesh.getPnt(NV[0]);
      EP[1]=mesh.getPnt(NV[1]);
      EP[2]=mesh.getPnt(NV[2]);
      EP[3]=mesh.getPnt(NV[3]);
      std::vector<MyPoint> GP=tmpEle.getGaussPnt();

	  std::vector<double> GW=tmpEle.getGaussWeight();
	  std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	  std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	  std::vector<std::vector<double> > basis_value;
      basisValue(q_point, EP, basis_value);
	  
	  for(int i = 0; i < q_point.size(); i++){
	      double Jxw = GW[i]*jacobian[i];
	      double f_value = g(q_point[i]);
              double u_h_val=0.0;
              for(int m=0; m<basis_value.size(); m++){
                  u_h_val+=f(ele_dof*j+m)*basis_value[m][i];
              }
              double df_val=fabs(f_value-u_h_val);
              err+=Jxw*df_val;
	  }
   }
   return err;
}

double BH2D::ComputeH1Error(Vector<double> &f, std::vector<double> (*g)(const double *))
{
//	int ele_dof=9;
    double err=0.0;
    int n_ele=mesh.n_element();
    for(int j=0; j<n_ele; j++){
        std::vector<int> NV=mesh.getEleVtx(j);
        std::vector<MyPoint> EP(4);
        EP[0]=mesh.getPnt(NV[0]);
        EP[1]=mesh.getPnt(NV[1]);
        EP[2]=mesh.getPnt(NV[2]);
        EP[3]=mesh.getPnt(NV[3]);
        std::vector<MyPoint> GP=tmpEle.getGaussPnt();
	    std::vector<double> GW=tmpEle.getGaussWeight();
  	    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	    std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	    std::vector<std::vector<std::vector<double> > > basis_grad;
    	basisGrad(q_point, EP, basis_grad);	  
	    for(int i = 0; i < q_point.size(); i++){
	        double Jxw = GW[i]*jacobian[i];
	        std::vector<double> f_value = g(q_point[i]);
            double u_h_grad1=0.0;
            double u_h_grad2=0.0;
            for(int m=0; m<basis_grad.size(); m++){
                u_h_grad1+=f(ele_dof*j+m)*basis_grad[m][i][0];
                u_h_grad2+=f(ele_dof*j+m)*basis_grad[m][i][1];
            }
            double df_val1=f_value[0]-u_h_grad1;
            double df_val2=f_value[1]-u_h_grad2;
            err+=Jxw*df_val1*df_val1;
            err+=Jxw*df_val2*df_val2;
	    }
    }
    err=sqrt(err);
    return err;
}

double BH2D::ComputeLINFError(Vector<double> &f, double (*g)(const double *))
{
//	int ele_dof=9;
    double err=0.0;
    int n_ele=mesh.n_element();
    for(int j=0; j<n_ele; j++){
        std::vector<int> NV=mesh.getEleVtx(j);
        std::vector<MyPoint> EP=mesh.getPnt(NV);
        std::vector<MyPoint> GP=tmpEle.getGaussPnt();
        std::vector<double> GW=tmpEle.getGaussWeight();
	    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
     	std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	    std::vector<std::vector<double> > basis_value;
        basisValue(q_point, EP, basis_value);
	    for(int i = 0; i < q_point.size(); i++){
	        double f_value = g(q_point[i]);
	        double u_h_val=0.0;
	        for(int m=0; m<basis_value.size(); m++){
                  u_h_val+=f(ele_dof*j+m)*basis_value[m][i];
            }
            double df_val=fabs(f_value-u_h_val);
            if(df_val>err)
                err=df_val;
	    }
	    
	    std::vector<std::vector<double> > basis_valueInt;
	    basisValue(EP, EP, basis_valueInt);
	    for(int i = 0; i < EP.size(); i++){
	        double f_value = g(EP[i]);
	        double u_h_val=0.0;
	        for(int m=0; m<basis_valueInt.size(); m++){
                  u_h_val+=f(ele_dof*j+m)*basis_valueInt[m][i];
            }
            double df_val=fabs(f_value-u_h_val);
            if(df_val>err)
            {
                err=df_val;
            }
	    }
    }
    return err;
}

double BH2D::ComputeSolAvg(Vector<double> &f)
{
//	int ele_dof=9;
	   double Avg=0.0;
	   int n_ele=mesh.n_element();
   for(int j=0; j<n_ele; j++){
      std::vector<int> NV=mesh.getEleVtx(j);
      std::vector<MyPoint> EP(4);
      EP[0]=mesh.getPnt(NV[0]);
      EP[1]=mesh.getPnt(NV[1]);
      EP[2]=mesh.getPnt(NV[2]);
      EP[3]=mesh.getPnt(NV[3]);
      std::vector<MyPoint> GP=tmpEle.getGaussPnt();

	  std::vector<double> GW=tmpEle.getGaussWeight();
	  std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	  std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	  std::vector<std::vector<double> > basis_value;
      basisValue(q_point, EP, basis_value);
	  
	  for(int i = 0; i < q_point.size(); i++){
	      double Jxw = GW[i]*jacobian[i];
//	      double f_value = g(t,q_point[i]);
//	      std::cout<<"f_value= "<<f_value<<std::endl;
              double u_h_val=0.0;
              for(int m=0; m<basis_value.size(); m++){
                  u_h_val+=f(ele_dof*j+m)*basis_value[m][i];
              }
              double df_val=u_h_val;
              Avg+=Jxw*df_val;
	  }
   }
//   Norm[0]=sqrt(Norm[0]);
//   Norm[1]=sqrt(Norm[1]);
   Avg = Avg/power(4.0*3.14159265358979323846,2);
   return Avg;
}

std::vector<double> BH2D::ComputeL2Norm(double t, Vector<double> &f, double (*g)(const double ,const double *))
{
//	int ele_dof=9;
   std::vector<double> Norm;
   Norm.resize(2,0.0);
   int n_ele=mesh.n_element();
   for(int j=0; j<n_ele; j++){
      std::vector<int> NV=mesh.getEleVtx(j);
      std::vector<MyPoint> EP(4);
      EP[0]=mesh.getPnt(NV[0]);
      EP[1]=mesh.getPnt(NV[1]);
      EP[2]=mesh.getPnt(NV[2]);
      EP[3]=mesh.getPnt(NV[3]);
      std::vector<MyPoint> GP=tmpEle.getGaussPnt();

	  std::vector<double> GW=tmpEle.getGaussWeight();
	  std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	  std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	  std::vector<std::vector<double> > basis_value;
      basisValue(q_point, EP, basis_value);
	  
	  for(int i = 0; i < q_point.size(); i++){
	      double Jxw = GW[i]*jacobian[i];
	      double f_value = g(t,q_point[i]);
//	      std::cout<<"f_value= "<<f_value<<std::endl;
              double u_h_val=0.0;
              for(int m=0; m<basis_value.size(); m++){
                  u_h_val+=f(ele_dof*j+m)*basis_value[m][i];
              }
              double df_val=u_h_val;
              Norm[0]+=Jxw*df_val;
              Norm[1]+=Jxw*f_value;
	  }
   }
//   Norm[0]=sqrt(Norm[0]);
//   Norm[1]=sqrt(Norm[1]);
   return Norm;
}

double BH2D::ComeputeMass(Vector<double> &f)
{
    double mass=0.0;
//    int n_ele=mesh.n_element();
    for(int i=0; i<n_ele; i++){
        std::vector<int> NV=mesh.getEleVtx(i);
        std::vector<MyPoint> EP(4);
        EP[0]=mesh.getPnt(NV[0]);
        EP[1]=mesh.getPnt(NV[1]);
        EP[2]=mesh.getPnt(NV[2]);
        EP[3]=mesh.getPnt(NV[3]);
        std::vector<MyPoint> GP=tmpEle.getGaussPnt();

	    std::vector<double> GW=tmpEle.getGaussWeight();
	    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	    std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	    std::vector<std::vector<double> > basis_value;
        basisValue(q_point, EP, basis_value);
	  
	  for(int j = 0; j < q_point.size(); j++){
	      double Jxw = GW[j]*jacobian[j];
		  
		  for(int l = 0; l < basis_value.size(); l++){
			  mass += Jxw*f(i*ele_dof+l)*basis_value[l][j];
		  }
      }
	  
    }
    return mass;
}

double BH2D::ComeputeEnergy()
{
    double energy=0.0;
//    int n_ele=mesh.n_element();
    for(int i=0; i<n_ele; i++){
        std::vector<int> NV=mesh.getEleVtx(i);
        std::vector<MyPoint> EP(4);
        EP[0]=mesh.getPnt(NV[0]);
        EP[1]=mesh.getPnt(NV[1]);
        EP[2]=mesh.getPnt(NV[2]);
        EP[3]=mesh.getPnt(NV[3]);
        std::vector<MyPoint> GP=tmpEle.getGaussPnt();

	    std::vector<double> GW=tmpEle.getGaussWeight();
	    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	    std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	    std::vector<std::vector<double> > basis_value;
        basisValue(q_point, EP, basis_value);
	  
	  for(int j = 0; j < q_point.size(); j++){
	      double Jxw = GW[j]*jacobian[j];
		  double phi_val=0.0;
		  double g1_val=0.0;
		  double g2_val=0.0;
		  
		  for(int l = 0; l < basis_value.size(); l++){
			  phi_val += u_h(i*ele_dof+l)*basis_value[l][j];
			  g1_val += g1_h(i*ele_dof+l)*basis_value[l][j];
			  g2_val += g2_h(i*ele_dof+l)*basis_value[l][j];
		  }
		  double M1_val = exp(-q_i[0]*phi_val);
		  double M2_val = exp(-q_i[1]*phi_val);
		  
		  // energy += Jxw*(0.5*q_h_val*q_h_val+U_DG[GaussPntNum*i+j]*U_DG[GaussPntNum*i+j]);
		  energy += Jxw*(g1_val*M1_val*mylog(g1_val*M1_val)+g2_val*M2_val*mylog(g2_val*M2_val));
		  energy += 0.5*Jxw*( q_i[0]*g1_val*M1_val+q_i[1]*g2_val*M2_val+_rho_(t-dt,q_point[j]) )*phi_val;
     }
	 
	 energy += getEnergy_g(i);
	 
    }
	
    return energy;
}

double BH2D::getEnergy_g(int n)
{
	double energyloc=0.0;
    std::vector<int> EV=mesh.getEleVtx(n);
    std::vector<MyPoint> EP=mesh.getPnt(EV);
    std::vector<MyPoint> GP=tmpEle.getGaussPnt();

    std::vector<double> GW=tmpEle.getGaussWeight();
    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
    std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
    std::vector<std::vector<double> > basis_value;
    basisValue(q_point, EP, basis_value);
	std::vector<std::vector<std::vector<double> > > basis_grad;
    basisGrad(q_point, EP, basis_grad);
		
//    std::vector<double> u_coef;
//    Read_u_Coef(n, u_coef);
//    Read_q_Coef(n, q_coef);
    
    double h=EP[2][0]-EP[0][0];
	
   	std::vector<double> ec(2);
   	ec[0]=(EP[0][0]+EP[1][0]+EP[2][0]+EP[3][0])/4.0;
   	ec[1]=(EP[0][1]+EP[1][1]+EP[2][1]+EP[3][1])/4.0;
    //     /*
    std::vector<int> eleEdg=mesh.getEleEdg(n);
    for(int i=0; i<eleEdg.size(); i++)
	{
		std::vector<int> edgVtx=mesh.getEdgVtx(eleEdg[i]);
// /*
		std::vector<MyPoint> edgi;
   		edgi=mesh.getPnt(edgVtx);
//		std::cout<<edgi[1][0]<<"   "<<edgi[1][1]<<std::endl;
// */
//  
/*
        std::vector<MyPoint> edgi(2);
   		edgi[0]=mesh.getPnt(edgVtx[0]);
   		edgi[1]=mesh.getPnt(edgVtx[1]);
//		std::cout<<edgi[1][0]<<"   "<<edgi[1][1]<<std::endl;
// */
        MyPoint midEdg=midpoint(edgi[0],edgi[1]);
		
//		std::cout<<ec[0]<<"   "<<ec[1]<<std::endl;

		int EdgLocation=1;
        if(fabs(midEdg[0]-ec[0])>1.0e-08){//该边平行y轴
			if(midEdg[0]-ec[0]>1.0e-08)
				EdgLocation=1;//该边是单元的右边
			else
				EdgLocation=3;//该边是单元的左边
		}
		else{//该边平行x轴
			if(midEdg[1]-ec[1]>1.0e-08)
				EdgLocation=2;//该边是单元的上边
			else
				EdgLocation=4;//该边是单元的下边
		}
//		int ik=element_patch[n][i];


        if(edge_patch[eleEdg[i]].size()==1)        // if the edge has one neighbor.
		{

		switch(EdgLocation){
			case 1: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();
	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double phi_val= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						phi_val += u_h(n*ele_dof+j)*basis_valueI[j][ie];
	
					}
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);								
							
//                  			Ele_Rhs[0][m]+=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);

//							Ele_Rhs[1][m]+=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);

//							Ele_Rhs[0][m]-=Jxw*beta0x/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0x/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                            Ele_Rhs[1][m]+=Jxw*beta0x/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0x/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);							
						}
//						std::cout<<"_u_grad_(t,q_pointE[ie])[0]= "<<_u_grad_(t,q_pointE[ie])[0]<<", _u_grad_(t,q_pointE[ie])[1]= "<<_u_grad_(t,q_pointE[ie])[1]<<std::endl;
//						Ele_Rhs[0][m]+=Jxw*beta_0/h*(_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);

//						Ele_Rhs[1][m]+=Jxw*zeroflux*M2_val*(exp(q_i[1]*_u_(t,q_pointE[ie]))*( _c2_grad_(t,q_pointE[ie])[0] + _c2_(t,q_pointE[ie])*q_i[1]*_u_grad_(t,q_pointE[ie])[0] )*basis_valueI[m][ie]);
						energyloc+=-0.5*Jxw*_u_grad_(t,q_pointE[ie])[0]*u_h(n*ele_dof+m)*basis_valueI[m][ie];


//						Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[0][m]-=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][0]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
						
//						Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]-=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][0]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
									
//						Ele_Rhs[1][m]+=Jxw*theta*(_u_grad_(t,q_pointE[ie])[0]*basis_valueI[m][ie]- _u_(t,q_pointE[ie])*basis_gradI[m][ie][0])
//						               +Jxw*antitheta*(_u_grad_(t-dt,q_pointE[ie])[0]*basis_valueI[m][ie]-_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
			
					}
					
					
          		}
				break;	
			}
			case 2: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double phi_val= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						phi_val += u_h(n*ele_dof+j)*basis_valueI[j][ie];
					}
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
//                  			Ele_Rhs[0][m]+=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_Rhs[1][m]+=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);

//							Ele_Rhs[0][m]-=Jxw*beta0y/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0y/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                           Ele_Rhs[1][m]+=Jxw*beta0y/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0y/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
				
//						Ele_Rhs[0][m]+=Jxw*zeroflux*M1_val*(exp(q_i[0]*_u_(t,q_pointE[ie]))*( _c1_grad_(t,q_pointE[ie])[1] + _c1_(t,q_pointE[ie])*q_i[0]*_u_grad_(t,q_pointE[ie])[1] )*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]+=Jxw*zeroflux*M2_val*(exp(q_i[1]*_u_(t,q_pointE[ie]))*( _c2_grad_(t,q_pointE[ie])[1] + _c2_(t,q_pointE[ie])*q_i[1]*_u_grad_(t,q_pointE[ie])[1] )*basis_valueI[m][ie]);
						energyloc+=-0.5*Jxw*_u_grad_(t,q_pointE[ie])[1]*u_h(n*ele_dof+m)*basis_valueI[m][ie];
					    
//					    Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//					    Ele_Rhs[0][m]-=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][1]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
					    
//					    Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//					    Ele_Rhs[1][m]-=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][1]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);


//						Ele_Rhs[1][m]+=Jxw*theta*(_u_grad_(t,q_pointE[ie])[1]*basis_valueI[m][ie]-_u_(t,q_pointE[ie])*basis_gradI[m][ie][1])
//						               +Jxw*antitheta*(_u_grad_(t-dt,q_pointE[ie])[1]*basis_valueI[m][ie]-_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
					
					}
//					energyloc+=0.5*Jxw*(_u_grad_(t,q_pointE[ie])[1]*phi_val);
          		}
				break;	
			}
			case 3: 
			{
				std::vector<double> edge_pnty(2, 0.0);
          
          		if(edgi[0][1]>edgi[1][1]){
                 	edge_pnty[0]=edgi[1][1];
                 	edge_pnty[1]=edgi[0][1];
             	}
             	else{                
                 	edge_pnty[0]=edgi[0][1];
                 	edge_pnty[1]=edgi[1][1];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pnty);
          		std::vector<double> GPy=tmpEdge.Local_to_Global(GPE, edge_pnty);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][0]=edgi[0][0];
             		q_pointE[ie][1]=GPy[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double phi_val= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						phi_val += u_h(n*ele_dof+j)*basis_valueI[j][ie];
					}
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
//                  			Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
											
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][0]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][0]);
//							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]-=Jxw*beta0x/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0x/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

 //                           Ele_Rhs[1][m]+=Jxw*beta0x/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0x/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
						
//						Ele_Rhs[0][m]-=Jxw*zeroflux*M1_val*(exp(q_i[0]*_u_(t,q_pointE[ie]))*( _c1_grad_(t,q_pointE[ie])[0] + _c1_(t,q_pointE[ie])*q_i[0]*_u_grad_(t,q_pointE[ie])[0] )*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]-=Jxw*zeroflux*M2_val*(exp(q_i[1]*_u_(t,q_pointE[ie]))*( _c2_grad_(t,q_pointE[ie])[0] + _c2_(t,q_pointE[ie])*q_i[1]*_u_grad_(t,q_pointE[ie])[0] )*basis_valueI[m][ie]);
						energyloc+=0.5*Jxw*_u_grad_(t,q_pointE[ie])[0]*u_h(n*ele_dof+m)*basis_valueI[m][ie];
						
//						Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[0][m]+=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][0]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
						
//						Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]+=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][0]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);


//						Ele_Rhs[1][m]+=Jxw*theta*(-_u_grad_(t,q_pointE[ie])[0]*basis_valueI[m][ie]+_u_(t,q_pointE[ie])*basis_gradI[m][ie][0])
//						               +Jxw*antitheta*(-_u_grad_(t-dt,q_pointE[ie])[0]*basis_valueI[m][ie]+_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][0]);
					}
//					energyloc+=0.5*Jxw*(_u_grad_(t,q_pointE[ie])[0]*phi_val);
          		}
				break;	
			}
			case 4: 
			{
				std::vector<double> edge_pntx(2, 0.0);
          
          		if(edgi[0][0]>edgi[1][0]){
                 	edge_pntx[0]=edgi[1][0];
                 	edge_pntx[1]=edgi[0][0];
             	}
             	else{                
                 	edge_pntx[0]=edgi[0][0];
                 	edge_pntx[1]=edgi[1][0];
             	}

          		std::vector<double> GPE=tmpEdge.getGaussPnt();

	  			std::vector<double> GWE=tmpEdge.getGaussWeight();
	  			std::vector<double> jacobianE = tmpEdge.Local_to_Global_jacobian(GPE, edge_pntx);
          		std::vector<double> GPx=tmpEdge.Local_to_Global(GPE, edge_pntx);

          		std::vector<MyPoint> q_pointE(GPE.size()); 
          		for(int ie=0; ie<GPE.size(); ie++){
             		q_pointE[ie][1]=edgi[0][1];
             		q_pointE[ie][0]=GPx[ie];
          		}

	  			std::vector<std::vector<double> > basis_valueI;
    	  		basisValue(q_pointE, EP, basis_valueI);
	  			std::vector<std::vector<std::vector<double> > > basis_gradI;
    	  		basisGrad(q_pointE, EP, basis_gradI);

	  			for(int ie = 0; ie < q_pointE.size(); ie++){
	      			double Jxw = GWE[ie]*jacobianE[ie];
					
					double phi_val= 0.0;
					for(int j=0; j<basis_valueI.size(); j++)
					{
						phi_val += u_h(n*ele_dof+j)*basis_valueI[j][ie];
					}
					
					double M1_val=exp(-q_i[0]*phi_val);
					double M2_val=exp(-q_i[1]*phi_val);
					
              		for(int m=0; m<basis_valueI.size(); m++){
              			for(int j=0; j<basis_valueI.size(); j++){
//							Ele_Rhs[0][m]-=Jxw*antitheta*beta_0/h3*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta_0/h3*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);
							
//                  			Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
											
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*(basis_gradI[j][ie][1]*basis_valueI[m][ie]+basis_valueI[j][ie]*basis_gradI[m][ie][1]);
							
//							Ele_Rhs[0][m]-=Jxw*antitheta*q_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
							
//							Ele_Rhs[1][m]-=0.0*Jxw*antitheta*u_coef[j]*beta0/h*(basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]-=Jxw*beta0y/h*antitheta*(q_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[0][m]+=Jxw*beta0y/h*antitheta*(q_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);

//                            Ele_Rhs[1][m]+=Jxw*beta0y/h*antitheta*(u_coef[j]*basis_valueI[j][ie]*basis_valueI[m][ie]);
											
//							Ele_Rhs[1][m]-=Jxw*beta0y/h*antitheta*(u_coef0[j]*basis_valueE[j][ie]*basis_valueI[m][ie]);											
						}
//						std::cout<<"t= "<<t<<", dt= "<<dt<<", q_pointE[ie]= "<<q_pointE[ie][0]<<" , q_pointE[ie][1]= "<<q_pointE[ie][1]<<",  _u_(t,q_pointE(ie))= "<<_u_(t,q_pointE[ie])<<std::endl;

//						Ele_Rhs[0][m]-=Jxw*zeroflux*M1_val*(exp(q_i[0]*_u_(t,q_pointE[ie]))*( _c1_grad_(t,q_pointE[ie])[1] + _c1_(t,q_pointE[ie])*q_i[0]*_u_grad_(t,q_pointE[ie])[1] )*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]-=Jxw*zeroflux*M2_val*(exp(q_i[1]*_u_(t,q_pointE[ie]))*( _c2_grad_(t,q_pointE[ie])[1] + _c2_(t,q_pointE[ie])*q_i[1]*_u_grad_(t,q_pointE[ie])[1] )*basis_valueI[m][ie]);
						energyloc+=0.5*Jxw*_u_grad_(t,q_pointE[ie])[1]*u_h(n*ele_dof+m)*basis_valueI[m][ie];
						
//						Ele_Rhs[0][m]+=Jxw*beta0/h*(theta*_q_(t,q_pointE[ie])*basis_valueI[m][ie]+antitheta*_q_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[0][m]+=Jxw*(theta*_q_(t,q_pointE[ie])*basis_gradI[m][ie][1]+antitheta*_q_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
//						
//						Ele_Rhs[1][m]+=Jxw*beta0/h*(theta*_u_(t,q_pointE[ie])*basis_valueI[m][ie]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_valueI[m][ie]);
//						Ele_Rhs[1][m]+=Jxw*(theta*_u_(t,q_pointE[ie])*basis_gradI[m][ie][1]+0.0*antitheta*_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
						
//						Ele_Rhs[1][m]+=Jxw*theta*(-_u_grad_(t,q_pointE[ie])[1]*basis_valueI[m][ie]+_u_(t,q_pointE[ie])*basis_gradI[m][ie][1])
//						               +Jxw*antitheta*(-_u_grad_(t-dt,q_pointE[ie])[1]*basis_valueI[m][ie]+_u_(t-dt,q_pointE[ie])*basis_gradI[m][ie][1]);
						               
					}
					
					energyloc+=0.5*Jxw*(_u_grad_(t,q_pointE[ie])[1]*phi_val);
          		}
				break;
			}
			default:
			{
				std::cout<<"Some problem in location the edg, please chech carefully!"<<std::endl;
				break;	
			}
		}
				
		} // end of the one neighbor case.
		
			
		
		
	}
	
	return energyloc;
	
}

double BH2D::mylog(double x)
{
	double y=0.0;
	if (x>0.0)
		y = log(x);
	if (x<=0.0)
	    std::cout<<"g is a negative value"<<std::endl;
		
	return y;
}

void BH2D::Compute_U_Real(Vector<double> &f, std::vector<MyPoint> &x_gnt, std::vector<double> &U_Gnt)
{
    int n_ele=mesh.n_element();
    std::vector<MyPoint> MP=tmpEle.getMeshPnt();
    int MeshPntNum = MP.size();
    U_Gnt.resize(MeshPntNum*n_ele);
    x_gnt.resize(MeshPntNum*n_ele);
	
   for(int j=0; j<n_ele; j++){
      std::vector<int> NV=mesh.getEleVtx(j);
      std::vector<MyPoint> EP(4);
      EP[0]=mesh.getPnt(NV[0]);
      EP[1]=mesh.getPnt(NV[1]);
      EP[2]=mesh.getPnt(NV[2]);
      EP[3]=mesh.getPnt(NV[3]);
      std::vector<MyPoint> GP=tmpEle.getGaussPnt();
	  std::vector<double> GW=tmpEle.getGaussWeight();
	  std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	  std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	  std::vector<std::vector<double> > basis_value;
      basisValue(q_point, EP, basis_value);
	  
	  std::vector<MyPoint> MP=tmpEle.getMeshPnt();
	  std::vector<MyPoint> Mq_point = tmpEle.Local_to_Global(MP, EP);
	  
	  std::vector<std::vector<double> > basis_valueM;
      basisValue(Mq_point, EP, basis_valueM);

      	for(int i = 0; i < Mq_point.size(); i++)
		{
			double u_h_val=0.0;
			for(int m=0; m<basis_valueM.size(); m++)
			{
				u_h_val+=f(ele_dof*j+m)*basis_valueM[m][i];
			}
            x_gnt[MeshPntNum*j+i] = Mq_point[i];
            U_Gnt[MeshPntNum*j+i] = u_h_val;	
		}
		
   }	
}

double BH2D::ComputeL2Error(double t, Vector<double> &f, double (*g)(const double ,const double *))
{
//	int ele_dof=9;
   double err=0.0;
   int n_ele=mesh.n_element();
   for(int j=0; j<n_ele; j++){
      std::vector<int> NV=mesh.getEleVtx(j);
      std::vector<MyPoint> EP(4);
      EP[0]=mesh.getPnt(NV[0]);
      EP[1]=mesh.getPnt(NV[1]);
      EP[2]=mesh.getPnt(NV[2]);
      EP[3]=mesh.getPnt(NV[3]);
      std::vector<MyPoint> GP=tmpEle.getGaussPnt();

	  std::vector<double> GW=tmpEle.getGaussWeight();
	  std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	  std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	  std::vector<std::vector<double> > basis_value;
      basisValue(q_point, EP, basis_value);
	  
	  for(int i = 0; i < q_point.size(); i++){
	      double Jxw = GW[i]*jacobian[i];
	      double f_value = g(t,q_point[i]);
              double u_h_val=0.0;
              for(int m=0; m<basis_value.size(); m++){
                  u_h_val+=f(ele_dof*j+m)*basis_value[m][i];
              }
              double df_val=f_value-u_h_val;
              err+=Jxw*df_val*df_val;
	  }
   }
   err=sqrt(err);
   return err;
}

double BH2D::ComputeL1Error(double t, Vector<double> &f, double (*g)(const double ,const double *))
{
//	int ele_dof=9;
   double err=0.0;
   int n_ele=mesh.n_element();
   for(int j=0; j<n_ele; j++){
      std::vector<int> NV=mesh.getEleVtx(j);
      std::vector<MyPoint> EP(4);
      EP[0]=mesh.getPnt(NV[0]);
      EP[1]=mesh.getPnt(NV[1]);
      EP[2]=mesh.getPnt(NV[2]);
      EP[3]=mesh.getPnt(NV[3]);
      std::vector<MyPoint> GP=tmpEle.getGaussPnt();

	  std::vector<double> GW=tmpEle.getGaussWeight();
	  std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	  std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	  std::vector<std::vector<double> > basis_value;
      basisValue(q_point, EP, basis_value);
	  
	  for(int i = 0; i < q_point.size(); i++){
	      double Jxw = GW[i]*jacobian[i];
	      double f_value = g(t,q_point[i]);
              double u_h_val=0.0;
              for(int m=0; m<basis_value.size(); m++){
                  u_h_val+=f(ele_dof*j+m)*basis_value[m][i];
              }
              double df_val=fabs(f_value-u_h_val);
              err+=Jxw*df_val;
	  }
   }
   return err;
}

double BH2D::ComputeH1Error(double t, Vector<double> &f, std::vector<double> (*g)(const double ,const double *))
{
//	int ele_dof=9;
    double err=0.0;
    int n_ele=mesh.n_element();
    for(int j=0; j<n_ele; j++){
        std::vector<int> NV=mesh.getEleVtx(j);
        std::vector<MyPoint> EP(4);
        EP[0]=mesh.getPnt(NV[0]);
        EP[1]=mesh.getPnt(NV[1]);
        EP[2]=mesh.getPnt(NV[2]);
        EP[3]=mesh.getPnt(NV[3]);
        std::vector<MyPoint> GP=tmpEle.getGaussPnt();
	    std::vector<double> GW=tmpEle.getGaussWeight();
  	    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	    std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	    std::vector<std::vector<std::vector<double> > > basis_grad;
    	basisGrad(q_point, EP, basis_grad);	  
	    for(int i = 0; i < q_point.size(); i++){
	        double Jxw = GW[i]*jacobian[i];
	        std::vector<double> f_value = g(t,q_point[i]);
            double u_h_grad1=0.0;
            double u_h_grad2=0.0;
            for(int m=0; m<basis_grad.size(); m++){
                u_h_grad1+=f(ele_dof*j+m)*basis_grad[m][i][0];
                u_h_grad2+=f(ele_dof*j+m)*basis_grad[m][i][1];
            }
            double df_val1=f_value[0]-u_h_grad1;
            double df_val2=f_value[1]-u_h_grad2;
            err+=Jxw*df_val1*df_val1;
            err+=Jxw*df_val2*df_val2;
	    }
    }
    err=sqrt(err);
    return err;
}

double BH2D::ComputeLINFError(double t, Vector<double> &f, double (*g)(const double ,const double *))
{
//	int ele_dof=9;
    double err=0.0;
    int n_ele=mesh.n_element();
    for(int j=0; j<n_ele; j++){
        std::vector<int> NV=mesh.getEleVtx(j);
        std::vector<MyPoint> EP=mesh.getPnt(NV);
        std::vector<MyPoint> GP=tmpEle.getGaussPnt();
        std::vector<double> GW=tmpEle.getGaussWeight();
	    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
     	std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	    std::vector<std::vector<double> > basis_value;
        basisValue(q_point, EP, basis_value);
	    for(int i = 0; i < q_point.size(); i++){
	        double f_value = g(t,q_point[i]);
	        double u_h_val=0.0;
	        for(int m=0; m<basis_value.size(); m++){
                  u_h_val+=f(ele_dof*j+m)*basis_value[m][i];
//                  std::cout<<"f(ele_dof*j+m)= "<<f(ele_dof*j+m)<<"   basis_valueInt[m][i]= "<<basis_value[m][i]<<std::endl;
            }
            double df_val=fabs(f_value-u_h_val);
            if(df_val>err)
            {
//		buildSparseMatrix();
                err=df_val;
			}
	    }
	    
	    std::vector<std::vector<double> > basis_valueInt;
	    basisValue(EP, EP, basis_valueInt);
	    for(int i = 0; i < EP.size(); i++){
	        double f_value = g(t,EP[i]);
	        double u_h_val=0.0;
	        for(int m=0; m<basis_valueInt.size(); m++){
                  u_h_val+=f(ele_dof*j+m)*basis_valueInt[m][i];
            }
            double df_val=fabs(f_value-u_h_val);
            if(fabs(df_val)>err)
            {
                err=df_val;
            }
	    }
    }
    return err;
}

void BH2D::L2Proj(double (*Init1)(const double *), double (*Init2)(const double *))
{
//	int ele_dof=9;
    int n_ele=mesh.n_element();
    c1_h.reinit(ele_dof*n_ele);
	c2_h.reinit(ele_dof*n_ele);
//    ui_h.reinit(ele_dof*n_ele);

    for(int i=0; i<n_ele; i++){

		std::vector<std::vector<double> > A_matrix, B_matrix;
        A_matrix.resize(ele_dof);
		B_matrix.resize(ele_dof);
        for(int j=0; j<ele_dof; j++){
            A_matrix[j].resize(ele_dof);
			B_matrix[j].resize(ele_dof);
        }
        for(int ai=0; ai<A_matrix.size(); ai++){
            for(int aj=0; aj<A_matrix[ai].size(); aj++){
                A_matrix[ai][aj]=0.0;
				B_matrix[ai][aj]=0.0;
             }
        }

         std::vector<double> LRHS(ele_dof), LRHS2(ele_dof);
         std::vector<double> LU(ele_dof), LU2(ele_dof);
		 
         std::vector<int> EV=mesh.getEleVtx(i);
         std::vector<MyPoint> EP=mesh.getPnt(EV);
         std::vector<MyPoint> GP=tmpEle.getGaussPnt();
	     std::vector<double> GW=tmpEle.getGaussWeight();
	     std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	     std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	     std::vector<std::vector<double> > basis_value;
	     
	     basisValue(q_point, EP, basis_value);
	  
	  for(int j = 0; j < q_point.size(); j++){
	      double Jxw = GW[j]*jacobian[j];
		  
	      double f1_value = Init1(q_point[j]);
		  double f2_value = Init2(q_point[j]);
		  
	      for(int k = 0; k < basis_value.size(); k++){
	         for(int m = 0; m < basis_value.size(); m++){
	             A_matrix[k][m] += Jxw*basis_value[m][j]*basis_value[k][j];
				 B_matrix[k][m] += Jxw*basis_value[m][j]*basis_value[k][j];
             }
              LRHS[k]+=Jxw*f1_value*basis_value[k][j];
			  LRHS2[k]+=Jxw*f2_value*basis_value[k][j];
	      }
	  }
	  
      GaussElimination(A_matrix, LU, LRHS);
	  GaussElimination(B_matrix, LU2, LRHS2);
      for(int j=0; j<ele_dof_Max; j++){
      	      c1_h(ele_dof*i+j)=LU[j];
			  c2_h(ele_dof*i+j)=LU2[j];
	  }
   }
}

void BH2D::Solg()
{

    for(int i=0; i<n_ele; i++){

		std::vector<std::vector<double> > A_matrix, B_matrix;
        A_matrix.resize(ele_dof);
		B_matrix.resize(ele_dof);
        for(int j=0; j<ele_dof; j++){
            A_matrix[j].resize(ele_dof);
			B_matrix[j].resize(ele_dof);
        }
        for(int ai=0; ai<A_matrix.size(); ai++){
            for(int aj=0; aj<A_matrix[ai].size(); aj++){
                A_matrix[ai][aj]=0.0;
				B_matrix[ai][aj]=0.0;
             }
        }

         std::vector<double> LRHS(ele_dof), LRHS2(ele_dof);
         std::vector<double> LU(ele_dof), LU2(ele_dof);
		 
         std::vector<int> EV=mesh.getEleVtx(i);
         std::vector<MyPoint> EP=mesh.getPnt(EV);
         std::vector<MyPoint> GP=tmpEle.getGaussPnt();
	     std::vector<double> GW=tmpEle.getGaussWeight();
	     std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	     std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	     std::vector<std::vector<double> > basis_value;
	     
	     basisValue(q_point, EP, basis_value);
	  
		for(int j = 0; j < q_point.size(); j++){
	      double Jxw = GW[j]*jacobian[j];
		  
		double phi_val=0.0, c1_val = 0.0, c2_val=0.0;
		for(int k=0; k<basis_value.size(); k++)
		{
			phi_val += u_h(i*ele_dof+k)*basis_value[k][j];
			c1_val += c1_h(i*ele_dof+k)*basis_value[k][j];
			c2_val += c2_h(i*ele_dof+k)*basis_value[k][j];
		}
		
		double M1_val=exp(-q_i[0]*phi_val);
		double M2_val=exp(-q_i[1]*phi_val);
		  
	      for(int k = 0; k < basis_value.size(); k++){
	         for(int m = 0; m < basis_value.size(); m++){
	             A_matrix[k][m] += Jxw*M1_val*basis_value[m][j]*basis_value[k][j];
				 B_matrix[k][m] += Jxw*M2_val*basis_value[m][j]*basis_value[k][j];
             }
              LRHS[k]+=Jxw*c1_val*basis_value[k][j];
			  LRHS2[k]+=Jxw*c2_val*basis_value[k][j];
	      }
	  }
	  
      GaussElimination(A_matrix, LU, LRHS);
	  GaussElimination(B_matrix, LU2, LRHS2);
      for(int j=0; j<ele_dof_Max; j++){
      	      g1_h(ele_dof*i+j)=LU[j];
			  g2_h(ele_dof*i+j)=LU2[j];
	  }
   }
}

void BH2D::Solgtest()
{

    for(int i=0; i<n_ele; i++){

		std::vector<std::vector<double> > A_matrix, B_matrix;
        A_matrix.resize(ele_dof);
		B_matrix.resize(ele_dof);
        for(int j=0; j<ele_dof; j++){
            A_matrix[j].resize(ele_dof);
			B_matrix[j].resize(ele_dof);
        }
        for(int ai=0; ai<A_matrix.size(); ai++){
            for(int aj=0; aj<A_matrix[ai].size(); aj++){
                A_matrix[ai][aj]=0.0;
				B_matrix[ai][aj]=0.0;
             }
        }

         std::vector<double> LRHS(ele_dof), LRHS2(ele_dof);
         std::vector<double> LU(ele_dof), LU2(ele_dof);
		 
         std::vector<int> EV=mesh.getEleVtx(i);
         std::vector<MyPoint> EP=mesh.getPnt(EV);
         std::vector<MyPoint> GP=tmpEle.getGaussPnt();
	     std::vector<double> GW=tmpEle.getGaussWeight();
	     std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	     std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	     std::vector<std::vector<double> > basis_value;
	     
	     basisValue(q_point, EP, basis_value);
	  
		for(int j = 0; j < q_point.size(); j++){
	      double Jxw = GW[j]*jacobian[j];
		  
		double phi_val=0.0, c1_val = 0.0, c2_val=0.0;
		for(int k=0; k<basis_value.size(); k++)
		{
			phi_val += u_h(i*ele_dof+k)*basis_value[k][j];
			c1_val += c1_h(i*ele_dof+k)*basis_value[k][j];
			c2_val += c2_h(i*ele_dof+k)*basis_value[k][j];
		}
		
		double M1_val=exp(-q_point[j][0]*q_point[j][1]);
		double M2_val=exp(-q_i[1]*phi_val);
		  
	      for(int k = 0; k < basis_value.size(); k++){
	         for(int m = 0; m < basis_value.size(); m++){
	             A_matrix[k][m] += Jxw*M1_val*basis_value[m][j]*basis_value[k][j];
				 B_matrix[k][m] += Jxw*M2_val*basis_value[m][j]*basis_value[k][j];
             }
              LRHS[k]+=Jxw*basis_value[k][j];
			  LRHS2[k]+=Jxw*c2_val*basis_value[k][j];
	      }
	  }
	  
      GaussElimination(A_matrix, LU, LRHS);
	  GaussElimination(B_matrix, LU2, LRHS2);
      for(int j=0; j<ele_dof_Max; j++){
      	      g1_h(ele_dof*i+j)=LU[j];
			  g2_h(ele_dof*i+j)=LU2[j];
	  }
   }
}

// /*
void BH2D::Limiter()
{
	std::vector<double> locgMin, locgMinLim;
	locgMin.resize(2,1.0e2);
	locgMinLim.resize(2,1.0e2);
    for(int i=0; i<n_ele; i++){
		
         std::vector<int> EV=mesh.getEleVtx(i);
         std::vector<MyPoint> EP=mesh.getPnt(EV);
         std::vector<MyPoint> GP=tmpEle.getGaussPnt();
	     std::vector<double> GW=tmpEle.getGaussWeight();
	     std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	     std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	     std::vector<std::vector<double> > basis_value;
	     basisValue(q_point, EP, basis_value);
		 
	     std::vector<MyPoint> MP=tmpEle.getMeshPnt();
	  
		std::vector<MyPoint> SP1=buildSamPnt(i,q_i[0]);
		std::vector<MyPoint> SP2=buildSamPnt(i,q_i[1]);
	    std::vector<double> g_avg = getAvg(i);
		std::vector<double> g_min = getMin(i,SP1,SP2);
// /*
		std::vector<double> g_minG = getMin(i,GP,GP);
		std::vector<double> g_minM = getMin(i,MP,MP);
//		std::cout<<g_min[0]<<std::endl; 
		g_min[0] = std::min(g_min[0], g_minG[0]);
		g_min[1] = std::min(g_min[1], g_minG[1]);
		g_min[0] = std::min(g_min[0], g_minM[0]);
		g_min[1] = std::min(g_min[1], g_minM[1]);
// */		
		if(g_avg[0]<=0.0 || g_avg[1]<=0.0)
			std::cout<<"A mistake -------------------"<<std::endl;
		
//		double theta1 = std::min(1.0, (g_avg[0]-LimTol)/(g_avg[0]-g_min[0]));
//		double theta2 = std::min(1.0, (g_avg[1]-LimTol)/(g_avg[1]-g_min[1]));

// /*		
		double theta1=1.0;
		double theta2=1.0;
		if(g_min[0]>0.0)
			theta1 = std::min(1.0, g_avg[0]/(g_avg[0]-g_min[0]));
		else
			theta1 = std::min(1.0, (g_avg[0]-std::min(LimTol,g_avg[0]))/(g_avg[0]-g_min[0]));
			
		if(g_min[1]>0.0)
			theta2 = std::min(1.0, g_avg[1]/(g_avg[1]-g_min[1]));
		else
			theta2 = std::min(1.0, (g_avg[1]-std::min(LimTol,g_avg[1]))/(g_avg[1]-g_min[1]));
// */	
		if (theta1*theta2!=1.0)
			std::cout<<"theta1= "<<theta1<<", theta2= "<<theta2<<std::endl;
	

		for(int j=0; j<ele_dof; j++)
		{
			g1_h(i*ele_dof+j)=theta1*g1_h(i*ele_dof+j);
			g2_h(i*ele_dof+j)=theta2*g2_h(i*ele_dof+j);
		}
//		std::cout<<g_avg[0]<<"  "<<g_avg[1]<<std::endl;
		g1_h(i*ele_dof)+=(1.0-theta1)*g_avg[0];
		g2_h(i*ele_dof)+=(1.0-theta2)*g_avg[1];
		
		std::vector<double> lg_min = getMin(i,SP1,SP2);
		std::vector<double> lg_minG = getMin(i,GP,GP);
		std::vector<double> lg_minM = getMin(i,MP,MP);
		
		locgMin[0] = std::min(locgMin[0], g_min[0]);
		locgMin[1] = std::min(locgMin[1], g_min[1]);
// /*		
		locgMin[0] = std::min(locgMin[0], g_minG[0]);
		locgMin[0] = std::min(locgMin[0], g_minM[0]);
		
		locgMin[1] = std::min(locgMin[1], g_minG[1]);
		locgMin[1] = std::min(locgMin[1], g_minM[1]);
// */		
		
		locgMinLim[0] = std::min(locgMinLim[0], lg_min[0]);
		locgMinLim[1] = std::min(locgMinLim[1], lg_min[1]);
// /*		
		locgMinLim[0] = std::min(locgMinLim[0], lg_minG[0]);
		locgMinLim[0] = std::min(locgMinLim[0], lg_minM[0]);
		locgMinLim[1] = std::min(locgMinLim[1], lg_minG[1]);
		locgMinLim[1] = std::min(locgMinLim[1], lg_minM[1]);
// */		
    }
	
	gStepMin[0] = locgMin[0];
	gStepMin[1] = locgMin[1];
	gStepMin[2] = locgMinLim[0];
	gStepMin[3] = locgMinLim[1];
	
}
// */

std::vector<MyPoint> BH2D::buildSamPnt(int n, int q_val)
{
	std::vector<double> GPE=tmpEdge.getGaussPnt();
	std::vector<double> GWE=tmpEdge.getGaussWeight();
	std::vector<MyPoint>  SamPnt;
	int GPNO = GPE.size();
	int SPNO = 3*GPNO;
	SamPnt.resize(2*SPNO);
	for(int i=1; i<3; i++)
		for(int j=0; j<GPNO; j++)
		{
			SamPnt[i*GPNO+j][0] = GPE[j];
			SamPnt[i*GPNO+j][1] = 2*i-3;
			SamPnt[SPNO+i*GPNO+j][0] = 2*i-3;
			SamPnt[SPNO+i*GPNO+j][1] = GPE[j];
		}
	
	std::vector<int> EV=mesh.getEleVtx(n);
	std::vector<MyPoint> EP=mesh.getPnt(EV);
	std::vector<MyPoint> q_pointE(GPNO); 
	
	for(int i=0; i<GPNO; i++)    // finding gamma^y
	{
		double leftint[2],rightint[2]; 
		leftint[0]=0.0;
		leftint[1]=0.0;
		rightint[0]=0.0;
		rightint[1]=0.0;
		
		for(int j=0; j<GPNO; j++)
		{
			q_pointE[j][0]=GPE[i];
			q_pointE[j][1]=GPE[j];
		}
		
		std::vector<std::vector<double> > basis_valueI;
    	basisValue(q_pointE, EP, basis_valueI);
		
		for(int ie = 0; ie < q_pointE.size(); ie++){
	      	double Jxw = GWE[ie];
			
			double phi_val=0.0;
			for(int k=0; k<basis_valueI.size(); k++)
			{
				phi_val += u_h(n*ele_dof+k)*basis_valueI[k][ie];
			}
		
			double M_val=exp(-q_val*phi_val);
//			double M2_val=exp(-q_i[1]*phi_val);
			
            for(int m=0; m<basis_valueI.size(); m++){
				leftint[0] += 0.5*Jxw*(q_pointE[ie][1]-q_pointE[ie][1]*q_pointE[ie][1])*M_val; //top
				leftint[1] += 0.5*Jxw*(1.0-q_pointE[ie][1])*M_val;
				rightint[0] += 0.5*Jxw*(q_pointE[ie][1]+q_pointE[ie][1]*q_pointE[ie][1])*M_val; //top
				rightint[1] += 0.5*Jxw*(1.0+q_pointE[ie][1])*M_val;
			}
		}
		
		double aj = leftint[0]/leftint[1];
		double bj = rightint[0]/rightint[1];
//		std::cout<<"aj= "<<aj<<", bj= "<<bj<<", gammaj= "<<bj-aj<<std::endl;
		SamPnt[i][0] = GPE[i];
		SamPnt[i][1] = (aj+bj)/2.0;
		if(aj>=bj)
		{
			std::cout<<"There is a mistake in find gamma!"<<std::endl;
		}
	} 

	for(int i=0; i<GPNO; i++)    // finding gamma^x
	{
		double leftint[2],rightint[2]; 
		leftint[0]=0.0;
		leftint[1]=0.0;
		rightint[0]=0.0;
		rightint[1]=0.0;
		
		for(int j=0; j<GPNO; j++)
		{
			q_pointE[j][0]=GPE[j];
			q_pointE[j][1]=GPE[i];
		}
		
		std::vector<std::vector<double> > basis_valueI;
    	basisValue(q_pointE, EP, basis_valueI);
		
		for(int ie = 0; ie < q_pointE.size(); ie++){
	      	double Jxw = GWE[ie];
			
			double phi_val=0.0;
			for(int k=0; k<basis_valueI.size(); k++)
			{
				phi_val += u_h(n*ele_dof+k)*basis_valueI[k][ie];
			}
		
			double M_val=exp(-q_val*phi_val);
//			double M2_val=exp(-q_i[1]*phi_val);
			
            for(int m=0; m<basis_valueI.size(); m++){
				leftint[0] += 0.5*Jxw*(q_pointE[ie][0]-q_pointE[ie][0]*q_pointE[ie][0])*M_val; //top
				leftint[1] += 0.5*Jxw*(1.0-q_pointE[ie][0])*M_val;
				rightint[0] += 0.5*Jxw*(q_pointE[ie][0]+q_pointE[ie][0]*q_pointE[ie][0])*M_val; //top
				rightint[1] += 0.5*Jxw*(1.0+q_pointE[ie][0])*M_val;
			}
		}
		
		double aj = leftint[0]/leftint[1];
		double bj = rightint[0]/rightint[1];
//		std::cout<<"aj= "<<aj<<", bj= "<<bj<<", gammaj= "<<bj-aj<<std::endl;
		SamPnt[SPNO+i][0] = (aj+bj)/2.0;
		SamPnt[SPNO+i][1] = GPE[i];
		if(aj>=bj)
		{
			std::cout<<"There is a mistake in find gamma!"<<std::endl;
		}
	}

	return SamPnt;
}

std::vector<double> BH2D::getAvg(int n)
{
    std::vector<int> EV=mesh.getEleVtx(n);
    std::vector<MyPoint> EP=mesh.getPnt(EV);
    std::vector<MyPoint> GP=tmpEle.getGaussPnt();

    std::vector<double> GW=tmpEle.getGaussWeight();
    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
    std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
    std::vector<std::vector<double> > basis_value;
    basisValue(q_point, EP, basis_value);
	
	std::vector<double> g_avg, g1_avg, g2_avg;
	g_avg.resize(2,0.0);
	g1_avg.resize(2,0.0);
	g2_avg.resize(2,0.0);
	
	for(int j = 0; j < q_point.size(); j++){
	    double Jxw = GW[j]*jacobian[j];
		  
		double phi_val=0.0;
		for(int k=0; k<basis_value.size(); k++)
		{
			phi_val += u_h(n*ele_dof+k)*basis_value[k][j];
		}
		
		double M1_val=exp(-q_i[0]*phi_val);
		double M2_val=exp(-q_i[1]*phi_val);
		  
	    for(int k = 0; k < basis_value.size(); k++){
            g1_avg[0] += Jxw*M1_val*g1_h(n*ele_dof+k)*basis_value[k][j];
			g2_avg[0] += Jxw*M2_val*g2_h(n*ele_dof+k)*basis_value[k][j];
	    }
		g1_avg[1] += Jxw*M1_val;
		g2_avg[1] += Jxw*M2_val;
	}
	
	g_avg[0] = g1_avg[0]/g1_avg[1];
	g_avg[1] = g2_avg[0]/g2_avg[1];
	
	return g_avg;
}

std::vector<double> BH2D::getMin(int n, std::vector<MyPoint> SP1, std::vector<MyPoint> SP2)
{
    std::vector<int> EV=mesh.getEleVtx(n);
    std::vector<MyPoint> EP=mesh.getPnt(EV);
//    std::vector<MyPoint> GP=tmpEle.getGaussPnt();

//    std::vector<double> GW=tmpEle.getGaussWeight();
    std::vector<MyPoint> q_point1 = tmpEle.Local_to_Global(SP1, EP); 
	std::vector<MyPoint> q_point2 = tmpEle.Local_to_Global(SP2, EP); 
//    std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(SP, EP);
    std::vector<std::vector<double> > basis_value1, basis_value2;
    basisValue(q_point1, EP, basis_value1);
	basisValue(q_point2, EP, basis_value2);
	
	std::vector<double> g_min;
	g_min.resize(2,1.0e2);
	
	for(int j = 0; j < q_point1.size(); j++){
		double g1_val=0.0, g2_val=0.0;
	    for(int k = 0; k < basis_value1.size(); k++){
            g1_val += g1_h(n*ele_dof+k)*basis_value1[k][j];
			g2_val += g2_h(n*ele_dof+k)*basis_value2[k][j];
	    }
		
		g_min[0] = std::min(g_min[0], g1_val);
		g_min[1] = std::min(g_min[1], g2_val);
	}
	
	return g_min;
}

double BH2D::getMinAvg(Vector<double> &f)
{
	double fmin=1.0e2;
	for(int n=0; n<n_ele; n++)
	{
    std::vector<int> EV=mesh.getEleVtx(n);
    std::vector<MyPoint> EP=mesh.getPnt(EV);
    std::vector<MyPoint> GP=tmpEle.getGaussPnt();

    std::vector<double> GW=tmpEle.getGaussWeight();
    std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
    std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
    std::vector<std::vector<double> > basis_value;
    basisValue(q_point, EP, basis_value);
	
	double h=EP[2][0]-EP[0][0];
	
	double locmin=0.0;
	
	for(int j = 0; j < q_point.size(); j++){
	    double Jxw = GW[j]*jacobian[j];
		  		  
	    for(int k = 0; k < basis_value.size(); k++){
            locmin += Jxw*f(n*ele_dof+k)*basis_value[k][j];
	    }
	}
	locmin = locmin/(h*h);
	fmin = std::min(fmin, locmin);
	}
	return fmin;
}

void BH2D::UpdateC()
{
    int n_ele=mesh.n_element();
	buildRHS_C();
    for(int i=0; i<n_ele; i++){
        std::vector<std::vector<double> > A_matrix, B_matrix;
        A_matrix.resize(ele_dof);
        B_matrix.resize(ele_dof);
        for(int j=0; j<ele_dof; j++){
            A_matrix[j].resize(ele_dof);
            B_matrix[j].resize(ele_dof);
        }
        for(int ai=0; ai<A_matrix.size(); ai++){
            for(int aj=0; aj<A_matrix[ai].size(); aj++){
                A_matrix[ai][aj]=0.0;
                B_matrix[ai][aj]=0.0;
             }
         }
         std::vector<double> LRHS(ele_dof);
         std::vector<double> LRHS2(ele_dof);
         std::vector<double> LU(ele_dof);
         std::vector<double> LU2(ele_dof);
         std::vector<int> EV=mesh.getEleVtx(i);
         std::vector<MyPoint> EP=mesh.getPnt(EV);
         std::vector<MyPoint> GP=tmpEle.getGaussPnt();

	     std::vector<double> GW=tmpEle.getGaussWeight();
	     std::vector<MyPoint> q_point = tmpEle.Local_to_Global(GP, EP); 
	     std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	     std::vector<std::vector<double> > basis_value;
    	 basisValue(q_point, EP, basis_value);
	  
	  for(int j = 0; j < q_point.size(); j++){
	      double Jxw = GW[j]*jacobian[j];

	      for(int k = 0; k < basis_value.size(); k++){
	         for(int m = 0; m < basis_value.size(); m++){
	             A_matrix[k][m] += Jxw*basis_value[m][j]*basis_value[k][j];
	             B_matrix[k][m] += Jxw*basis_value[m][j]*basis_value[k][j];
             }
              LRHS[k]=RHS0(ele_dof*i+k);
              LRHS2[k]=RHS0(v_size+ele_dof*i+k);
	      }
	  }
      GaussElimination(A_matrix, LU, LRHS);
      GaussElimination(B_matrix, LU2, LRHS2);
      for(int j=0; j<ele_dof; j++){
      	  c1_h(ele_dof*i+j)=LU[j];
      	  c2_h(ele_dof*i+j)=LU2[j];
	  }
   }
}

void BH2D::writePoint()
{
   std::ofstream outdata;
//   outdata.open("x.txt");

    char FileName[64]="";
	if (ele_dof>6)
	{
		int fig=std::round(log(sqrt(n_ele)/4)/log(2));
		if (sqrt(n_ele)==64)
			sprintf(FileName, "xR.txt");
		else
			sprintf(FileName, "x%d.txt",fig);
	}
	else
	{
		int fig=std::round(log(sqrt(n_ele)/8)/log(2));
		if (sqrt(n_ele)==128)
			sprintf(FileName, "xR.txt");
		else
			sprintf(FileName, "x%d.txt",fig);
	}

	puts(FileName);
	outdata.open(FileName);
	
   int n_ele=mesh.n_element();
   int x_dim=x_Gnt.size();
   outdata<<x_dim<<"\t"<<0<<"\n";
   for(int j=0; j<x_dim; j++){
      outdata<<x_Gnt[j][0]<<"\t"<<x_Gnt[j][1]<<"\n";
   }
   outdata.close();
}

void BH2D::writeData()
{
    std::ofstream outdata;
    char FileName[64]="";
	if (ele_dof>6)
	{
		int fig=std::round(log(sqrt(n_ele)/4)/log(2));
		if (sqrt(n_ele)==64)
			sprintf(FileName, "uR.txt");
		else
			sprintf(FileName, "u%d.txt",fig);
	}
	else
	{
		int fig=std::round(log(sqrt(n_ele)/8)/log(2));
		if (sqrt(n_ele)==128)
			sprintf(FileName, "uR.txt");
		else
			sprintf(FileName, "u%d.txt",fig);
	}
		
//	sprintf(FileName, "u0.txt");
	puts(FileName);
	outdata.open(FileName);
	int u_dim=u_DG.size();
	outdata<<u_dim<<"\n";
	outdata.precision(16);
	for(int j=0; j<u_DG.size(); j++){
        outdata<<u_DG[j]<<"\n";
    }
    outdata.close();
}

void BH2D::writePointF()
{
   std::ofstream outdata;
   outdata.open("x.txt");
	
   int n_ele=mesh.n_element();
   int x_dim=x_Gnt.size();
   outdata<<x_dim<<"\t"<<0<<"\n";
   for(int j=0; j<x_dim; j++){
      outdata<<x_Gnt[j][0]<<"\t"<<x_Gnt[j][1]<<"\n";
   }
   outdata.close();
}

void BH2D::writeDataF(int &fig)
{
    std::ofstream outdata;
    char FileName[64]="";
	sprintf(FileName, "ue%d.txt",fig);
	puts(FileName);
	outdata.open(FileName);
	int u_dim=u_DG.size();
	outdata<<u_dim<<"\n";
	outdata.precision(16);
	for(int j=0; j<u_DG.size(); j++){
        outdata<<u_DG[j]<<"\n";
    }
    outdata.close();
}

void BH2D::writeDataC1(int &fig)
{
    std::ofstream outdata;
    char FileName[64]="";
	sprintf(FileName, "c1e%d.txt",fig);
	puts(FileName);
	outdata.open(FileName);
	int u_dim=c1_DG.size();
	outdata<<u_dim<<"\n";
	outdata.precision(16);
	for(int j=0; j<c1_DG.size(); j++){
        outdata<<c1_DG[j]<<"\n";
    }
    outdata.close();
}

void BH2D::writeDataC2(int &fig)
{
    std::ofstream outdata;
    char FileName[64]="";
	sprintf(FileName, "c2e%d.txt",fig);
	puts(FileName);
	outdata.open(FileName);
	int u_dim=c2_DG.size();
	outdata<<u_dim<<"\n";
	outdata.precision(16);
	for(int j=0; j<c2_DG.size(); j++){
        outdata<<c2_DG[j]<<"\n";
    }
    outdata.close();
}

void BH2D::writePointT()
{
   std::ofstream outdata;
//   outdata.open("x.txt");

    char FileName[64]="";
	sprintf(FileName, "xxR.txt");

	puts(FileName);
	outdata.open(FileName);
	
   int n_ele=mesh.n_element();
   int x_dim=x_Gnt.size();
   outdata<<x_dim<<"\t"<<0<<"\n";
   for(int j=0; j<x_dim; j++){
      outdata<<x_Gnt[j][0]<<"\t"<<x_Gnt[j][1]<<"\n";
   }
   outdata.close();
}

void BH2D::writeDataT()
{
    std::ofstream outdata;
    char FileName[64]="";
    int fig=t_step-1;
    if (fig==4)
		sprintf(FileName, "uuR.txt");
	else
	    sprintf(FileName, "uu%d.txt",fig);

	puts(FileName);
	outdata.open(FileName);
	int u_dim=u_DG.size();
	outdata<<u_dim<<"\n";
	outdata.precision(16);
	for(int j=0; j<u_DG.size(); j++){
        outdata<<u_DG[j]<<"\n";
    }
    outdata.close();
}

void BH2D::writeDataC()
{
    std::ofstream outdata;
    char FileName[64]="";
	sprintf(FileName, "uh.txt");

	puts(FileName);
	outdata.open(FileName);
	int u_dim=u_h.size();
	outdata<<u_dim<<"\n";
	outdata.precision(16);
	for(int j=0; j<u_h.size(); j++){
        outdata<<u_h[j]<<"\n";
    }
    outdata.close();
}

void BH2D::writeEnergy()
{
   std::ofstream outdata;
   outdata.open("energy.txt");
   int n_eng=Energy.size();
   outdata<<n_eng<<"\t"<<0<<"\n";
   for(int j=0; j<n_eng; j++){
      outdata<<Energy[j][0]<<"\t"<<Energy[j][1]<<"\n";
   }
   outdata.close();
}

void BH2D::writeMassC1()
{
   std::ofstream outdata;
   outdata.open("massc1.txt");
   int n_mass=Mass.size();
   outdata<<n_mass<<"\t"<<0<<"\n";
   for(int j=0; j<n_mass; j++){
      outdata<<Mass[j][0]<<"\t"<<Mass[j][1]<<"\n";
   }
   outdata.close();
}

void BH2D::writeMassC2()
{
   std::ofstream outdata;
   outdata.open("massc2.txt");
   int n_mass=Mass.size();
   outdata<<n_mass<<"\t"<<0<<"\n";
   for(int j=0; j<n_mass; j++){
      outdata<<Mass[j][0]<<"\t"<<Mass[j][2]<<"\n";
   }
   outdata.close();
}

void BH2D::writeAvgC1()
{
   std::ofstream outdata;
   outdata.open("avgc1.txt");
   int n_mass=cAvg.size();
   outdata<<n_mass<<"\t"<<0<<"\n";
   for(int j=0; j<n_mass; j++){
      outdata<<cAvg[j][0]<<"\t"<<cAvg[j][1]<<"\n";
   }
   outdata.close();
}

void BH2D::writeAvgC2()
{
   std::ofstream outdata;
   outdata.open("avgc2.txt");
   int n_mass=cAvg.size();
   outdata<<n_mass<<"\t"<<0<<"\n";
   for(int j=0; j<n_mass; j++){
      outdata<<cAvg[j][0]<<"\t"<<cAvg[j][2]<<"\n";
   }
   outdata.close();
}

void BH2D::writeg1()
{
   std::ofstream outdata;
   outdata.open("ming1.txt");
   int n_mass=gMin.size();
   outdata<<n_mass<<"\t"<<0<<"\n";
   for(int j=0; j<n_mass; j++){
      outdata<<gMin[j][0]<<"\t"<<gMin[j][1]<<"\n";
   }
   outdata.close();
}

void BH2D::writeg2()
{
   std::ofstream outdata;
   outdata.open("ming2.txt");
   int n_mass=gMin.size();
   outdata<<n_mass<<"\t"<<0<<"\n";
   for(int j=0; j<n_mass; j++){
      outdata<<gMin[j][0]<<"\t"<<gMin[j][2]<<"\n";
   }
   outdata.close();
}

void BH2D::writeg1Lim()
{
   std::ofstream outdata;
   outdata.open("ming1lim.txt");
   int n_mass=gMin.size();
   outdata<<n_mass<<"\t"<<0<<"\n";
   for(int j=0; j<n_mass; j++){
      outdata<<gMin[j][0]<<"\t"<<gMin[j][3]<<"\n";
   }
   outdata.close();
}

void BH2D::writeg2Lim()
{
   std::ofstream outdata;
   outdata.open("ming2lim.txt");
   int n_mass=gMin.size();
   outdata<<n_mass<<"\t"<<0<<"\n";
   for(int j=0; j<n_mass; j++){
      outdata<<gMin[j][0]<<"\t"<<gMin[j][4]<<"\n";
   }
   outdata.close();
}

void BH2D::run()
{

	poly=2;
    ele_dof_Max = 6;      // This should be 1, 4, 9, 16, 25 for tensor product mesh.
	ele_dof0=ele_dof_Max;           // Don't change this. It is >= 1.
	ele_dof=ele_dof_Max;
	
	
	n_ele=mesh.n_element();
	v_size=ele_dof*n_ele;
	u_h.reinit(v_size);
	
	init();
    beta_0=16.0;           // This is not penalty, the one specially for u.
	beta0=beta_0;            // This is penalty, the one for u and q.
	beta1 = 1.0/6.0;       // beta1 for Poisson equation
	
	zeroflux = 0.0;       // if zeroflux=0.0, otherwize zeroflux = 1.0;
	LimTol = 1.0e-10;

    t=0.0;
    switch(ele_dof)
    {
        case 3:
		{
			dt = 1e-3;
			break;
		}
		case 4:
		{
			dt = 1e-6;
			break;
		}
		case 6:
		{
			dt = 1e-6;
			break;
		}
		case 9:
		{
			dt = 1e-6;
			break;
		}
		case 10:
		{
	//		dt = 1e-4;
			dt = 2e-5;
			break;
		}
		case 16:
		{
			dt = 5e-5;
			break;
		}
		case 25:
		{
			dt = 1e-5;
			break;
		}
		default:
		{
			dt = 5e-4;
			std::cout<<"The default dt is used."<<std::endl;
			break;
		}
	}
	
//  
/*
	if(n_ele == 900)
		dt = 2e-6;
	else if(n_ele == 100)
		dt = 2e-5;
	else if(n_ele <= 1600)
		dt = (1600/n_ele)*1.0e-6;
	else
// */	
	
	dt = 1.0e-6;
	double tend=1e-4;

    t_step = 0;    
// time order in t 
/*    
//   	dt=0.001;

    t_step = 6;   // t_step = 1, 2, 3, 4, 5
    
    double t_init=1.0;
	dt=t_init*pow(0.5,t_step-1);
//	dt=1.0e-4;
// */

    NT=std::round(tend/dt);
//    std::cout<<"NT= "<<dt<<std::endl;
    
	gStepMin.resize(4, 0.0);
    int EngTotal =std::min(NT,1000);
	int MassTotal =std::min(NT,100);
	int MinTotal = std::min(NT,100);
    Energy.resize(EngTotal+1);
	Mass.resize(MassTotal+1);
	gMin.resize(MinTotal+1);
	cAvg.resize(MinTotal+1);
    for(int i=0; i<Energy.size();i++)
        Energy[i].resize(2,0.0);	
	for(int i=0; i<Mass.size();i++)
	{
        Mass[i].resize(3,0.0);
		gMin[i].resize(5,0.0);
		cAvg[i].resize(3,0.0);
	}
	int Eng=0;   //count the step of Energy;
    int Fig=0;
	int Mas=0;
	int MinNO=0;
    
    L2Proj(&_c1_0_, &_c2_0_);
	double L2Errorc10=ComputeL1Error(c1_h, &_c1_0_);
	double L2Errorc20=ComputeL1Error(c2_h, &_c2_0_);
	std::cout<<std::setiosflags(std::ios::scientific)<<std::setprecision(5);
	std::cout<<"L1Errorc10= "<<L2Errorc10<<", L1Errorc20= "<<L2Errorc20<<std::endl;
	std::cout.unsetf(std::ios::scientific);
	
	q_i.resize(2);
	q_i[0]=1.0;
	q_i[1]=-1.0;
	
	g1_h.reinit(ele_dof*n_ele);
	g2_h.reinit(ele_dof*n_ele);
	c1_pre.reinit(ele_dof*n_ele);
	c1_pre.reinit(ele_dof*n_ele);
	
//for(int i=0;i<101;i++)
//	std::cout<<_phi_(-PI+0.02*i*PI)<<std::endl;

//    L2Proj(&_f_);
//	int fig=0;
//	Compute_U_Real(u_h,x_Gnt,u_DG);
//	writePointF();
//	writeDataF(fig);

//    L2ProjQ();
//    InitU(&_U_0_);
//	L2ProjU(&_U_0_);
		
//  Printing the L^2 error for each time step /*
//    double L2Error0=ComputeL2Error(u_h, &_f_);
//    double H1Error0=ComputeH1Error(u_h, &_u_grad_0_);
//    double LinfError0=ComputeLINFError(u_h, &_f_);
//    double LinfError0=ComputeLINFError(t, u_h, &_u_);
//    std::vector<double> L2Norm0=ComputeL2Norm(0.0, u_h, &_u_);
//    std::cout<<"L2Error0= "<<L2Error0<<", LinfError0= "<<LinfError0<<std::endl;
	
//    double L2ErrorU0=ComputeL2Error(U_h, &_U_0_);
//    double H1ErrorU0=ComputeH1Error(U_h, &_u_grad_0_);
//    double LinfErrorU0=ComputeLINFError(U_h, &_U_0_);
//	std::cout<<"t= "<<t<<", L2ErrorU0= "<<L2ErrorU0<<", LinfErrorU0= "<<LinfErrorU0<<std::endl;
	
//    Energy=ComeputeEnergy();
//	std::cout<<"Energy0= "<<Energy<<std::endl;

// */
// 
/*
    int Eng=0;   //count the step of Energy;
    int Fig=0;
    if(t==0.0)
	{
		 std::cout<<"---------Begin---------"<<std::endl;
         Compute_U_Real(u_h,x_Gnt,u_DG);
         Energy[Eng][0]=t;
         Energy[Eng][1]=ComeputeEnergy();
         std::cout<<"Energy[0]= "<<Energy[Eng][1]<<std::endl;
         writePointF();
         writeDataF(Fig);
         Fig++;
         Eng++;	
	
	}
// */

//  Printing the L^2 error of q_h for each time step 
/*
    double QL2Error0=ComputeL2Error(q_c, &_q_0_);
    double QH1Error0=ComputeH1Error(q_c, &_q_grad_0_);
    double QLinfError0=ComputeLINFError(q_c, &_q_0_);
    std::cout<<"t= "<<t<<", QL2Error0= "<<QL2Error0<<", QH1Error0= "<<QH1Error0<<", QLinfError0= "<<QLinfError0<<std::endl;
// */
// /*
//		std::cout<<"I am here! "<<std::endl;
//std::cout<<std::atan2(-0.5,-0.5)<<std::endl;
	buildSparsityPattern();
	buildSparseMatrix();
    beta1 = 1.0/6.0;   // beta1 for time dependent equation
	for(int step=0; step<NT+1; step++)
    {
		c1_pre = c1_h;
		c2_pre = c2_h;
		
		t=dt*(step);
		buildRHS();
		Solve();
		Solg();	
		Limiter();
		UpdateC();
		t=dt*(step+1);
		
//  Check error/* printing the solution. 

    if ( (step)==std::round(1.0/EngTotal*Eng*NT) && step<=NT )
	{
         Energy[Eng][0]=t-dt;
         Energy[Eng][1]=ComeputeEnergy();
//		 std::cout<<"---------Energy---------"<<std::endl;
         std::cout<<"Energy["<<Eng<<"]= "<<Energy[Eng][1]<<std::endl;
         Eng++;
	}
// */	

//  Check error /* printing the solution. 

    if ( (step)==std::round(1.0/MassTotal*Mas*NT) && step<=NT )
	{
         Mass[Mas][0]=t-dt;
         Mass[Mas][1]=ComeputeMass(c1_pre);
         Mass[Mas][2]=ComeputeMass(c2_pre);
//		 std::cout<<"---------Mass---------"<<std::endl;
         std::cout<<"c1Mass["<<Mas<<"]= "<<Mass[Mas][1]<<", ";
		 std::cout<<"c2Mass["<<Mas<<"]= "<<Mass[Mas][2]<<std::endl;
         Mas++;
	}
// */

//  Check error /* printing the solution. 

    if ( (step)==std::round(1.0/MinTotal*MinNO*NT) && step<=NT )
	{
         gMin[MinNO][0]=t-dt;
		 for(int lmin=0; lmin<4; lmin++)
			gMin[MinNO][lmin+1]=gStepMin[lmin];
		 cAvg[MinNO][0]=t-dt;
		 cAvg[MinNO][1]=getMinAvg(c1_pre);
		 cAvg[MinNO][2]=getMinAvg(c2_pre);
//		 std::cout<<"---------Mass---------"<<std::endl;
         std::cout<<"gMim["<<MinNO<<"]= "<<gMin[MinNO][1]<<", "<<gMin[MinNO][3]<<", "<<gMin[MinNO][2]<<", "<<gMin[MinNO][4]<<std::endl;
         std::cout<<"cAvg["<<MinNO<<"]= "<<cAvg[MinNO][1]<<", "<<cAvg[MinNO][2]<<std::endl;
		 MinNO++;
	}
// */

//    if(t==1.0/10.0*Fig*tend)
    double FigTol = 1.0e-7;
	double pt = t-dt;
//    if ( fabs(pt-0.0)<FigTol || fabs(pt-0.01)<FigTol || fabs(pt-0.1)<FigTol || fabs(pt-0.2)<FigTol || fabs(pt-0.5)<FigTol || fabs(pt-1.0)<FigTol )
//   Check error /*
    if ( (step)==std::round(1.0/10.0*Fig*NT) && step<=NT )
	{
		 std::cout<<"---------Begin---------"<<std::endl;
         Compute_U_Real(c1_pre,x_Gnt,c1_DG);
		 Compute_U_Real(c2_pre,x_Gnt,c2_DG);
		 Compute_U_Real(u_h,x_Gnt,u_DG);
		 writePointF();
		 writeDataF(Fig);
         writeDataC1(Fig);
		 writeDataC2(Fig);
         Fig++;
    }
// */
		
// showing the executaion place 
/* 
    if((step+1)%1==0){
//     std::cout<<"Steps= "<<step<<std::endl;
     double L2Error=ComputeL2Error(t-dt, u_h, &_u_);
	 double H1Error=ComputeH1Error(t-dt, u_h, &_u_grad_);
     double LinfError=ComputeLINFError(t-dt, u_h, &_u_);
//     std::vector<double> L2Norm=ComputeL2Norm(t, u_h, &_u_);
     std::cout<<"t= "<<t<<", L2Error= "<<L2Error<<", H1Error= "<<H1Error<<", LinfError= "<<LinfError<<std::endl;
	 
	}
// */

// showing the executaion place,  Check error /* 
    if((step+1)%1==0){
//     std::cout<<"Steps= "<<step<<std::endl;
     double L1ErrorC1=ComputeL1Error(t-dt, c1_pre, &_c1_);
	 double L1ErrorC2=ComputeL1Error(t-dt, c2_pre, &_c2_);
	 double L1ErrorPhi=ComputeL1Error(t-dt, u_h, &_u_);
//	 double H1Error=ComputeH1Error(t, c1_h, &_c1_grad_);
//     double LinfError=ComputeLINFError(t, c1_h, &_c1_);
//     std::vector<double> L2Norm=ComputeL2Norm(t, u_h, &_u_);
     std::cout<<"t= "<<t-dt<<", L1ErrorC1= "<<L1ErrorC1<<", L1ErrorC2= "<<L1ErrorC2<<", L1ErrorPhi= "<<L1ErrorPhi<<std::endl;
	}
 
// */

// showing the executaion place 
/* 
    if((step+1)%1==0){
//     std::cout<<"Steps= "<<step<<std::endl;
     double L2Error=ComputeL2Error(t, c2_h, &_c2_);
	 double H1Error=ComputeH1Error(t, c2_h, &_c2_grad_);
     double LinfError=ComputeLINFError(t, c2_h, &_c2_);
//     std::vector<double> L2Norm=ComputeL2Norm(t, u_h, &_u_);
     std::cout<<"t= "<<t<<", L2Error= "<<L2Error<<", H1Error= "<<H1Error<<", LinfError= "<<LinfError<<std::endl;
	 
	}
 
// */

// showing the executaion place 
/* 
    if((step+1)%1==0){
//     std::cout<<"Steps= "<<step<<std::endl;
     double L2Error=ComputeL2Error(g1_h, &_g_);
     double LinfError=ComputeLINFError(g1_h, &_g_);
//     std::vector<double> L2Norm=ComputeL2Norm(t, u_h, &_u_);
     std::cout<<"t= "<<t<<", L2Error= "<<L2Error<<", LinfError= "<<LinfError<<std::endl;
	 
	}
 
// */
		
	}
//	buildRHS();
//	Solve();
// Check error /*
//	writeDataC();
//  Compute_U_Real(u_h,x_Gnt,u_DG);
//	fig=1;
//	writeDataF(fig);
//	writePoint();
	writeEnergy();
	writeMassC1();
	writeMassC2();
	writeAvgC1();
	writeAvgC2();
	writeg1();
	writeg2();
	writeg1Lim();
	writeg2Lim();
	
// */
// Printing the L^1 error for each time step /*
     double L1ErrorC1=ComputeL1Error(t-dt, c1_pre, &_c1_);
	 double L1ErrorC2=ComputeL1Error(t-dt, c2_pre, &_c2_);
	 double L1ErrorPhi=ComputeL1Error(t-dt, u_h, &_u_);
//	 double H1Error=ComputeH1Error(t, u_h, &_u_grad_);
//     double LinfError=ComputeLINFError(t, u_h, &_u_);

	 std::cout<<std::setiosflags(std::ios::scientific)<<std::setprecision(5);
     std::cout<<" "<<L1ErrorC1<<"  "<<L1ErrorC2<<"  "<<L1ErrorPhi<<std::endl;
	 std::cout.unsetf(std::ios::scientific);

// */      

//     double L2ErrorU=ComputeL2Error(t, U_h, &_U_);
 //    double LinfErrorU=ComputeLINFError(t, U_h, &_U_);
 //    std::cout<<L2ErrorU<<"  "<<LinfErrorU<<std::endl;
	   
}



void BH2D::GaussSeidel(std::vector<std::vector<double> > &A_matrix, std::vector<double> &u, std::vector<double> &rhs)
{
   std::vector<double> last_u(A_matrix.size());  
   do{
      last_u=u; 
      double error=0.0;     
      for(int i=0;i<A_matrix.size();i++){
         double temp=rhs[i];
         for(int j=0; j<A_matrix[i].size(); j++){
             if(j != i){
                temp-=A_matrix[i][j]*u[j];
             }   
         }
         u[i]=temp/A_matrix[i][i];
      }
      u_int Vsize=u.size();
      for(int k=0;k<Vsize;k++){
         error+=(last_u[k]-u[k])*(last_u[k]-u[k]);
      }
      error=sqrt(error);
//      std::cout<<error<<std::endl;
//      std::getchar();
      if(error<1.0e-10)
          break;
   }while(1); 
}

void BH2D::GaussElimination(std::vector<std::vector<double> > &A_matrix, std::vector<double> &u, std::vector<double> &rhs)
{
   for(int s=0; s<rhs.size()-1; s++){
      CheckIsHaveSolution(A_matrix, s);
      Pivot(A_matrix, rhs, s);
      for(int i=s+1; i<rhs.size(); i++){
         double c=A_matrix[i][s]/A_matrix[s][s];
         rhs[i]-=rhs[s]*c;
         for(int j=s+1; j<rhs.size(); j++){
             A_matrix[i][j]-=A_matrix[s][j]*c;
         }
      }
   }
////////////////////////////////back replace;
   u_int Vsize=u.size();
   u[Vsize-1]=rhs[Vsize-1]/A_matrix[Vsize-1][Vsize-1];
   for(int i=Vsize-2; i>=0; i--){
      double S=0.0;
      for(int j=i+1; j<Vsize; j++){
         S+=A_matrix[i][j]*u[j];
      }
      u[i]=(rhs[i]-S)/A_matrix[i][i];
   }
}


void BH2D::NumberSwap(double &a, double &b)
{
   double tmp=a;
   a=b;
   b=tmp;
}

void BH2D::CheckIsHaveSolution(std::vector<std::vector<double> > &A_matrix, int k)
{
   int i=0;
   for(i=k; i<A_matrix.size(); i++){
      if(A_matrix[i][k] != 0.0)
        break;
   }
   if(i==A_matrix.size()){
      std::cout<<"No solution or no unique solution!"<<std::endl;
      std::cout<<"Please input CONTRAL+C to stop!"<<std::endl;
      getchar();
   } 
}

void BH2D::Pivot(std::vector<std::vector<double> > &A_matrix, std::vector<double> &Rhs, int k)
{
   int t=k;
   for(int j=k; j<A_matrix.size(); j++){
      if(fabs(A_matrix[j][k])>fabs(A_matrix[t][k])){
         t=j;
      }
   }
   if(t!=k){
      NumberSwap(Rhs[k],Rhs[t]);
      for(int j=k; j<A_matrix.size(); j++){
         NumberSwap(A_matrix[k][j], A_matrix[t][j]);
      }
   }
}

void BH2D::writeData_RHS()
{
	std::ofstream outdataRHS;
	outdataRHS.open("RHS.txt", std::ios::app);
	outdataRHS<<2*v_size<<"\n"<<"\n";
	for(int j=0; j<2*v_size; j++){
		outdataRHS<<RHS(j)<<"\n";
	}
	outdataRHS.close();
}

void BH2D::writeData_Mat()
{
	int n_ele=mesh.n_element();
	std::ofstream outdataMat;
	outdataMat.open("StiffMatrix.txt", std::ios::app);
	outdataMat<<n_ele*ele_dof*ele_dof*12<<"\n"<<"\n";

	//For Matrix A
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
			    outdataMat<<ele_dof*i+j<<"    "<<ele_dof*i+k<<"    "<<StiffMatrix.el(ele_dof*i+j, ele_dof*i+k)<<"\n";
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
					if(beta_0!=0.0)
					    outdataMat<<ele_dof*i+j<<"    "<<ele_dof*im+k<<"    "<<StiffMatrix.el(ele_dof*i+j, ele_dof*im+k)<<"\n";
				}
			} 
		}
	}
// /*
	//For Matrix B
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
			    outdataMat<<ele_dof*i+j<<"    "<<v_size+ele_dof*i+k<<"    "<<StiffMatrix.el(ele_dof*i+j, v_size+ele_dof*i+k)<<"\n";
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
				    outdataMat<<ele_dof*i+j<<"    "<<v_size+ele_dof*im+k<<"    "<<StiffMatrix.el(ele_dof*i+j, v_size+ele_dof*im+k)<<"\n";
				}
			} 
		}
	}
	//For Matrix C
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
			    outdataMat<<v_size+ele_dof*i+j<<"    "<<ele_dof*i+k<<"    "<<StiffMatrix.el(v_size+ele_dof*i+j, ele_dof*i+k)<<"\n";
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
				    outdataMat<<v_size+ele_dof*i+j<<"    "<<ele_dof*im+k<<"    "<<StiffMatrix.el(v_size+ele_dof*i+j, ele_dof*im+k)<<"\n";
				}
			} 
		}
	}
	//For Matrix D
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
			    outdataMat<<v_size+ele_dof*i+j<<"    "<<v_size+ele_dof*i+k<<"    "<<StiffMatrix.el(v_size+ele_dof*i+j, v_size+ele_dof*i+k)<<"\n";
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
//				    sparsity_pattern.add(v_size+ele_dof*i+j, v_size+ele_dof*im+k);
				}
			} 
		}
	}
	outdataMat.close();
}

void BH2D::fpwriteData_RHS()
{
	FILE *fp = NULL;
    fp = fopen("RHS.txt", "w");
    fprintf(fp, "%d\n", 2*v_size);
//    fprintf(fp, "%d %d %d\n", num_rows, num_cols, num_nonzeros);
    for (int i = 0; i < 2*v_size; i ++)
    {
        fprintf(fp, "%.15le\n", RHS(i));
    }
    fclose(fp);
}

void BH2D::fpwriteData_Mat()
{
	int n_ele=mesh.n_element();
	FILE *fp = NULL;
	fp = fopen("StiffMatrix.txt", "w");
	fprintf(fp, "%d\n", n_ele*ele_dof*ele_dof*12);

	//For Matrix A
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
				fprintf(fp, "%d %d %.15le\n", ele_dof*i+j, ele_dof*i+k, StiffMatrix.el(ele_dof*i+j, ele_dof*i+k));
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
					if(beta_0!=0.0)
					    fprintf(fp, "%d %d %.15le\n", ele_dof*i+j, ele_dof*im+k, StiffMatrix.el(ele_dof*i+j, ele_dof*im+k));
				}
			} 
		}
	}
// /*
	//For Matrix B
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
			    fprintf(fp, "%d %d %.15le\n", ele_dof*i+j, v_size+ele_dof*i+k, StiffMatrix.el(ele_dof*i+j, v_size+ele_dof*i+k));
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
				    fprintf(fp, "%d %d %.15le\n", ele_dof*i+j, v_size+ele_dof*im+k, StiffMatrix.el(ele_dof*i+j, v_size+ele_dof*im+k));
				}
			} 
		}
	}
	//For Matrix C
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
			    fprintf(fp, "%d %d %.15le\n", v_size+ele_dof*i+j, ele_dof*i+k, StiffMatrix.el(v_size+ele_dof*i+j, ele_dof*i+k));
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
				    fprintf(fp, "%d %d %.15le\n", v_size+ele_dof*i+j, ele_dof*im+k, StiffMatrix.el(v_size+ele_dof*i+j, ele_dof*im+k));
				}
			} 
		}
	}
	//For Matrix D
	for(int i=0; i<n_ele; i++){
        for(int j=0; j<ele_dof; j++){
			for(int k=0; k<ele_dof; k++){
			    fprintf(fp, "%d %d %.15le\n", v_size+ele_dof*i+j, v_size+ele_dof*i+k, StiffMatrix.el(v_size+ele_dof*i+j, v_size+ele_dof*i+k));
			}
		}		
		for(int m=0; m<element_patch[i].size(); m++){
			int im=element_patch[i][m];
			for(int j=0; j<ele_dof; j++){
			    for(int k=0; k<ele_dof; k++){
//				    sparsity_pattern.add(v_size+ele_dof*i+j, v_size+ele_dof*im+k);
				}
			} 
		}
	}
	fclose(fp);
}
