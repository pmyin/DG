/***********************************
 **@file BH2D.cpp
 **@author Peimeng Yin
 **date 2017 4 3
 **
 **brief   This is an example in 2d.
 **
 **********************************/
#include <algorithm>
#include <fstream>
#include <iterator>
#include <numeric>
#include <set>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdarg.h>

#include <cmath>
#include <time.h>

#include <lac/sparse_decomposition.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <lac/precondition_selector.h>
#include <lac/sparse_ilu.h>
#include <lac/sparse_mic.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_gmres.h>
#include <lac/solver_cg.h>
#include <lac/full_matrix.h>
// #include <lac/compressed_sparsity_pattern.h>
#include <lac/solver_richardson.h>


using namespace dealii;


class MyPoint;

MyPoint midpoint(const MyPoint&, const MyPoint&);
double distance(const MyPoint&, const MyPoint&);
MyPoint barycenter(const std::vector<MyPoint>&, const double * = NULL);
MyPoint operator+(const MyPoint&, const MyPoint&);
MyPoint operator-(const MyPoint&, const MyPoint&);
std::istream& operator>>(std::istream&, MyPoint&);
std::ostream& operator<<(std::ostream&, const MyPoint&);

class MyPoint
{
 private:
  double x[2]; /**< Coordinate of the point. */
 public:
  MyPoint(); /**< Default constructor. */
  MyPoint(const double *); /**< Constructor with data from a double array. */
  MyPoint(const MyPoint&); /**< Copy constructor. */
  MyPoint(double, ...); /**< Constructor taking parameters as the entries of the coordinate. */
  ~MyPoint(); /**< Destructor. */
 public:
  MyPoint& operator=(const MyPoint&); /**< Copy a point */
  operator const double *() const; /**< Casting to double pointer. */
  operator double *(); /**< Casting to double pointer. */
  const double& operator[](int) const; /**< Acess to entry. */
  double& operator[](int); /**< Acess to entry. */
  double length() const; /**< Length of the vector from the origin to the point. */
  MyPoint& operator+=(const MyPoint&);
  MyPoint& operator-=(const MyPoint&);
  MyPoint& operator*=(const double&);
  MyPoint& operator/=(const double&);
 public:
  friend MyPoint midpoint(const MyPoint&, const MyPoint&); /**< The middle point of two points. */
  friend double distance(const MyPoint&, const MyPoint&); /**< The middle point of two points. */
  friend MyPoint barycenter(const std::vector<MyPoint>&, const double *);
  friend MyPoint operator+ (const MyPoint&, const MyPoint&); /**< Add the coordinate of two point together. */
  friend MyPoint operator- (const MyPoint&, const MyPoint&); /**< Minus the coordinate. */
  friend std::istream& operator>> (std::istream&, MyPoint&); /**< Stream input. */
  friend std::ostream& operator<< (std::ostream&, const MyPoint&); /**< Stream output. */
};

class TmpEle
{
 private:
         std::vector<MyPoint>  Pnt;
         std::vector<MyPoint>  GaussPnt;
         std::vector<double>  GaussWeight;
         
         std::vector<MyPoint>  MeshPnt;
         
 
 public:
         TmpEle(){};
         virtual ~TmpEle(){};
         
 public:
        void buildTE(int i); 

	MyPoint getPnt(int );

	double getVolume();

	std::vector<MyPoint> getGaussPnt();

	std::vector<double> getGaussWeight();
	std::vector<MyPoint> getMeshPnt();
        
	MyPoint Local_to_Global(const MyPoint& , const std::vector<MyPoint> &) const; 
	MyPoint Global_to_Local(const MyPoint& , const std::vector<MyPoint> &) const; 
	double Local_to_Global_jacobian(const MyPoint& , const std::vector<MyPoint> &) const; 
	double Global_to_Local_jacobian(const MyPoint& , const std::vector<MyPoint> &) const; 
	std::vector<MyPoint> Local_to_Global(const std::vector<MyPoint>& , const std::vector<MyPoint> &) const; 
	std::vector<MyPoint> Global_to_Local(const std::vector<MyPoint>& , const std::vector<MyPoint> &) const;
	std::vector<double> Local_to_Global_jacobian(const std::vector<MyPoint>& , const std::vector<MyPoint> &) const;
	std::vector<double> Global_to_Local_jacobian(const std::vector<MyPoint>& , const std::vector<MyPoint> &) const;
};


class TmpEdge
{//参考单元
 private:
         std::vector<double>  Pnt;//参考单元的顶点
         std::vector<double>  GaussPnt;//参考单元上的高斯点
         std::vector<double>  GaussWeight;//参考单元上的高斯点对应的权系数
 
 public:
         TmpEdge(){};
         virtual ~TmpEdge(){};
         
 public:
        void buildTE(int i); //建立参考单元

	double getPnt(int );

	double getVolume();

	std::vector<double> getGaussPnt();

	std::vector<double> getGaussWeight();
        
	double Local_to_Global(const double& , const std::vector<double> &) const; //将参考单元上的点映射到实际单元
	double Global_to_Local(const double& , const std::vector<double> &) const; //将实际单元上的点映射到参考单元
	double Local_to_Global_jacobian(const double& , const std::vector<double> &) const; //从参考单元到实际单元的jacobian行列式
	double Global_to_Local_jacobian(const double& , const std::vector<double> &) const; //从实际单元到参考单元的jacobian行列式
	std::vector<double> Local_to_Global(const std::vector<double>& , const std::vector<double> &) const; 
	std::vector<double> Global_to_Local(const std::vector<double>& , const std::vector<double> &) const;
	std::vector<double> Local_to_Global_jacobian(const std::vector<double>& , const std::vector<double> &) const;
	std::vector<double> Global_to_Local_jacobian(const std::vector<double>& , const std::vector<double> &) const;
};

class Mesh
{
 private:
         std::vector<MyPoint>  Pnt;
         std::vector<std::vector<int> >  Edg;
         std::vector<std::vector<int> >  Ele;
 
 public:
         Mesh(){};
         virtual ~Mesh(){};
         
 public:
        void readData(const std::string& f); 

        int getEleVtx(int i, int j);
        std::vector<int> getEleVtx(int i);

        int getEleEdg(int i, int j);
        std::vector<int> getEleEdg(int i);

        int getEdgVtx(int i, int j);
        std::vector<int> getEdgVtx(int i);

	    MyPoint getPnt(int i);
    
        int n_edge();

        int n_element();
    
        int n_point();

	    std::vector<MyPoint> getPnt(std::vector<int> i);

		int getEleBndInfo(int i);

		int getEdgBndInfo(int i);
};

class BH2D
{//求解
 private:
  Mesh mesh;//网格
  TmpEle  tmpEle;//参考单元
  TmpEdge tmpEdge;
  Vector<double> ur_h;//有限元解函数
  Vector<double> ui_h;//有限元解函数
  Vector<double> u_h;//有限元解函数
  Vector<double> u_c;//有限元解函数
  Vector<double> q_c;//有限元解函数
  Vector<double> u_pre;
  Vector<double> u_pre_1;
  Vector<double> U_h;
  Vector<double> U_pre;
  Vector<double> c1_h;
  Vector<double> c2_h;
  Vector<double> c1_pre;
  Vector<double> c2_pre;
  Vector<double> c1_mid;
  Vector<double> c2_mid;
  Vector<double> g1_h;
  Vector<double> g2_h;
  Vector<double> g1_pre;
  Vector<double> g2_pre;
  Vector<double> g1_mid;
  Vector<double> g2_mid;
  
  std::vector<double> q_i;
	
	std::vector<std::vector<double> > u_star_R;
	std::vector<std::vector<double> > u_star_I;
	std::vector<std::vector<double> > u_dstar_R;
	std::vector<std::vector<double> > u_dstar_I;
	
	std::vector<std::vector<double> > u_theta_R;
	std::vector<std::vector<double> > u_theta_I;
	std::vector<std::vector<double> > u_n_R;
	std::vector<std::vector<double> > u_n_I;
	
	std::vector<std::vector<double> > u_starBnd_R;
	std::vector<std::vector<double> > u_starBnd_I;
	std::vector<std::vector<double> > u_dstarBnd_R;
	std::vector<std::vector<double> > u_dstarBnd_I;
	
	std::vector<std::vector<double> > u_thetaBnd_R;
	std::vector<std::vector<double> > u_thetaBnd_I;
	std::vector<std::vector<double> > u_nBnd_R;
	std::vector<std::vector<double> > u_nBnd_I;

  SparsityPattern  sparsity_pattern;
  SparseMatrix<double>   StiffMatrix;
  Vector<double> RHS;
  Vector<double> RHS0;
  
  std::vector<std::vector<int> > edge_patch;
  std::vector<std::vector<int> > NPedge_patch;
  std::vector<std::vector<int> > element_patch;
  std::vector<std::vector<int> > element_patch_edge;
  std::vector<int> edge_patch_edge;

  double t;
  double dt;
  double dtSys;
  int NT;
  int ele_dof;
  int ele_dof0;
  int ele_dof_Max;
  int v_size;
  double theta;
  double antitheta;
  double beta_0;
  double beta0;
  double beta0x;
  double beta0y;
  double beta1;
  int GaussPntNum;
  int gaussBndPntNum;
  int poly;
  double zeroflux;
  double LimTol;
  
  std::vector<std::vector<double>> Energy;
  std::vector<std::vector<double>> Mass, cAvg, gMin;
  std::vector<double> gStepMin;
  
  std::vector<MyPoint> x_Gnt;
  std::vector<double> u_DG;
  std::vector<double> c1_DG;  
  std::vector<double> c2_DG;    
    
  int n_ele;
  int t_step;


 public:
  BH2D(const std::string& f);
  virtual ~BH2D(){};

 public:
  void init();//初始化
  
  double power(double &f, int &n);
  double power(double f, int n);

  void basisValue(MyPoint &p, std::vector<MyPoint> &vtx, std::vector<double> & val);

  void basisGrad(MyPoint &p, std::vector<MyPoint> &vtx, std::vector<std::vector<double> > & val);
  
  void basis2Grad(MyPoint &p, std::vector<MyPoint> &v, std::vector<std::vector<double> > & val);

  void basisValue(std::vector<MyPoint> &p, std::vector<MyPoint> &vtx, std::vector<std::vector<double> > & val);

  void basisGrad(std::vector<MyPoint> &p, std::vector<MyPoint> &vtx, std::vector<std::vector<std::vector<double> > >& val);
  
  void basis2Grad(std::vector<MyPoint> &p, std::vector<MyPoint> &v, std::vector<std::vector<std::vector<double> > > & val);

  void Solve();

  void run();
  
  void step_Forward();

  void L2Proj(double (*gr)(const double *));
  void L2Proj(double (*Init1)(const double *), double (*Init2)(const double *));
  void Solg();
  void Solgtest();
  void L2ProjQ();
  void L2ProjU(double (*gr)(const double *));
  void L2ProjU();
  void InitU(double (*gr)(const double *));
  void UpdateC();
  void UpdateC2();
  void Limiter();
  std::vector<MyPoint> buildSamPnt(int i, int q_val);
  std::vector<double> getAvg(int n);
  std::vector<double> getMin(int n, std::vector<MyPoint> SP1, std::vector<MyPoint> SP2);
  double getMinAvg(Vector<double> &f);
  double getMin(int n, std::vector<MyPoint> Pnt, Vector<double> &f);
  
  void GaussSeidel(std::vector<std::vector<double> > &A_matrix, std::vector<double> &u, std::vector<double> &rhs);

  void GaussSeidelFS(SparseMatrix<double> &A_matrix, Vector<double> &u, Vector<double> &rhs);
  
  void GaussElimination(std::vector<std::vector<double> > &A_matrix, std::vector<double> &u, std::vector<double> &rhs);

  void NumberSwap(double &a, double &b);

  void CheckIsHaveSolution(std::vector<std::vector<double> > &A_matrix, int k);

  void Pivot(std::vector<std::vector<double> > &A_matrix, std::vector<double> &Rhs, int k);

//////////////////////////////////////////////////////////////////////////////////////////////
  void build_Edge_Patch();
  void build_Element_Patch();

  void buildSparsityPattern();
  void buildSparseMatrix();
  void buildSparseMatrixRight();
  void buildSparseMatrixUp();
  void buildSparseMatrixLeft();
  void buildSparseMatrixDown();
  void buildRHS();
  void buildRHS_C();
  void buildRHS_C2();
  void buildRHSRight();
  void buildRHSUp();
  void buildRHSLeft();
  void buildRHSDown();
  void getF();
  double mylog(double x);


	void getElementMatrix_A(int n, std::vector<std::vector<double> > &Ele_matrix);
	
	void getElementMatrix(int n, std::vector<std::vector< std::vector<double> > > &Ele_matrix);

    void getElementRHS(int n, std::vector< std::vector<double> > &Ele_Rhs);
	void getElementRHS(int n, std::vector<double> &Ele_Rhs);
	void getElementRHS_C(int n, std::vector< std::vector<double> > &Ele_Rhs);
	void getElementRHS_C2(int n, std::vector< std::vector<double> > &Ele_Rhs);
	void getElementRHS_g(int n, std::vector< std::vector<double> > &Ele_Rhs);
	
	void Read_u_Coef(int i, std::vector<double> &u_coef);
	void Read_q_Coef(int i, std::vector<double> &q_coef);

    void getElementF_I(int n, std::vector<double> &Ele_Rhs);

  double FEMFunctionValue(Vector<double> u_tmp, MyPoint p, int n);

  std::vector<double>  FEMFunctionValue(Vector<double> u_tmp, std::vector<MyPoint> p, int n);

//////////////////////////////////////////////////////////////////////////////////////////////

  double ComputerL2Error(double t, std::vector<std::vector<double> > &f, double (*g)(const double ,const double *));//计算L2范数误差

  double ComputerL2Error(std::vector<std::vector<double> > &f, double (*g)(const double *));//计算L2范数误差
  
  double ComputeL2Error(Vector<double> &f, double (*g)(const double *)); //计算L2范数误差
  double ComputeL1Error(Vector<double> &f, double (*g)(const double *)); //计算L1范数误差
  double ComputeH1Error(Vector<double> &f, std::vector<double> (*g)(const double *));  //计算H1范数误差
  double ComputeLINFError(Vector<double> &f, double (*g)(const double *)); //计算L infty范数误差
  
  double ComputeL2Error(double t, Vector<double> &f, double (*g)(const double ,const double *)); //计算L2范数误差
  double ComputeL1Error(double t, Vector<double> &f, double (*g)(const double ,const double *)); //计算L1范数误差
  double ComputeH1Error(double t, Vector<double> &f, std::vector<double> (*g)(const double ,const double *));  //计算H1范数误差
  double ComputeLINFError(double t, Vector<double> &f, double (*g)(const double ,const double *)); //计算L infty范数误差
  
  double ComputeSolAvg(Vector<double> &f);
  std::vector<double> ComputeL2Norm(double t, Vector<double> &f, double (*g)(const double ,const double *)); //计算L2范数
  
  double ComeputeEnergy();
  double getEnergy_g(int n);
  double ComeputeMass(Vector<double> &f);
  
  void Compute_U_Real(Vector<double> &f, std::vector<MyPoint> &x_gnt, std::vector<double> &U_Gnt);

  double ComputerLINFError(double t, std::vector<std::vector<double> > &f, double (*g)(const double ,const double *));//计算L infty范数误差

  double ComputerLINFError(std::vector<std::vector<double> > &f, double (*g)(const double *));//计算L infty范数误差

  double ComputerL2Norm(double t, double (*g)(const double ,const double *));//计算L2范数

  double ComputerL2Norm(double (*g)(const double *));//计算L2范数

  double ComputerL2Norm(std::vector<std::vector<double> > &f);//计算L2范数

/////////////////////////////////////////////////////数据输出

    void writeData_RHS();
    void writeData_Mat();
    void fpwriteData_RHS();
    void fpwriteData_Mat();
	
	void writeError();
	
	void writeData_u(double st);
	
	void writeData_U(double st);
	
	void writeMod_u(double st);
	
	void writeMod_U(double st);
	
	
	void writePoint();
	void writeData();
	void writeEnergy();
	void writeMassC1();
	void writeMassC2();
	void writeAvgC1();
	void writeAvgC2();
	void writeg1();
	void writeg2();
	void writeg1Lim();
	void writeg2Lim();
	
	void writePointT();
	void writeDataT();
	
	void writePointF();
	void writeDataF(int & fig);
	void writeDataC1(int & fig);
	void writeDataC2(int & fig);
	
	void writeDataC();
	

//////////////////////////////////////////////////////////////////////////////////////////////

    double calMod(std::vector<std::vector<double> > &uR_tmp, std::vector<std::vector<double> > &uI_tmp);
	
	void updata_u_star();

	void updata_u_dstar();

	void updata_u_theta();

	void updata_u_n();
	
	void updata_uBnd_star();

	void updata_uBnd_dstar();

	void updata_uBnd_theta();

	void updata_uBnd_n();
};

