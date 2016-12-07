#include <cmath>
#include "Tools.h"
#include "ImplicitQRSVD.h"
#include "SymmetricTridiagonal.h"
#include "SimulationDriver.h"
#include "EnergyTests.h"
#include "ThreeDOutput.h"
#include "KrylovSolver.h"

void EnergyTest(){
  typedef double T;
  typedef Eigen::Matrix<T,Eigen::Dynamic,1> TVect;
  int N=5;
  T a=(T)0,b=(T)1;
  T dX=(b-a)/(T)(N-1);
    T tb= 1;
    bool fixed_a=true;
  JIXIE::NeoHookean<T> nh((T)1);
  JIXIE::LinearElasticity<T> le((T)1);
  JIXIE::FEMHyperelasticity<T> fem(a,dX,N,nh);
  TVect x(N);
  for(int i=0;i<N;i++) x(i)=(T)1*(a+dX*(T)i);

  JIXIE::EnergyTest<T> et("output",fem,10, tb, fixed_a);
  et.RefinementTest(x);

}

void ElasticitySimulation(){


  typedef double T;
  typedef Eigen::Matrix<T,Eigen::Dynamic,1> TVect;

  JIXIE::ElasticityParameters<T> parameters;
  parameters.N=20;
  parameters.a=(T)0;
  T b=(T)1;
  parameters.dX=(b-parameters.a)/(T)(parameters.N-1);
  parameters.dt=(T).01;
  parameters.output_dir=std::string("output");
  parameters.rho=(T)1;
  parameters.k=(T)1;
  parameters.Newton_tol=(T)1e-8;
  parameters.max_newton_it=40;
  parameters.final_time=(T)10;
  parameters.frames_per_second=40;
    parameters.tb = 1;
    parameters.fixed_a = true;
  JIXIE::ElasticityDriver<T> driver(parameters);
  bool verbose=true;
  driver.RunSimulation(verbose);
    driver.writeDX();
}

void ConvertBinaryToDat(){
  typedef double T;
  typedef Eigen::Matrix<T,Eigen::Dynamic,1> TVect;
  TVect x,v; int N=0,frame=0;

  std::string data_dir("output");
  std::string output_dat_dir("output/matlab");

  while(JIXIE::ElasticityDriver<T>::Read_State(x,v,N,data_dir,frame)){
    char str[12];
    sprintf(str, "%d", frame++);
    std::string frame_name(str);
    std::string positions_string(std::string("particle_x_")+frame_name);
    std::string velocities_string(std::string("particle_v_")+frame_name);
    FILE_IO::Write_DAT_File(std::string(output_dat_dir+std::string("/")+positions_string+std::string(".dat")),x);
    FILE_IO::Write_DAT_File(std::string(output_dat_dir+std::string("/")+velocities_string+std::string(".dat")),x);
  }
}


void ConvertBinaryToObj(){
    typedef double T;
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> TVect;
    TVect x,v, F; int N=0,frame=0;T dX;T l=0.5, h=0.5;
    
    std::string data_dir("output");
    std::string output_dat_dir("output/matlab");
    
    JIXIE::ElasticityDriver<T>::readDX(dX, data_dir);
    std::cout<<"dX is "<<dX<<std::endl;
    while(JIXIE::ElasticityDriver<T>::Read_State(x,v,N,data_dir,frame)){
        char str[12];
        sprintf(str, "%d", frame++);
        std::string frame_name(str);
        std::string positions_string(std::string("particle_x_")+frame_name);
        FILE_IO::Write_DAT_File(std::string(output_dat_dir+std::string("/")+positions_string+std::string(".dat")),x);
        F.resize(N-1);
        JIXIE::ObjBody<T> objBody;
        for(int i=0;i<x.size()-1;i++){
            F(i) = (x(i+1)-x(i))/dX;
            T a = JIXIE::MATH_TOOLS::rsqrt(F(i));
            JIXIE::Cube<T> unitCube((x(i)+x(i+1))/2,0,0,1,x(i+1)-x(i),l*a,h*a);
            objBody.addCube(unitCube);
        }
        objBody.writeOutput("outputObject_"+frame_name);
    }
}




void testThreeDOutput(){
    typedef double T;
    JIXIE::Cube<T> unitCube1(-0.5,0,0,1,1,3,2);
    JIXIE::Cube<T> unitCube2(0.5,0,0,1,1,3,2);
    JIXIE::ObjBody<T> objBody;
//    JIXIE::Shape<T> shape(5.0,5.0,5.0,.5);
//    shape.insertObj(objBody);
    objBody.addCube(unitCube1);
    objBody.addCube(unitCube2);
    objBody.writeOutput("outputObject");
}


void testQR(){
    typedef double T;
    Eigen::Matrix<T,2,2> A;
    typedef Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> TVect;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TMat;
    TVect Q, R;
    TMat m(2,2);
    JIXIE::SymmetricTridiagonal<T> ST(A);
    ST.SetToZero();
    ST(1,1) = 1;
    ST(0,0) = 1;
    ST(0,1) = 2;
    ST.QR();
    ST.Set_Q(Q);
    ST.Set_R(R);
//    m.setZero();
//    m.resize(3,3);
//    m.row(2).setZero();
//    m.col(2).setZero();
    m<<1,2,2,1;
    Eigen::HouseholderQR<TMat> qr(m);
    TMat QQ(m);
    TMat RR = qr.matrixQR().triangularView<Eigen::Upper>();
    std::cout<<"Q is "<<Q;
    std::cout<<"R is "<<R;
    std::cout<<"m is "<<m;
    std::cout<<"Q is "<<QQ;
    std::cout<<"R is "<<RR;
}

int main()
{
//  EnergyTest();
  ElasticitySimulation();
  ConvertBinaryToDat();
//    ConvertBinaryToObj();
//    testThreeDOutput();
}
