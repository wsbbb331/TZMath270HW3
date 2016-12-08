//
//  KrylovSolver.h
//  Math270AHW2
//
//  Created by Tianyang Zhang on 12/6/16.
//  Copyright © 2016 UCLA. All rights reserved.
//

#ifndef KrylovSolver_h
#define KrylovSolver_h
#include <Eigen/Dense>
#include <Eigen/Core>

namespace JIXIE {
    template <class T>
    class KrylovSolver{
    public:
        typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TVect;
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TMat;
        inline void arnoldiIteration(const TVect b, const int n, TMat& Q, const TMat& A, TMat& h){
            int m = b.size();
//            std::cout<<"Q is "<<std::endl<<Q<<std::endl;
            TVect v=A * Q.col(n);
            
            for (int j=0;j<=n;j++){
                h(j,n) = Q.col(j).transpose()*v;
                v = v - h(j,n) * Q.col(j);
            }
            h(n+1,n)=v.norm();
//            std::cout<<"Q.col(n+1) is "<<std::endl<<Q.col(n+1)<<std::endl;
            Q.col(n+1)=v/h(n+1,n);
//            for (int i=0;i<m;i++){
//                std::cout<<"v(i) is "<<std::endl<<v(i)<<std::endl;
//                std::cout<<"h(n+1,n) is "<<std::endl<<h(n+1,n)<<std::endl;
//                std::cout<<"Q(i,n+1) is "<<std::endl<<Q(i,n+1)<<std::endl;
//                std::cout<<"Q.cols() is "<<std::endl<<Q.cols()<<std::endl;
//                std::cout<<"Q.rows() is "<<std::endl<<Q.rows()<<std::endl;
//                std::cout<<"n is "<<std::endl<<n<<std::endl;
//                std::cout<<"i is "<<std::endl<<i<<std::endl;
//                Q(i,n+1) = v(i)/h(n+1,n);
//            }
//            std::cout<<"Q is "<<std::endl<<Q<<std::endl;
        }
        
        inline void computeGMRES(const TVect b, const TMat& A, TVect& x){
            TMat Q,h;
//            std::cout<<"A is "<<std::endl<<A<<std::endl;
            int m = b.size();
            Q.resize(m,m+1);
            h.resize(m+1,m);
            Q.setZero();
            h.setZero();
            for(int i=0;i<m;i++){
                Q(i,0) = b(i)/(b.norm());
            }
//            std::cout<<"b norm is "<<std::endl<<b.norm()<<std::endl;
//            std::cout<<"b is "<<std::endl<<b<<std::endl;
//            std::cout<<"Q is "<<std::endl<<Q<<std::endl;
            for(int i=0;i<m;i++){
                arnoldiIteration(b, i, Q, A, h);
                TMat hTemp = h.block(0,0,i+2,i+1);
                TMat QTemp = Q.block(0,0,m,i+1);
                TVect newb(i+2);
                newb.setZero();
                newb(0) = b.norm();
//                std::cout<<"hTemp.transpose*htemp is "<<std::endl<<hTemp.transpose()*hTemp<<std::endl;
//                std::cout<<"hTemp.transpose*newb is "<<std::endl<<hTemp.transpose()*newb<<std::endl;
//                TVect y = hTemp.colPivHouseholderQr().solve(newb);
//                TVect y = hTemp.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(newb);
                TVect y = (hTemp.transpose() * hTemp).ldlt().solve(hTemp.transpose() * newb);
                
//                std::cout<<"hTemp is "<<std::endl<<hTemp<<std::endl;
//                std::cout<<"y is "<<std::endl<<y<<std::endl;
//                Eigen::HouseholderQR<TMat> qr(hTemp);
//                TMat hQ0 = qr.householderQ();
//                TMat hQ = hQ0.block(0,0,i+2,i+1);
//                std::cout<<"hQ is "<<std::endl<<hQ<<std::endl;
//                TMat R = qr.matrixQR();
//                std::cout<<"R is "<<std::endl<<R<<std::endl;
//                getUpperTriangle(R);
//                std::cout<<"R is "<<std::endl<<R<<std::endl;

//                std::cout<<"newb is "<<std::endl<<newb<<std::endl;
//                TVect d=hQ.transpose()*newb;
//                std::cout<<"d is "<<std::endl<<d<<std::endl;
//                TVect y(d);
//                backSubstitution(R, d, y);
//                std::cout<<"y is "<<std::endl<<y<<std::endl;
                x = QTemp*y;
//                std::cout<<"x is "<<std::endl<<x<<std::endl;
            }
        }
        
        inline void backSubstitution(const TMat& A, const TVect& b, TVect& x){
            int n = A.rows();
            for (int i = n-1; i >= 0; --i) {
                T s = 0;
                for (int j = i + 1; j < A.cols(); ++j) s = s + A(i,j)*x(j);
                x(i) = (b(i)-s)/A(i,i);
            }
        }
        inline void getUpperTriangle(TMat& R){
            int n=R.cols();
            for(int j=0;j<R.cols();j++){
                for(int i=j+1;i<R.rows();i++){
                    R(i,j) = 0;
                }
            }
            TMat RTemp = R.block(0,0,n,n);
            R.resize(n,n);
            R = RTemp;
        }
    };
}

#endif /* KrylovSolver_h */
