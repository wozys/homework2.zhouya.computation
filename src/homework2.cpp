#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

inline arma::mat soft_threshD(const arma::mat &X,const double &tao){
  arma::mat U,U1;
  arma::vec s;
  arma::mat V,V1;
  arma::mat D;
  arma::mat Dtao;
  svd(U1,s,V1,X);
  int nx=X.n_rows;
  int mx=X.n_cols;
  int n=s.n_elem;
  U=U1.submat(0,0,nx-1,n-1);
  V=V1.submat(0,0,mx-1,n-1);
  for(int i=0;i<=n-1;i++){
    if(s(i)-tao>=0)
      s(i)=s(i)-tao;
    else
      s(i)=0;
  }
  D=diagmat(s);
  Dtao=U*D*V.t();
  return Dtao;
}

inline arma::mat Pomega(const arma::mat &Omega,const arma::mat &X){
  int n=Omega.n_rows;
  int nx=X.n_rows,mx=X.n_cols;
  arma::mat Y=zeros<mat>(nx,mx);
  for (int i=0;i<=n-1;i++){
    Y(Omega(i,0)-1,Omega(i,1)-1)=X(Omega(i,0)-1,Omega(i,1)-1);
  }
  return Y;
}

// [[Rcpp::export]]
arma::mat homework2(arma::mat X,arma::mat Omega, double lambda){
  arma::mat Xk,Yk,oldXk,Dinitial;
  arma::mat U1,U,V1,V;
  arma::vec s;
  int n=X.n_rows;
  int m=X.n_cols;
  int deltak=1;
  svd(U1,s,V1,X);
  int l=s.n_elem;
  U=U1.submat(0,0,n-1,l-1);
  V=V1.submat(0,0,m-1,l-1);
  
  for(int i=0;i<=l-1;i++){
    if(s(i)-2*lambda*deltak>=0)
      s(i)=s(i)-2*lambda*deltak;
    else
      s(i)=0;
  }
  Dinitial=diagmat(s);
  //迭代初值
  Xk=U*Dinitial*V.t();

  int k=0;
  double tk=1,tk1=1;
  arma::mat Zk=Xk;
  //最大迭代次数默认为1000
  for(int i=1;i<=1000;i++){
    k++;
    oldXk=Xk;
    Yk=Zk+deltak*Pomega(Omega,X-Zk);
    Xk=soft_threshD(Yk,2*lambda*deltak);
    tk1=0.5*(1+sqrt(1+4*tk*tk));
    Zk=Xk+((tk-1)/tk1)*(Xk-oldXk);
    tk=tk1;
    if(norm(Xk-oldXk,"fro")<=0.1)
      break;
  }
  //把迭代后的结果处理，使得满足基本条件
  for(int j=0;j<=m-1;j++){
    for(int i=0;i<=n-1;i++){
      if(Xk(i,j)<=0)
        Xk(i,j)=0;
      else{
        if(Xk(i,j)>1){
          Xk(i,j)=1;
        }
      }
    }
  }
  //显示迭代次数
  std::cout<<k<<std::endl;
  return Xk;
}
