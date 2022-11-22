#include <iostream>
#include <fstream>
#include "cuba.h"
#include <cmath>
#include <Eigen/Dense>
#include "spline.h"
#include "LHAPDF/PDF.h"

using namespace std;
using namespace Eigen;
using namespace LHAPDF;
/*****************************begin const****************************************/
#define M_PI 3.14159265358979323846
double sq(double a){return pow(a,2);}
double tao,Q2,W,y,qmax,maxr,factor,alphaem,Nc,sumeq2,CA,CF,A,Sperp,x0,lambda,Qs2,qp,forL,forT,xB,Lambd;
int nplot;
double taomin,taobin;
double a[2][500];
tk::spline Fmv;
const PDF* FrF= mkPDF("JAM20-SIDIS_FF_hadron_nlo");
void constants() {
    x0=2.47e-5;
    lambda=0.282;
    Lambd=0.241;
    Sperp=23.58;
    CA=3.0;
    CF=4.0/3;
    Nc=3.0;
    alphaem=1.0/127;
    Q2=1.0;
    W=90;
    xB=Q2/sq(W);
    nplot=100;
    y= 0.8;
    qmax=100.0;
    maxr=100.0;
    taomin=-2.0;
    taobin=4.0/(nplot-1);
    sumeq2=2.0/3;
    Qs2=pow(x0/xB,lambda)*A;
    qp=sqrt(Q2/2/xB);
    factor=2*Nc*alphaem*sumeq2;//2*Nc*alphaem*sumeq*Sperp;
    forL=alphaem/xB/y/M_PI*(1.0-y);
    forT=alphaem/xB/y/M_PI*(1.0+sq(1.0-y))/2.0;
}

void dushujv(){
    ifstream inFile;
    inFile.open("/home/matata/CLionProjects/GTMD/MVn.txt");
    for(int i=0;i<2;i++){
        for(int j=0;j<500;j++){
            inFile >> a[i][j];
        }
    }
}
void Fmvcs(){
    vector<double> X(500), Y(500);
    int good=0;
    for(int i=0;i<500;i++)
    {
        if (good==0) {
            X[i]=a[0][i];
            Y[i] = a[1][i];
            if (Y[i]<1e-7) {
                Y[i] = 0.0;
                good=1;
            }
        }
        else{
            X[i]=a[0][i];
            Y[i] = 0.0;
        }
    }
    Fmv.set_points(X,Y);
    Fmv.make_monotonic();
}
/*****************************end const****************************************/

/*****************************hard function****************************************/

double HT(double k, double q, double xi){
    double epsilon=xi*(1.0-xi)*Q2,ans;
    double k2=sq(k);
    double l2=sq(q);
    ans=k2/sq(k2+epsilon)-(1.0+(k2-l2-epsilon)/sqrt(sq(k2+epsilon+l2)-4.0*k2*l2))/(k2+epsilon)+
        ((k2+l2)*(k2+l2+epsilon)-4.0*k2*l2)/pow(sq(k2+epsilon+l2)-4.0*k2*l2,3.0/2.0);
    return (sq(xi)+sq(1.0-xi))*ans;
}
double HL(double k, double q, double xi){
    double epsilon=xi*(1.0-xi)*Q2,ans;
    double k2=sq(k);
    double l2=sq(q);
    ans=1.0/sq(k2+epsilon)-2.0/(k2+epsilon)/sqrt(sq(k2+epsilon+l2)-4*k2*l2)+(k2+l2+epsilon)/pow(sq(k2+epsilon+l2)-4.0*k2*l2,3.0/2.0);
    return 4.0*sq(xi)*sq(1.0-xi)*ans*Q2;
}

/*****************************end hard function****************************************/

/*****************************general pdf****************************************/
double F1(Vector2d q){
    return exp(-q.dot(q)/Qs2)/M_PI/Qs2;
}

/*****************************end general pdf****************************************/
/***************************** FF ****************************************/
double FF(double z){
    return ((FrF->xfxQ2(1, z,Q2)+FrF->xfxQ2(-1, z,Q2))*8/9
            +(FrF->xfxQ2(2, z,Q2)+FrF->xfxQ2(-2, z,Q2))/9)/pow(z,3);
}
/*****************************end FF ****************************************/
double integrandL(double xi,double k,double q){
    return Fmv(q)*HL(k,q,xi);
}

double integrandT(double xi,double k,double q){
    return Fmv(q)*HT(k,q,xi);
}

static int dsigmaL(const int *ndim, const cubareal xx[],
                   const int *ncomp, cubareal ff[], void *userdata) {
    double k,q;
    k=tao;
    q=qmax*xx[0];
    ff[0]=sq(qmax)*xx[0]*integrandL(0.5,k,q);
    return 0;
}


static int dsigmaT(const int *ndim, const cubareal xx[],
                   const int *ncomp, cubareal ff[], void *userdata) {
    double k,q;
    k=tao;
    q=qmax*xx[0];
    ff[0]=sq(qmax)*xx[0]*integrandT(0.5,k,q);
    return 0;
}

static int dsigmaLT(const int *ndim, const cubareal xx[],
                    const int *ncomp, cubareal ff[], void *userdata) {
    //xi=xx[0], r=xx[1]*maxr
    ff[0]=4*sq(xx[0]*(1.0-xx[0]))*Q2*cyl_bessel_k(0.0,maxr*xx[1]*sqrt(Q2*xx[0]*(1.0-xx[0])))
          *cyl_bessel_k(0.0,maxr*xx[1]*sqrt(Q2*xx[0]*(1.0-xx[0])))*(1-exp(-CA/4.0/CF*Qs2*sq(xx[1]*maxr)*log(1.0/Lambd/maxr/xx[1]+exp(1))));
    return 0;
}
static int dsigmaTT(const int *ndim, const cubareal xx[],
                    const int *ncomp, cubareal ff[], void *userdata) {
    //xi=xx[0], r=xx[1]*maxr
    ff[0]=(sq(xx[0])+sq(1.0-xx[0]))*Q2*xx[0]*(1.0-xx[0])*cyl_bessel_k(1.0,maxr*xx[1]*sqrt(Q2*xx[0]*(1.0-xx[0])))
          *cyl_bessel_k(1.0,maxr*xx[1]*sqrt(Q2*xx[0]*(1.0-xx[0])))*(1-exp(-CA/4.0/CF*Qs2*sq(xx[1]*maxr)*log(1.0/Lambd/maxr/xx[1]+exp(1))));
    return 0;
}
int main() {
    ofstream mycout("/home/matata/CLionProjects/GTMD/SIDISEEC2NOEEC.txt");
    constants();
    dushujv();
    Fmvcs();
    double xT[nplot],yT[nplot];
    double xL[nplot],yL[nplot];
    double normT[nplot],normL[nplot];
    double sigmaT,sigmaL;
    int neval, fail;
    cubareal integral[1], error[1], prob[1];
    Vegas(2, 1, dsigmaTT, nullptr, 1,
          1e-3, 1e-12, 1, 0,
          1e2, 1e8, 1e4, 1e4, 1000,
          0, nullptr, nullptr,
          &neval, &fail, integral, error, prob);
    sigmaT=integral[0];
    Vegas(2, 1, dsigmaLT, nullptr, 1,
          1e-3, 1e-12, 1, 0,
          1e2, 1e8, 1e4, 1e4, 1000,
          0, nullptr, nullptr,
          &neval, &fail, integral, error, prob);
    sigmaL=integral[0];
    for (int i=0;i<nplot;i+=1) {
        tao = pow(10,taomin+taobin*i);
        Vegas(4, 1, dsigmaT, nullptr, 1,
              1e-2, 1e-12, 1, 0,
              1e2, 1e8, 1e4, 1e4, 1000,
              0, nullptr, nullptr,
              &neval, &fail, integral, error, prob);
        xT[i] = factor*integral[0]/2;
        yT[i] = factor*error[0]/2;
        normT[i] = integral[0]/sigmaT;
    }
    for (int i=0;i<nplot;i+=1) {
        tao = pow(10,taomin+taobin*i);
        Vegas(4, 1, dsigmaL, nullptr, 1,
              1e-2, 1e-12, 1, 0,
              1e2, 1e8, 1e4, 1e4, 1000,
              0, nullptr, nullptr,
              &neval, &fail, integral, error, prob);
        xL[i] = factor*integral[0]/2;
        yL[i] = factor*error[0]/2;
        normL[i] = integral[0]/sigmaL;
    }
    mycout<<"tao2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<pow(10,taomin+taobin*i)<<',';
    }
    mycout<<']'<<endl;
    mycout<<"xT2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<xT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"yT2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<yT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"xL2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<xL[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"yL2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<yL[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"tot2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<forL*xL[i]+forT*xT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"nxT2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<normT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"nxL2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<normL[i]<<',';
    }
    mycout<<']'<<endl;

    mycout<<"ntot2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<(forL*xL[i]+forT*xT[i])/(forL*sigmaL+forT*sigmaT)<<',';
    }
    mycout<<']'<<endl;
    mycout.close();
    cout<<sigmaT<<endl;
    cout<<sigmaL<<endl;
    cout<<Qs2;
    return 0;
}
//
// Created by matata on 3/4/22.
//

