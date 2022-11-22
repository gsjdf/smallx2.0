#include <iostream>
#include <fstream>
#include "cuba.h"
#include <cmath>


using namespace std;

/*****************************begin const****************************************/
double sq(double a){return pow(a,2);}
double tao,Q2,W,maxl,factor,alphaem,Nc,sumeq2,Sperp,qp,xB,mq2;
int nplot;
double taomin,taobin;
double ap[4100][3],axp[100][2],an[4100][3],axn[100][2];
void constants() {
    mq2=sq(0.14);
    Sperp=23.58;
    Nc=3.0;
    alphaem=1.0/137;
    Q2=10.0;
    W=90.0;
    xB=(Q2+mq2)/(sq(W));
    nplot=100;
    maxl=1000.0;
    taomin=-4.0;
    taobin=4.0/(nplot-1);
    sumeq2=2.0/3;
    qp=sqrt(Q2/2.0/xB);
    factor=Nc*alphaem*sumeq2;//2*Nc*alphaem*sumeq*Sperp;
}


void readlinep(){
    //read line
    ifstream inFile;
    inFile.open("/home/matata/CLionProjects/GTMD/MVP.txt");
    for(int i=0;i<4100;i++){
        for(int j=0;j<3;j++){
            inFile >> ap[i][j];
        }
    }
    //choose x
    for(int i=0;i<100;i++){
        axp[i][0]=log(ap[i][1]);//ap[i][1];
        double Y=log(0.01/xB);
        int nY=int(Y/0.2)*100;
        axp[i][1]=ap[nY+i][2]+(ap[nY+100+i][2]-ap[nY+i][2])*(Y-ap[nY+i][0])/(ap[nY+100+i][0]-ap[nY+i][0]);
    }
}
void readlinen(){
    //read line
    ifstream inFile;
    inFile.open("/home/matata/CLionProjects/GTMD/MVN.txt");
    for(int i=0;i<4100;i++){
        for(int j=0;j<3;j++){
            inFile >> an[i][j];
        }
    }
    //choose x
    for(int i=0;i<100;i++){
        axn[i][0]=log(an[i][1]);//an[i][1];
        double Y=log(0.01/xB);
        int nY=int(Y/0.2)*100;
        axn[i][1]=an[nY+i][2]+(an[nY+100+i][2]-an[nY+i][2])*(Y-an[nY+i][0])/(an[nY+100+i][0]-an[nY+i][0]);
    }
}
double Fmvp(double p){
    int ip = int(log(p/0.0001)/log(1.176811952435));
    if(ip<0)
    {return 0;}
    if(ip>=99)
    {return 0;}
    return axp[ip][1]+(axp[ip+1][1]-axp[ip][1])*(log(p)-axp[ip][0])/(axp[ip+1][0]-axp[ip][0]);
    //return axp[ip][1]+(axp[ip+1][1]-axp[ip][1])*(p-axp[ip][0])/(axp[ip+1][0]-axp[ip][0]);
}
double Fmvn(double p){
    int ip = int(log(p/0.0001)/log(1.176811952435));
    if(ip<0)
    {return 0;}
    if(ip>=99)
    {return 0;}
    return axn[ip][1]+(axn[ip+1][1]-axn[ip][1])*(log(p)-axn[ip][0])/(axn[ip+1][0]-axn[ip][0]);
    //return axn[ip][1]+(axn[ip+1][1]-axn[ip][1])*(p-axn[ip][0])/(axn[ip+1][0]-axn[ip][0]);
}

/*****************************end const****************************************/

/*****************************hard function****************************************/

double HT(double k, double q, double xi){
    double epsilon=xi*(1.0-xi)*Q2+mq2,ans;
    double k2=sq(k);
    double l2=sq(q);
    ans=k2/sq(k2+epsilon)-(1.0+(k2-l2-epsilon)/sqrt(sq(k2+epsilon+l2)-4.0*k2*l2))/(k2+epsilon)+
        ((k2+l2)*(k2+l2+epsilon)-4.0*k2*l2)/pow(sq(k2+epsilon+l2)-4.0*k2*l2,3.0/2.0);
    return (sq(xi)+sq(1.0-xi))*ans;
}
double HL(double k, double q, double xi){
    double epsilon=xi*(1.0-xi)*Q2+mq2,ans;
    double k2=sq(k);
    double l2=sq(q);
    ans=1.0/sq(k2+epsilon)-2.0/(k2+epsilon)/sqrt(sq(k2+epsilon+l2)-4*k2*l2)+(k2+l2+epsilon)/pow(sq(k2+epsilon+l2)-4.0*k2*l2,3.0/2.0);
    return 4.0*sq(xi)*sq(1.0-xi)*ans*Q2;
}

/*****************************end hard function****************************************/

static int dsigmaL(const int *ndim, const cubareal xx[],
                   const int *ncomp, cubareal ff[], void *userdata) {
    //k is k ,l is xx[0]*maxl, xi is xx[1],
    double k,l;
    k=xx[1]*qp*sqrt(tao);
    l=maxl*xx[0];
    double temhl=HL(k,l,xx[1]);
    ff[0]=maxl*l*pow(xx[1],3)*Fmvp(l)*temhl;
    ff[1]=maxl*l*pow(xx[1],3)*Fmvn(l)*temhl;
    return 0;
}


static int dsigmaT(const int *ndim, const cubareal xx[],
                   const int *ncomp, cubareal ff[], void *userdata) {
    double k,l;
    k=xx[1]*qp*sqrt(tao);
    l=maxl*xx[0];
    double temht=HT(k,l,xx[1]);
    ff[0]=maxl*l*pow(xx[1],3)*Fmvp(l)*temht;
    ff[1]=maxl*l*pow(xx[1],3)*Fmvn(l)*temht;
    return 0;
}


int main() {
    ofstream mycout("/home/matata/CLionProjects/GTMD/SIDISEEC1.txt");
    constants();
    readlinep();
    readlinen();
    double xT[nplot],yT[nplot],xTn[nplot],yTn[nplot];
    double xL[nplot],yL[nplot],xLn[nplot],yLn[nplot];
    int neval, fail;
    cubareal integral[2], error[2], prob[2];
    for (int i=0;i<nplot;i+=1) {
        tao = pow(10,taomin+taobin*i);
        Vegas(2, 2, dsigmaT, nullptr, 1,
              1e-2, 1e-12, 1, 0,
              1e2, 1e8, 1e4, 1e4, 1000,
              0, nullptr, nullptr,
              &neval, &fail, integral, error, prob);
        xT[i] = sq(qp)*factor*integral[0];
        yT[i] = sq(qp)*factor*error[0];
        xTn[i] = sq(qp)*factor*integral[1];
        yTn[i] = sq(qp)*factor*error[1];
    }
    for (int i=0;i<nplot;i+=1) {
        tao = pow(10,taomin+taobin*i);
        Vegas(2, 2, dsigmaL, nullptr, 1,
              1e-2, 1e-12, 1, 0,
              1e2, 1e8, 1e4, 1e4, 1000,
              0, nullptr, nullptr,
              &neval, &fail, integral, error, prob);
        xL[i] = sq(qp)*factor*integral[0];
        yL[i] = sq(qp)*factor*error[0];
        xLn[i] = sq(qp)*factor*integral[1];
        yLn[i] = sq(qp)*factor*error[1];
    }
    mycout<<"tao=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<pow(10,taomin+taobin*i)<<',';
    }
    mycout<<']'<<endl;
    mycout<<"xT=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<xT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"yT=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<yT[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"xL=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<xL[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"yL=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<yL[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"xT2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<xTn[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"yT2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<yTn[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"xL2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<xLn[i]<<',';
    }
    mycout<<']'<<endl;
    mycout<<"yL2=[";
    for (int i=0;i<nplot;i+=1){
        mycout<<yLn[i]<<',';
    }
    mycout<<']'<<endl;
    mycout.close();
    return 0;
}
//
// Created by matata on 3/4/22.
//

