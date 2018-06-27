//  main.cpp
//  Asian_Call_Control_m
//  Created by Dongliang Yi on 4/23/16.
//  Copyright © 2016 Dongliang Yi. All rights reserved.

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

double se;
double expectation;
double _K=2;
double _T=2.00;
double _r=0.05;
double _q=0;
double _sigma=0.5;
double _S0=2;
int _m ; // 50 stock price averaging
double inv_m;
double _deltaT;
double value[1000], control[1000];
double _b,_pho;
double Expectation_Geometric_Call;

// New methods on generating Unit RV!
#define m1 2147483647
#define m2 2145483479
#define a12 63308
#define a13 -183326
#define a21 86098
#define a23 -539608
#define q12 33921
#define q13 11714
#define q21 24919
#define q23 3976
#define r12 12979
#define r13 2883
#define r21 7417
#define r23 2071
#define Invmp1 4.656612873077393e-10
int x10,x11,x12,x20,x21,x22;

int Random()
{
    int h,p12,p13,p21,p23;
    h=x10/q13;p13=-a13*(x10-h*q13)-h*r13;
    h=x11/q12;p12=a12*(x11-h*q12)-h*r12;
    if (p13<0) {
        p13=p13+m1;
    }
    if (p12<0) {
        p12=p12+m1;
    }
    x10=x11;x11=x12;x12=p12-p13;
    if (x12<0) {
        x12=x12+m1;
    }
    
    h=x20/q23;p23=-a23*(x20-h*q23)-h*r23;
    h=x22/q21;p21=a21*(x22-h*q21)-h*r21;
    if (p23<0) {
        p23=p23+m2;
    }
    if (p21<0) {
        p21=p21+m2;
    }
    
    x20 =x21; x21 =x22; x22 =p21 - p23; if(x22 <0) x22 =x22 +m2;
    
    if (x12<x22)
    {
        return (x12-x22+m1);
    }
    else
        return (x12-x22);
}

double Uniform01()
{
    int Z;
    Z=Random();if(Z==0) Z=m1; return (Z*Invmp1);
}

double V[15]={1.253314137315500,0.6556795424187985,0.4213692292880545,0.3045902987101033,0.2366523829135607,0.1928081047153158,0.1623776608968675,0.1401041834530502,0.1231319632579329,0.1097872825783083,0.09902859647173193,0.09017567550106468,0.08276628650136917,0.0764757610162485,0.07106958053885211};
double c=0.918938533204672;

double Min(double a, double b) 
{
    return (b < a )? b:a;
}

double get_cdf(double x)
{
    if (x>15) {
        return 1;
    }
    if (x<-15) {
        return 0;
    }
    
    double y,a,b,q,s,h;
    int j,z;
    j=floor(Min((abs(x)+0.5), 14));
    z=j;
    h=abs(x)-z;
    a=V[j];
    b=z*a-1;
    q=1;
    s=a+h*b;
    for (int i=2; i<24-j; i=i+2) {
        a=(a+z*b)/i;
        b=(b+z*a)/(i+1);
        q=q*h*h;
        s=s+q*(a+h*b);
    }
    y=s*exp(-0.5*x*x-c);
    if (x>0) {
        y=1-y;
    }
    return y;
}

//Inverse Transform Method
double a[]={2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};
double b[]={-8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833};
double d[]={0.3374754822726147, 0.9761690190917186, 0.1607979714918209, 0.0276438810333863, 0.0038405729373609, 0.0003951896511919, 0.0000321767881768, 0.0000002888167364, 0.0000003960315187};

double Beasley_method(double u)
{
    double y = u-0.5;
    double r;
    double x;
    if (abs(y)<0.42)
    {
        r=y*y;
        x=y*(((a[3]*r+a[2])*r+a[1])*r+a[0])/((((b[3]*r+b[2])*r+b[1])*r+b[0])*r+1);
    }
    else
    {
        r=u;
        if (y>0) {
            r=1-u;
        }
        r=log(-log(r));
        x=d[0]+r*(d[1]+r*(d[2]+r*(d[3]+r*(d[4]+r*(d[5]+r*(d[6]+r*(d[7]+r*d[8])))))));
        
        if (y<0) {
            x=-x;
        }
    }
    return x;
}

double get_gaussian_inverse()
{
    double U=Uniform01();
    double x0=Beasley_method(U);
    double cdf = get_cdf(x0);
    double x1 = x0 - (cdf-U)/(exp(-pow(x0, 2)/2)/sqrt(2*3.141592654));
    return x1;
}

double max(double a, double b) {
    return (b < a )? a:b;
}

int monte_carlo_initial()
{
    double Price[_m];
    double RV = get_gaussian_inverse();
    Price[0]=_S0*exp((_r-_q-pow(_sigma,2)*0.5)*_deltaT*(1)+_sigma*sqrt(_deltaT*(1))*RV);
    for (int M = 1; M<_m; M++) {
        RV = get_gaussian_inverse();
        Price[M] = Price[M-1]*exp((_r-_q-pow(_sigma,2)*0.5)*_deltaT*(1)+_sigma*sqrt(_deltaT*(1))*RV);
    }
    
    double Inter_Price = pow(Price[0], inv_m);
    for (int M = 1; M<_m; M++) {
        Inter_Price = Inter_Price * pow(Price[M], inv_m);
    }
    Inter_Price = max(Inter_Price-_K, 0);
    control[0]=Inter_Price;
    
    double x = 0;
    for (int j = 0; j<_m; j++) {
        x += Price[j];    }
    x = x/_m;
    x = max(x-_K,0);
    value[0]=x;
    
    for (int i=1; i<1000; i++) {
        RV = get_gaussian_inverse();
        Price[0]=_S0*exp((_r-_q-pow(_sigma,2)*0.5)*_deltaT*(1)+_sigma*sqrt(_deltaT*(1))*RV);
        for (int M = 1; M<_m; M++) {
            RV = get_gaussian_inverse();
            Price[M] = Price[M-1]*exp((_r-_q-pow(_sigma,2)*0.5)*_deltaT*(1)+_sigma*sqrt(_deltaT*(1))*RV);
        }
        
        double Inter_Price = pow(Price[0],inv_m);
        for (int M = 1; M<_m; M++) {
            Inter_Price = Inter_Price * pow(Price[M],inv_m);
        }
        Inter_Price = max(Inter_Price - _K, 0);
        control[i]=Inter_Price;
        
        long double x = 0;
        for (int j = 0; j<_m; j++) {
            x += Price[j];    }
        x = x/_m;
        x = max(x-_K,0);
        value[i]=x;
    }
    return 1;
}

int b_Cal()
{
    long double X = 0;
    long double Y = 0;
    for (int i=0; i<1000; i++) {
        X+=control[i];
        Y+=value[i];
    }
    
    X=X/1000;
    Y=Y/1000;
    long double xx=0.0;
    long double xy=0.0;
    long double yy=0.0;
    
    for (int i=0; i<1000; i++) {
        long double temp=control[i]-X;
        long double temp2=value[i]-Y;
        xx+=pow(temp,2);
        xy+=temp2*temp;
        yy+=pow(temp2,2);
        
    }
    _b=xy/xx;
    _pho=xy/sqrt(xx*yy);
    return 1;
}

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility
                                       const double& time,   // time to maturity
                                       const double& q)
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+(r-q)*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return S*exp(-q*time)*get_cdf(d1) - K*exp(-r*time)*get_cdf(d2);
};

int Expectation_Geometric_Call_cal()
{
    double _r_ = _r;
    double _q_ = _q + 0.5*_sigma*_sigma-0.5*(2*_m+1)*_sigma*_sigma/(3*_m);
    double _sigma_ = sqrt((2*_m+1)*_sigma*_sigma/(3*_m));
    double _K_ = _K;
    double _T_ = 0.5*(_T+ _deltaT);
    double _S0_ =_S0;
    
    Expectation_Geometric_Call = exp(_r_*_T_)*option_price_call_black_scholes( _S0_, _K_, _r_, _sigma_, _T_, _q_);
    //cout << "Geometric: " << Expectation_Geometric_Call<<endl;
    //cout <<";;;;;;;;;;;;;;;;;;;;;;" <<endl;
    return 1;
}

int monte_carlo_inverse(int no_of_trials)
{
    double x,x2;
    double Price[_m];
    double RV = get_gaussian_inverse();
    Price[0]=_S0*exp((_r-_q-pow(_sigma,2)*0.5)*_deltaT*(1)+_sigma*sqrt(_deltaT*(1))*RV);
    for (int M = 1; M<_m; M++) {
        RV = get_gaussian_inverse();
        Price[M] = Price[M-1]*exp((_r-_q-pow(_sigma,2)*0.5)*_deltaT*(1)+_sigma*sqrt(_deltaT*(1))*RV);
    }
    
    double Inter_Price = pow(Price[0],inv_m);
    for (int M = 1; M<_m; M++) {
        Inter_Price = Inter_Price * pow(Price[M],inv_m);
    }
    Inter_Price = max(Inter_Price - _K, 0);
    
    x = 0;
    for (int j = 0; j<_m; j++) {
        x += Price[j];
        //cout << "Price " << j+1 << " : " << Price[j]<<endl;
    }
    x = x/_m;
    x = max(x-_K,0) + _b*(Expectation_Geometric_Call-Inter_Price);
    
    //cout << "first Price:" << x<<endl;
    x2 = pow(x, 2);
    
    for (int i=1; i<no_of_trials; i++) {
        RV = get_gaussian_inverse();
        Price[0]=_S0*exp((_r-_q-pow(_sigma,2)*0.5)*_deltaT*(1)+_sigma*sqrt(_deltaT*(1))*RV);
        for (int M = 1; M<_m; M++) {
            double RV = get_gaussian_inverse();
            Price[M] = Price[M-1]*exp((_r-_q-pow(_sigma,2)*0.5)*_deltaT*(1)+_sigma*sqrt(_deltaT*(1))*RV);
        }
        double Inter_Price = pow(Price[0],inv_m);
        for (int M = 1; M<_m; M++) {
            Inter_Price = Inter_Price * pow(Price[M],inv_m);
        }
        Inter_Price = max(Inter_Price - _K, 0);
        long double Temp_x = 0;
        for (int j = 0; j<_m; j++) {
            Temp_x += Price[j];
        }
        Temp_x = Temp_x/_m;
        Temp_x = max(Temp_x-_K, 0) + _b*(Expectation_Geometric_Call-Inter_Price);
        double Temp_x2 = pow(Temp_x, 2);
        x += Temp_x;
        x2 += Temp_x2;
    }
    expectation=x/no_of_trials;
    se=sqrt(((x2/no_of_trials)-pow(expectation, 2))/(no_of_trials-1));
    return 1;
}

int main(int argc, const char * argv[])
{
    // insert code here...
    srand(time(NULL));
    x10= random() % 10000000;
    x11= random() % 20000000;
    x12= random() % 30000000;
    x20= random() % 40000000;
    x21= random() % 50000000;
    x22= random() % 60000000;
    
    int no_of_trials =10000;
    sscanf (argv[1], "%d", &_m);
    inv_m = 1/double(_m);
    _deltaT = _T/double(_m);
    clock_t start, finish;
    double duration;
    monte_carlo_initial();
    //monte_carlo_inverse_2();
    b_Cal();
    Expectation_Geometric_Call_cal();
    //monte_carlo_inverse1(no_of_trials);
    
    start = clock();
    monte_carlo_inverse(no_of_trials);
    finish = clock();
    
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    double discount=exp(-_r*_T);
    cout << "Number of trials: " << no_of_trials<<endl;
    cout << "Number of dimensions: " << _m<<endl;
    cout << "-------------------------------" <<endl;
    cout << "Asian Options with Control Variable:" <<endl;
    cout << "The price is: " << discount*expectation<<endl;
    cout << "b: " <<_b<<endl;
    cout << "pho: " <<_pho<<endl;
    cout << "The Standard Error is: "<< discount*se<<endl;//here I have discount on the se
    cout << "The Computational Time(s): "<<duration<<endl;
    cout << "-------------------------------" <<endl;
}

