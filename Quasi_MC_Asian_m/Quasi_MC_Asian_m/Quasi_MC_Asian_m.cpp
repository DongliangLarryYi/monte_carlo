// Pricing asian call option with quasi monte carlo method
// Inverse method used to generate gaussian random variables
// Sobol' sequence is used here
//  Created by Dongliang Yi on 4/23/16.
//  Copyright © 2016 Dongliang Yi. All rights reserved.

#include "sobol.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

double se;
double expectation;
double _K=2.0;
double _T=2.00;
double _r=0.05;
double _q=0;
double _sigma=0.5;
double _S0=2;
int _m; // 50 stock price averaging
double inv_m ;
double _deltaT ;
double _alpha=1.96; //5% significancedouble control1[1000],control2[1000];

// Generating Unit RV with L'Ecuyer algorithm on p.52 (Fig 2.3 ) [Paul_Glasserman]_Monte_Carlo_Methods_in_Financial
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
    for (int i=2; i<24-j; i=i+2) 
    {
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

// Beasley method
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

double get_gaussian_determin(double s)
{
    double U=s;
    double x0=Beasley_method(U);
    double cdf = get_cdf(x0);
    double x1 = x0 - (cdf-U)/(exp(-pow(x0, 2)/2)/sqrt(2*3.141592654));
    return x1;
}

double max(double a, double b) {
    return (b < a )? a:b;
}

// quasi monte carlo method with sobol sequence
double Quasi_Monte_Carlo_Sobol(int number, int batch)
{
    double y, y2;
    y=0;
    y2=0;
    srand((int)time(0));
    for (int j = 0 ; j<batch; j++) 
    {
        double x;
        x =0;
        double U[_m];
        for (int i =0 ; i<_m; i++) 
        {
            srand(((int)time(0)*10000*(j+1))*batch);
            U[i]=Uniform01();
        }
        for (unsigned long long i = 0; i < number; ++i)
        {
            // Print a few dimensions of each point.
            double Price[_m];
            const double s = sobol::sample(i+1, 0);
            double inter = s+U[0] - floor(s+U[0]);
            const double RV = get_gaussian_determin(inter);
            Price[0]=_S0*exp((_r-_q-pow(_sigma,2)*0.5)*_deltaT*(1)+_sigma*sqrt(_deltaT*(1))*RV);
            for (unsigned d = 1; d < _m; ++d)
            {
                const double s = sobol::sample(i+1, d);
                double inter = s+U[d] - floor(s+U[d]);
                const double RV = get_gaussian_determin(inter);
                Price[d]=Price[d-1]*exp((_r-_q-pow(_sigma,2)*0.5)*_deltaT*(1)+_sigma*sqrt(_deltaT*(1))*RV);
            }
            double Temp_x = 0;
            for (int j = 0; j<_m; j++) {
                Temp_x += Price[j];
            }
            Temp_x = Temp_x/_m;
            Temp_x = max(Temp_x-_K, 0);
            x += Temp_x;
        }
        x = x/number;
        double x2= pow(x, 2);
        y += x;
        y2 +=x2;
    }
    expectation = y/batch;
    se = sqrt(((y2/batch)-pow(expectation, 2))/(batch-1));
    return expectation;
}

int main(int argc, const char * argv[])
{
    srand(time(NULL));
    x10= random() % 10000000;
    x11= random() % 20000000;
    x12= random() % 30000000;
    x20= random() % 40000000;
    x21= random() % 50000000;
    x22= random() % 60000000;
    int no_of_trials =500;
    sscanf (argv[1], "%d", &_m);
    inv_m = 1/double(_m);
    _deltaT = _T/double(_m);
    sscanf (argv[1], "%d", &_m);
    int no_of_batchs =20;
    clock_t start, finish;
    double duration;
    
    // Iterate over points.
    start = clock();
    double Asian_Price;
    Asian_Price = Quasi_Monte_Carlo_Sobol(no_of_trials, no_of_batchs);
    double discount=exp(-_r*_T);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout << "The dimension is: " << _m<<endl;
    cout << "-------------------------------" <<endl;
    cout << "Total trials: " <<no_of_trials<<" * "<<no_of_batchs<< " = "<<no_of_batchs*no_of_trials<<endl<<endl;
    cout << "The price is: " << discount*expectation<<endl;
    cout << "The Standard Error is: "<< discount*se<<endl;//here I have discount on the se
    cout << "The Computational Time(s): "<<duration<<endl;
    cout << "-------------------------------" <<endl;
}

