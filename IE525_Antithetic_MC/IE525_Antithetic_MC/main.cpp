//
//  main.cpp
//  IE525_Antithetic_MC
//
//  Created by Dongliang Yi on 2/10/16.
//  Copyright © 2016 Dongliang Yi. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

double se;
double expectation;
double _K=1870;
double _T=double(1)/double(52);
double _r=0.003866;
double _q=0.0232;
double _sigma=0.2979;
double _S0=1868.99;
double _alpha=1.96; //5% significance


double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
}

double get_gaussian()
{
    return (sqrt(-2.0*log(get_uniform()))*cos(2*3.141592654*get_uniform()));
}

double max(double a, double b) {
    return (b < a )? a:b;
}

int monte_carlo(int no_of_trials)
{
    double x,x2;
    double temp=(_r-_q-pow(_sigma,2)*0.5)*_T;
    double temp2=_sigma*sqrt(_T);
    //cout << _T <<endl;
    double temp3,temp4;
    double RV = get_gaussian();
    x = (max(_S0*exp(temp+RV*temp2)-_K,0)+max(_S0*exp(temp-RV*temp2)-_K,0))/2;
    //cout<<x<<endl;
    x2 = pow(x, 2);
    for (int i=2; i<=no_of_trials; i++) {
        RV = get_gaussian();
        temp3= (max(_S0*exp(temp+RV*temp2)-_K,0)+max(_S0*exp(temp-RV*temp2)-_K,0))/2;
        temp4=pow(temp3, 2);
        x+=temp3;
        x2+=temp4;
        //x=(1-1/i)*x+temp3/i;
        //x2=(1-1/i)*x2+temp4/i;
    }
    //se= sqrt((x2-pow(x, 2))/(no_of_trials-1));
    //expectation=x;
    expectation=x/no_of_trials;
    se=sqrt(((x2/no_of_trials)-pow(expectation, 2))/(no_of_trials-1));
    return 1;
}

int main(int argc, const char * argv[])
{
    // insert code here...
    
    int no_of_trials;
    sscanf (argv[1], "%d", &no_of_trials);
    
    monte_carlo(no_of_trials);
    double discount=exp(-_r*_T);
    
    cout << "Number of trials: " << no_of_trials<<endl;
    //cout << expectation <<endl;
    cout << "The call price is: " << discount*expectation<<endl;
    cout << "The Standard Error is: "<< discount*se<<endl;//here I have discount on the se
    cout << "The 95% CI: ["<<discount*expectation-discount*se*_alpha<<", "<<discount*expectation+discount*se*_alpha<<"]"<<endl;
    cout << "The CI Interval length(cents): "<<200*discount*se*_alpha<<endl;
    
    
}

