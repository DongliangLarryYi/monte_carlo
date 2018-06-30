//  Monte carlo simulation with antithetic approach for variance reduction
//  Created by Dongliang Yi on 2/10/16.
//  Copyright © 2016 Dongliang Yi. All rights reserved.

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

// generate unifor Random variables
double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
}

// generate gaussian random variables
double get_gaussian()
{
    return (sqrt(-2.0*log(get_uniform()))*cos(2*3.141592654*get_uniform()));
}

// return bigger number
double max(double a, double b) {
    return (b < a )? a:b;
}

// montel carlo simulation with antithetic approach
int monte_carlo(int no_of_trials)
{
    double x,x2;
    double temp=(_r-_q-pow(_sigma,2)*0.5)*_T;
    double temp2=_sigma*sqrt(_T);
    double temp3,temp4;
    double RV = get_gaussian();
    // antithetic method
    x = (max(_S0*exp(temp+RV*temp2)-_K,0)+max(_S0*exp(temp-RV*temp2)-_K,0))/2;
    x2 = pow(x, 2);
    for (int i=2; i<=no_of_trials; i++) {
        RV = get_gaussian();
        temp3= (max(_S0*exp(temp+RV*temp2)-_K,0)+max(_S0*exp(temp-RV*temp2)-_K,0))/2;
        temp4=pow(temp3, 2);
        x+=temp3;
        x2+=temp4;
    }
    // expected value and standard deviation
    expectation=x/no_of_trials;
    se=sqrt(((x2/no_of_trials)-pow(expectation, 2))/(no_of_trials-1));
    return 1;
}

int main(int argc, const char * argv[])
{
    int no_of_trials;
    sscanf (argv[1], "%d", &no_of_trials);
    monte_carlo(no_of_trials);
    double discount=exp(-_r*_T);
    
    cout << "Number of trials: " << no_of_trials<<endl;
    cout << "The call price is: " << discount*expectation<<endl;
    cout << "The Standard Error is: "<< discount*se<<endl;//here I have discount on the se
    cout << "The 95% CI: ["<<discount*expectation-discount*se*_alpha<<", "<<discount*expectation+discount*se*_alpha<<"]"<<endl;
    cout << "The CI Interval length(cents): "<<200*discount*se*_alpha<<endl;
}

