#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include <constants.h>
#include <complex>
#include <cmath>
#include <future>
#include <option_pricing.h>
#include <sim_engine.h>

using namespace std;

// double du = 0.01;
std::complex<double> phi_1(complex<double> u){
    
    complex<double> i(0.0, 1.0);

    double b1 = kappa + lambda - (rho * sigma);
    complex<double> a = kappa * theta;

    complex<double> beta = b1 - (rho * u * i * sigma);

    complex<double> d = sqrt( (beta * beta) - (sigma * sigma * ( (u * i) - (u * u) )));
    complex<double> g = (beta + d) / (beta - d);

    complex<double> C = (r * u * i * tau) + ((a/(sigma * sigma)) * (((beta + d) * tau) - (2.0 * log ((1.0 - (g * exp(d*tau))) / (1.0 - g)))));
    complex<double> D = ((beta + d)/(sigma * sigma)) * (( 1.0 - exp( d * tau )) / ( 1.0 - (g * exp(d * tau))));
    std::complex<double> v = std::exp(C + D * v0 + i * u * std::log(S0));
    std::cout<<"Value phi1 "<<v<<std::endl;           
    cout<<"b1 "<<b1<<endl<<"a "<<a<<endl<<"beta "<<beta<<endl<<"d "<<d<<endl<<"g "<<g<<endl<<"C "<<C<<endl<<"D "<<D<<endl;

    return exp(C + D * v0 + i * u * std::log(S0));
}

std::complex<double> characteristicFunction(std::complex<double> u, int j) {
    std::complex<double> i(0.0, 1.0);
    
    // Parameters for P1 and P2
    double b = (j == 1) ? kappa - rho * sigma : kappa;
    std::complex<double> a = kappa * theta;
    
    std::complex<double> alpha = -0.5 * u * u - 0.5 * i * u + (double)j * i * u;
    std::complex<double> beta = b - rho * sigma * i * u;
    std::complex<double> gamma = 0.5 * sigma * sigma;
    
    std::complex<double> d = std::sqrt(beta * beta - 4.0 * alpha * gamma);
    std::complex<double> rplus = (beta + d) / (2.0 * gamma);
    std::complex<double> rminus = (beta - d) / (2.0 * gamma);
    
    std::complex<double> g = rplus / rminus;
    
    std::complex<double> exp_dT = std::exp(d * T);
    std::complex<double> D = rplus * (1.0 - exp_dT) / (1.0 - g * exp_dT);
    std::complex<double> C = (r * u * i * T) + (a * (rplus * T - 2.0 / (sigma * sigma) * 
                                 std::log((1.0 - g * exp_dT) / (1.0 - g))));

    std::complex<double> v = std::exp(C + D * v0 + i * u * std::log(S0));
    std::cout<<"Value "<<v<<std::endl;                             
    cout<<"b1 "<<b<<endl<<"a "<<a<<endl<<"beta "<<beta<<endl<<"d "<<d<<endl<<"g "<<g<<endl<<"C "<<C<<endl<<"D "<<D<<endl;

    return std::exp(C + D * v0 + i * u * std::log(S0));
}

std::complex<double> phi_2(complex<double> u){
    
    complex<double> i(0.0, 1.0);

    double b2 = kappa + lambda;
    double a = kappa * theta;

    complex<double> beta = b2 - (rho * u * i * sigma);

    complex<double> d = sqrt( (beta * beta) + (sigma * sigma * ( (u * i) + (u * u) )));
    complex<double> g = (beta + d) / (beta - d);

    complex<double> C = (r * u * i * tau) + ((a/(sigma * sigma)) * (((beta + d) * tau) - (2.0 * log ((1.0 - (g * exp(d*tau))) / (1.0 - g)))));
    complex<double> D = ((beta + d)/(sigma * sigma)) * (( 1.0 - exp( d * tau )) / ( 1.0 - (g * exp(d * tau))));

    std::complex<double> v = std::exp(C + D * v0 + i * u * std::log(S0));
    std::cout<<"Value "<<v<<std::endl;  
    return exp(C + D * v0 + i * u * std::log(S0));
}

int main(){

    std::complex<double> p1 =phi_1(0.022);
    std::complex<double> p2 = characteristicFunction(0.022, 1);
    std::complex<double> p3 = phi_2(0.022);
    std::complex<double> p4 = characteristicFunction(0.022, 0);


    cout<<p1<<endl<<p2<<endl<<p3<<endl<<p4;
}