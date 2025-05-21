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
const std::complex<double> z1(0.0, 1.0);
std::complex<double> b1 = kappa + lambda - (rho * sigma);
std::complex<double> a = kappa * theta;
int phi_1(std::complex<double> u){
    
    std::complex<double> d1 = sqrt(pow(rho * sigma * z1 * u - b1, 2) - sigma * sigma * (z1 * u - u * u));
    std::complex<double> g1 = (b1 - rho * sigma * z1 * u + d1) / (b1 - rho * sigma * z1 * u - d1);
    std::complex<double> C1 = (r * z1 * u * tau) + ((a / pow(sigma, 2)) * ((b1 - (rho * sigma * z1 * u) + d1)) * tau - 2.0 * log((1.0 - g1 * exp(d1 * tau)) / (1.0 - g1)));
    std::complex<double> D1 = (b1 - rho * sigma * z1 * u + d1) / pow(sigma, 2) * ((1.0 - exp(d1 * tau)) / (1.0 - g1 * exp(d1 * tau)));

    cout<< d1 << endl;
    cout<< g1 << endl;
    cout<< C1 << endl;
    cout<< D1 << endl;
    cout<< log(1.0 - g1 * exp(d1 * tau)) - log(1.0 - g1) <<endl;

    cout << exp(C1 + D1 * v0 + z1 * u * log(S0)) <<endl;
    cout<< (std::exp(-z1 * u * std::log(K)) /  (z1 * u)) <<endl;

    return 0;
}

int main(){
    phi_1(0.022);
}