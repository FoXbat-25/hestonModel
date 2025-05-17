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

using namespace std;

double du = 0.01;

double heston_charteristic_func(u){
    std::complex<double> z1(0.0, 1.0);
    std::complex<double> xi = kappa - rho * sigma * z1 * u;
    std::complex<double> d = sqrt(pow(rho * sigma * z1 * u - xi, 2) - pow(sigma, 2) * (-u*z1 - pow(u,2)));
    std::complex<double> g = (xi - rho * sigma * z1 * u - d) / (xi - rho * sigma * z1 * u + d);
    std::complex<double> C = r * z1 * u * T + (kappa * theta) / pow(sigma, 2) * ((xi - rho * sigma * z1 * u - d) * T - 2.0 * log((1.0 - g * exp(-d * T)) / (1.0 - g)));
    std::complex<double> D = (xi - rho * sigma * z1 * u - d) / pow(sigma, 2) * ((1.0 - exp(-d * T)) / (1.0 - g * exp(-d * T)));

    return exp(C + D * v0 + z1 * u * log(S0));
}

double heston_call_prices(){

    double integrand = exp(z1 * u * log(K)).real() / (z1 * u) * heston_charteristic_func(u-z1)
}
