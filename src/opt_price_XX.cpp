#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <functional>
#include <constants.h>

    // Characteristic function for Heston model
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
        
        std::complex<double> g = rminus / rplus;
        
        std::complex<double> exp_dT = std::exp(d * T);
        std::complex<double> D = rminus * (1.0 - exp_dT) / (1.0 - g * exp_dT);
        std::complex<double> C = a * (rminus * T - 2.0 / (sigma * sigma) * 
                                     std::log((1.0 - g * exp_dT) / (1.0 - g)));
        
        return std::exp(C + D * v0 + i * u * std::log(S0));
    }

    // Probability calculation using integration
    double calculateProbability(int j) {
        const int N = 10000;  // Number of integration points
        const double upper = 100.0;  // Upper integration limit
        const double du = upper / N;
        
        double integral = 0.0;
        std::complex<double> i(0.0, 1.0);
        
        for (int n = 1; n <= N; ++n) {
            double u = n * du;
            std::complex<double> integrand = std::exp(-i * u * std::log(K)) * 
                                           characteristicFunction(u, j) / (i * u);
            integral += integrand.real() * du;
        }
        
        return 0.5 + integral / M_PI;
    }

int main(){
    std::cout << "=== Heston Model Call Option Pricing ===" << std::endl;
    std::cout << std::endl;

    // Example parameters
    double S0 = 100.0;      // Current stock price
    double K = 100.0;       // Strike price
    double T = 1.0;         // Time to expiration (1 year)
    double r = 0.05;        // Risk-free rate (5%)
    double v0 = 0.04;       // Initial variance (20% vol)
    double kappa = 2.0;     // Mean reversion speed
    double theta = 0.04;    // Long-term variance
    double sigma = 0.3;     // Vol of vol
    double rho = -0.7;      // Correlation


    std::cout<<calculateProbability(0);
}