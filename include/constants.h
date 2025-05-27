#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <vector>
#include <sim_engine.h>
#include <complex>

extern const std::vector<double> t;
extern const std::complex<double> i;
constexpr double S0 = 100.0;
constexpr double r = 0.05; //risk_free_rate
constexpr int num_paths = 1;
constexpr double T = 1; //total time
constexpr int N = 1000; // time steps
constexpr double dt = T/N;
constexpr double mu = 0.1; //expected return
constexpr double v0 = 0.04;//Initial volatility
constexpr double kappa = 2.0; //Mean reversion pullback force
constexpr double theta = 0.04; //Long term volatility
constexpr double sigma = 0.1; //Volatility of volatility
constexpr double K = 100.0; //Strike price
constexpr double rho = -0.5; //Coorelation coefficient
constexpr double lambda = 0.0; //lamda
constexpr double tau = T;
constexpr double alpha = 2.0;

#endif