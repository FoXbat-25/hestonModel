#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <vector>
#include <sim_engine.h>

extern const std::vector<double> t;
constexpr double RISK_FREE_RATE = 0.05;
constexpr int num_paths = 100;
constexpr float T = 1.0; //total time
constexpr int N = 1000; // time steps
constexpr double dt = T/N;
constexpr float mu = 0.1; //expected return
constexpr float v0 = 0.1;//Initial volatility
constexpr float kappa = 3.0; //Mean reversion pullback force
constexpr float theta = 0.15; //Long term volatility
constexpr float sigma = 0.1; //Volatility
constexpr double K = 110.0; //Price above ?

#endif