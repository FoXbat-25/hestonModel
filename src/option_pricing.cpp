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

double du = 0.01;
const std::complex<double> z1(0.0, 1.0);

std::vector<double> laguerre_nodes = {
    0.044489365833267, 0.234526109519619, 0.576884629301887, 1.07244875381782,
    1.72240877644465, 2.52833670642579, 3.49221327302199, 4.61645676974977,
    5.90395850417424, 7.35812673318624, 8.9829409242126, 10.78301863254,
    12.7636979867427, 14.9311397555226, 17.2924543367153, 19.8558609403361,
    22.6308890131968, 25.6286360224592, 28.8621018163235, 32.3466291539647,
    36.100494805752, 40.1457197715394, 44.5092079957549, 49.2243949873,
    54.3347290344528, 59.8972646396137, 65.9913482991345, 72.7320382460062,
    80.2769241042407, 88.8461168783834, 98.8295428682839, 111.751398097938
};

std::vector<double> laguerre_weights = {
    0.109218341952385, 0.210443107938813, 0.235213229669848, 0.195903335972881,
    0.129983786286071, 0.070578623865717, 0.031760912509175, 0.011918214834838,
    0.003738816294611, 0.000980803306614, 0.000214864918801, 3.920341969175e-05,
    5.934541612868e-06, 7.416404578667e-07, 7.604567879121e-08, 6.350602226625e-09,
    4.281382971040e-10, 2.231172482824e-11, 8.876106310440e-13, 2.646771450114e-14,
    5.527781380316e-16, 7.534216236920e-18, 6.409941131294e-20, 3.016329866515e-22,
    7.145698456577e-25, 6.975560190105e-28, 2.385330845472e-31, 2.380756066018e-35,
    4.119758119200e-40, 5.002592957451e-46, 5.090229640621e-54, 1.544155350067e-65
};

std::complex<double> heston_characteristic_func(std::complex<double> u){
    
    std::complex<double> xi = kappa - rho * sigma * z1 * u;
    std::complex<double> d = sqrt(pow(rho * sigma * z1 * u - xi, 2) - pow(sigma, 2) * (-u*z1 - pow(u,2)));
    std::complex<double> g = (xi - rho * sigma * z1 * u - d) / (xi - rho * sigma * z1 * u + d);
    std::complex<double> C = r * z1 * u * T + (kappa * theta) / pow(sigma, 2) * ((xi - rho * sigma * z1 * u - d) * T - 2.0 * log((1.0 - g * exp(-d * T)) / (1.0 - g)));
    std::complex<double> D = (xi - rho * sigma * z1 * u - d) / pow(sigma, 2) * ((1.0 - exp(-d * T)) / (1.0 - g * exp(-d * T)));

    return exp(C + D * v0 + z1 * u * log(S0));
}

prices heston_prices_parellel(){

    auto integrand = [&](double u) -> double {
        std::complex<double> value = ((std::exp(z1 * u * std::log(K)) /  (z1 * u)) *  heston_characteristic_func(u - z1));
            return value.real();
    };    

    size_t N = laguerre_nodes.size();
    size_t num_threads = std::thread::hardware_concurrency();
    size_t chunk = (N + num_threads - 1) / num_threads;

    std::vector<std::future<double>> futures;

    for (size_t t = 0; t < num_threads; ++t) {
        futures.push_back(std::async(std::launch::async, [=, &integrand]() -> double {
            double sum = 0.0;
            size_t start = t * chunk;
            size_t end = std::min(start + chunk, laguerre_nodes.size());
            for (size_t i = start; i < end; ++i) {
                sum += laguerre_weights[i] * integrand(laguerre_nodes[i]);
            }
            return sum;
        }));
    }

    double integral = 0.0;
    for (auto& f : futures) integral += f.get();

    double call_opt = std::exp(-r * T) * (0.5 * S0 - integral / M_PI);
    double put_opt = std::exp(-r * T) / M_PI * integral - S0 + K * exp(-r * T); 

    prices opt_price;
    opt_price.call = call_opt;
    opt_price.put = put_opt;
    return opt_price; 
}

// double heston_put_prices_parellel(){
    
//     auto integrand = [&](double u) -> double {
//         std::complex<double> value = ((std::exp(z1 * u * std::log(K)) /  (z1 * u)) *  heston_characteristic_func(u - z1));
//             return value.real();
//     };

//     size_t N = laguerre_nodes.size();
//     size_t num_threads = std::thread::hardware_concurrency();
//     size_t chunk = (N + num_threads - 1) / num_threads;

//     std::vector<std::future<double>> futures;
// }

// double heston_call_prices(){
    
//     double integral = trap_integration(integrand, 1e-8, 100.0, 100000);

//     return exp(-r * T) * 0.5 * S0 - exp(-r * T) / M_PI * integral;
// }