#include <iostream>
#include <stdio.h>
#include <sim_engine.h>
#include <constants.h>

using namespace std;

int main() {

    auto start = std::chrono::high_resolution_clock::now();
    cout<< "Initiating heston-model"<<endl;
    all_S_v sim = gbm_sim();
    plottMultiple(t, sim.S, "Monte Carlo Paths: Stock Price");
    plottMultiple(t, sim.v, "Monte Carlo Paths: Volatility");

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Main execution time: " << duration.count() << " ms" << std::endl;

    return 0;
}