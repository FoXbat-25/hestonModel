#include <iostream>
#include <stdio.h>
#include <sim_engine.h>
#include <constants.h>
#include <option_pricing.h>

using namespace std;

int main() {

    auto start = std::chrono::high_resolution_clock::now();
    cout<< "Initiating heston-model"<<endl;
    all_S_v sim = gbm_sim();
    plottMultiple(t, sim.S, "Monte Carlo Paths: Stock Price");
    plottMultiple(t, sim.v, "Monte Carlo Paths: Volatility");

    prices opt_price = heston_prices_parellel();
    double call_option = opt_price.call;
    double put_option = opt_price.put;

    std::cout << "Call option is priced at " << call_option << endl;
    std::cout << "Put option is priced at " << put_option << endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Main execution time: " << duration.count() << " ms" << std::endl;

    return 0;
}