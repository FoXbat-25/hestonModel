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
// const std::complex<double> z1(0.0, 1.0);

std::vector<double> laguerre_nodes = { 0.022415874146706646, 0.11812251209676998, 0.2903657440180373, 0.5392862212279776, 0.8650370046481125, 1.2678140407752432, 1.7478596260594357, 2.3054637393075073, 2.9409651567252513, 3.654752650207291, 4.447266343313096, 5.318999254496392, 6.270499046923655, 7.302370002587394, 8.415275239483027, 9.609939192796109, 10.887150383886373, 12.2477645042443, 13.692707845547504, 15.222981111524728, 16.839663652648735, 18.54391817085919, 20.336995948730234, 22.220242665950874, 24.195104875933254, 26.263137227118484, 28.426010527501028, 30.685520767525972, 33.04359923643783, 35.50232389114121, 38.06393216564647, 40.73083544445863, 43.50563546642153, 46.391142978616195, 49.39039902562469, 52.506699341346305, 55.74362241327838, 59.10506191901711, 62.59526440015139, 66.21887325124756, 69.98098037714684, 73.88718723248296, 77.94367743446313, 82.1573037783193, 86.53569334945652, 91.08737561313309, 95.82194001552074, 100.75023196951398, 105.88459946879995, 111.23920752443958, 116.8304450513065, 122.67746026853858, 128.80287876923768, 135.23378794952583, 142.00312148993152, 149.15166590004938, 156.73107513267115, 164.8086026551505, 173.47494683642427, 182.85820469143147, 193.15113603707292, 204.67202848505946, 218.03185193532852, 234.80957917132616 };
std::vector<double> laguerre_weights = { 0.056252842339248155, 0.11902398731245409, 0.15749640386206856, 0.1675470504157258, 0.15335285577918614, 0.12422105360931551, 0.0903423009864627, 0.05947775576833932, 0.03562751890402631, 0.019480410431161187, 0.00974359489937938, 0.00446431036416511, 0.0018753595813226526, 0.0007226469815748107, 0.00025548753283343134, 8.287143534394733e-05, 2.4656863967879446e-05, 6.726713878827902e-06, 1.681785369963643e-06, 3.850812981545725e-07, 8.068728040988378e-08, 1.5457237067572897e-08, 2.7044801476167667e-09, 4.3167754754259984e-10, 6.277752541759839e-11, 8.3063173762867e-12, 9.984031787217486e-13, 1.0883538871163783e-13, 1.0740174034413078e-14, 9.575737231571912e-16, 7.69702802364661e-17, 5.564881137452472e-18, 3.609756409009412e-19, 2.0950953695484468e-20, 1.0847933010972538e-21, 4.994699486362548e-23, 2.037836974598286e-24, 7.339537564276772e-26, 2.3237830821981357e-27, 6.438234706906655e-29, 1.5531210957879239e-30, 3.2442500920186983e-32, 5.832386267834391e-34, 8.963254833100631e-36, 1.1687039895504068e-37, 1.2820559843596414e-39, 1.1720949374046017e-41, 8.835339672325698e-44, 5.424955590304308e-46, 2.6755426666782054e-48, 1.0429170314111074e-50, 3.152902351956925e-53, 7.229541910645761e-56, 1.2242353012297833e-58, 1.4821685049012373e-61, 1.232519348814245e-64, 6.6914990045686734e-68, 2.2204659418499165e-71, 4.120946094738085e-75, 3.7743990618958994e-79, 1.4141150529168788e-83, 1.59183306404108e-88, 2.9894843488597464e-94, 2.0890635084364083e-101 };

std::complex<double> heston_characteristic_func(std::complex<double> u){
    
    std::complex<double> xi = kappa - rho * sigma * z1 * u;
    std::complex<double> d = sqrt(pow(rho * sigma * z1 * u - xi, 2) + sigma * sigma * (z1 * u + u * u));

    std::complex<double> g = (xi - rho * sigma * z1 * u - d) / (xi - rho * sigma * z1 * u + d);
    std::complex<double> C = (r * z1 * u * T) + ((kappa * theta) / pow(sigma, 2)) * ((xi - rho * sigma * z1 * u - d) * T - 2.0 * log((1.0 - g * exp(-d * T)) / (1.0 - g)));
    std::complex<double> D = (xi - rho * sigma * z1 * u - d) / pow(sigma, 2) * ((1.0 - exp(-d * T)) / (1.0 - g * exp(-d * T)));

    return exp(C + D * v0 + z1 * u * log(S0));
}

prices heston_prices_parellel(){

    auto integrand1 = [&](double u) -> double {
        std::complex<double> value = (std::exp(-z1 * u * std::log(K)) /  (z1 * u)) *  heston_characteristic_func(u - z1);
            return value.real();
    };

    auto integrand2 = [&](double u) -> double {
        std::complex<double> value = (std::exp(-z1 * u * std::log(K)) /  (z1 * u)) *  heston_characteristic_func(u);
            return value.real();
    };

    size_t N = laguerre_nodes.size();
    size_t num_threads = std::thread::hardware_concurrency();
    size_t chunk = (N + num_threads - 1) / num_threads;

    std::vector<std::future<std::pair<double, double>>> futures;

    for (size_t t = 0; t < num_threads; ++t) {
        futures.push_back(std::async(std::launch::async, [=, &integrand1, &integrand2]() -> std::pair<double, double>  {
            double sum1 = 0.0;
            double sum2 = 0.0;
            size_t start = t * chunk;
            size_t end = std::min(start + chunk, laguerre_nodes.size());
            for (size_t i = start; i < end; ++i) {
                sum1 += laguerre_weights[i] * integrand1(laguerre_nodes[i]);
                sum2 += laguerre_weights[i] * integrand2(laguerre_nodes[i]);
            }
            return std::make_pair(sum1, sum2);
        }));
    }

    double integral1 = 0.0, integral2 = 0.0;
    for (auto& f : futures) {
        auto [part1, part2] = f.get();
        integral1 += part1;
        integral2 += part2;
    }

    // Step 4: Compute P1, P2
    double P1 = 0.5 + (integral1 / M_PI);
    double P2 = 0.5 + (integral2 / M_PI);

    double call_price = S0 * P1 - K * std::exp(-r * T) * P2;
    double put_price = call_price - S0 + K * std::exp(-r * T);
    prices opt_price;
    opt_price.call = call_price;
    opt_price.put = put_price;
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

// std::vector<double> laguerre_nodes = {
//     0.044489365833267, 0.234526109519619, 0.576884629301887, 1.07244875381782,
//     1.72240877644465, 2.52833670642579, 3.49221327302199, 4.61645676974977,
//     5.90395850417424, 7.35812673318624, 8.9829409242126, 10.78301863254,
//     12.7636979867427, 14.9311397555226, 17.2924543367153, 19.8558609403361,
//     22.6308890131968, 25.6286360224592, 28.8621018163235, 32.3466291539647,
//     36.100494805752, 40.1457197715394, 44.5092079957549, 49.2243949873,
//     54.3347290344528, 59.8972646396137, 65.9913482991345, 72.7320382460062,
//     80.2769241042407, 88.8461168783834, 98.8295428682839, 111.751398097938
// };

// std::vector<double> laguerre_weights = {
//     0.109218341952385, 0.210443107938813, 0.235213229669848, 0.195903335972881,
//     0.129983786286071, 0.070578623865717, 0.031760912509175, 0.011918214834838,
//     0.003738816294611, 0.000980803306614, 0.000214864918801, 3.920341969175e-05,
//     5.934541612868e-06, 7.416404578667e-07, 7.604567879121e-08, 6.350602226625e-09,
//     4.281382971040e-10, 2.231172482824e-11, 8.876106310440e-13, 2.646771450114e-14,
//     5.527781380316e-16, 7.534216236920e-18, 6.409941131294e-20, 3.016329866515e-22,
//     7.145698456577e-25, 6.975560190105e-28, 2.385330845472e-31, 2.380756066018e-35,
//     4.119758119200e-40, 5.002592957451e-46, 5.090229640621e-54, 1.544155350067e-65
// };