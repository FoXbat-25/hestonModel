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
#include <fftw3.h>

using namespace std;

const std::complex<double> i(0.0, 1.0);

std::vector<double> laguerre_nodes = { 0.022415874146706646, 0.11812251209676998, 0.2903657440180373, 0.5392862212279776, 0.8650370046481125, 1.2678140407752432, 1.7478596260594357, 2.3054637393075073, 2.9409651567252513, 3.654752650207291, 4.447266343313096, 5.318999254496392, 6.270499046923655, 7.302370002587394, 8.415275239483027, 9.609939192796109, 10.887150383886373, 12.2477645042443, 13.692707845547504, 15.222981111524728, 16.839663652648735, 18.54391817085919, 20.336995948730234, 22.220242665950874, 24.195104875933254, 26.263137227118484, 28.426010527501028, 30.685520767525972, 33.04359923643783, 35.50232389114121, 38.06393216564647, 40.73083544445863, 43.50563546642153, 46.391142978616195, 49.39039902562469, 52.506699341346305, 55.74362241327838, 59.10506191901711, 62.59526440015139, 66.21887325124756, 69.98098037714684, 73.88718723248296, 77.94367743446313, 82.1573037783193, 86.53569334945652, 91.08737561313309, 95.82194001552074, 100.75023196951398, 105.88459946879995, 111.23920752443958, 116.8304450513065, 122.67746026853858, 128.80287876923768, 135.23378794952583, 142.00312148993152, 149.15166590004938, 156.73107513267115, 164.8086026551505, 173.47494683642427, 182.85820469143147, 193.15113603707292, 204.67202848505946, 218.03185193532852, 234.80957917132616 };
std::vector<double> laguerre_weights = { 0.056252842339248155, 0.11902398731245409, 0.15749640386206856, 0.1675470504157258, 0.15335285577918614, 0.12422105360931551, 0.0903423009864627, 0.05947775576833932, 0.03562751890402631, 0.019480410431161187, 0.00974359489937938, 0.00446431036416511, 0.0018753595813226526, 0.0007226469815748107, 0.00025548753283343134, 8.287143534394733e-05, 2.4656863967879446e-05, 6.726713878827902e-06, 1.681785369963643e-06, 3.850812981545725e-07, 8.068728040988378e-08, 1.5457237067572897e-08, 2.7044801476167667e-09, 4.3167754754259984e-10, 6.277752541759839e-11, 8.3063173762867e-12, 9.984031787217486e-13, 1.0883538871163783e-13, 1.0740174034413078e-14, 9.575737231571912e-16, 7.69702802364661e-17, 5.564881137452472e-18, 3.609756409009412e-19, 2.0950953695484468e-20, 1.0847933010972538e-21, 4.994699486362548e-23, 2.037836974598286e-24, 7.339537564276772e-26, 2.3237830821981357e-27, 6.438234706906655e-29, 1.5531210957879239e-30, 3.2442500920186983e-32, 5.832386267834391e-34, 8.963254833100631e-36, 1.1687039895504068e-37, 1.2820559843596414e-39, 1.1720949374046017e-41, 8.835339672325698e-44, 5.424955590304308e-46, 2.6755426666782054e-48, 1.0429170314111074e-50, 3.152902351956925e-53, 7.229541910645761e-56, 1.2242353012297833e-58, 1.4821685049012373e-61, 1.232519348814245e-64, 6.6914990045686734e-68, 2.2204659418499165e-71, 4.120946094738085e-75, 3.7743990618958994e-79, 1.4141150529168788e-83, 1.59183306404108e-88, 2.9894843488597464e-94, 2.0890635084364083e-101 };

std::complex<double> phi_1(complex<double> u){
    
    complex<double> i(0.0, 1.0);

    double b1 = kappa + lambda - (rho * sigma);
    double a = kappa * theta;

    complex<double> beta = b1 - (rho * u * i * sigma);

    complex<double> d = sqrt( (beta * beta) - (sigma * sigma * ( (u * i) - (u * u) )));
    if (d.real() < 0) d = -d;

    complex<double> g = (beta + d) / (beta - d);

    complex<double> C = (r * u * i * tau) + ((a/(sigma * sigma)) * (((beta + d) * tau) - (2.0 * log ((1.0 - (g * exp(d*tau))) / (1.0 - g)))));
    complex<double> D = ((beta + d)/(sigma * sigma)) * (( 1.0 - exp( d * tau )) / ( 1.0 - (g * exp(d * tau))));

    return exp(C + D * v0 + i * u * std::log(S0));
}

std::complex<double> phi_2(complex<double> u){
    
    complex<double> i(0.0, 1.0);

    double b2 = kappa + lambda;
    double a = kappa * theta;

    complex<double> beta = b2 - (rho * u * i * sigma);

    complex<double> d = sqrt( (beta * beta) + (sigma * sigma * ( (u * i) + (u * u) )));
    complex<double> g = (beta - d) / (beta + d);

    // if (std::abs(g) > 1.0)
    // g *= 0.99999999 / std::abs(g); 

    complex<double> C = (r * u * i * tau) + ((a/(sigma * sigma)) * (((beta - d) * tau) - (2.0 * log ((1.0 - (g * exp(-d*tau))) / (1.0 - g)))));
    complex<double> D = ((beta - d)/(sigma * sigma)) * (( 1.0 - exp( -d * tau )) / ( 1.0 - (g * exp(-d * tau))));

    return exp(C + D * v0 + i * u * std::log(S0));
}

pairr gaussLaguerreProb() {
    
    double sum1 = 0.0;
    double sum2 = 0.0;
    
    double integral = 0.0;
    std::complex<double> i(0.0, 1.0);

    auto integrand1 = [&](double u) -> double{
        std::complex<double> q = std::exp(u) * std::exp(-i * u * std::log(K)) * phi_1(u) / (i * u);
        return q.real();
    };

    for (size_t p = 0; p < laguerre_nodes.size(); ++p) {
        sum1+= laguerre_weights[p] * integrand1(laguerre_nodes[p]);
    }

    auto integrand2 = [&](double u) -> double{
        std::complex<double> t = std::exp(u) * std::exp(-i * u * std::log(K)) * phi_2(u) / (i * u);
        return t.real();
    };

    for (size_t p = 0; p < laguerre_nodes.size(); ++p) {
        sum2 += laguerre_weights[p] * integrand2(laguerre_nodes[p]);
    }


    double P1 = 0.5 + sum1 / M_PI;
    cout<<P1<<endl;
    double P2 = 0.5 + sum2 / M_PI;
    cout<<P2<<endl;

    double call_price = S0 * P1 - K * std::exp(-r * T) * P2;
    double put_price = call_price - S0 + K * std::exp(-r * T);
    pairr opt_price;
    opt_price.first = call_price;
    opt_price.second = put_price;
    return opt_price;
    
}

double callFFT(double K_target){
    const int N = 4096; 
    const double alpha = 1.5; 
    const double eta = 0.25; 
    const double l = 2.0 * M_PI / (N * eta);
    double b = M_PI / eta;
    int m = N;
    
    std::complex<double> i(0.0, 1.0);
    std::vector<std::complex<double>> fft_input(m);
    std::vector<std::complex<double>> psii(m);
    std::vector<std::complex<double>> integrandd(m);

    for (int j = 0; j<m; ++j){
        double v = eta * j;
        std::complex<double> u = v - i * (alpha + 1.0);
        std::complex<double> psi = std::exp(-r * T) * phi_2(u);
        std::complex<double> integrand = (std::exp(i * v * b) * psi )/ 
                                        ((alpha * alpha) + alpha - (v * v) + (i * (2.0 * alpha + 1.0) * v));

        double simpson = (j == 0 || j == m ) ? 1.0 : ((j % 2 == 0) ? 2.0 : 4.0);
        fft_input[j] = integrand * eta * simpson/3.0;  
        psii[j] = psi.real();
        integrandd[j] = integrand;
    }

    std::vector<std::complex<double>> fft_output(m);
    fftw_complex* in = fftw_alloc_complex(m);
    fftw_complex* out = fftw_alloc_complex(m);

    for (int j = 0; j < m; ++j) {
        in[j][0] = fft_input[j].real();
        in[j][1] = fft_input[j].imag();
    }

    fftw_plan plan = fftw_plan_dft_1d(m, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int j = 0; j < m; ++j) {
        fft_output[j] = std::complex<double>(out[j][0], out[j][1]) / static_cast<double>(m);
    }

    fftw_free(in);
    fftw_free(out);

    std::vector<double> strikes(m);
    std::vector<double> call_prices(m);
    std::vector<double> indxx(m);

    for (int j = 0; j < m; ++j) {
        double k = -b + j * l; 
        double kk = std::exp(k + log(S0));   

        //double k = -b + j * lambda;
        std::complex<double> damping = std::exp(-alpha * k);
        double call_price = (damping * fft_output[j]).real() / M_PI;

        //indxx[j] = j;
        strikes[j] = k;
        call_prices[j] = call_price;
    }

    double log_K_target = std::log(K_target);
    double float_idx = (log_K_target + b) / l;
    int idx = static_cast<int>(std::floor(float_idx));
    if (idx < 0) idx = 0;
    if (idx >= m - 1) idx = m - 2;
    double w = float_idx - idx;
    double price = (1 - w) * call_prices[idx] + w * call_prices[idx + 1];

    std::vector<double> fft_output_real(m);
    std::vector<double> fft_input_real(m);
    std::vector<double> intgrnd(m);
    
    for (int j = 0; j < m; ++j) {
        fft_output_real[j] = fft_output[j].real();
        fft_input_real[j] = fft_input[j].real();
        indxx[j] = j;
        intgrnd[j] = integrandd[j].imag();
    }

    plott(indxx, strikes, std::to_string(eta));
    plott(indxx, fft_output_real, "FFT output real");
    plott(indxx, fft_input_real, "FFT input real");
    plott(indxx, intgrnd, "Integrand real");

    return price;
}

// First, let's add a debug version to check intermediate values
std::complex<double> phi_2_debug(std::complex<double> u) {
    std::complex<double> i(0.0, 1.0);

    double b2 = kappa + lambda;
    double a = kappa * theta;

    std::complex<double> beta = b2 - (rho * u * i * sigma);
    std::complex<double> d = sqrt((beta * beta) + (sigma * sigma * ((u * i) + (u * u))));
    std::complex<double> g = (beta + d) / (beta - d);

    // Check for numerical issues
    std::complex<double> exp_dt = exp(d * tau);
    std::complex<double> denom = 1.0 - (g * exp_dt);
    
    if (std::abs(denom) < 1e-12) {
        // Handle near-zero denominator
        denom = std::complex<double>(1e-12, 0.0);
    }

    std::complex<double> C = (r * u * i * tau) + ((a / (sigma * sigma)) * 
        (((beta + d) * tau) - (2.0 * log((1.0 - (g * exp_dt)) / (1.0 - g)))));
    
    std::complex<double> D = ((beta + d) / (sigma * sigma)) * 
        ((1.0 - exp_dt) / denom);

    std::complex<double> result = exp(C + D * v0 + i * u * std::log(S0));
    
    // Debug output for problematic values
    if (!std::isfinite(result.real()) || !std::isfinite(result.imag())) {
        std::cout << "Warning: Non-finite characteristic function at u = " << u << std::endl;
        return std::complex<double>(1.0, 0.0); // Return neutral value
    }
    
    return result;
}

double callFFT_corrected(double K_target) {
    const int N = 4096;
    const double alpha = 1.2;  // Must be > 1 for call options
    const double eta = 0.1;
    const double lambda = 2.0 * M_PI / (N * eta);  // Using lambda instead of l
    const double b = (N * lambda) / 2.0;
    
    std::complex<double> i(0.0, 1.0);
    std::vector<std::complex<double>> fft_input(N);

    // Debug: Check parameter ranges
    double log_strike_min = -b;
    double log_strike_max = b;
    double strike_min = exp(log_strike_min);
    double strike_max = exp(log_strike_max);
    
    std::cout << "Strike range: [" << strike_min << ", " << strike_max << "]" << std::endl;
    std::cout << "Target strike: " << K_target << std::endl;

    for (int j = 0; j < N; ++j) {
        double v = eta * j;  // This is the integration variable
        
        // The key correction: u should be the Fourier variable
        std::complex<double> u = v - i * (alpha + 1.0);
        
        // Get characteristic function
        std::complex<double> phi_u = phi_2_debug(u);
        
        // Apply discount factor
        std::complex<double> psi = std::exp(-r * tau) * phi_u;
        
        // Carr-Madan integrand
        std::complex<double> denominator = alpha * alpha + alpha - v * v + i * (2.0 * alpha + 1.0) * v;
        
        // Check for near-zero denominator
        if (std::abs(denominator) < 1e-12) {
            fft_input[j] = std::complex<double>(0.0, 0.0);
            continue;
        }
        
        std::complex<double> integrand = (std::exp(i * v * b) * psi) / denominator;

        // Simpson's rule weights
        double simpson = (j == 0 || j == N - 1) ? 1.0 : ((j % 2 == 0) ? 2.0 : 4.0);
        
        fft_input[j] = integrand * eta * simpson / 3.0;
        
        // Debug first few values
        if (j < 5) {
            std::cout << "j=" << j << ", v=" << v << ", u=" << u << ", phi=" << phi_u 
                      << ", integrand=" << integrand << std::endl;
        }
    }

    // FFT computation
    fftw_complex* in = fftw_alloc_complex(N);
    fftw_complex* out = fftw_alloc_complex(N);

    for (int j = 0; j < N; ++j) {
        in[j][0] = fft_input[j].real();
        in[j][1] = fft_input[j].imag();
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    std::vector<std::complex<double>> fft_output(N);
    for (int j = 0; j < N; ++j) {
        fft_output[j] = std::complex<double>(out[j][0], out[j][1]) / static_cast<double>(N);
    }

    fftw_free(in);
    fftw_free(out);

    // Extract option prices
    std::vector<double> strikes(N);
    std::vector<double> call_prices(N);

    for (int j = 0; j < N; ++j) {
        double k = -b + j * lambda;  // log-strike
        double K = std::exp(k);      // strike price

        double Re_fft = fft_output[j].real();
        double price = std::exp(-alpha * k) * Re_fft / M_PI;
        
        // Ensure non-negative prices
        price = std::max(0.0, price);

        strikes[j] = K;
        call_prices[j] = price;
        
        // Debug some prices
        if (j % 500 == 0) {
            std::cout << "Strike=" << K << ", Price=" << price << ", FFT_real=" << Re_fft << std::endl;
        }
    }

    // Check if all prices are the same (straight line issue)
    bool all_same = true;
    for (int j = 1; j < N; ++j) {
        if (std::abs(call_prices[j] - call_prices[0]) > 1e-10) {
            all_same = false;
            break;
        }
    }
    
    if (all_same) {
        std::cout << "ERROR: All prices are the same! Check characteristic function and parameters." << std::endl;
    }

    // Linear interpolation
    double log_K_target = std::log(K_target);
    double float_idx = (log_K_target + b) / lambda;
    
    if (float_idx <= 0) {
        return call_prices[0];
    }
    if (float_idx >= N - 1) {
        return call_prices[N - 1];
    }
    
    int idx = static_cast<int>(std::floor(float_idx));
    double w = float_idx - idx;
    double price = (1.0 - w) * call_prices[idx] + w * call_prices[idx + 1];

    std::cout << "Final interpolated price: " << price << std::endl;
    
    // Call your plotting function
    plott(call_prices, strikes, "");

    return std::max(0.0, price);
}

// Alternative implementation with different approach
double callFFT_alternative(double K_target) {
    const int N = 4096;
    const double alpha = 1.5;
    const double eta = 0.25;
    const double lambda = 2.0 * M_PI / (N * eta);
    const double b = N * lambda / 2.0;
    
    std::complex<double> i(0.0, 1.0);
    
    // Create the integrand array
    std::vector<std::complex<double>> x(N);
    
    for (int j = 0; j < N; ++j) {
        double v_j = eta * j;
        
        // Fourier transform variable
        std::complex<double> u = v_j - i * (alpha + 1.0);
        
        // Characteristic function
        std::complex<double> phi_u = phi_2_debug(u);
        
        // Modified characteristic function
        std::complex<double> psi = std::exp(-r * tau) * phi_u;
        
        // Denominator
        std::complex<double> denom = std::pow(alpha + i * v_j, 2.0) + alpha;
        
        if (std::abs(denom) < 1e-12) {
            x[j] = std::complex<double>(0.0, 0.0);
        } else {
            x[j] = std::exp(i * b * v_j) * psi / denom;
        }
        
        // Apply quadrature weight
        x[j] *= eta;
    }
    
    // Apply FFT
    fftw_complex* in = fftw_alloc_complex(N);
    fftw_complex* out = fftw_alloc_complex(N);
    
    for (int j = 0; j < N; ++j) {
        in[j][0] = x[j].real();
        in[j][1] = x[j].imag();
    }
    
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    // Extract prices
    std::vector<double> call_prices(N);
    std::vector<double> strikes(N);
    
    for (int j = 0; j < N; ++j) {
        double k_j = -b + lambda * j;
        double K_j = std::exp(k_j);
        
        double fft_real = out[j][0] / N;  // Normalize
        double price = std::exp(-alpha * k_j) * fft_real / M_PI;
        
        call_prices[j] = std::max(0.0, price);
        strikes[j] = K_j;
    }
    
    fftw_free(in);
    fftw_free(out);
    
    // Find target price
    double log_K = std::log(K_target);
    double idx_float = (log_K + b) / lambda;
    
    if (idx_float <= 0) return call_prices[0];
    if (idx_float >= N-1) return call_prices[N-1];
    
    int idx = (int)std::floor(idx_float);
    double w = idx_float - idx;
    
    return (1.0 - w) * call_prices[idx] + w * call_prices[idx + 1];
}

int main(){
    pairr prob = gaussLaguerreProb();
    cout<<prob.first<<endl;
    cout<<prob.second<<endl;

    double c = callFFT(100);
    cout<<c;
    // const int N = 4096;
    // std::vector<double> s(50);
    // std::vector<double> t(50);
    // const double alpha = 1.5;
    // const double eta = 0.25;
    // const double lambda = 2.0 * M_PI / (N * eta);
    // const double b = N * lambda / 2.0;
    // std::complex<double> im(0.0, 1.0);  
    // for (int p=0; p<50; ++p){
    //     double v = eta * p;
    //     std::complex<double> u = v - i * (alpha + 1.0);
    //     std::complex<double> psi = std::exp(-r * T) * phi_2(u);
    //     std::complex<double> integrand = (std::exp(i * v * b) * psi )/ 
    //                                     ((alpha * alpha) + alpha - (v * v) + (i * (2.0 * alpha + 1.0) * v));
    //     s[p] = integrand.real();
    //     t[p] = p;
    // }

    // plott(s,t, "");
    
}


