#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <functional>
#include <sim_engine.h>

// For optimization - you may need to install a library like NLopt or implement simple optimization
#include <limits>

class HestonCarrMadan {
private:
    // Model parameters
    double S0, r, q, kappa, theta, sigma, rho, v0;
    
    // FFT implementation using Cooley-Tukey algorithm
    void fft(std::vector<std::complex<double>>& x, bool invert = false) {
        int n = x.size();
        
        // Bit-reversal permutation
        for (int i = 1, j = 0; i < n; i++) {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1) {
                j ^= bit;
            }
            j ^= bit;
            if (i < j) {
                std::swap(x[i], x[j]);
            }
        }
        
        // FFT computation
        for (int len = 2; len <= n; len <<= 1) {
            double ang = 2 * M_PI / len * (invert ? -1 : 1);
            std::complex<double> wlen(cos(ang), sin(ang));
            
            for (int i = 0; i < n; i += len) {
                std::complex<double> w(1);
                for (int j = 0; j < len / 2; j++) {
                    std::complex<double> u = x[i + j];
                    std::complex<double> v = x[i + j + len / 2] * w;
                    x[i + j] = u + v;
                    x[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }
        
        if (invert) {
            for (auto& xi : x) {
                xi /= n;
            }
        }
    }
    
    // Simpson's rule weights
    std::vector<double> getSimpsonWeights(int n) {
        std::vector<double> weights(n, 1.0);
        for (int i = 1; i < n - 1; i += 2) {
            weights[i] = 4.0;
        }
        for (int i = 2; i < n - 1; i += 2) {
            weights[i] = 2.0;
        }
        weights[0] = 1.0;
        weights[n-1] = 1.0;
        
        for (auto& w : weights) {
            w /= 3.0;
        }
        
        return weights;
    }
    
    // Simple linear interpolation
    double interpolate(const std::vector<double>& x, const std::vector<double>& y, double xi) {
        if (xi <= x[0]) return y[0];
        if (xi >= x.back()) return y.back();
        
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        int idx = std::distance(x.begin(), it);
        
        if (idx == 0) return y[0];
        if (idx >= x.size()) return y.back();
        
        double x1 = x[idx-1], x2 = x[idx];
        double y1 = y[idx-1], y2 = y[idx];
        
        return y1 + (y2 - y1) * (xi - x1) / (x2 - x1);
    }
    
public:
    // Constructor
    HestonCarrMadan(double S0 = 100.0, double r = 0.05, double q = 0.0,
                    double kappa = 2.0, double theta = 0.04, double sigma = 0.3,
                    double rho = -0.7, double v0 = 0.04)
        : S0(S0), r(r), q(q), kappa(kappa), theta(theta), 
          sigma(sigma), rho(rho), v0(v0) {}
    
    // Getters and setters
    void setParameters(double S0_, double r_, double q_, double kappa_, 
                      double theta_, double sigma_, double rho_, double v0_) {
        S0 = S0_; r = r_; q = q_; kappa = kappa_;
        theta = theta_; sigma = sigma_; rho = rho_; v0 = v0_;
    }
    
    std::vector<double> getParameters() const {
        return {S0, r, q, kappa, theta, sigma, rho, v0};
    }
    
    // Heston characteristic function
    std::complex<double> characteristicFunction(std::complex<double> u, double T) {
        std::complex<double> i(0.0, 1.0);
        
        // Complex coefficients
        std::complex<double> d = sqrt(pow(rho * sigma * i * u - kappa, 2.0) - 
                                     sigma * sigma * (-i * u - u * u));
        
        std::complex<double> g = (kappa - rho * sigma * i * u - d) / 
                                (kappa - rho * sigma * i * u + d);
        
        std::complex<double> exp_dT = exp(-d * T);
        
        std::complex<double> C = (r - q) * i * u * T + 
                                (kappa * theta / (sigma * sigma)) * 
                                ((kappa - rho * sigma * i * u - d) * T - 
                                 2.0 * log((1.0 - g * exp_dT) / (1.0 - g)));
        
        std::complex<double> D = ((kappa - rho * sigma * i * u - d) / (sigma * sigma)) * 
                                ((1.0 - exp_dT) / (1.0 - g * exp_dT));
        
        return exp(C + D * v0 + i * u * log(S0));
    }
    
    // Carr-Madan FFT option pricing
    std::pair<std::vector<double>, std::vector<double>> 
    carrMadanFFT(double T, double alpha = 1.1, int N = 4096, double eta = 0.25) {
        // Ensure N is power of 2
        int n = 1;
        while (n < N) n *= 2;
        N = n;
        
        // Calculate lambda for strike spacing
        double lambda = 2.0 * M_PI / (N * eta);
        
        // Log strike range
        double b = N * lambda / 2.0;
        std::vector<double> k(N);
        for (int i = 0; i < N; i++) {
            k[i] = -b + i * lambda;
        }
        
        // Integration points
        std::vector<double> v(N);
        for (int i = 0; i < N; i++) {
            v[i] = i * eta;
        }
        
        // Modified characteristic function for Carr-Madan
        auto psi = [this, alpha, T](double v_j) -> std::complex<double> {
            std::complex<double> i(0.0, 1.0);
            std::complex<double> u = v_j - (alpha + 1.0) * i;
            return exp(-r * T) * characteristicFunction(u, T) / 
                   (alpha * alpha + alpha - v_j * v_j + i * (2.0 * alpha + 1.0) * v_j);
        };
        
        // Calculate the integrand
        std::vector<std::complex<double>> x(N);
        for (int i = 0; i < N; i++) {
            if (i == 0) {
                x[i] = 0.5 * eta * psi(v[i]);
            } else {
                x[i] = eta * psi(v[i]);
            }
        }
        
        std::vector<double> indxx(N);
        // Apply Simpson's rule weights
        std::vector<double> weights = getSimpsonWeights(N);
        for (int i = 0; i < N; i++) {
            x[i] *= weights[i];
            indxx[i] = i;

        }
        
        // Apply FFT
        fft(x);
        
        // Extract call prices and strikes
        std::vector<double> strikes, call_prices;
        strikes.reserve(N);
        call_prices.reserve(N);
        
        for (int i = 0; i < N; i++) {
            double strike = exp(k[i] + log(S0));
            double call_price = std::real(exp(-alpha * k[i]) * x[i] / M_PI);
            
            // Filter out negative prices and extreme strikes
            if (call_price > 0.001 && strike > 0.1 * S0 && strike < 5.0 * S0) {
                strikes.push_back(strike);
                call_prices.push_back(call_price);
            }
        }

        std::vector<double> fft_output_real(N);
        std::vector<double> fft_input_real(N);

        for (int j = 0; j < N; ++j) {
            fft_output_real[j] = x[j].real();
        }

        plott(indxx, fft_output_real, "");
        return std::make_pair(strikes, call_prices);
    }
    
    // Get single option price using FFT and interpolation
    double getOptionPrice(double K, double T, const std::string& option_type = "call") {
        auto [strikes, call_prices] = carrMadanFFT(T);
        
        double call_price;
        if (K < strikes.front() || K > strikes.back()) {
            // Boundary extrapolation
            if (K < strikes.front()) {
                call_price = std::max(0.0, S0 * exp(-q * T) - K * exp(-r * T));
            } else {
                call_price = 0.001;
            }
        } else {
            call_price = interpolate(strikes, call_prices, K);
        }
        
        if (option_type == "call") {
            return std::max(call_price, 0.0);
        } else if (option_type == "put") {
            // Put-call parity
            double put_price = call_price - S0 * exp(-q * T) + K * exp(-r * T);
            return std::max(put_price, 0.0);
        } else {
            throw std::invalid_argument("option_type must be 'call' or 'put'");
        }
    }
    
    // Market data structure
    struct MarketData {
        double K;
        double T;
        double market_price;
        std::string option_type;
    };
    
    // Simple optimization using grid search (for demonstration)
    // In practice, you'd use a proper optimization library
    struct CalibrationResult {
        std::vector<double> params;
        bool success;
        double rmse;
        int iterations;
    };
    
    CalibrationResult calibrate(const std::vector<MarketData>& market_data,
                               const std::vector<double>& initial_guess = {2.0, 0.04, 0.3, -0.7, 0.04},
                               int max_iterations = 100) {
        
        double best_error = std::numeric_limits<double>::max();
        std::vector<double> best_params = initial_guess;
        
        // Simple grid search around initial guess (basic optimization)
        std::vector<double> scales = {0.5, 0.8, 1.0, 1.2, 1.5};
        
        for (double scale : scales) {
            std::vector<double> test_params = initial_guess;
            for (auto& param : test_params) {
                param *= scale;
            }
            
            // Ensure valid parameters
            test_params[0] = std::max(0.1, std::min(10.0, test_params[0]));  // kappa
            test_params[1] = std::max(0.01, std::min(1.0, test_params[1])); // theta
            test_params[2] = std::max(0.01, std::min(2.0, test_params[2])); // sigma
            test_params[3] = std::max(-0.99, std::min(0.99, test_params[3])); // rho
            test_params[4] = std::max(0.001, std::min(1.0, test_params[4])); // v0
            
            // Check Feller condition
            if (2 * test_params[0] * test_params[1] <= test_params[2] * test_params[2]) {
                continue;
            }
            
            // Store original parameters
            auto original_params = getParameters();
            
            // Set test parameters
            kappa = test_params[0];
            theta = test_params[1];
            sigma = test_params[2];
            rho = test_params[3];
            v0 = test_params[4];
            
            double error = 0.0;
            bool valid = true;
            
            for (const auto& data_point : market_data) {
                try {
                    double model_price = getOptionPrice(data_point.K, data_point.T, 
                                                       data_point.option_type);
                    error += pow(model_price - data_point.market_price, 2.0);
                } catch (...) {
                    valid = false;
                    break;
                }
            }
            
            if (valid) {
                error = sqrt(error / market_data.size());
                if (error < best_error) {
                    best_error = error;
                    best_params = test_params;
                }
            }
            
            // Restore original parameters
            S0 = original_params[0];
            r = original_params[1];
            q = original_params[2];
            kappa = original_params[3];
            theta = original_params[4];
            sigma = original_params[5];
            rho = original_params[6];
            v0 = original_params[7];
        }
        
        // Set best parameters
        kappa = best_params[0];
        theta = best_params[1];
        sigma = best_params[2];
        rho = best_params[3];
        v0 = best_params[4];
        
        return {best_params, best_error < 1e6, best_error, (int)scales.size()};
    }
    
    // Black-Scholes formula for implied volatility calculation
    double blackScholesCall(double S, double K, double T, double r, double q, double sigma) {
        double d1 = (log(S/K) + (r - q + 0.5*sigma*sigma)*T) / (sigma*sqrt(T));
        double d2 = d1 - sigma*sqrt(T);
        
        // Approximation of cumulative normal distribution
        auto normcdf = [](double x) {
            return 0.5 * (1.0 + erf(x / sqrt(2.0)));
        };
        
        return S*exp(-q*T)*normcdf(d1) - K*exp(-r*T)*normcdf(d2);
    }
    
    // Simple implied volatility calculation using binary search
    double impliedVolatility(double option_price, double S, double K, double T, double r, double q) {
        double low = 0.001, high = 5.0;
        double epsilon = 1e-6;
        
        for (int i = 0; i < 100; i++) {
            double mid = (low + high) / 2.0;
            double bs_price = blackScholesCall(S, K, T, r, q, mid);
            
            if (abs(bs_price - option_price) < epsilon) {
                return mid;
            }
            
            if (bs_price < option_price) {
                low = mid;
            } else {
                high = mid;
            }
        }
        
        return (low + high) / 2.0;
    }
    
    // Generate implied volatility surface
    std::vector<std::vector<double>> 
    impliedVolatilitySurface(const std::vector<double>& strikes, 
                            const std::vector<double>& maturities) {
        std::vector<std::vector<double>> iv_surface(maturities.size(), 
                                                   std::vector<double>(strikes.size()));
        
        for (size_t i = 0; i < maturities.size(); i++) {
            for (size_t j = 0; j < strikes.size(); j++) {
                double option_price = getOptionPrice(strikes[j], maturities[i], "call");
                iv_surface[i][j] = impliedVolatility(option_price, S0, strikes[j], 
                                                   maturities[i], r, q);
            }
        }
        
        return iv_surface;
    }
    
    // Print results to console
    void printOptionPrices(double T, const std::string& filename = "") {
        auto [strikes, call_prices] = carrMadanFFT(T);
        
        std::ostream* out = &std::cout;
        std::ofstream file;
        
        if (!filename.empty()) {
            file.open(filename);
            out = &file;
        }
        
        *out << "Strike,Call Price,Put Price" << std::endl;
        *out << std::fixed << std::setprecision(4);
        
        for (size_t i = 0; i < strikes.size(); i++) {
            double put_price = call_prices[i] - S0 * exp(-q * T) + strikes[i] * exp(-r * T);
            *out << strikes[i] << "," << call_prices[i] << "," << put_price << std::endl;
        }
        
        if (!filename.empty()) {
            file.close();
        }
    }
    
    // Print model parameters
    void printParameters() {
        std::cout << "Heston Model Parameters:" << std::endl;
        std::cout << "S0 (Initial Stock Price): " << S0 << std::endl;
        std::cout << "r (Risk-free Rate): " << r << std::endl;
        std::cout << "q (Dividend Yield): " << q << std::endl;
        std::cout << "kappa (Mean Reversion Speed): " << kappa << std::endl;
        std::cout << "theta (Long-term Variance): " << theta << std::endl;
        std::cout << "sigma (Vol of Vol): " << sigma << std::endl;
        std::cout << "rho (Correlation): " << rho << std::endl;
        std::cout << "v0 (Initial Variance): " << v0 << std::endl;
    }
};

// Example usage and testing
int main() {
    std::cout << "Heston Model with Carr-Madan FFT - C++ Implementation" << std::endl;
    std::cout << "======================================================" << std::endl;
    
    // Initialize Heston model
    HestonCarrMadan heston(100.0, 0.05, 0.0, 2.0, 0.04, 0.3, -0.7, 0.04);
    
    heston.printParameters();
    std::cout << std::endl;
    
    // Test single option pricing
    double K = 100.0, T = 0.25;
    double call_price = heston.getOptionPrice(K, T, "call");
    double put_price = heston.getOptionPrice(K, T, "put");
    
    std::cout << "Single Option Pricing (K=" << K << ", T=" << T << "):" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Call Price: $" << call_price << std::endl;
    std::cout << "Put Price: $" << put_price << std::endl;
    std::cout << std::endl;
    
    // Test FFT pricing for multiple strikes
    auto [strikes, call_prices] = heston.carrMadanFFT(0.25);
    
    std::cout << "FFT Pricing Results (T=0.25):" << std::endl;
    std::cout << "Number of strikes: " << strikes.size() << std::endl;
    std::cout << "Strike range: $" << *std::min_element(strikes.begin(), strikes.end()) 
              << " - $" << *std::max_element(strikes.begin(), strikes.end()) << std::endl;
    std::cout << std::endl;
    
    // Print first few option prices
    std::cout << "Sample Option Prices:" << std::endl;
    std::cout << "Strike\tCall\tPut" << std::endl;
    for (size_t i = 0; i < std::min(size_t(10), strikes.size()); i++) {
        double put = call_prices[i] - heston.getParameters()[0] * exp(-heston.getParameters()[2] * T) + 
                    strikes[i] * exp(-heston.getParameters()[1] * T);
        std::cout << strikes[i] << "\t" << call_prices[i] << "\t" << put << std::endl;
    }
    std::cout << std::endl;
    
    // Example calibration
    std::cout << "Example Model Calibration:" << std::endl;
    std::vector<HestonCarrMadan::MarketData> market_data = {
        {90.0, 0.25, 12.5, "call"},
        {100.0, 0.25, 5.8, "call"},
        {110.0, 0.25, 1.9, "call"},
        {95.0, 0.5, 9.2, "call"},
        {105.0, 0.5, 4.1, "call"}
    };
    
    auto calibration_result = heston.calibrate(market_data);
    std::cout << "Calibration successful: " << (calibration_result.success ? "Yes" : "No") << std::endl;
    std::cout << "RMSE: " << std::setprecision(6) << calibration_result.rmse << std::endl;
    std::cout << "Calibrated parameters:" << std::endl;
    std::cout << "  kappa: " << std::setprecision(4) << calibration_result.params[0] << std::endl;
    std::cout << "  theta: " << calibration_result.params[1] << std::endl;
    std::cout << "  sigma: " << calibration_result.params[2] << std::endl;
    std::cout << "  rho: " << calibration_result.params[3] << std::endl;
    std::cout << "  v0: " << calibration_result.params[4] << std::endl;
    std::cout << std::endl;
    
    // Generate implied volatility surface
    std::cout << "Generating implied volatility surface..." << std::endl;
    std::vector<double> strikes_iv = {80, 90, 100, 110, 120};
    std::vector<double> maturities_iv = {0.1, 0.25, 0.5, 1.0};
    
    auto iv_surface = heston.impliedVolatilitySurface(strikes_iv, maturities_iv);
    
    std::cout << "Implied Volatility Surface:" << std::endl;
    std::cout << "Maturity\\Strike\t";
    for (double K : strikes_iv) {
        std::cout << K << "\t";
    }
    std::cout << std::endl;
    
    for (size_t i = 0; i < maturities_iv.size(); i++) {
        std::cout << maturities_iv[i] << "\t\t";
        for (size_t j = 0; j < strikes_iv.size(); j++) {
            std::cout << std::setprecision(3) << iv_surface[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    
    // Save option prices to CSV file
    std::cout << "\nSaving option prices to 'heston_prices.csv'..." << std::endl;
    heston.printOptionPrices(0.25, "heston_prices.csv");
    
    return 0;
}