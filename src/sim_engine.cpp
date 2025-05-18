#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include <constants.h>

using namespace std;

vector<double> rng(double mean, double dt, int N){
    
    std::random_device rd;     // Only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // Random-number engine used (Mersenne-Twister in this case)
    std::normal_distribution<> dist(mean ,sqrt(dt)); // Guaranteed unbiased

    std :: vector <double> dW(N);

    for (int i = 0; i<N; ++i){
        dW[i] = dist(rng);
    }
    return dW;
}

vector <double> linspace(double start, double end, int num) {
    
    vector <double> result;

    if (num == 0){
        return result;
    }

    if (num == 1){
        result.push_back(start);
        return result;
    }

    double step = (end - start)/(num -1);
    for (int i = 0; i < num; ++i){
        result.push_back(start + i*step);
    }

    return result;
}

const std::vector <double> t = linspace(0.0, T, N+1);

void plott(const vector<double>& x, const vector<double>& y, const string& title) {
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");

    if (!gnuplotPipe) {
        cerr<< "Couldn't open gnuplot. \n";
        return;
    }

    fprintf(gnuplotPipe, "set title '%s'\n", title.c_str());
    fprintf(gnuplotPipe, "plot '-' with lines title '%s'\n", title.c_str());

    for (size_t i = 0; i < x.size(); ++i) {
        fprintf(gnuplotPipe, "%f %f\n", x[i], y[i]);
    }

    fprintf(gnuplotPipe, "e\n");
    fflush(gnuplotPipe);
    pclose(gnuplotPipe);
}

void plottMultiple(const vector<double>& x, const vector<vector<double>>& y_series, const string& title) {
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");

    if (!gnuplotPipe) {
        cerr << "Couldn't open gnuplot.\n";
        return;
    }

    fprintf(gnuplotPipe, "set title '%s'\n", title.c_str());
    fprintf(gnuplotPipe, "plot ");

    for (size_t i = 0; i < y_series.size(); ++i) {
        if (i > 0) fprintf(gnuplotPipe, ", ");
        fprintf(gnuplotPipe, "'-' with lines notitle"); //title 'Path %lu'", i+1
    }
    fprintf(gnuplotPipe, "\n");

    for (const auto& y : y_series) {
        for (size_t j = 0; j < x.size(); ++j) {
            fprintf(gnuplotPipe, "%f %f\n", x[j], y[j]);
        }
        fprintf(gnuplotPipe, "e\n");
    }

    fflush(gnuplotPipe);
    pclose(gnuplotPipe);
}

struct PriceStats {
    double avg;
    std::vector<double> prices;
};

PriceStats calc_fmetrics(int num_paths, int N, vector<vector<double>> all_S){
    double avg_final_price = 0;
    vector <double> prices;
    for (int p = 0; p < num_paths; ++p){
        avg_final_price += all_S[p][N];
        prices.push_back(all_S[p][N]);
    }
    avg_final_price /= num_paths;

    PriceStats result;
    result.avg = avg_final_price;
    result.prices = prices;
    return {avg_final_price, prices}; 
}

double price_Kplus_prob(const std::vector<double>& prices, double K){
    
    int count = 0;
    for (int i =0; i<prices.size(); ++i){
        if (prices[i] > K) count ++;
    }

    return static_cast<double> (count)/prices.size();
}

double trap_integration(std::function<double(double)> f, double a, double b, int n){
    double h = (b - a)/n;
    double summ_ht = 0.5 * (f(a) + f(b)); //Area under the curve
    for (int i=1; i<=n; ++i){
        summ_ht += f(a + i *h);
    }
    double auc = summ_ht * h;
    return auc;
}

all_S_v gbm_sim(){

    auto start = std::chrono::high_resolution_clock::now();
    
    vector<vector<double>> all_S(num_paths, vector<double>(N+1));
    vector<vector<double>> all_v(num_paths, vector<double>(N+1));

    for (int p = 0; p < num_paths; ++p) {
        vector<double> dw1 = rng(0.0, dt, N);
        vector<double> dw2 = rng(0.0, dt, N);

        all_S[p][0] = S0;
        all_v[p][0] = v0;

        for (int i = 1; i <= N; ++i) {
            all_v[p][i] = all_v[p][i-1] + kappa * (theta - all_v[p][i-1]) * dt + sigma * sqrt(all_v[p][i-1]) * dw1[i-1];
            all_S[p][i] = all_S[p][i-1] * (1 + mu * dt + sqrt(all_v[p][i-1]) * dw2[i-1]);
        }
    }

    PriceStats stats = calc_fmetrics(num_paths, N, all_S);
    std::cout << "Average final price: " << stats.avg << std::endl;

    cout<<fixed<<setprecision(4);
    double prob = price_Kplus_prob(stats.prices, K);
    std::cout << "Probability of price above " << K << " is " << prob <<endl;


    //Timer ends
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "GBM execution time: " << duration.count() << " ms" << std::endl;

    all_S_v gbm_stats;
    gbm_stats.S = all_S;
    gbm_stats.v = all_v;
    return gbm_stats;
}