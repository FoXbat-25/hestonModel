#ifndef SIM_ENGINE_H
#define SIM_ENGINE_H

#include <vector>
#include <string>


void plott(const std::vector<double>& x, const std::vector<double>& y, const std::string& title);
void plottMultiple(const std::vector<double>& x, const std::vector<std::vector<double>>& y_series, const std::string& title);
struct all_S_v{
    std::vector<std::vector<double>> S;
    std::vector<std::vector<double>> v;
};
all_S_v gbm_sim();


#endif // SIM_ENGINE_H