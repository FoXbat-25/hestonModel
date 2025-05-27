#ifndef OPTION_PRICING_H
#define OPTION_PRICING_H

struct pairr{
    double first;
    double second;
};

pairr heston_prices_parellel();
double heston_call_price();
pairr price();
#endif