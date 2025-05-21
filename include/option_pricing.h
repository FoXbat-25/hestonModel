#ifndef OPTION_PRICING_H
#define OPTION_PRICING_H

struct prices{
    double call;
    double put;
};

prices heston_prices_parellel();
double heston_call_price();
prices price();
#endif