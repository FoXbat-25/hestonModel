# hestonModel

The model assumes that volatility follows a mean-reverting process, which means it tends to revert to a long-term average over time.

This feature of the model allows it to capture the volatility smile observed in the market, where options with different strike prices but the same maturity may have different implied volatilities.

--

## Heston model parameters

The Heston model has several parameters that describe the dynamics of the underlying asset's price and volatility.

The main parameters of the Heston model are:

Initial asset price (S0): The current price of the underlying asset.
Mean reversion rate (κ): The speed at which the volatility reverts to its long-term average.
Long-term average volatility (θ): The long-term average volatility level to which the volatility reverts.
Volatility of volatility (ν): The volatility of the volatility process. It determines the amplitude of volatility fluctuations.
Correlation between asset price and volatility (ρ): The correlation between the asset price and its volatility process. This parameter determines how changes in the asset price affect its volatility.
Risk-free interest rate (r): The risk-free interest rate.
Time to maturity (T): The time until the option's expiration.
Strike price (K): The price at which the option holder has the right to buy or sell the underlying asset.

--

## Heston model stochastic equations

1. The first equation represents Stock Price Dynamics, and is as follows:

dS(t) = µS(t)dt + √v(t)S(t)dW₁(t)

This equation describes the logarithmic price movement of the underlying asset. Here's a breakdown of the terms:

dS(t): Infinitesimal change in the price of the underlying asset at time t.
µ: Coefficient, representing the expected return of the asset per unit time.
S(t): Price of the underlying asset at time t.
dt: Infinitesimal time change.
√v(t): Volatility factor, representing the standard deviation of the asset's return per unit time. This term incorporates the concept of stochastic volatility, meaning the volatility can fluctuate over time.
dW₁(t): Infinitesimal Wiener process (Brownian motion), representing a random shock or innovation term.

This equation captures the price movement of the asset, considering both the expected return and random fluctuations influenced by volatility.

2. The second equation represents Volatility Dynamics and is as follows:

dv(t) = κ(θ - v(t))dt + συ(t)dW2(t)

Here's a breakdown of the terms:

dv(t): Infinitesimal change in the volatility of the asset at time t.
κ: Mean reversion rate, representing the speed at which volatility tends to move back towards its long-term average (θ).
θ: Long-term average volatility.
v(t): Volatility of the asset at time t.
dt: Infinitesimal time change.
σ: Volatility of volatility, representing the standard deviation of the volatility process.
dW2(t): Infinitesimal Wiener process (Brownian motion), independent of the first Wiener process (dW₁(t)), representing a random shock affecting the volatility.
This equation models the volatility itself as a separate process with its own mean reversion and random fluctuations.

--

## Assumptions while using the Heston model

When using the Heston model, several assumptions are made:

Continuous trading: Trading occurs continuously, and the market is frictionless.
No dividends: The underlying asset does not pay any dividends during the option's life.
No transaction costs: There are no transaction costs associated with trading the option or the underlying asset.
No arbitrage opportunities: There are no risk-free arbitrage opportunities in the market.
Constant risk-free rate: The risk-free interest rate is constant and known.
Constant parameters: The parameters of the model (such as mean reversion rate, long-term average volatility, volatility of volatility, and correlation) are assumed to be constant over time.
However, it's important to note that in reality, some of these assumptions may not hold true, and adjustments may be necessary when applying the model to real-world situations.