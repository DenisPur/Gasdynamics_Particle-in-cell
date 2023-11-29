#pragma once
#include <functional>

std::function<double(double, double)> van_der_waals_gas (double a, double b) {
    auto gas_state = [a, b](double eps, double rho)
    {
        return rho * eps / (3 * (1 - b * rho)) - a*a * rho*rho;
    };
    return gas_state;
}

