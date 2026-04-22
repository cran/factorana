#ifndef GAUSS_HERMITE_H
#define GAUSS_HERMITE_H

#include <vector>

// Compute Gauss-Hermite quadrature nodes and weights
// for integrating functions against the weight function W(x) = exp(-x^2)
//
// Parameters:
//   n - number of quadrature points (1 <= n <= 190)
//   x - output vector of nodes (size n)
//   w - output vector of weights (size n)
//
// Note: For integrating against the standard normal density N(0,1),
// you need to rescale: x_scaled = sqrt(2) * x, w_scaled = w / sqrt(pi)

void calcgausshermitequadrature(int n, std::vector<double>& x, std::vector<double>& w);

// Convenience function that returns a pair of vectors
std::pair<std::vector<double>, std::vector<double>> gauss_hermite_nodes_weights(int n);

#endif // GAUSS_HERMITE_H
