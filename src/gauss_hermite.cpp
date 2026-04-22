#include "gauss_hermite.h"
#include <cmath>
#include <utility>

/*************************************************************************
Computation of nodes and weights for a Gauss-Hermite quadrature formula

The algorithm calculates the nodes and weights of the Gauss-Hermite
quadrature formula on domain (-infinity, +infinity) with weight function
W(x)=Exp(-x*x).

Input parameters:
    n   -   a required number of nodes.
            1 <= n <= 190.

Output parameters:
    x   -   vector of nodes (resized to n).
    w   -   vector of weighting coefficients (resized to n).

The algorithm was designed by using information from the QUADRULE library.
*************************************************************************/

void calcgausshermitequadrature(int n, std::vector<double>& x, std::vector<double>& w)
{
    // Resize output vectors
    x.resize(n);
    w.resize(n);

    int i;
    int j;
    double r = 0;
    double r1;
    double p1;
    double p2;
    double p3;
    double dp3;
    double pipm4;
    const double pi = 3.14159265358979323846;

    pipm4 = std::pow(pi, -0.25);

    for(i = 0; i <= (n+1)/2-1; i++)
    {
        // Initial guess for root
        if( i==0 )
        {
            r = std::sqrt(double(2*n+1)) - 1.85575 * std::pow(double(2*n+1), -double(1)/double(6));
        }
        else
        {
            if( i==1 )
            {
                r = r - 1.14 * std::pow(double(n), 0.426) / r;
            }
            else
            {
                if( i==2 )
                {
                    r = 1.86*r - 0.86*x[0];
                }
                else
                {
                    if( i==3 )
                    {
                        r = 1.91*r - 0.91*x[1];
                    }
                    else
                    {
                        r = 2*r - x[i-2];
                    }
                }
            }
        }

        // Newton iteration to refine root
        do
        {
            p2 = 0;
            p3 = pipm4;
            for(j = 0; j <= n-1; j++)
            {
                p1 = p2;
                p2 = p3;
                p3 = p2*r*std::sqrt(double(2)/double(j+1)) - p1*std::sqrt(double(j)/double(j+1));
            }
            dp3 = std::sqrt(double(2*j)) * p2;
            r1 = r;
            r = r - p3/dp3;
        }
        while((std::fabs(r-r1) >= std::fabs(r)*1e-15) && (std::fabs(r-r1) >= 1e-100));

        // Store node and weight, using symmetry
        x[i] = r;
        w[i] = 2 / (dp3*dp3);
        x[n-1-i] = -x[i];
        w[n-1-i] = w[i];
    }
}

std::pair<std::vector<double>, std::vector<double>> gauss_hermite_nodes_weights(int n)
{
    std::vector<double> x, w;
    calcgausshermitequadrature(n, x, w);
    return std::make_pair(x, w);
}
