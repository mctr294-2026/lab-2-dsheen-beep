#include "roots.hpp"
#include <cmath>     // for std::abs
#include <limits>    // for std::numeric_limits

// Required tolerance and iteration limit
static const double TOL = 1e-6; // Tolerance for convergence
static const int MAX_ITER = 1000000; // Maximum iterations

// Bisection Method
// Guaranteed to converge if f(a) and f(b) have opposite signs.

bool bisection(std::function<double(double)> f, 
               double a, double b, // interval [a,b]
               double *root) // root is a pointer where we store the answer
{
    double fa = f(a); // f at the endpoints
    double fb = f(b); // f at the endpoints

    // Root must be bracketed, the condition f(a)*f(b) < 0
    if (fa * fb > 0)
        return false;

    for (int i = 0; i < MAX_ITER; i++) // limit iterations until convergence
    {
        double mid = (a + b) / 2.0; // compute midpoint
        double fmid = f(mid); // evaluate function at midpoint

        // Check convergence, within the defined tolerance
        if (std::abs(fmid) < TOL || std::abs(b - a) < TOL) // Stop if f(mid) is close enough to zero OR the interval is extremely small
        {
            *root = mid; // Save the root approximation (midpoint) into the output pointer
            return true; // Root found successfully
        }

        // Decide which half contains the root
        if (fa * fmid < 0) // If f(a) and f(mid) have opposite signs, the root lies in the left half
        {
            b = mid; // Update b to mid
            fb = fmid; // Update f(b) since b has changed
        }
        else
        {
            a = mid; // Update a to mid
            fa = fmid; // Update f(a) since a has changed
        }
    }

    return false; // If we reach here, maximum iterations were exceeded without convergence
}


// Regula Falsi (False Position)
// Uses linear interpolation instead of midpoint.

bool regula_falsi(std::function<double(double)> f, // function pointer
                  double a, double b, // interval [a,b]
                  double *root) // root is a pointer where we store the answer
{
    double fa = f(a); // f at the endpoints
    double fb = f(b); // f at the endpoints

    if (fa * fb > 0) // Root must be bracketed
        return false; // No root in [a,b]

    for (int i = 0; i < MAX_ITER; i++) // limit iterations until convergence
    {
        // False-position formula
        double c = (a * fb - b * fa) / (fb - fa); // Compute the point where the line through (a,f(a)) and (b,f(b)) crosses the x-axis
        double fc = f(c); // Evaluate function at c

        if (std::abs(fc) < TOL) // Check convergence
        {
            *root = c; // Save the root approximation into the output pointer
            return true; // Root found successfully
        }

        // Update bracket
        if (fa * fc < 0) // Root lies in [a,c]
        {
            b = c; // Update b to c
            fb = fc; // Update f(b) since b has changed
        }
        else
        {
            a = c; // Update a to c
            fa = fc; // Update f(a) since a has changed
        }

        if (std::abs(b - a) < TOL) // Interval is sufficiently small
        {
            *root = (a + b) / 2.0; // Save midpoint as root approximation
            return true; // Root found successfully
        }
    }

    return false; // Maximum iterations exceeded without convergence
}


// Newton-Raphson Method
// Requires derivative g(x) and starts at initial guess c

bool newton_raphson(std::function<double(double)> f, // function pointer
                    std::function<double(double)> g, // derivative pointer
                    double a, double b, double c, // interval [a,b] and initial guess c
                    double *root) // root is a pointer where we store the answer
{
    double x = c; // initial guess c

    for (int i = 0; i < MAX_ITER; i++) // limit iterations until convergence
    {
        double fx = f(x); // evaluate function at current guess
        double gx = g(x); // evaluate derivative at current guess

        // Converged
        if (std::abs(fx) < TOL) // Check convergence with the given tolerance
        {
            *root = x; // Save the root approximation into the output pointer
            return true; // Root found successfully
        }

        // Derivative cannot be zero
        if (std::abs(gx) < 1e-12) // Prevent division by zero error
            return false; // Failure due to zero derivative

        // Newton step
        double x_next = x - fx / gx; // Update guess using Newton-Raphson formula

        // Must stay inside [a,b]
        if (x_next < a || x_next > b) // Check if the next guess is within the interval, otherwise fail
            return false; // Failure due to out-of-bounds of a and b

        x = x_next; // Move to next guess
    }

    return false; // Maximum iterations exceeded without convergence
}

// Secant Method
// Approximates derivative using two points and starts with x0 = c and x1 = (c + small offset)
 
bool secant(std::function<double(double)> f, // function pointer
            double a, double b, double c, // interval [a,b] and initial guess c
            double *root) // root is a pointer where we store the answer
{
    double x0 = c; // the initial guess c
    double x1 = c + 1e-3; // small offset vlalue we use of 0.001

    // Ensure x1 stays inside interval
    if (x1 > b) // If offset goes beyond b
        x1 = c - 1e-3; // Move offset to the left side

    for (int i = 0; i < MAX_ITER; i++) // limit iterations until convergence
    {
        double f0 = f(x0); // f at previous point
        double f1 = f(x1); // f at current point

        if (std::abs(f1) < TOL) // Check convergence at given tolerance
        {
            *root = x1; // Save the root approximation into the output pointer
            return true; // Root found successfully
        }

        // Prevent division by zero
        double denom = (f1 - f0); // Denominator of the secant formula
        if (std::abs(denom) < 1e-12) // If denominator is too small
            return false; // Failure due to potential division by zero

        // Secant formula
        double x2 = x1 - f1 * (x1 - x0) / denom; // Compute next approximation using secant method

        // Must remain in interval
        if (x2 < a || x2 > b) // Check if the next guess is within the interval
            return false; // Failure due to out-of-bounds of a and b

        // Shift forward
        x0 = x1; // Move previous point to current
        x1 = x2; // Move current point to next
    } 

    return false; // Maximum iterations exceeded without convergence
}