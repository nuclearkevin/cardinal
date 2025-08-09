#include "RootFinding.h"

namespace SecantMethod
{
  Real
  root(std::function<Real(Real)> const & f, Real x1, Real x2, Real tol, unsigned int max_itr)
  {
    // The two initial guesses must be different when using secant method.
    if (x1 == x2)
      throw MooseException("x1 and x2 are the same in SecantMethod::root.");

    Real x_a = x1;
    Real x_b = x2;
    Real x_c;

    Real f_a = f(x_a);
    Real f_b = f(x_b);

    for (unsigned int i = 1; i <= max_itr; ++i)
    {
      x_c = x_b - f_b * (x_b - x_a) / (f_b - f_a);

      x_a = x_b;
      f_a = f_b;
      x_b = x_c;
      f_b = f(x_c);

      // Return the root as we've hit the convergence criteria.
      if (std::abs(x_c - x_b) < tol)
        return x_c;
    }

    throw MooseException("Maximum number of iterations exceeded in SecantMethod::root.");

    return 0.0;
  }
} // namespace SecantMethod
