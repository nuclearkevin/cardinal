// MOOSE includes
#include "Moose.h"
#include "MooseTypes.h"

// C++ includes
#include <functional>

namespace SecantMethod
{
  /**
   * Compute the roots of a function with the secant method.
   * @param[in] f the function to compute a root of
   * @param[in] x1 the first guess for the root of the function
   * @param[in] x2 the second guess for the root of the function
   * @param[in] tol a convergence tolerance with a default of 1e-12
   * @param[in] max_itr the maximum number of allowed iterations with a default of 100
   */
  Real root(std::function<Real(Real)> const & f, Real x1, Real x2,
            Real tol = 1.0e-12, unsigned int max_itr = 100);
} // namespace SecantMethod
