/********************************************************************/
/*                  SOFTWARE COPYRIGHT NOTIFICATION                 */
/*                             Cardinal                             */
/*                                                                  */
/*                  (c) 2021 UChicago Argonne, LLC                  */
/*                        ALL RIGHTS RESERVED                       */
/*                                                                  */
/*                 Prepared by UChicago Argonne, LLC                */
/*               Under Contract No. DE-AC02-06CH11357               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*             Prepared by Battelle Energy Alliance, LLC            */
/*               Under Contract No. DE-AC07-05ID14517               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*                 See LICENSE for full restrictions                */
/********************************************************************/

#pragma once

#include "ElementVectorPostprocessor.h"

/**
 * A vector post-processor which computes a histogram of the number of
 * elements which contain a variable value within the edges of a bin.
 * This vector post-processor only works with CONSTANT MONOMIAL field
 * variables, as variables with higher order basis functions could be
 * double-counted.
 */
class VariableValueCountHistogram : public ElementVectorPostprocessor
{
public:
  static InputParameters validParams();

  VariableValueCountHistogram(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  /// Number of histogram bins.
  const unsigned int _nbins;

  /// Minimum variable value.
  const Real & _min_value;

  /// Maximum variable value.
  const Real & _max_value;

  /// Bin width.
  Real _delta_var;

  /// Coupled variable that is being binned.
  const VariableValue & _value;

  /// Whether the bins should be log-spaced or not.
  const bool _log;

  /// The left value of the bins.
  std::vector<Real> _left_edge;

  /// The right value of the bins.
  std::vector<Real> _right_edge;

  /// Value mid point of the bin.
  VectorPostprocessorValue & _bin_center;

  /// Aggregated counts for the given bin.
  VectorPostprocessorValue & _counts;
};
