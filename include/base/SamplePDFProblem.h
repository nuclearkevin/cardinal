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

#include <random>

#include "CardinalProblem.h"

class SamplePDFProblem : public CardinalProblem
{
public:
  static InputParameters validParams();

  SamplePDFProblem(const InputParameters & parameters);

  /**
   * TODO
   */
  virtual void externalSolve() override;

  /**
   * TODO
   */
  virtual void syncSolutions(Direction direction) override;

  /**
   * TODO
   */
  virtual void addExternalVariables() override;

  /**
   * TODO
   */
  virtual bool solverSystemConverged(const unsigned int) override;

protected:
  /**
   * TODO
   */
  Point rejectSampleCoordPDF(THREAD_ID tid);

  /**
   * TODO
   */
  Real sampleNumber(THREAD_ID tid);

  /// The number of Monte Carlo samples required when sampling the given CDF.
  const unsigned int & _samples;

  /// The number of retries to use when rejection sampling.
  const unsigned int & _max_num_attempts;

  /// The x position CDF.
  const Function * _x_pdf;

  /// The x-coordinate blanket constant for rejection sampling.
  const Real & _x_blanket;

  /// The y position CDF.
  const Function * _y_pdf;

  /// The y-coordinate blanket constant for rejection sampling.
  Real _y_blanket;

  /// The z position CDF.
  const Function * _z_pdf;

  /// The z-coordinate blanket constant for rejection sampling.
  Real _z_blanket;

  /// Minimum bounding box corner for rejection sampling.
  const Point & _bb_min;

  /// Maximum bounding box corner for rejection sampling.
  const Point & _bb_max;

  /// The extents of the bounding box.
  const Point _bb_extents;

  /// The results variable name.
  const std::string & _result_var_name;

  /// Storage for the intermediate "tally" results.
  std::vector<Real> _pseudo_tally_sum;
  std::vector<Real> _pseudo_tally_sum_sq;

  /// A mapping between tally bins (active elements) and all elements.
  std::vector<dof_id_type> _bin_to_elem_map;

  /// A mapping between elements and tally bins (active elements).
  std::vector<int> _elem_to_bin_map;

  /// The random number generators.
  std::vector<std::default_random_engine> _rng;

  /// Point locators.
  std::vector<std::unique_ptr<libMesh::PointLocatorBase>> _pl;

  /// The number associated with the mean variable added.
  unsigned int _mean_var_number;

  /// The number associated with the standard deviation variable added.
  unsigned int _std_dev_var_number;

  /// The number associated with the relative error variable added.
  unsigned int _rel_var_number;

  /**
   * Time spent computing and storing the solution. This is accumulated
   * over all time / adaptivity steps.
   */
  Real _time;
};
