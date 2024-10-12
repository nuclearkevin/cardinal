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

#include "Indicator.h"

#include "OpenMCCellAverageProblem.h"

class OpenMCIndicator : public Indicator
{
public:
  static InputParameters validParams();

  OpenMCIndicator(const InputParameters & parameters);

protected:
  /// The OpenMCCellAverageProblem associated with this indicator.
  const OpenMCCellAverageProblem * _openmc_problem;

  /// The field variable holding the results of this indicator.
  MooseVariable & _field_var;

  /// The current element.
  const Elem * const & _current_elem;
};
