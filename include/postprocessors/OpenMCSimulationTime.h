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

#include "GeneralPostprocessor.h"

#include "OpenMCBase.h"

/**
 * Get the total simulation time reported by OpenMC. This can be one of the following:
 * - Time spent in inactive batches;
 * - Time spent in active batches;
 * - Total time spent in both active and inactive batches
 */
class OpenMCSimulationTime : public GeneralPostprocessor, public OpenMCBase
{
public:
  static InputParameters validParams();

  OpenMCSimulationTime(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}

  virtual Real getValue() const override;

protected:
  /// An enum to select the time to pluck from OpenMC and report.
  const enum class TimeType
  {
    Inactive = 0,
    Active   = 1,
    Total    = 2
  } _reported_time;
};
