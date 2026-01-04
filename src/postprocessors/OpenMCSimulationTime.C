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

#ifdef ENABLE_OPENMC_COUPLING

#include "OpenMCSimulationTime.h"

#include "openmc/timer.h"

registerMooseObject("CardinalApp", OpenMCSimulationTime);

InputParameters
OpenMCSimulationTime::validParams()
{
  auto params = GeneralPostprocessor::validParams();
  params += OpenMCBase::validParams();
  params.addClassDescription(
      "A class which reports the time spent in an OpenMC simulation (in seconds).");
  params.addParam<MooseEnum>(
      "reported_time",
      MooseEnum("inactive active total", "total"),
      "The time to report. Options are the time spent in inactive batches, "
      "the time spent in active batches, or the total OpenMC simulation time.");

  return params;
}

OpenMCSimulationTime::OpenMCSimulationTime(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    OpenMCBase(this, parameters),
    _reported_time(getParam<MooseEnum>("reported_time").getEnum<TimeType>())
{ }

Real
OpenMCSimulationTime::getValue() const
{
  switch (_reported_time)
  {
    case TimeType::Inactive:
      return openmc::simulation::time_inactive.elapsed();
    case TimeType::Active:
      return openmc::simulation::time_active.elapsed();
    case TimeType::Total:
      return openmc::simulation::time_inactive.elapsed() + openmc::simulation::time_active.elapsed();
    default:
    {
      mooseError("Unhandled TimeType enum in OpenMCSimulationTime!");
      return 0.0;
    }
  }
}

#endif
