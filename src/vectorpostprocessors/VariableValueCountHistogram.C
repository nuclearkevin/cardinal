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

#include "VariableValueCountHistogram.h"

#include "MooseVariableFE.h"

registerMooseObject("CardinalApp", VariableValueCountHistogram);

InputParameters
VariableValueCountHistogram::validParams()
{
  auto params = ElementVectorPostprocessor::validParams();
  params.addClassDescription(
      "Compute a histogram of element counts binned according to variable values.");
  params.addCoupledVar("variable", "Variable to bin element counts by. Must be of type CONSTANT MONOMIAL.");
  params.addParam<unsigned int>("bin_number", 50, "Number of histogram bins");
  params.addRequiredParam<PostprocessorName>("min_value", "A post-processor computing the minimum variable value.");
  params.addRequiredParam<PostprocessorName>("max_value", "A post-processor computing the maximum variable value.");
  params.addParam<bool>("log_space", false, "Whether the histogram should be log spaced or not.");

  return params;
}

VariableValueCountHistogram::VariableValueCountHistogram(const InputParameters & parameters)
  : ElementVectorPostprocessor(parameters),
    _nbins(getParam<unsigned int>("bin_number")),
    _min_value(getPostprocessorValue("min_value")),
    _max_value(getPostprocessorValue("max_value")),
    _value(coupledValue("variable")),
    _log(getParam<bool>("log_space")),
    _bin_center(declareVector(coupledName("variable"))),
    _counts(declareVector("counts"))
{
  if (coupledComponents("variable") != 1)
    paramError("variable", "VariableValueCountHistogram only supports variables with a single component!");

  if (getFieldVar("variable", 0)->feType() != FEType(libMesh::CONSTANT, libMesh::MONOMIAL))
    paramError("variable",
      "VariableValueCountHistogram only supports CONSTANT MONOMIAL shape functions. Please "
      "ensure that 'variable' is of type MONOMIAL and order CONSTANT.");

  // initialize the bin center value vector
  _bin_center.resize(_nbins);
  _left_edge.resize(_nbins);
  _right_edge.resize(_nbins);
}

void
VariableValueCountHistogram::initialize()
{
  // Reset the histogram.
  _counts.assign(_nbins, 0.0);

  if (_log)
  {
    const Real log_max = std::log10(_max_value);
    const Real log_min = std::log10(_min_value);
    const Real delta = (log_max - log_min) / static_cast<Real>(_nbins);
    for (unsigned int i = 0; i < _nbins; ++i)
    {
      _left_edge[i] = std::pow(10.0, i * delta + log_min);
      _bin_center[i] = std::pow(10.0, (i + 0.5) * delta + log_min);
      _right_edge[i] = std::pow(10.0, (i + 1.0) * delta + log_min);
    }
  }
  else
  {
    const Real delta = (_max_value - _min_value) / static_cast<Real>(_nbins);
    for (unsigned int i = 0; i < _nbins; ++i)
    {
      _left_edge[i] = i * delta + _min_value;
      _bin_center[i] = (i + 0.5) * delta + _min_value;
      _right_edge[i] = (i + 1.0) * delta + _min_value;
    }
  }
}

void
VariableValueCountHistogram::execute()
{
  // Linear search for the bin.
  // TODO: binary search
  for (unsigned int bin = 0; bin < _nbins; ++bin)
  {
    if (_value[0] >= _left_edge[bin] && _value[0] < _right_edge[bin])
    {
      _counts[bin] += 1.0;
      break;
    }
  }
}

void
VariableValueCountHistogram::finalize()
{
  gatherSum(_counts);
}

void
VariableValueCountHistogram::threadJoin(const UserObject & y)
{
  const auto & uo = static_cast<const VariableValueCountHistogram &>(y);
  mooseAssert(uo._counts.size() == _counts.size(),
              "Inconsistent counts vector lengths across threads.");

  for (unsigned int i = 0; i < _counts.size(); ++i)
    _counts[i] += uo._counts[i];
}
