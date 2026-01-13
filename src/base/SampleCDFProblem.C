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

#include "SampleCDFProblem.h"

#include "Function.h"

#include <omp.h>

registerMooseObject("CardinalApp", SampleCDFProblem);

InputParameters
SampleCDFProblem::validParams()
{
  auto params = CardinalProblem::validParams();
  params.addClassDescription("A problem which samples a given cumulative density function with "
                             "Monte Carlo integration.");
  params.addRequiredParam<std::string>(
      "result_var_name",
      "The variable name to use for the mean value. '_std_dev' is appended to this to store the "
      "standard deviation, and '_rel_error' is appended to store the relative error.");
  params.addRequiredParam<FunctionName>("x_coord_cdf",
                                        "The cumulative density function for the x position.");
  params.addParam<FunctionName>("y_coord_cdf",
                                "The cumulative density function for the y position.");
  params.addParam<FunctionName>("z_coord_cdf",
                                "The cumulative density function for the z position.");
  params.addRequiredParam<unsigned int>(
      "samples",
      "The number of samples to use when sampling the given cumulative density function.");

  return params;
}

SampleCDFProblem::SampleCDFProblem(const InputParameters & parameters)
  : CardinalProblem(parameters),
    _samples(getParam<unsigned int>("samples")),
    _x_cdf(getFunction("x_coord_cdf")),
    _result_var_name(getParam<std::string>("result_var_name"))
{
  if (mesh().spatialDimension() > 2)
  {
    if (!isParamValid("y_coord_cdf"))
      paramError("y_coord_cdf", "In 2D or 3D, y_coord_cdf is required!");
    else
      _y_cdf = &getFunction("y_coord_cdf");

    if (mesh().spatialDimension() > 3)
    {
      if (!isParamValid("z_coord_cdf"))
        paramError("z_coord_cdf", "In 3D, z_coord_cdf is required!");
      else
        _z_cdf = &getFunction("z_coord_cdf");
    }
  }

  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    _rng.emplace_back(tid);
    _pl.emplace_back(mesh().getMeshPtr()->sub_point_locator());
    _pl.back()->set_contains_point_tol(1e-8);
    _pl.back()->enable_out_of_mesh_mode();
  }
}

void
SampleCDFProblem::externalSolve()
{
  auto msh = mesh().getMeshPtr();
  const auto num_active_elem = mesh().nActiveElem();

  // Reset "tally" data structures.
  _pseudo_tally_sum.clear();
  _pseudo_tally_sum.resize(num_active_elem, 0.0);

  // Reset mapping data structures.
  _bin_to_elem_map.clear();
  _bin_to_elem_map.reserve(num_active_elem);
  _elem_to_bin_map.clear();
  _elem_to_bin_map.resize(mesh().nElem(), -1);
  for (const auto & elem : libMesh::as_range(msh->active_elements_begin(), msh->active_elements_end()))
  {
    _bin_to_elem_map.push_back(elem->id());
    _elem_to_bin_map[elem->id()] = _bin_to_elem_map.size() - 1;
  }

  // CDF sampling loop.
  #pragma omp parallel for
  for (unsigned int i = 0; i < _samples; ++i)
  {
    const auto tid = omp_get_thread_num();
    // Treat time (t) as the uniform random number.
    const Real x = _x_cdf.value(sampleNumber(tid), Point(0.0, 0.0, 0.0));
    const Real y = _y_cdf ? _y_cdf->value(sampleNumber(tid), Point(0.0, 0.0, 0.0)) : 0.0;
    const Real z = _z_cdf ? _z_cdf->value(sampleNumber(tid), Point(0.0, 0.0, 0.0)) : 0.0;

    // Sample the mesh to find a bin.
    const auto elem = (*_pl[tid])(Point(x, y, z));
    if (!elem)
      continue;

    // Store it in a "tally" bin.
    const auto bin = _elem_to_bin_map[elem->id()];
    #pragma omp atomic
    _pseudo_tally_sum[bin] += 1.0;
  }

  // Normalize the results by the number of samples.
  #pragma omp parallel for
  for (unsigned int bin = 0; bin < _pseudo_tally_sum.size(); ++bin)
    _pseudo_tally_sum[bin] /= static_cast<Real>(_samples);
}

void
SampleCDFProblem::syncSolutions(Direction direction)
{
  _aux->serializeSolution();

  switch (direction)
  {
    case ExternalProblem::Direction::TO_EXTERNAL_APP:
      break;
    case ExternalProblem::Direction::FROM_EXTERNAL_APP:
      #pragma omp parallel for
      for (unsigned int bin = 0; bin < _pseudo_tally_sum.size(); ++bin)
      {
        const auto elem = _mesh.queryElemPtr(_bin_to_elem_map[bin]);
        if (!elem)
          continue;

        auto dof_idx = elem->dof_number(_aux->number(), _var_number, 0);
        _aux->solution().set(dof_idx, _pseudo_tally_sum[bin]);
      }
      break;
  }

  _aux->solution().close();
  _aux->system().update();
}

void
SampleCDFProblem::addExternalVariables()
{
  auto var_params = _factory.getValidParams("MooseVariable");
  var_params.set<MooseEnum>("family") = "MONOMIAL";
  var_params.set<MooseEnum>("order") = "CONSTANT";

  checkDuplicateVariableName(_result_var_name, "test");
  addAuxVariable("MooseVariable", _result_var_name, var_params);

  _var_number = _aux->getFieldVariable<Real>(0, _result_var_name).number();
}

Real
SampleCDFProblem::sampleNumber(THREAD_ID tid)
{
  auto dis = std::uniform_real_distribution<Real>(0.0, 1.0);
  return dis(_rng[tid]);
}
