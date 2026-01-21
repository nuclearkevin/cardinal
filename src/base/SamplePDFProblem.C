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

#include "SamplePDFProblem.h"

#include "Function.h"

#include <chrono>
#include <omp.h>

registerMooseObject("CardinalApp", SamplePDFProblem);

InputParameters
SamplePDFProblem::validParams()
{
  auto params = CardinalProblem::validParams();
  params.addClassDescription("A problem which samples a given cumulative density function with "
                             "Monte Carlo integration.");
  params.addRequiredParam<std::string>(
      "result_var_name",
      "The variable name to use for the mean value. '_std_dev' is appended to this to store the "
      "standard deviation, and '_rel_error' is appended to store the relative error.");
  params.addRequiredParam<FunctionName>("x_coord_pdf",
                                        "The cumulative density function for the x position.");
  params.addParam<FunctionName>("y_coord_pdf",
                                "The cumulative density function for the y position.");
  params.addParam<FunctionName>("z_coord_pdf",
                                "The cumulative density function for the z position.");

  params.addRequiredParam<Real>("x_blanket", "The envelope around the target PDF in x to use for rejection sampling.");
  params.addParam<Real>("y_blanket", "The envelope around the target PDF in y to use for rejection sampling.");
  params.addParam<Real>("z_blanket", "The envelope around the target PDF in z to use for rejection sampling.");
  params.addRequiredParam<Point>("bb_min", "The minimum point to use for computing the rejection sampling bounding box.");
  params.addRequiredParam<Point>("bb_max", "The maximum point to use for computing the rejection sampling bounding box.");

  params.addRequiredParam<unsigned int>(
      "samples",
      "The number of samples to use when sampling the given cumulative density function.");

  params.addParam<unsigned int>(
      "max_attempts",
      100,
      "The maximum number of attempts to use when rejection sampling the given PDFs. Defaults to 100.");

  return params;
}

SamplePDFProblem::SamplePDFProblem(const InputParameters & parameters)
  : CardinalProblem(parameters),
    _samples(getParam<unsigned int>("samples")),
    _max_num_attempts(getParam<unsigned int>("max_attempts")),
    _x_blanket(getParam<Real>("x_blanket")),
    _bb_min(getParam<Point>("bb_min")),
    _bb_max(getParam<Point>("bb_max")),
    _bb_extents(_bb_max - _bb_min),
    _result_var_name(getParam<std::string>("result_var_name")),
    _time(0.0)
{
  if (n_processors() > 1)
    mooseError("SamplePDFProblem cannot be executed with MPI!");
}

void
SamplePDFProblem::externalSolve()
{
  auto t_start = std::chrono::high_resolution_clock::now();

  auto msh = mesh().getMeshPtr();
  const auto num_active_elem = mesh().nActiveElem();

  // Reset RNG and point locators. Needed to ensure the mesh changes are taken into
  // account when "tallying", and that each adaptivity step sees the same RNG stream.
  _rng.clear();
  _pl.clear();
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    _rng.emplace_back(tid);
    _pl.emplace_back(mesh().getMeshPtr()->sub_point_locator());
    _pl.back()->set_contains_point_tol(1e-8);
    _pl.back()->enable_out_of_mesh_mode();
  }

  // Reset "tally" data structures.
  _pseudo_tally_sum.clear();
  _pseudo_tally_sum.resize(num_active_elem, 0.0);
  _pseudo_tally_sum_sq.clear();
  _pseudo_tally_sum_sq.resize(num_active_elem, 0.0);

  // Reset mapping data structures.
  _bin_to_elem_map.clear();
  _bin_to_elem_map.reserve(num_active_elem);
  _elem_to_bin_map.clear();
  _elem_to_bin_map.resize(mesh().nElem(), -1);
  for (const auto & elem :
       libMesh::as_range(msh->active_elements_begin(), msh->active_elements_end()))
  {
    _bin_to_elem_map.push_back(elem->id());
    _elem_to_bin_map[elem->id()] = _bin_to_elem_map.size() - 1;
  }

// CDF sampling loop.
#pragma omp parallel for
  for (unsigned int i = 0; i < _samples; ++i)
  {
    const auto tid = omp_get_thread_num();
    const auto p = sampleNumber(tid);
    const auto pp = p * p;

    // Use rejection sampling (with a uniform distribution) to sample the given PDFs.
    const auto sample_point = rejectSampleCoordPDF(tid);

    // Sample the mesh to find a bin.
    const auto elem = (*_pl[tid])(sample_point);
    if (!elem)
      continue;

    // Store it in a "tally" bin.
    const auto bin = _elem_to_bin_map[elem->id()];
#pragma omp atomic
    _pseudo_tally_sum[bin] += 1.0;

#pragma omp atomic
    _pseudo_tally_sum_sq[bin] += 1.0;
  }

  auto t_end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count();
  _time += 1e-6 * static_cast<Real>(duration);
}

void
SamplePDFProblem::syncSolutions(Direction direction)
{
  switch (direction)
  {
    case ExternalProblem::Direction::TO_EXTERNAL_APP:
      break;
    case ExternalProblem::Direction::FROM_EXTERNAL_APP:
    {
      auto t_start = std::chrono::high_resolution_clock::now();
      _aux->serializeSolution();

      const auto realizations = static_cast<Real>(_samples);

#pragma omp parallel for
      for (unsigned int bin = 0; bin < _pseudo_tally_sum.size(); ++bin)
      {
        const auto elem = _mesh.queryElemPtr(_bin_to_elem_map[bin]);
        if (!elem)
          continue;

        const auto mean = _pseudo_tally_sum[bin] / realizations;
        const auto sum_sq = _pseudo_tally_sum_sq[bin];
        const auto std_dev =
            std::sqrt(std::max(0.0, (sum_sq / realizations - mean * mean) / (realizations - 1)));

        auto mean_dof_idx = elem->dof_number(_aux->number(), _mean_var_number, 0);
        _aux->solution().set(mean_dof_idx, mean / elem->volume());

        auto std_dev_dof_idx = elem->dof_number(_aux->number(), _std_dev_var_number, 0);
        _aux->solution().set(std_dev_dof_idx, std_dev / elem->volume());

        auto rel_dof_idx = elem->dof_number(_aux->number(), _rel_var_number, 0);
        _aux->solution().set(rel_dof_idx, std_dev / mean);
      }

      _aux->solution().close();
      _aux->system().update();

      auto t_end = std::chrono::high_resolution_clock::now();
      auto duration =
          std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count();
      _time += 1e-6 * static_cast<Real>(duration);

      setPostprocessorValueByName("cdf_sampler_runtime", _time);
    }
    break;
  }
}

void
SamplePDFProblem::addExternalVariables()
{
  // Fetch the inverted CDFs to sample positions.
  _x_pdf = &getFunction(getParam<FunctionName>("x_coord_pdf"));
  if (mesh().spatialDimension() > 2)
  {
    if (!isParamValid("y_coord_pdf"))
      paramError("y_coord_pdf", "In 2D or 3D, y_coord_pdf is required!");
    else
      _y_pdf = &getFunction(getParam<FunctionName>("y_coord_pdf"));

    if (!isParamValid("y_blanket"))
      paramError("y_blanket", "In 2D or 3D, y_blanket is required!");
    else
      _y_blanket = getParam<Real>("y_blanket");

    if (mesh().spatialDimension() > 3)
    {
      if (!isParamValid("z_coord_pdf"))
        paramError("z_coord_pdf", "In 3D, z_coord_pdf is required!");
      else
        _z_pdf = &getFunction(getParam<FunctionName>("z_coord_pdf"));

      if (!isParamValid("z_blanket"))
        paramError("z_blanket", "In 2D or 3D, z_blanket is required!");
      else
        _z_blanket = getParam<Real>("z_blanket");
    }
  }

  // Add a post-processor to store the "simulation" time
  auto pp_params = _factory.getValidParams("Receiver");
  addPostprocessor("Receiver", "cdf_sampler_runtime", pp_params);

  // Add a variable to store the volumetric mean value.
  auto var_params = _factory.getValidParams("MooseVariable");
  var_params.set<MooseEnum>("family") = "MONOMIAL";
  var_params.set<MooseEnum>("order") = "CONSTANT";

  checkDuplicateVariableName(_result_var_name + "_mean", "test");
  addAuxVariable("MooseVariable", _result_var_name + "_mean", var_params);
  _mean_var_number = _aux->getFieldVariable<Real>(0, _result_var_name + "_mean").number();

  // Add a variable to store the volumetric standard deviation value.
  checkDuplicateVariableName(_result_var_name + "_std_dev", "test");
  addAuxVariable("MooseVariable", _result_var_name + "_std_dev", var_params);
  _std_dev_var_number = _aux->getFieldVariable<Real>(0, _result_var_name + "_std_dev").number();

  // Add a variable to store the relative error value.
  checkDuplicateVariableName(_result_var_name + "_rel_error", "test");
  addAuxVariable("MooseVariable", _result_var_name + "_rel_error", var_params);
  _rel_var_number = _aux->getFieldVariable<Real>(0, _result_var_name + "_rel_error").number();
}

// TODO: allow users to improve the sampling efficiency by providing proposal distributions.
Point
SamplePDFProblem::rejectSampleCoordPDF(THREAD_ID tid)
{
  Point final_coord(0.0, 0.0, 0.0);
  // x-coord is always required.
  unsigned int draws;
  for (draws = 0; draws < _max_num_attempts; ++draws)
  {
    const auto x_p_1 = sampleNumber(tid);
    const auto x_p_2 = sampleNumber(tid);
    const auto x_coord = _bb_min(0) + x_p_2 * _bb_extents(0);

    if (x_p_1 <= _x_pdf->value(x_coord, Point(0.0, 0.0, 0.0)) / (_x_blanket / _bb_extents(0)))
    {
      final_coord(0) = x_coord;
      break;
    }
  }
  if (draws == _max_num_attempts)
    paramError("max_attempts", "Rejection sampling failed! Please increase the maximum number of attempts for resampling!");

  // Sample the the y-coord if required.
  if (_y_pdf)
  {
    for (draws = 0; draws < _max_num_attempts; ++draws)
    {
      const auto y_p_1 = sampleNumber(tid);
      const auto y_p_2 = sampleNumber(tid);
      const auto y_coord = _bb_min(1) + y_p_2 * _bb_extents(1);

      if (y_p_1 <= _y_pdf->value(y_coord, Point(0.0, 0.0, 0.0)) / (_y_blanket / _bb_extents(1)))
      {
        final_coord(1) = y_coord;
        break;
      }
    }
  }
  if (draws == _max_num_attempts)
    paramError("max_attempts", "Rejection sampling failed! Please increase the maximum number of attempts for resampling!");

  // Sample the the z-coord if required.
  if (_z_pdf)
  {
    for (draws = 0; draws < _max_num_attempts; ++draws)
    {
      const auto z_p_1 = sampleNumber(tid);
      const auto z_p_2 = sampleNumber(tid);
      const auto z_coord = _bb_min(2) + z_p_2 * _bb_extents(2);

      if (z_p_1 <= _z_pdf->value(z_coord, Point(0.0, 0.0, 0.0)) / (_z_blanket / _bb_extents(2)))
      {
        final_coord(2) = z_coord;
        break;
      }
    }
  }
  if (draws == _max_num_attempts)
    paramError("max_attempts", "Rejection sampling failed! Please increase the maximum number of attempts for resampling!");

  return final_coord;
}

Real
SamplePDFProblem::sampleNumber(THREAD_ID tid)
{
  auto dis = std::uniform_real_distribution<Real>(0.0, 1.0);
  return dis(_rng[tid]);
}

bool
SamplePDFProblem::solverSystemConverged(const unsigned int)
{
  return true;
}
