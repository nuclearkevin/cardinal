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

#include "TallyFoMAux.h"

registerMooseObject("CardinalApp", TallyFoMAux);

InputParameters
TallyFoMAux::validParams()
{
  auto params = OpenMCAuxKernel::validParams();
  params.addClassDescription(
      "An auxkernel which a Figure of Merit (FoM) for a tally variable.");
  params.addCoupledVar(
      "tally_value",
      "The variable containing the value of the tally. TallyFoMAux "
      "assumes this is a volumetric quantity. This is required for the discrepancy "
      "figures of merit.");
  params.addRequiredCoupledVar(
      "tally_std_dev",
      "The variable containing the standard deviation of the tally. TallyFoMAux "
      "assumes this is a volumetric quantity.");
  params.addCoupledVar(
      "ref_value",
      "The variable containing the reference solution. TallyFoMAux "
      "assumes this is a volumetric quantity. This is required for the discrepancy "
      "figures of merit.");

  params.addRequiredParam<PostprocessorName>("sim_time", "The OpenMC simulation time.");

  params.addRequiredParam<MooseEnum>("fom_type",
                                     MooseEnum("var_red rel_dis abs_dis"),
                                     "The type of Figure of Merit (FoM) to compute. Options are the standard "
                                     "variance reduction FoM (var_red), the relative discrepancy FoM "
                                     "(rel_dis), or the absolute discrepancy FoM (abs_dis).");

  params.addParam<bool>(
      "divide_by_volume",
      false,
      "Whether the FoM should be divided by the element volume or not.");

  return params;
}

TallyFoMAux::TallyFoMAux(const InputParameters & parameters)
  : OpenMCAuxKernel(parameters),
    _tally_val(isCoupled("tally_value") ? &coupledValue("tally_value") : nullptr),
    _tally_std_dev(coupledValue("tally_std_dev")),
    _reference_val(isCoupled("ref_value") ? &coupledValue("ref_value") : nullptr),
    _sim_time(getPostprocessorValue("sim_time")),
    _fom_type(getParam<MooseEnum>("fom_type").getEnum<FoMType>()),
    _divide_by_volume(getParam<bool>("divide_by_volume"))
{
  if (_var.feType() != FEType(libMesh::CONSTANT, libMesh::MONOMIAL))
    paramError("variable",
               "FDTallyGradAux only supports CONSTANT MONOMIAL shape functions. Please "
               "ensure that 'variable' is of type MONOMIAL and order CONSTANT.");

  if (_tally_val)
    if (getFieldVar("tally_value", 0)->feType() != FEType(libMesh::CONSTANT, libMesh::MONOMIAL))
      paramError("tally_value",
        "TallyFoMAux only supports CONSTANT MONOMIAL shape functions. Please "
        "ensure that 'tally_value' is of type MONOMIAL and order CONSTANT.");

  if (getFieldVar("tally_std_dev", 0)->feType() != FEType(libMesh::CONSTANT, libMesh::MONOMIAL))
    paramError("tally_std_dev",
      "TallyFoMAux only supports CONSTANT MONOMIAL shape functions. Please "
      "ensure that 'tally_std_dev' is of type MONOMIAL and order CONSTANT.");

  if (_reference_val)
    if (getFieldVar("ref_value", 0)->feType() != FEType(libMesh::CONSTANT, libMesh::MONOMIAL))
      paramError("ref_value",
        "TallyFoMAux only supports CONSTANT MONOMIAL shape functions. Please "
        "ensure that 'ref_value' is of type MONOMIAL and order CONSTANT.");

  // Error-check the different FoM options.
  switch (_fom_type)
  {
    case FoMType::VarRed:
      break;
    case FoMType::RelDis:
    case FoMType::AbsDis:
    {
      if (!_tally_val)
        paramError("tally_value",
          "A tally variable must be provided when computing using the following figure of merit: ",
          getParam<MooseEnum>("fom_type"));

      if (!_reference_val)
        paramError("ref_value",
          "A reference solution must be provided when computing using the following figure of merit: ",
          getParam<MooseEnum>("fom_type"));
      break;
    }
    default:
    {
      mooseError("Unhandled FoMType enum in TallyFoMAux!");
      break;
    }
  }
}

Real
TallyFoMAux::computeValue()
{
  // Undo the divide-by-volume (since we're assuming the variables provided are volumetric).
  const auto std_2 = std::pow(_tally_std_dev[0] * _current_elem->volume(), 2.0);

  // Every FoM uses the variance reduction FoM as a base.
  auto fom = 1.0 / (_sim_time * std_2);
  switch (_fom_type)
  {
    case FoMType::RelDis:
    {
      fom /= (std::abs((*_tally_val)[0] - (*_reference_val)[0]) / (*_reference_val)[0]);
      break;
    }
    case FoMType::AbsDis:
    {
      // Undo the divide-by-volume (since we're assuming the variables provided are volumetric).
      fom /= (_current_elem->volume() * std::abs((*_tally_val)[0] - (*_reference_val)[0]));
      break;
    }
    default:
      break;
  }

  return _divide_by_volume ? fom / _current_elem->volume() : fom;
}

#endif
