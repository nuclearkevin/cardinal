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

#include "OpticalDepthMarker.h"

registerMooseObject("CardinalApp", OpticalDepthMarker);

InputParameters
OpticalDepthMarker::validParams()
{
  auto params = Marker::validParams();
  params.addClassDescription(
    "A class which marks elements for refinement/coarsening based on the optical depth "
    "'seen' by a neutron/photon within an element. If the optical depth is greater than "
    "'refine' the element is refined. If the optical depth is less than 'coarsen', the "
    "element is coarsened.");
  params.addRequiredCoupledVar("reaction_rate", "The particle reaction rate to use.");
  params.addRequiredCoupledVar("scalar_flux", "The scalar particle flux.");
  params.addRequiredParam<Real>("refine", "The refinement threshold.");
  params.addRequiredParam<Real>("coarsen", "The coarsening threshold.");

  return params;
}

OpticalDepthMarker::OpticalDepthMarker(const InputParameters & parameters)
  : Marker(parameters),
    Coupleable(this, false),
    _refine(getParam<Real>("refine")),
    _coarsen(getParam<Real>("coarsen")),
    _rxn_rate(coupledValue("reaction_rate")),
    _scalar_flux(coupledValue("scalar_flux"))
{
  if (getCoupledVars().at("reaction_rate")[0]->feType() != FEType(CONSTANT, MONOMIAL))
    paramError("reaction_rate", "At present, only CONSTANT MONOMIAL variables are supported.");
  if (getCoupledVars().at("scalar_flux")[0]->feType() != FEType(CONSTANT, MONOMIAL))
    paramError("scalar_flux", "At present, only CONSTANT MONOMIAL variables are supported.");
}

Marker::MarkerValue
OpticalDepthMarker::computeElementMarker()
{
  auto od = _scalar_flux[0] < libMesh::TOLERANCE ? 0.0 : _rxn_rate[0] / _scalar_flux[0] * _current_elem->hmax();

  if (od > _refine)
    return Marker::MarkerValue::REFINE;
  else if (od < _coarsen)
    return Marker::MarkerValue::COARSEN;
  else
    return Marker::MarkerValue::DO_NOTHING;
}
