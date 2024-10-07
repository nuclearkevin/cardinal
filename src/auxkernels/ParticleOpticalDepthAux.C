#include "ParticleOpticalDepthAux.h"

registerMooseObject("CardinalApp", ParticleOpticalDepthAux);

InputParameters
ParticleOpticalDepthAux::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("A class which computes an estimate of the optical depth 'seen' by a particle within an element.");
  params.addRequiredCoupledVar("reaction_rate", "The particle reaction rate to use.");
  params.addRequiredCoupledVar("scalar_flux", "The scalar particle flux.");

  return params;
}

ParticleOpticalDepthAux::ParticleOpticalDepthAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _rxn_rate(coupledValue("reaction_rate")),
    _scalar_flux(coupledValue("scalar_flux"))
{
  if (this->mooseVariable()->feType() != FEType(CONSTANT, MONOMIAL))
    paramError("variable", "At present, only CONSTANT MONOMIAL variables are supported.");
}

Real
ParticleOpticalDepthAux::computeValue()
{
  // Prevent a divide by zero when the scalar flux is small / has bad statistics.
  return _scalar_flux[_qp] < libMesh::TOLERANCE ? 0.0 : _rxn_rate[_qp] / _scalar_flux[_qp] * _current_elem->hmax();
}
