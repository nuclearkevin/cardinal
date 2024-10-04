#include "ParticleOpticalDepthAux.h"

registerMooseObject("CardinalApp", ParticleOpticalDepthAux);

InputParameters
ParticleOpticalDepthAux::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("A class which computes an estimate of the particle mean free path in an element.");
  params.addRequiredCoupledVar("total_reaction_rate", "The total particle reaction rate.");
  params.addRequiredCoupledVar("scalar_flux", "The scalar particle flux.");

  return params;
}

ParticleOpticalDepthAux::ParticleOpticalDepthAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _total_rxn_rate(coupledValue("total_reaction_rate")),
    _scalar_flux(coupledValue("scalar_flux"))
{
  if (this->mooseVariable()->feType() != FEType(CONSTANT, MONOMIAL))
    paramError("variable", "At present, only CONSTANT MONOMIAL variables are supported.");
}

Real
ParticleOpticalDepthAux::computeValue()
{
  // Prevent a divide by zero when the scalar flux is small / has bad statistics.
  return _scalar_flux[_qp] < libMesh::TOLERANCE ? 0.0 : _total_rxn_rate[_qp] / _scalar_flux[_qp] * _current_elem->hmax();
}
