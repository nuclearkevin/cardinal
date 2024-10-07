#pragma once

#include "AuxKernel.h"

/**
 * A AuxKernel which computes an estimate of the particle mean free path in an element given a variable
 * containing the total reaction rate and a variable containing the scalar flux.
 */
class ParticleOpticalDepthAux : public AuxKernel
{
public:
  static InputParameters validParams();

  ParticleOpticalDepthAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// The variable containing the total reaction rate.
  const VariableValue & _rxn_rate;

  /// The variable containing the scalar flux.
  const VariableValue & _scalar_flux;
};
