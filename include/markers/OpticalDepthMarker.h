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

#include "Marker.h"
#include "Coupleable.h"
#include "MaterialPropertyInterface.h"

/**
 * A Marker which marks elements for refinement/coarsening based on an estimate of the optical depth of a
 * neutron/photon which traverses the element.
 */
class OpticalDepthMarker : public Marker,
                           public Coupleable
{
public:
  static InputParameters validParams();

  OpticalDepthMarker(const InputParameters & parameters);

protected:
  virtual MarkerValue computeElementMarker() override;

  /// The refinement threshold.
  const Real & _refine;

  /// The coarsening threshold.
  const Real & _coarsen;

  /// The variable containing the total reaction rate.
  const VariableValue & _rxn_rate;

  /// The variable containing the scalar flux.
  const VariableValue & _scalar_flux;
};
