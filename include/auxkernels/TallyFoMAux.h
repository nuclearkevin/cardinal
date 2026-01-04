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

#include "OpenMCAuxKernel.h"

/**
 * A class which computes a Figure of Merit (FoM) for a tally quantity.
 * This can be one of several candidate FoMs:
 * - The variance reduction FoM:   FoM = 1 / (T * sigma^2)
 * - The relative discrepancy FoM: FoM = ref / (T * sigma^2 * |var - ref|)
 * - The absolute discrepancy FoM: FoM = 1 / (T * sigma^2 * |var - ref|)
 * Optionally, the FoM may be divided by the element volume.
 */
class TallyFoMAux : public OpenMCAuxKernel
{
public:
  static InputParameters validParams();

  TallyFoMAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// The variable containing the tally score.
  const VariableValue * _tally_val;

  /// The variables containing the tally score standard deviation.
  const VariableValue & _tally_std_dev;

  /// The variable containing the "reference" solution.
  const VariableValue * _reference_val;

  /// The simulation time.
  const Real & _sim_time;

  /// The type of FoM to compute.
  const enum class FoMType
  {
    VarRed = 0,
    RelDis = 1,
    AbsDis = 2
  } _fom_type;

  /// Whether the FoM should be divided by volume or not.
  const bool _divide_by_volume;
};
