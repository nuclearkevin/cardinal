#include "ThresholdHeuristicsUserObject.h"

registerMooseObject("CardinalApp", ThresholdHeuristicsUserObject);

InputParameters
ThresholdHeuristicsUserObject::validParams()
{

  InputParameters params = ClusteringUserObject::validParams();
  params.addRequiredParam<double>("threshold", " The value against which the clustering process is compared.");
  params.addParam<bool>(
      "value_crosses_threshold", true," Return true if the value is more than the threshold");
  params.addClassDescription("A special type of ClusterUserObject that applies threshold"
                        "based heuristics on the element pairs.");

  return params;
}

ThresholdHeuristicsUserObject::ThresholdHeuristicsUserObject(const InputParameters & parameters)
  : ClusteringUserObject(parameters),
    _threshold(getParam<double>("threshold")),
    _value_crosses_threshold(getParam<bool>("value_crosses_threshold"))
{}

bool
ThresholdHeuristicsUserObject::belongsToCluster(libMesh::Elem * elem, libMesh::Elem * neighbor_elem)
{

  return _value_crosses_threshold
             ? ((getMetricData(elem) > _threshold && getMetricData(neighbor_elem) > _threshold))
             : ((getMetricData(elem) < _threshold && getMetricData(neighbor_elem) < _threshold));
}


