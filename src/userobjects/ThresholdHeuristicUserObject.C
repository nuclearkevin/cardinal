#include "ThresholdHeuristicsUserObject.h"

registerMooseObject("CardinalApp",ThresholdHeuristicsUserObject);

InputParameters
ThresholdHeuristicsUserObject::validParams(){

    InputParameters params = ClusteringUserObject::validParams();
    params.addRequiredParam<double>("threshold"," threshold ");
    params.addParam<bool>("value_crosses_threshold",true," return true of the value is more than the threshold");

    return params;
}

ThresholdHeuristicsUserObject::ThresholdHeuristicsUserObject(const InputParameters & parameters):
                                ClusteringUserObject(parameters),
                                _threshold(getParam<double>("threshold")),
                                _value_crosses_threshold(getParam<bool>("value_crosses_threshold"))
{
}

bool
ThresholdHeuristicsUserObject::belongsToCluster(libMesh::Elem *elem, libMesh::Elem *neighbor_elem){

    return _value_crosses_threshold
    ? ((getMetricData(elem) > _threshold &&
        getMetricData(neighbor_elem) > _threshold))
    : ((getMetricData(elem) < _threshold &&
        getMetricData(neighbor_elem) < _threshold));

}
