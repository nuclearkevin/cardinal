#include "ValueRangeHeuristicUserObject.h"

InputParameters
ValueRangeHeuristicUserObject::validParams(){

    InputParameters params == ClusteringUserObject::validParams();
    params.addRequiredParam<Real>("lower_limit","Lower limit of the range");
    params.addRequiredParam<Real>("upper_limit","Upper limit of the range");
    params.addClassDescription("Mimics the ValueRangeMarker in moose");

    return params;
}

ValueRangeHeuristicUserObject::ValueRangeHeuristicUserObject(const InputParameters & parameters):
                                ClusteringUserObject(parameters),
                                _lower_limit(getParam<Real>("lower_limit")),
                                _upper_limit(getParam<Real>("upper_limit"))
{}

ValueRangeHeuristicUserObject::belongsToCluster(libMesh::Elem * base_element, libMesh::Elem* neighbor_element){

    return _lower_limit<getMetricData(base_element)<_upper_limit && _lower_limit<getMetricData(neighbor_element)<_upper_limit;
}