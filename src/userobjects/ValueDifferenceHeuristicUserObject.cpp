#include "ValueDifferenceHeuristicUserObject.h"

InputParameters
ValueDifferenceHeuristicUserObject::validParams(){

    InputParameters params == ClusteringUserObject::validParams();
    params.addRequiredParam<Real>("ValueDifferenceHeuristicUserObject","Lower limit of the range");
    params.addClassDescription("equal neighbor heuristic");

    return params;
}

ValueDifferenceHeuristicUserObject::ValueDifferenceHeuristicUserObject(const InputParameters & parameters):
        ClusteringUserObject(parameters),
        _tolerance(getParam<Real>("tolerance"))
{}

bool
ValueDifferenceHeuristicUserObject::belongsToCluster(libMesh::Elem * base_element, libMesh::Elem* neighbor_element){
    return std::abs(getMetricData(base_element) - getMetricData(neighbor_element))<_tolerance;
}