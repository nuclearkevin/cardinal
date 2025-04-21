#include "ThresholdHeuristicsUserObject.h"

registerMooseObject("CardinalApp",ThresholdHeuristicsUserObject);
ThresholdHeuristicsUserObject::validParams(){

    InputParameters & params = ClusteringUserObject::validParams();
    params.addParam<double>("tolerance",0.001," Tolerance ");
    return params;
}

ThresholdHeuristicsUserObject::ThresholdHeuristicsUserObject(const InputParameters & params):
                                ClusteringUserObject(params),
                                _tolerance(getParam<double>("tolerance"))
{
}

bool ThresholdHeuristicsUserObject::belongsToCluster(libMesh::Elem *base_elem, libMesh::Elem *neighbor_elem){

    return std::abs(getMetric(base_elem) - getMetric (neighbor_elem))<=_tolerance;

}
