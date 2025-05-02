#pragma once

#include "ClusteringUserObject.h"

class ValueRangeHeuristicUserObject: public ClusteringUserObject{

public:
    static InputParameters validParams();
    ValueRangeHeuristicUserObject(const InputParameters params);
    virtual bool belongsToCluster(libMesh::Elem* base_element,libMesh::Elem* neighbor_elem) override;
protected:
    Real _lower_limit;
    Real _upper_limit;
};