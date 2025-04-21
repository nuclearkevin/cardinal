#pragma once

#include "ClusteringUserObject.h"

class ThresholdHeuristicsUserObject: public ClusteringUserObject{

public:

    static InputParameters validParams();
    ThresholdHeuristicsUserObject(const InputParameters & parameters);

protected:
    virtual bool belongsToCluster(libMesh::Elem *base_element, libMesh::Elem *neighbor_elem) override;

private:
    double _tolerance;

};