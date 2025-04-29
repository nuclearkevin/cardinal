#pragma once

#include "ClusteringUserObject.h"

class ErrorFractionHeuristicUserObject:public ClusteringUserObject{

public:
    static InputParameters validParams();
    ErrorFractionHeuristicUserObject(const InputParameters& params);
protected:
    virtual bool belongsToCluster(libMesh::Elem* base_element,libMesh::Elem*) override;
    virtual void execute() override;
    void extremesFinder();

    Real _max  =0; //maybe not a good presumption
    Real _min  =0;

    Real _upper_fraction;
    Real _lower_fraction;

    Real _upper_cut_off;
    Real _lower_cut_off;

};