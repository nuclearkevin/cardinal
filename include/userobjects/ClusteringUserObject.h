#pragma once

#include "GeneralUserObject.h"
#include "NonlinearSystemBase.h"

class ClusteringUserObject: public GeneralUserObject{

public:

    static InputParameters validParams();
    ClusteringUserObject(const InputParameters & parameters);

    virtual void execute() override;
    virtual void initialize() override{};
    virtual void finalize() override {};

    //void findCluster();
    Real getMetricData(const Elem * elem);
    //virtual bool belongsToCluster(libMesh::Elem *base_element, libMesh::Elem *neighbor_elem);

    const ExtraElementIDName _id_name;
    const unsigned int _extra_integer_index=0;
    libMesh::MeshBase& _mesh;
    const VariableName _metric_variable_name;
    MooseVariableFEBase & _metric_variable;
    NonlinearSystemBase & _nl ;
    libMesh::DofMap & _dof_map;

    const unsigned int _metric_variable_index;


};
