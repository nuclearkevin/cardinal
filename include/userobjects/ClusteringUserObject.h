#pragma once

#include "GeneralUserObject.h"

class AuxiliarySystem;

class ClusteringUserObject: public GeneralUserObject{

public:

    static InputParameters validParams();
    ClusteringUserObject(const InputParameters & parameters);

    virtual void execute() override;
    virtual void initialize() override{};
    virtual void finalize() override {};

    void findCluster();
    Real getMetricData(const Elem * elem);
    void applyNoClusteringInitialCondition();
    virtual bool belongsToCluster(libMesh::Elem *base_element, libMesh::Elem *neighbor_elem) =0;

    const ExtraElementIDName _id_name;
    libMesh::MeshBase & _mesh;
    const AuxVariableName _metric_variable_name;
    MooseVariableFEBase & _metric_variable;
    AuxiliarySystem & _auxiliary_system ;
    libMesh::DofMap & _dof_map;

    unsigned int _metric_variable_index;
    unsigned int _extra_integer_index=0;
    int not_visited = -1;


};