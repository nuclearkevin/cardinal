#include "ClusteringUserObject.h"
#include "AuxiliarySystem.h"


registerMooseObject("CardinalApp", ClusteringUserObject);

InputParameters
ClusteringUserObject::validParams(){
    InputParameters params = GeneralUserObject::validParams();
    params.addRequiredParam<ExtraElementIDName>("id_name","initial value of the clustering ID");
    params.addRequiredParam<AuxVariableName>("metric_variable_name",
                                          "The name of the variable based on which clustering will be done");
    params.addClassDescription(" clustering object base that for amalgamation post processor"
                               " or maybe another user obejct. I have to talk to Prof Wilson");

    return params;
}


ClusteringUserObject::ClusteringUserObject(const InputParameters & parameters)
        :GeneralUserObject(parameters),
         _id_name(getParam<ExtraElementIDName>("id_name")),
         _mesh(_fe_problem.mesh()),
         _metric_variable_name(getParam<AuxVariableName>("metric_variable_name")),
         _metric_variable(_fe_problem.getVariable(_tid,_metric_variable_name)),
         _auxiliary_system(_fe_problem.getAuxiliarySystem()),
         _dof_map(_auxiliary_system.dofMap()),
         _metric_variable_index(_auxiliary_system.getVariable(_tid, _metric_variable_name).number())

{
    //have to move to aux kernel system
    //non linear system is only related to
    //when some sort of physics is solved
}


void ClusteringUserObject::execute(){
    std::cout<<_id_name<<"\n";

    applyNoClusteringInitialCondition();

    for (auto & elem:_mesh.element_ptr_range()){

        libMesh::Point centroid = elem->vertex_average();
        //std::cout<<elem->get_extra_integer(0)<<" ";
        //elem->set_extra_integer(0, 10);
        //std::cout<<elem->get_extra_integer(0)<<"\n";

        std::cout<<getMetricData(elem);

    }


}


Real
ClusteringUserObject::getMetricData(const libMesh::Elem *elem) {

    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<double> solution_value(1);
    _dof_map.dof_indices(elem, dof_indices, _metric_variable_index);
    _auxiliary_system.solution().get(dof_indices, solution_value);

    return static_cast<Real>(solution_value[0]);
}

void ClusteringUserObject::applyNoClusteringInitialCondition(){
    for (auto &elem:_mesh.element_ptr_range()){
        if (elem){
            elem->set_extra_integer(_extra_integer_index,-1);
        }
    }
}
