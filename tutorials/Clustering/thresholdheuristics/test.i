[Mesh]
    [mesh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 20
        ny = 20
        extra_element_integers = 'cluster_id'
    []
[]

[AuxVariables]
    [metric_var]
        order = FIRST
        family = LAGRANGE
    []
    [cluster_id_aux]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[AuxKernels]
    [create_metric]
        type = FunctionAux
        variable = metric_var
        function = 'sqrt(x*x + y*y)'
        execute_on = 'TIMESTEP_BEGIN'
    []
    [store_element_id]
        type=ExtraElementIDAux
        extra_id_name ="cluster_id"
        variable=cluster_id_aux
    []
[]

[UserObjects]
    [clustering]
        type = ThresholdHeuristicsUserObject
        execute_on = 'TIMESTEP_END'

        id_name = 'cluster_id'
        metric_variable_name = 'metric_var'
        threshold = 1.1
    []
[]

[Problem]
    type = FEProblem
    solve = false
[]

[Executioner]
    type = Transient
    solve = false
    dt = 0.1
    num_steps = 2
[]

[Outputs]
    exodus = true
[]
