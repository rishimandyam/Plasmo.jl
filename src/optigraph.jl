const OptiElement = Union{OptiNode,OptiEdge}

mutable struct OptiGraph <: AbstractOptiGraph
    label::Symbol

    # topology: TODO: OrderedSets
    optinodes::OrderedSet{OptiNode}       # local optinodes
    optiedges::OrderedSet{OptiEdge}       # local optiedges
    subgraphs::OrderedSet{OptiGraph}      # local subgraphs
    optiedge_map::OrderedDict{Set{OptiNode},OptiEdge}

    # subgraphs keep a reference to their parent
    parent_graph::Union{Nothing,OptiGraph}

    # it is possible nodes and edges may use a parent graph as their model backend
    # this is the case if constructing an optigraph from subgraphs
    optimizer_graph::OptiGraph

    # track node membership in other graphs; nodes use this to query different backends
    node_to_graphs::OrderedDict{OptiNode,Vector{OptiGraph}}
    edge_to_graphs::OrderedDict{OptiEdge,Vector{OptiGraph}}

    # intermediate moi backend that maps graph elements to an MOI model
    backend::MOI.ModelLike

    node_obj_dict::OrderedDict{Tuple{OptiNode,Symbol},Any} # object dictionary for nodes
    edge_obj_dict::OrderedDict{Tuple{OptiEdge,Symbol},Any} # object dictionary for edges
    obj_dict::Dict{Symbol,Any}
    ext::Dict{Symbol,Any}

    bridge_types::Set{Any}
    is_model_dirty::Bool

    #Constructor
    function OptiGraph(;name::Symbol=Symbol(:g,gensym()))
        optigraph = new()
        optigraph.optinodes = OrderedSet{OptiNode}()
        optigraph.optiedges = OrderedSet{OptiEdge}()
        optigraph.subgraphs = OrderedSet{OptiGraph}()
        optigraph.optiedge_map = OrderedDict{Set{OptiNode},OptiEdge}()
        optigraph.parent_graph = nothing
        optigraph.optimizer_graph = optigraph

        optigraph.node_to_graphs = OrderedDict{OptiNode,Vector{OptiGraph}}()
        optigraph.node_obj_dict = OrderedDict{Tuple{OptiNode,Symbol},Any}()
        optigraph.edge_to_graphs = OrderedDict{OptiEdge,Vector{OptiGraph}}()
        optigraph.edge_obj_dict = OrderedDict{Tuple{OptiEdge,Symbol},Any}()

        optigraph.backend = GraphMOIBackend(optigraph)
        optigraph.obj_dict = Dict{Symbol,Any}()
        optigraph.ext = Dict{Symbol,Any}()
        optigraph.label = name

        optigraph.bridge_types = Set{Any}()
        optigraph.is_model_dirty = false 
        return optigraph
    end
end

function Base.string(graph::OptiGraph)
    return "OptiGraph" * " " * string(graph.label)
end
Base.print(io::IO, graph::OptiGraph) = Base.print(io, Base.string(graph))
Base.show(io::IO, graph::OptiGraph) = Base.print(io, graph)

# TODO: numerical precision like JuMP Models do
# JuMP.value_type(::Type{OptiGraph{T}}) where {T} = T

function graph_backend(graph::OptiGraph)
    return graph.optimizer_graph.backend
end

### Add subgraph

function add_subgraph(
    graph::OptiGraph; 
    optimizer_graph=nothing,
    name::Symbol=Symbol(:sg,gensym())
)
    subgraph = OptiGraph(; name=name)
    subgraph.parent_graph=graph
    if optimizer_graph != nothing
        if optimizer_graph in traverse_parents(subgraph)
            subgraph.optimizer_graph = optimizer_graph
        else
            error("Invalid optigraph passed as `optimizer_graph`")
        end
    else
        subgraph.optimizer_graph = subgraph
    end
    push!(graph.subgraphs, subgraph)
    return subgraph
end

function traverse_parents(graph::OptiGraph)
    parents = OptiGraph[]
    if graph.parent_graph != nothing
        push!(parents, graph.parent_graph)
        append!(parents, traverse_parents(graph.parent_graph))
    end
    return parents
end

### Add nodes

function add_node(
    graph::OptiGraph; 
    label=Symbol(graph.label,Symbol(".n"),length(graph.optinodes)+1)
)
    node_index = NodeIndex(length(graph.optinodes)+1)
    optinode = OptiNode{OptiGraph}(graph, node_index, label)
    push!(graph.optinodes, optinode)
    add_node(graph_backend(graph), optinode)
    return optinode
end

"""
    all_nodes(graph::OptiGraph)::Vector{OptiNode}

Recursively collect all optinodes in `graph` by traversing each of its subgraphs.
"""
function all_nodes(graph::OptiGraph)
    nodes = graph.optinodes
    for subgraph in graph.subgraphs
        nodes = [nodes; all_nodes(subgraph)]
    end
    return nodes
end

### Add edges

function get_edge(graph::OptiGraph, nodes::Set{OptiNode{OptiGraph}})
    return graph.optiedge_map[nodes]
end

function has_edge(graph::OptiGraph, nodes::Set{OptiNode{OptiGraph}})
    if haskey(graph.optiedge_map, nodes)
        return true
    else
        return false
    end
end

function add_edge(
    graph::OptiGraph,
    nodes::OptiNode...;
    label=Symbol(graph.label, Symbol(".e"), length(graph.optiedges)+1)
)
    if has_edge(graph, Set(nodes))
        edge = get_edge(graph, Set(nodes))
    else
        edge = OptiEdge{OptiGraph}(graph, label, OrderedSet(collect(nodes)))
        push!(graph.optiedges, edge)
        graph.optiedge_map[Set(collect(nodes))] = edge
        add_edge(graph_backend(graph), edge)
    end
    return edge
end

"""
    all_edges(graph::OptiGraph)::Vector{OptiNode}

Recursively collect all optiedges in `graph` by traversing each of its subgraphs.
"""
function all_edges(graph::OptiGraph)
    edges = graph.optiedges
    for subgraph in graph.subgraphs
        edges = [edges; all_edges(subgraph)]
    end
    return edges
end

### Add subgraphs

"""
    get_subgraphs(graph::OptiGraph)::Vector{OptiGraph}

Retrieve the local subgraphs of `graph`.
"""
function get_subgraphs(optigraph::OptiGraph)
    return optigraph.subgraphs
end

### MOI Extension

function MOI.get(graph::OptiGraph, attr::MOI.AnyAttribute)
    MOI.get(graph_backend(graph), attr)
end

function MOI.set(graph::OptiGraph, attr::MOI.AnyAttribute, args...)
    MOI.set(graph_backend(graph), attr, args...)
end

### JuMP Extension

function JuMP.index(graph::OptiGraph, vref::NodeVariableRef)
    gb = graph_backend(graph)
    return gb.element_to_graph_map[vref]
end

function JuMP.value(graph::OptiGraph, nvref::NodeVariableRef; result::Int = 1)
    return MOI.get(graph_backend(graph), MOI.VariablePrimal(result), nvref)
end

function JuMP.backend(graph::OptiGraph)
    return graph_backend(graph).moi_backend
end

function JuMP.object_dictionary(graph::OptiGraph)
    return graph.obj_dict
end

function JuMP.add_nonlinear_operator(
    graph::OptiGraph,
    dim::Int,
    f::Function,
    args::Vararg{Function,N};
    name::Symbol = Symbol(f),
) where {N}
    nargs = 1 + N
    if !(1 <= nargs <= 3)
        error(
            "Unable to add operator $name: invalid number of functions " *
            "provided. Got $nargs, but expected 1 (if function only), 2 (if " *
            "function and gradient), or 3 (if function, gradient, and " *
            "hesssian provided)",
        )
    end
    MOI.set(graph, MOI.UserDefinedFunction(name, dim), tuple(f, args...))
    return JuMP.NonlinearOperator(f, name)
end

### Link Constraints

function JuMP.add_constraint(
    graph::OptiGraph, con::JuMP.AbstractConstraint, name::String=""
)
    nodes = _collect_nodes(JuMP.jump_function(con))
    @assert length(nodes) > 0
    length(nodes) > 1 || error("Cannot create a linking constraint on a single node")
    edge = add_edge(graph, nodes...)
    con = JuMP.model_convert(edge, con)
    cref = _moi_add_edge_constraint(edge, con)
    return cref
end

function _collect_nodes(
    jump_func::Union{
        NodeVariableRef, 
        JuMP.GenericAffExpr, 
        JuMP.GenericQuadExpr,
        JuMP.GenericNonlinearExpr
    }
)
    vars = _extract_variables(jump_func)
    nodes = JuMP.owner_model.(vars)
    return collect(nodes)
end

### Objective Function

function JuMP.objective_function(
    graph::OptiGraph,
    ::Type{F},
) where {F<:MOI.AbstractFunction}
    func = MOI.get(JuMP.backend(graph), MOI.ObjectiveFunction{F}())::F
    return JuMP.jump_function(graph, func)
end

function JuMP.objective_function(graph::OptiGraph, ::Type{T}) where {T}
    return JuMP.objective_function(graph, JuMP.moi_function_type(T))
end

function JuMP.objective_function(graph::OptiGraph)
    F = MOI.get(JuMP.backend(graph), MOI.ObjectiveFunctionType())
    return JuMP.objective_function(graph, F)
end

function JuMP.set_objective(
    graph::OptiGraph, sense::MOI.OptimizationSense, func::JuMP.AbstractJuMPScalar
)
    JuMP.set_objective_sense(graph, sense)
    JuMP.set_objective_function(graph, func)
    return
end

function JuMP.set_objective_sense(graph::OptiGraph, sense::MOI.OptimizationSense)
    MOI.set(graph_backend(graph), MOI.ObjectiveSense(), sense)
    return
end

function JuMP.set_objective_function(
    graph::OptiGraph, 
    expr::JuMP.GenericAffExpr{C,NodeVariableRef}
) where C <: Real
    _moi_set_objective_function(graph, expr)
    return
end

function JuMP.set_objective_function(
    graph::OptiGraph, 
    expr::JuMP.GenericQuadExpr{C,NodeVariableRef}
) where C <: Real
    _moi_set_objective_function(graph, expr)
    return
end

function JuMP.set_objective_function(
    graph::OptiGraph, 
    expr::JuMP.GenericNonlinearExpr{NodeVariableRef}
)
    _moi_set_objective_function(graph, expr)
    return
end

function _moi_set_objective_function(
    graph::OptiGraph, 
    expr::JuMP.GenericAffExpr{C,NodeVariableRef}
) where C <: Real
    moi_func = JuMP.moi_function(expr)
    
    # add variables to backend if using subgraphs
    _add_backend_variables(graph_backend(graph), expr)

    # update the moi function variable indices
    graph_moi_func = _create_graph_moi_func(graph_backend(graph), moi_func, expr)

    MOI.set(
        graph_backend(graph),
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{C}}(),
        graph_moi_func,
    )
    return
end

function _moi_set_objective_function(
    graph::OptiGraph, 
    expr::JuMP.GenericQuadExpr{C,NodeVariableRef}
) where C <: Real
    moi_func = JuMP.moi_function(expr)
    
    # add variables to backend if using subgraphs
    _add_backend_variables(graph_backend(graph), expr)

    # update the moi function variable indices
    graph_moi_func = _create_graph_moi_func(graph_backend(graph), moi_func, expr)

    MOI.set(
        graph_backend(graph),
        MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{C}}(),
        graph_moi_func,
    )
    return
end

function _moi_set_objective_function(
    graph::OptiGraph, 
    expr::JuMP.GenericNonlinearExpr{NodeVariableRef}
)
    moi_func = JuMP.moi_function(expr)
    
    # add variables to backend if using subgraphs
    _add_backend_variables(graph_backend(graph), expr)

    # update the moi function variable indices
    graph_moi_func = _create_graph_moi_func(graph_backend(graph), moi_func, expr)

    MOI.set(
        graph_backend(graph),
        MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
        graph_moi_func,
    )
    return
end

"""
    JuMP.objective_value(graph::OptiGraph)

Retrieve the current objective value on optigraph `graph`.
"""
function JuMP.objective_value(graph::OptiGraph)
    return MOI.get(backend(graph), MOI.ObjectiveValue())
end

