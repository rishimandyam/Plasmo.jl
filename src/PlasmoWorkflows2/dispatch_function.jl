mutable struct DispatchFunction
    func::Function                   #the function to call
    args::Vector{Any}                #the function args
    kwargs::Dict{Any,Any}
    result::Nullable{Any}            #the result after calling the event
end
DispatchFunction() = DispatchFunction(() -> nothing,[],Dict(),nothing)
DispatchFunction(func::Function) =  DispatchFunction(func,[],Dict(),nothing)
run!(dfunc::DispatchFunction) = dfunc.result = dfunc.func(dfunc.args...,dfunc.kwargs...)
getresult(dfunc::DispatchFunction) = get(dfunc.result)
