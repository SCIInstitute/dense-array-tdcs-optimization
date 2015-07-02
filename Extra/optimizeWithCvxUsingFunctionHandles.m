function result = optimizeWithCvxUsingFunctionHandles(objectiveFunction, ...
    constraintFunctionArray, constraintBoundArray, variableSize)

nConstraint = numel(constraintFunctionArray);

assert(isa(objectiveFunction,'function_handle'),['Objective function'...
    'is not type function_handle.']);

for iConstraint = 1:nConstraint
    assert(isa(constraintFunctionArray{iConstraint},'function_handle'),...
        [num2str(iConstraint) 'th constraint function is not type function_handle.']);
end

assert(nConstraint==numel(constraintBoundArray),['Number of constraints'...
    'dont match the number of constraint bounds']);


cvx_begin
variable x(variableSize)
minimize objectiveFunction(x)
subject to
for iConstraint = 1:nConstraint
    constraintFunctionArray{iConstraint}(x) <= constraintBoundArray(iConstraint);
end
cvx_end

result = x;
