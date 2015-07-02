function [electrodeCurrents,fval,dualVariables] = optimizationUsingCvxToolbox(w,Q,tot,ind,pmax)
%% FINDS THE BEST ELECTRODE CURRENT STIMULUS CONFIGURATION USING CVX TOOLKIT
% IN ADDITION IT FINDS THE FUNCTION VALUE AT SOLUTION AND DUAL VARIABLES (FOR ACTIVE SET 
% INVESTIGATIN PURPOSES) IF REQUIRED.
%
%Written by: Seyhmus Guler, 4/29/14
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
%w: Linear weights for the objective function.
%sqrtQ: Square roots of matrices for quadratic constraint.
%tot: The total current bound
%ind: Individual electrode current bound.
%pmax: The bound on the electrical power in avoidance regions
%
%OUTPUTS:
%electrodeCurrents: best solution for electrode current array.
%fval: Objective function, the current in the ROI along the desired
    %direction.
%DualVariables: dual variables for the constraints at solution. The
    % variables corresponding to first total bound, second individual 
    % lower bound, third individual upper bound and lastly power bound 
    % is stored in an array of size #constraints x 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
L = numel(w); %Number of electrodes
pp = numel(Q); % Number of power constraints

sqrtQ = cell(1,pp);
for i=1:pp
    %Cholesky factor for reasons explained above
    [sqrtQ{i},p] = chol(Q{i});
    if p>0
        warning('quadratic matrix is not positive definite.');
        minEig = min(eig(Q{i}));
        sqrtQ{i} = chol(Q{i}+(-minEig+eps)*eye(size(Q{i})));
    end
end
if size(ind,1) == L %In case reference electrode bound is not defined
    ind(L+1,:) = ind(end,:); %Make ref elec have the same const as last one
end

if size(ind,2) == 1 %Lower bound = - Upper bound
    ind = [-ind ind];
end

cvx_begin quiet
cvx_precision high
cvx_solver sedumi

%optimization variable
variable x(L)

%expression of electrode current array, reference added
expression y(L+1)
y = [x; -sum(x)]; %reference elec current = - (sum of other electrodes)

%dual variables to check active/passive constraints
dual variable totConstV
dual variable indConstLB
dual variable indConstUB
dual variable powConstV{pp}

%maximize current density in the ROI along the desired direction
maximize w*x

%subject to total and ind electrode currents, and the power constraints
subject to
totConstV : norm(y,1) <= 2*tot; %#ok<*VUNUS> %total current bound
indConstLB : ind(:,1) <= y; %lower bound on individual currents
indConstUB : y <= ind(:,2); %#ok<*BDSCI> %upper bound on individual currents
%power in the avoid regions is bounded:
for i = 1:pp
    powConstV{i} : norm(sqrtQ{i}*x) <= sqrt(pmax(i)); 
end
cvx_end

if ~strcmp(cvx_status,'Solved')
    fprintf('%s\n','No solution with high precision, trying lower precision settings');
    cvx_begin quiet
    cvx_precision medium
    cvx_solver sedumi
    
    variable x(L)
   
    expression y(L+1)
    y = [x; -sum(x)];
    
    dual variable totConstV
    dual variable indConstLB
    dual variable indConstUB
    dual variable powConstV{pp}
    
    %maximize current density in the ROI along the desired direction
    maximize w*x
    
    %subject to total curr., ind elec currents and the power constraints
    subject to
    totConstV : norm(y,1) <= 2*tot; %total current bound
    indConstLB : ind(:,1) <= y; %lower bound on individual currents
    indConstUB : y <= ind(:,2); %upper bound on individual currents
    for i = 1:pp
        powConstV{i} : norm(sqrtQ{i} * x) <= sqrt(pmax(i)); %power in the avoid regions is bounded
    end
    cvx_end
end

if ~strcmp(cvx_status,'Solved')
    fprintf('%s\n','Still no solution, trying SDPT3 as solver');
    
    cvx_begin quiet
    cvx_solver SDPT3
    variable x(L)
   
    expression y(L+1)
    y = [x; -sum(x)];
    
    dual variable totConstV
    dual variable indConstLB
    dual variable indConstUB
    dual variable powConstV{pp}
    
    %maximize current density in the ROI along the desired direction
    maximize w*x
    
    %subject to total curr., ind elec currents and the power constraints
    subject to
    totConstV : norm(y,1) <= 2*tot; %total current bound
    indConstLB : ind(:,1) <= y; %lower bound on individual currents
    indConstUB : y <= ind(:,2); %upper bound on individual currents
    for i = 1:pp
        powConstV{i} : norm(sqrtQ{i}*x) <= sqrt(pmax(i)); %power in the avoid regions is bounded
    end
    cvx_end
end

electrodeCurrents = x;
fval = cvx_optval;
dualVariables.tot = totConstV;
dualVariables.indLB = indConstLB;
dualVariables.indUB = indConstUB; 
dualVariables.pow = cell2mat(powConstV); 

fprintf('%s%f%s\n','Electrode current optimization is solved in ',toc,' seconds.');
end
