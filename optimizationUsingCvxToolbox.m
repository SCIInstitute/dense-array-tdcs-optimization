function [I,fval,actCons] = optimizationUsingCvxToolbox(w,Q,tot,ind,pmax,actTh)
%FINDS THE BEST ELECTRODE CURRENT STIMULUS CONFIGURATION USING CVX
%IN ADDITION IT FINDS THE FUNCTION VALUE AT SOLUTION AND THE SET OF ACTIVE
%CONSTRAINTS IF REQUIRED.
%
%Written by: Seyhmus Guler, Revisited: Moritz Dannhauer
%Last edit: 3/9/14 by Guler,S
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
    %w: Linear weights for the objective function.
    %Q: Quadratic constraint matrices for each avoid region
    %tot: the total current bound
    %ind: individual electrode current bounds
    %pmax: the bound on the electrical power in avoidance regions
    %actTh: A constraint is considered active if the value is at most 
    % actTh far from the constraint boundary. 
%OUTPUTS:
    %I: best solution for electrode current array
    %fval: objective function, the current in the ROI along the dd
    %actCons: a logical array showing the active constraints at solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
L = numel(w); %Number of electrodes

if size(ind,1) == L %In case reference electrode bound is not defined
    ind(L+1,:) = ind(end,:);
end

%Finding best electrode configuration using CVX toolkit
cvx_begin quiet
cvx_precision high
cvx_solver sedumi
variable x(L)
expression y(L+1) 
y = [x; -sum(x)]; %expression with reference electrode added.
maximize w*x %maximize CD in the ROI along the desired direction
subject to
-ind <= y <= ind; %Individual electrode current bound
norm(y,1) <= 2*tot; %total current bound
for i = 1:numel(Q)
quad_form(x,Q{i}) <= pmax(i); %power in the avoid regions is bounded
end
cvx_end
I = x;
fval = cvx_optval;

if (nargout >=3 && nargin >= 6)
    threshold = actTh;
    actCons = false(1+2*(L+1)+numel(Q),1); %total # constraints
    actCons(1) = norm([x;-sum(x)],1)-2*tot <= threshold;
    actCons(2:L+2) = abs([x;-sum(x)]+ind) <= threshold;
    actCons(L+3:2*L+3) = abs([x;-sum(x)]-ind) <= threshold;
    for i=1:numel(Q)
        actCons(2*L+3+i) = abs(x' * Q{i} * x - pmax(i)) <= threshold;
    end
end
fprintf('%s%f%s\n','Electrode current optimization is solved in ',toc,' seconds.');
end
