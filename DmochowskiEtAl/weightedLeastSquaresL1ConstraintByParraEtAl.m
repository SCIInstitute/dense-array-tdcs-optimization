function [I, fVal, dV] = weightedLeastSquaresL1ConstraintByParraEtAl(sqrtQ_weighted, v_weighted, ind)
% Weighted least squares solution to hd_tDCS electrode current stimulus 
% optimization with L1 constraint.
%
% Synopsis: [I, fval, dV] = weightedLeastSquaresL1ConstraintByParraEtAl(...
%           sqrtQ_weighted, v_weighted, ind)
%
% Input:    sqrtQ_weighted  = chol factor of weighted quadratic matrix.
%           v_weighted      = weighted right hand side vector.
%           ind             = individual electrode current bound.
%
% Output:   I       =   array of electrode currents.
%           fVal    =   least squares error for best solution. 
%           dV      =   dual variables for the constraints.

% Notes:    1. Use the equation in section 3.2. of " Optimized multi-electrode 
%           stimulation increases focality and intensity at the target.",
%           Jacek P Dmochowski, et al., Journal of neural engineering 
%           8.4 (2011): 046011.
%
%           2. We convert quadratic objective function to norm 2 in order
%           to speed up the convergence. This is the reason behind using
%           Cholesky factor of the quadratic matrix instead of the matrix
%           itself. minimize(x' * Q * x - 2*v*x) is equivalent to
%           minimize( norm(chol(Q) * (x - m), 2 ) ) where m = Q\v.

tic;
L = size(sqrtQ_weighted,2); %number of electrodes
Q = sqrtQ_weighted' * sqrtQ_weighted;
m = Q\v_weighted;

if size(ind,1) == L
    ind(L+1,:) = ind(end,:);
end

cvx_begin quiet
cvx_solver sedumi
cvx_precision high
variable x(L)
expression y(L+1)
y = [x; -sum(x)];
dual variable indConstLB
dual variable indConstUB
minimize norm(sqrtQ_weighted * (x-m));
subject to
indConstLB : -ind <= y;
indConstUB : y <= ind;
cvx_end

if ~strcmp(cvx_status,'Solved')
    fprintf('%s\n','No solution with high precision, trying low precision.');
    cvx_begin quiet
    cvx_solver sedumi
    cvx_precision low
    variable x(L)
    expression y(L+1)
    y = [x; -sum(x)];
    dual variable indConstLB
    dual variable indConstUB
    minimize norm(sqrtQ_weighted * (x-m));
    subject to
    indConstLB : -ind <= y;
    indConstUB : y <= ind;
    cvx_end
end

I = x;
fVal = cvx_optval;
dV = [indConstLB; indConstUB];
fprintf('%s%f%s\n', 'Weighted least squares solution is found in ', toc, ...
    ' seconds.');
end
