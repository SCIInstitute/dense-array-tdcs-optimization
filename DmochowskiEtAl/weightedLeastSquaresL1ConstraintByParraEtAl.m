function [I, fVal, dV] = weightedLeastSquaresL1ConstraintByParraEtAl(Q_weighted, v_weighted, ind)
% Weighted least squares solution to dense array tDCS stimulus pattern
% optimization with L1 (individual electrode) constraint.
%
% Synopsis: [I, fval, dV] = weightedLeastSquaresL1ConstraintByParraEtAl(...
%           Q_weighted, v_weighted, ind)
%
% Input:    Q_weighted      = Weighted quadratic matrix.
%           v_weighted      = weighted right hand side vector.
%           ind             = individual electrode current bound.
%
% Output:   I       =   array of electrode currents.
%           fVal    =   least squares error for best solution.
%           dV      =   dual variables for the constraints.

% Notes:    1. The implementation is based on the equation (11) in section 3.2 
%           of the article: Dmochowski et al., "Optimized multi-electrode 
%           stimulation increases focality and intensity at the target," 
%           Journal of neural engineering, 8.4 : 046011, 2011.
%
%           2. We convert quadratic objective function to norm 2 in order
%           to speed up the convergence. This is the reason behind using
%           Cholesky factor of the quadratic matrix instead of the matrix
%           itself. minimize(x' * Q * x - 2*v*x) is equivalent to
%           minimize( norm(chol(Q) * (x - m), 2 ) ) where m = Q\v.
%
%           3. To calculate the matrix Q_weighted and the vector
%           v_weighted, the functions calcQuadraticLinearTerms() and 
%           findWeightedQuadratic() can be used.  

tic;

%Cholesky factor
[sqrtQ_weighted,p] = chol(Q_weighted);
if p>0
    warning('quadratic matrix is not positive definite.');
    minEig = min(eig(Q_weighted));
    sqrtQ_weighted = chol(Q_weighted + (-minEig+eps)*eye(size(Q_weighted)));
end

L = size(sqrtQ_weighted,2); %number of electrodes (reference excluded)
m = Q_weighted\v_weighted;

if size(ind,1) == L
    ind(L+1,:) = ind(end,:); %reference electrode bound
end

cvx_begin quiet
cvx_solver sedumi
cvx_precision high
variable x(L) % electrode current array (ref excluuded)
expression y(L+1) % electrode current array (ref included)
y = [x; -sum(x)]; % reference electrode current = - sum(the rest)
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
dV.indLB = indConstLB;
dV.indUB = indConstUB;
fprintf('%s%f%s\n', 'Weighted least squares solution is found in ', toc, ...
    ' seconds.');
end
