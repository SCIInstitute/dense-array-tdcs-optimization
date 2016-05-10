function [I, fVal, dV] = weightedLeastSquaresConstrainedByParraEtAl(Q_weighted, v_weighted, tot)
% Weighted least squares solution to dense array tDCS stimulus pattern
% optimization with total injected current constraint
%
% Synopsis: [I, fVal, dV] = weightedLeastSquaresConstrainedByParraEtAl( ...
%           Q_weighted, v_weighted, tot)
%
% Input:    Q_weighted      =   weighted quadratic matrix.
%           v_weighted      =   weighted linear term.
%           tot             =   total current bound.
%
% Output:   I       =   array of electrode currents.
%           fVal    =   least squares error for best solution. 
%           dV      =   dual variables for the constraints.

% Notes:    1. The implementation is based on the equation (10) in section 3.1 
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
%           3. The conversion from original paper to the formula in the
%           script:
%           minimize_over_I ( || sqrt(W) (Ed - A*I) || ) ^2
%           = (Ed - A*I)' * sqrt(W)' * sqrt(W)* (Ed -A*I)
%           = (Ed - A*I)' * W * (Ed - A*I)
%           = Ed'*W*Ed + I'*A'*W*A*I - I'*A'*W*Ed - Ed'*W*A*I 
%           = Ed'*W*Ed + I'*(A'*W*A)*I - 2*(Ed'*W*A)*I
%           = constant + I'*(  Q   )*I - 2*(   v   )*I  where
%           Q := A' * W * A, v := Ed' * W * A.
%
%           4. Minimizing above expression is equivalent to minimizing
%           (according to note 2 and assuming sqrtQ'*sqrtQ = Q) :
%           norm ( sqrtQ * (x - inv(Q) * v)).
%
%
%           5. To calculate the matrix Q_weighted and the vector
%           v_weighted, the functions calcQuadraticLinearTerms() and 
%           findWeightedQuadratic() can be used.  

tic;
[sqrtQ_weighted,p] = chol(Q_weighted);
if p>0
    sqrtQ_weighted = chol(Q_weighted + 1e-9*eye(size(Q_weighted,1)));
end
L = size(sqrtQ_weighted,2); %number of electrodes
m = Q_weighted\v_weighted;

cvx_begin quiet
cvx_solver sedumi
cvx_precision high
variable x(L)
expression y(L+1)
y = [x; -sum(x)];
dual variable totConst
minimize norm(sqrtQ_weighted*(x-m));
subject to
totConst : norm(y,1) <= 2*tot;
cvx_end
if ~strcmp(cvx_status,'Solved')
    fprintf('%s\n','No solution with high precision, trying low precision.');
    cvx_begin quiet
    cvx_solver sedumi
    cvx_precision low
    variable x(L)
    expression y(L+1)
    y = [x; -sum(x)];
    dual variable totConst
    minimize norm(sqrtQ_weighted*(x-m));
    subject to
    totConst : norm(y,1) <= 2*tot;
    cvx_end
end

I = x;
fVal = cvx_optval;
dV = totConst;
fprintf('%s%f%s\n', 'Weighted least squares solution is found in ', toc, ...
    ' seconds.');
end




