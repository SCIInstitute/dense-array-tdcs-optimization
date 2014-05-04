function [I, fVal, dV] = weightedLeastSquaresConstrainedByParraEtAl(sqrtQ_weighted, v_weighted, tot)
% Weighted least squares solution to hd_tDCS electrode current stimulus 
% optimization with total current constraint.
%
% Synopsis: [I, fVal, dV] = weightedLeastSquaresConstrainedByParraEtAl( ...
%           sqrtQ_weighted, v_weighted, tot)
%
% Input:    sqrtQ_weighted  = chol factor of weighted quadratic matrix.
%           v_weighted      = weighted linear term.
%           tot             = total current bound.
%
% Output:   I       =   array of electrode currents.
%           fVal    =   least squares error for best solution. 
%           dV      =   dual variables for the constraints.

% Notes:    1. Use the equation in section 3.1. of " Optimized multi-electrode 
%           stimulation increases focality and intensity at the target.",
%           Jacek P Dmochowski, et al., Journal of neural engineering 
%           8.4 (2011): 046011.
%
%           2. We convert quadratic objective function to norm 2 in order
%           to speed up the convergence. This is the reason behind using
%           Cholesky factor of the quadratic matrix instead of the matrix
%           itself. minimize(x' * Q * x - 2*v*x) is equivalent to
%           minimize( norm(chol(Q) * (x - m), 2 ) ) where m = Q\v.
%           
%           3. The conversion from paper equation to the formula in the
%           script:
%           minimize_over_I ( || sqrt(W) (Jd - A*I) || ) ^2
%           = (Jd - A*I)' * sqrt(W)' * sqrt(W)* (Jd -A*I)
%           = (Jd - A*I)' * W * (Jd - A*I)
%           = Jd'*W*Jd + I'*A'*W*A*I - I'*A'*W*Jd - Jd'*W*A*I 
%           = Jd'*W*Jd + I'*(A'*W*A)*I - 2*(Jd'*W*A)*I
%           = constant + I'*(  Q   )*I - 2*(   v   )*I  where
%           Q := A' * W * A, v := Jd' * W * A.
%
%           Minimizing above expression is equivalent to minimizing
%           (according to note 2. and assuming sqrtQ'*sqrtQ = Q) :
%           norm ( sqrtQ * (x - inv(Q) * v)).

tic;
L = size(sqrtQ_weighted,2); %number of electrodes
Q = sqrtQ_weighted' * sqrtQ_weighted;
m = Q\v_weighted;

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




