function [I, fVal, dV] = weightedLeastSquaresTotalIndividual(Q, b, tot, ind)
% Finds hd-tDCS electrode current stimulus pattern using weighted least squares
% with total and individual electrode currents.
%
%
% Synopsis: [I, fval, dV] = weightedLeastSquaresTotalIndividual( ...
%           sqrtQ_weighted, v_weighted, tot, ind)
%
% Input:    sqrtQ_weighted  = chol factor of weighted quadratic matrix.
%           v_weighted      = weighted right hand side vector.
%           tot             = total applied current bound.
%           ind             = individual electrode current bound.
%
% Output:   I       =   array of electrode currents.
%           fVal    =   least squares error for best solution. 
%           dV      =   dual variables for the constraints.

% Notes:    1. This implementation is to mimic the results in the paper 
%              'Optimization of multifocal transcranial current stimulation
%              for weighted cortical pattern targeting from realistic 
%              modeling of electric fields' by Ruffini et al., 2014.
%
%           2. We convert quadratic objective function to norm 2 in order
%           to speed up the convergence. This is the reason behind using
%           Cholesky factor of the quadratic matrix instead of the matrix
%           itself. minimize(x' * Q * x + b' * x) is equivalent to
%           minimize( norm(chol(Q) * (x + m), 2 ) ) where m = 2*Q'\v.

tic;
L = size(Q,2); %number of electrodes
[sqrtQ,p] = chol(Q);
if p>0
    warning('quadratic matrix is not positive definite.');
    minEig = min(eig(Q));
    sqrtQ = chol(Q + (-minEig+eps)*eye(size(Q)));
end
m = 2*Q'\b;

if size(ind,1) == L
    ind(L+1,:) = ind(end,:);
end

cvx_begin quiet
cvx_solver sedumi
cvx_precision high
variable x(L)
expression y(L+1)
y = [x; -sum(x)];

dual variable totConstV
dual variable indConstLB
dual variable indConstUB

minimize norm(sqrtQ * (x+m));
subject to
totConstV : norm(y,1) <= 2*tot;
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
    
    dual variable totConstV
    dual variable indConstLB
    dual variable indConstUB
    
    minimize norm(sqrtQ * (x+m));
    subject to
    totConstV : norm(y,1) <= 2*tot;
    indConstLB : -ind <= y;
    indConstUB : y <= ind;
    cvx_end
end

if ~strcmp(cvx_status,'Solved')
    fprintf('%s\n','No solution with sedumi, trying sdpt3.');
    
    cvx_begin quiet
    cvx_solver SDPT3
    variable x(L)
    expression y(L+1)
    y = [x; -sum(x)];
    
    dual variable totConstV
    dual variable indConstLB
    dual variable indConstUB
    
    minimize norm(sqrtQ * (x+m));
    subject to
    totConstV : norm(y,1) <= 2*tot;
    indConstLB : -ind <= y;
    indConstUB : y <= ind;
    cvx_end
end

I = x;
fVal = cvx_optval;
dV = [totConstV; indConstLB; indConstUB];
fprintf('%s%f%s\n', 'Weighted least squares solution is found in ', toc, ...
    ' seconds.');
end
