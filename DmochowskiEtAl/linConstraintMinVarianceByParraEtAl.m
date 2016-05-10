function [I, fVal, dV] = linConstraintMinVarianceByParraEtAl(Q,C,J0,Smax,opt)
% Linearly constrained minimum variance solution to dense  array tDCS 
% stimulus pattern optimization with total injected or individual electrode 
% constraint.
%
% Synopsis: [I, fVal, dV] = linConstraintMinVarianceByParraEtAl( ...
%           Q, C, J0, Smax, opt)
%
% Input:    Q       =   quadratic matrix.
%           C       =   matrix linking electrode currents to electric
%                       field at the target node(s).
%           EO      =   electric field at the target.
%           Smax    =   Constraint bound on either total or individual
%                       electrode currents
%           opt     =   '1' if individual constraint, '0' if total constraint. 
%
% Output:   I       =   array of electrode currents.
%           fVal    =   least squares error for the best solution. 
%           dV      =   dual variables for the constraints.

% Notes:    1. The implementation is based on the equation (12) in section 3.2 
%           of the article: Dmochowski et al., "Optimized multi-electrode 
%           stimulation increases focality and intensity at the target," 
%           Journal of neural engineering, 8.4 : 046011, 2011.
%   
%           2. The solution is not based on the equation (14) in section
%           3.3 of the original paper. Instead, we use CVX to solve the problem
%           in (12).

tic;
L = size(Q,2); %number of electrodes


if (size(J0,1) ~= size(C,1))
    error('Size mismatch between J0 and C');
end

%Cholesky factor
[sqrtQ,p] = chol(Q);
if p>0
    warning('quadratic matrix is not positive definite.');
    minEig = min(eig(Q));
    sqrtQ = chol(Q + (-minEig+eps)*eye(size(Q)));
end

if opt == 1 %individual electrode constraint
    if size(Smax,1) == L %reference electrode bound is not defined
        Smax(L+1,:) = Smax(1,:);
    end
    
    cvx_begin quiet
    cvx_solver sedumi
    cvx_precision high
    
    variable x(L);
    expression y(L+1);
    y = [x; -sum(x)];
    dual variable linConst
    dual variable indConstLB
    dual variable indConstUB
    
    minimize norm(sqrtQ*x);
    subject to
    linConst : C * x == J0; %#ok<EQEFF>
    indConstLB : -Smax <= y; %#ok<VUNUS>
    indConstUB : y <= Smax; %#ok<VUNUS>
    cvx_end
    
    if ~strcmp(cvx_status,'Solved')
        fprintf('%s\n','No solution with high precision, trying lower precision.');
        cvx_begin quiet
        cvx_solver sedumi
        cvx_precision low
        
        variable x(L);
        expression y(L+1);
        y = [x; -sum(x)];
        dual variable linConst
        dual variable indConstLB
        dual variable indConstUB
        
        minimize norm(sqrtQ*x);
        subject to
        linConst : C * x == J0; %#ok<EQEFF>
        indConstLB : -Smax <= y; %#ok<VUNUS>
        indConstUB : y <= Smax; %#ok<VUNUS>
        cvx_end
    end
    
    I = x;
    fVal = cvx_optval;
    dV = [linConst; indConstLB; indConstUB];
    
    
elseif opt == 0 %total current constraint
    
    cvx_begin quiet
    cvx_solver sedumi
    cvx_precision high
    
    variable x(L)
    expression y(L+1);
    y = [x; -sum(x)];
    dual variable linConst
    dual variable totConst
    
    minimize norm(sqrtQ*x);
    subject to
    linConst : C * x == J0; %#ok<EQEFF>
    totConst : norm(y,1) <= 2 * Smax; %#ok<VUNUS>
    cvx_end
    
    if ~strcmp(cvx_status,'Solved')
        fprintf('%s\n','No solution with high precision, trying lower precision.');
        cvx_begin quiet
        cvx_solver sedumi
        cvx_precision low
        
        variable x(L)
        expression y(L+1);
        y = [x; -sum(x)];
        dual variable linConst
        dual variable totConst
        
        minimize norm(sqrtQ*x);
        subject to
        linConst : C * x == J0; %#ok<EQEFF>
        totConst : norm(y,1) <= 2 * Smax; %#ok<VUNUS>
        cvx_end
    end
    I = x;
    fVal = cvx_optval;
    dV = [linConst; totConst];
end

fprintf('%s%f%s\n', 'Linearly constrained minimum variance solution is found in ', toc, ...
    ' seconds.');

end



