function [I, fVal, dV] = linConstraintMinVarianceByParraEtAl(sqrtQ,C,J0,Smax,opt)
% Linearly constraint minimum variance solution to hd_tDCS electrode 
% current stimulus optimization with total or individual current constraint.
%
% Synopsis: [I, fVal, dV] = linConstraintMinVarianceByParraEtAl( ...
%           sqrtQ, C, J0, Smax, opt)
%
% Input:    sqrtQ   =   chol factor of quadratic matrix.
%           C       =   matrix linking electrode currents to current
%                       density at the target node(s).
%           JO      =   current density at the target.
%           Smax    =   constraint. Either tot or ind.  
%           opt     =   1 if individual constraint, 0 if total constraint. 
%
% Output:   I       =   array of electrode currents.
%           fVal    =   least squares error for best solution. 
%           dV      =   dual variables for the constraints.

% Notes:    1. Use the equation in section 3.3. of " Optimized multi-electrode 
%           stimulation increases focality and intensity at the target.",
%           Jacek P Dmochowski, et al., Journal of neural engineering 
%           8.4 (2011): 046011.
%

tic;
L = size(sqrtQ,2); %number of electrodes
if (size(J0,1) ~= size(C,1))
    error('Size mismatch between J0 and C');
end

if opt == 1 %individual electrode constraint
    if size(Smax,1) == L %In case reference electrode bound is not defined
        Smax(L+1,:) = Smax(end,:); %Make ref elec have the same const as last one
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
    indConstLB : -Smax <= y;
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
        indConstLB : -Smax <= y;
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



