function [I,dV] = linConstraintMinVarianceByParraEtAl(sqrtQ,C,J0,Smax,opt)
%Section 3.3 Linearly constrained minimum variance
%Based on  the equations in the paper with the title "Optimized
%Multi-electrode stimulation increases focality and intensity at the target
%by Dmochowski et al, 2011.
%%INPUTS:
%sqrtQ:
%C:
%J0:
%Smax:
%opt:
%%OUTPUTS:
%I: current Array

%I = min_overI ||A*I||^2 subject to C*I = J0. Hard linear constraint
% A = LFM * T and ||A*I||^2 = I' * A' * A * I = I' * Q * I where
% A' * A = T' * LFM' * LFM * T = Q

% if (opt = 1) Individual electrode constraint
% else (opt = 0) Total current constraint

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
        fprintf('%s\n','No solution with high resolution, trying lower precision.');
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
        fprintf('%s\n','No solution with high resolution, trying lower precision.');
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
    dV = [linConst; totConst];
end
end



