function [ I, dV ] = optimizeIntensityByParraEtAl(C, J0, tot)
% Optimizes the electrode currents to yield maximum current intensity
% at target points.
%
% Synopsis: [ I, dV ] = optimizeIntensityByParraEtAl(C, J0, tot)
%
% Input:    C       =   matrix linking electrode currents to current
%                       density at the target node(s).
%           JO      =   desired current density orientation at the target.
%           tot     =   total current constraint
%
% Output:   I       =   array of electrode currents.
%           dV      =   dual variable(s) corresponding to constraints.

% Notes:    1. Use the equation in section 3.4. of " Optimized multi-electrode 
%           stimulation increases focality and intensity at the target.",
%           Jacek P Dmochowski, et al., Journal of neural engineering 
%           8.4 (2011): 046011.
%

v = J0' * C; % Linear coefficients for objective function
L = size(C,2); % Number of electrodes

cvx_begin quiet
cvx_solver sedumi
cvx_precision high
variable x(L);
expression y(L+1);
y = [x;-sum(x)];
dual variable totConst

maximize v * x %maximize the current projected onto the desired direction
subject to
totConst : norm(y,1) <= 2*tot;
cvx_end

if ~strcmp(cvx_status,'Solved')
    fprintf('%s\n','No solution with high precision, trying low precision.');
    cvx_begin quiet
    cvx_solver sedumi
    cvx_precision low
    variable x(L);
    expression y(L+1);
    y = [x;-sum(x)];
    dual variable totConst
    
    maximize v * x
    subject to
    totConst : norm(y,1) <= 2*tot;
    cvx_end
end

I = x;
dV = totConst;

fprintf('%s%f%s\n', 'Optimization for intensity is finished in ', toc, ...
    ' seconds.');
end

