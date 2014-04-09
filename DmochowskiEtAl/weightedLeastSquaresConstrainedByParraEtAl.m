function [I,dV] = weightedLeastSquaresConstrainedByParraEtAl(sqrtQ_weighted,v_weighted,Smax)
%Finds electrode currents based on least squares with weights.
%Written by: Seyhmus Guler, 3/8/14

%Section 3.1 b) Weighted least squares with norm1 constraint.
%Based on  the equations in the paper with the title "Optimized
%Multi-electrode stimulation increases focality and intensity at the target
%by Dmochowski et al, 2011.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%INPUTS:
%Q_weighted : T' * Jmat' * W * Jmat * T. Weighted quadratic term. Size: LxL
%v_veighted : T' * Jmat' * W * Jd. Weighted linear term. Size: Lx1

%I = min_overI || sqrt(W) (Jd - A*I)||^2 subject to total current constr.
% (Jd - A*I)' * sqrt(W)' * sqrt(W)* (Jd -A*I)
% = (Jd - A*I)' * W * (Jd - A*I) = (Jd'*W*Jd + I' * A'*W*A * I - I'*A'*W*Jd
% - Jd'*W*A*I = Jd'*W*Jd + I'*(A'*W*A)*I - 2 (Jd'*W*A)*I

% minimizing above expression is equivalent to minimizing:
% I'* (A'*W*A)* I - 2(Jd'*W*A) * I =~ I' * Q * I - 2 * v' * I

% A = LFM * T =>
% Q = T' * LFM' * W * LFM * T
% v' = Jd' * W * LFM * T

% I'QI - 2v'I = (I - inv(Q)v)' Q (I - inv(Q)v) + constant

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
totConst : norm(y,1) <= 2*Smax;
cvx_end
if ~strcmp(cvx_status,'Solved')
    fprintf('%s\n','No solution with high precision, trying low precision.');
    cvx_begin quiet
    cvx_solver sedumi
    cvx_precision high
    variable x(L)
    expression y(L+1)
    y = [x; -sum(x)];
    dual variable totConst
    minimize norm(sqrtQ_weighted*(x-m));
    subject to
    totConst : norm(y,1) <= 2*Smax;
    cvx_end
end
I = x;
dV = totConst;
end




