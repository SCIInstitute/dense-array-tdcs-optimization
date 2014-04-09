function [I,fval,dV] = weightedLeastSquaresL1ConstraintByParraEtAl(sqrtQ_weighted,v_weighted,Smax)
%Section 3.2 Weighted least squares with L1 constraint.
%Based on  the equations in the paper with the title "Optimized
%Multi-electrode stimulation increases focality and intensity at the target
%by Dmochowski et al, 2011.

%I = min_overI || sqrt(W) (Jd - A*I)||^2 subject to L1 current constr.
%See weightedLeastSquaresConstrainedByParraEtAl for more detail.

L = size(sqrtQ_weighted,2); %number of electrodes
Q = sqrtQ_weighted' * sqrtQ_weighted;
m = Q\v_weighted;

if size(Smax,1) == L
    Smax(L+1,:) = Smax(end,:);
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
indConstLB : -Smax <= y;
indConstUB : y <= Smax;
cvx_end
if ~strcmp(cvx_status,'Solved')
    fprintf('%s\n','No solution with high precision, trying low precision.');
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
    indConstLB : -Smax <= y;
    indConstUB : y <= Smax;
    cvx_end
end
I = x;
fval = cvx_optval;
dV = [indConstLB; indConstUB];
end
