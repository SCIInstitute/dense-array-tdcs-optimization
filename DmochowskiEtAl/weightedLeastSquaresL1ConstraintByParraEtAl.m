function I = weightedLeastSquaresL1ConstraintByParraEtAl(Q_weighted,v_weighted,Smax)
%Section 3.2 Weighted least squares with L1 constraint.
%Based on  the equations in the paper with the title "Optimized 
%Multi-electrode stimulation increases focality and intensity at the target 
%by Dmochowski et al, 2011.

%I = min_overI || sqrt(W) (Jd - A*I)||^2 subject to L1 current constr.
%See weightedLeastSquaresConstrainedByParraEtAl for more detail.

L = size(Q_weighted,2); %number of electrodes
m = Q_weighted\v_weighted;

cvx_begin
variable x(L)
expression y(L+1)
y = [x; -sum(x)];
minimize quad_form(x-m,Q_weighted);
subject to
-Smax <= y;
y <= Smax;
cvx_end
I = x;
end
