function I = linConstraintMinVarianceByParraEtAl(Q,C,J0,Smax,opt)
%Section 3.3 Linearly constrained minimum variance
%Based on  the equations in the paper with the title "Optimized
%Multi-electrode stimulation increases focality and intensity at the target
%by Dmochowski et al, 2011.
%%INPUTS:
    %Q;
%%OUTPUTS:
    %

%I = min_overI ||A*I||^2 subject to C*I = J0. Hard linear constraint
% A = LFM * T and ||A*I||^2 = I' * A' * A * I = I' * Q * I where
% A' * A = T' * LFM' * LFM * T = Q

% if (opt = 1) Individual electrode constraint
% else (opt = 0) Total current constraint

L = size(Q,2); %number of electrodes

if opt == 1 %individual electrode constraint
    cvx_begin
    variable x(L);
    expression y(L+1);
    y = [x; -sum(x)];
    minimize quad_form(x,Q);
    subject to
    C * x == J0; %#ok<EQEFF>
    -Smax <= y;
    y <= Smax; %#ok<VUNUS>
    cvx_end
    I = x;
    
elseif opt == 0 %total current constraint
    cvx_begin
    variable x(L)
    expression y(L+1);
    y = [x; -sum(x)];
    minimize quad_form(x,Q);
    subject to 
    C * x == J0; %#ok<EQEFF>
    norm(y,1) <= 2 *Smax; %#ok<VUNUS>
    cvx_end
    I = x;
end

        
        
