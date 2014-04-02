function I = optimizeIntensityByParraEtAl(C, J0, Smax)
%Optimizes the electrode currents to yield maximum current intensity 
%at target points.
%Written by: Seyhmus Guler, 3/8/14

%Section 3.4 Optimizing for intensity
%Based on  the equations in the paper with the title "Optimized
%Multi-electrode stimulation increases focality and intensity at the target
%by Dmochowski et al, 2011.

%%%%%%%%%%%%%%%%%
%INPUTS:
    %C: The matrix linking electrode currents to the current density at 
    %   target regions
    %J0: The desired CD at target area
    %tot: total current allowed
%OUTPUTS:
    %I: Electrode currents
%%%%%%%%%%%%%%%%%

v = J0' * C; %Linear coefficients for objective function
L = size(C,2); %Number of electrodes

cvx_begin
variable x(L);
expression y(L+1);
y = [x;-sum(x)];
maximize v * x
subject to
norm(y,1) <= 2*Smax;
cvx_end
I = x;
end

