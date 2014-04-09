function I = leastSquaresbyParraEtAl(sqrtQ,v)
%Finds the electrode current stimulus pattern based on least squares
%Written by: Seyhmus Guler, 3/8/14

%Section 3.1 a) Least squares with no constraints.
%Based on  the equations in the paper with the title "Optimized 
%Multi-electrode stimulation increases focality and intensity at the target 
%by Dmochowski et al, 2011.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
    % Q: A' * A = T' * Jmat' * Jmat * T. Size: LxL, L being # electrodes
    % v: A' * Jd = T' * Jmat' * Jd. Size: Lx1
%OUTPUTS:
    % I: Electrode currents Size: Lx1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %minimize_overI ||desiredJ - matrixA * I ||^2
% %solution to this least squares problem is 
% % I = inv(A'* A)*A' * desiredJ

% A = LeadFieldMatrix * T; 3MxL = (3MxN) x (NxL)
    %Current density everywhere = Lead field matrix * transfer matrix *
    %electrode currents
    
% A' = T' * LFM' 
% A = LFM' * T => A' * A = T' * (LFM' * LFM) * T
% and (A'*A) \ A' * Jd = Q \ v where
% Q = T' * (LFM' * LFM) * T and v = (T' * LFM') * Jd = T' * LFM' * Jd

Q = sqrtQ' * sqrtQ;

if size(Q,1) == size(Q,2) && size(Q,1) == size(v,1)
    I = Q \ v;
else
    error('mismatch in matrix sizes');
end

    
    

    
    




