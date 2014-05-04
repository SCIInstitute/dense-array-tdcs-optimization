function I = leastSquaresbyParraEtAl(sqrtQ,v)
% Least squares solution to hd_tDCS electrode current stimulus optimization.
%
% Synopsis: I = leastSquare(sqrtQ,v)
%
% Input:    sqrtQ   =   square root (Chol factor) of left hand side matrix.
%           v       =   right hand side vector.
%
% Output:   I       =   array of electrode currents.

% Notes:    1. Use the equation in section 3.1. of " Optimized multi-electrode 
%           stimulation increases focality and intensity at the target.",
%           Jacek P Dmochowski, et al., Journal of neural engineering 
%           8.4 (2011): 046011.
%
%           2. Below is the relationship between matrices used in the paper 
%           and formula used in this script:
%           A = G * T => A' * A = T' * (G' * G) * T
%           and (A'*A) \ A' * Jd = Q \ v where
%           Q = T' * (LFM' * LFM) * T and 
%           v = (T' * LFM') * Jd = T' * LFM' * Jd

Q = sqrtQ' * sqrtQ;

if size(Q,1) == size(Q,2) && size(Q,1) == size(v,1)
    I = Q \ v;
else
    error('Mismatch in matrix sizes.');
end

    
    

    
    




