function I = leastSquaresbyParraEtAl(Q,v)
% Least squares solution to dense tDCS optimization with no constraints
%
% Synopsis: I = leastSquare(sqrtQ,v)
%
% Input:    Q       =   quadratic matrix for the objective funtion.
%           v       =   right hand side vector.
%
% Output:   I       =   array of electrode currents.

% Notes:    1. The implementation is based on the equation (6) in section 3.1 
%           of the article: Dmochowski et al., "Optimized multi-electrode 
%           stimulation increases focality and intensity at the target," 
%           Journal of neural engineering, 8.4 : 046011, 2011.
%
%           2. Below is the relationship between the matrices used in the paper 
%           and formula used in this script:
%           A = G * T => A' * A = T' * (G' * G) * T
%           and (A'*A) \ A' * Ed = Q \ v where
%           Q = A' * A = T' * (G' * G) * T and 
%           v = (T' * G') * Ed = T' * (G' * Ed)

if size(Q,1) == size(Q,2) && size(Q,1) == size(v,1)
    I = Q \ v;
else
    error('Mismatch in matrix sizes.');
end

    
    

    
    




