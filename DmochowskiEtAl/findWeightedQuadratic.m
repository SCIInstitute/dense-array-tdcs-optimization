function [Q_weighted,v_weighted] = findWeightedQuadratic(Q_ROI,Q_nonROI,v,k,numOfROINodes,numOfnonROINodes)
% Weights the contribution of ROI and nonROI regions on the objective.
%
% Synopsis: [sqrtQ_weighted, v_weighted] = findWeightedQuadratic(Q_ROI, ...
%           Q_nonROI, v, k, numOfROINodes, numOfnonROINodes)
%
% Input:    Q_ROI       =   matrix representing ROI 
%           Q_nonROI    =   matrix representing non ROI 
%           v           = right hand side vector
%           k           = Weight coefficient
%           numOfROINodes = # ROI nodes
%           numOfnonROINodes = # non ROI nodes
%
% Output:   sqrtQ_weighted  =   Chol factor of weighted audratic matrix.
%           v_weighted      = weighted right hand side. 

% Notes:    1. Use the equation in section 3.1. of " Optimized multi-electrode 
%           stimulation increases focality and intensity at the target.",
%           Jacek P Dmochowski, et al., Journal of neural engineering 
%           8.4 (2011): 046011.
%

w = k*numOfnonROINodes/numOfROINodes;
wC = 1/(w+1);
w = w/(w+1);

Q_weighted = w * Q_ROI + wC * Q_nonROI;
v_weighted = w * v;