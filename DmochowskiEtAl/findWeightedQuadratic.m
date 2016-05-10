function [Q_weighted,v_weighted] = findWeightedQuadratic(Q_Roi,Q_nonRoi,v,k,numOfRoiNodes,numOfNonRoiNodes)
% Weights the contribution of ROI and nonROI regions on the objective.
%
% Synopsis: [Q_weighted, v_weighted] = findWeightedQuadratic(Q_ROI, ...
%           Q_nonROI, v, k, numOfROINodes, numOfnonROINodes)
%
% Input:    Q_Roi            =  matrix representing ROI 
%           Q_nonRoi         =  matrix representing nonROI 
%           v                =  right hand side vector
%           k                =  weight coefficient
%           numOfRoiNodes    =  # ROI nodes
%           numOfNonRoiNodes =  # nonROI nodes
%
% Output:   Q_weighted      =   weighted quadratic matrix.
%           v_weighted      =   weighted right hand side vector. 

% Notes:    1. The implementation is based on the equation (8) in section 3.1 
%           of the article: Dmochowski et al., "Optimized multi-electrode 
%           stimulation increases focality and intensity at the target," 
%           Journal of neural engineering, 8.4 : 046011, 2011.

w = k*numOfNonRoiNodes/numOfRoiNodes;
wC = 1/(w+1);
w = w/(w+1);

Q_weighted = w * Q_Roi + wC * Q_nonRoi;
v_weighted = w * v;