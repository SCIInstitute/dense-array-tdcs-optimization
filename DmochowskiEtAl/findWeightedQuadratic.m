function [sqrtQ_weighted,v_brain_weighted] = findWeightedQuadratic(Q_ROI,Q_nonROI,v_brain,k,roi,brainLabels)
%weights the quadratic and other Jd term such that the contribution of ROI
%is larger. 

%Written by Seyhmus Guler 4/8/14
nonroi = brainLabels;
nonroi(roi==1) = 0;
w = k*nnz(nonroi)/nnz(roi);
wC = 1/(w+1);
w = w/(w+1);

Q_brain_weighted = w * Q_ROI + wC * Q_nonROI;
sqrtQ_weighted = chol(Q_brain_weighted);
v_brain_weighted = w * v_brain;