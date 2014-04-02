function variables = findRequiredVariables4ParraEtAlFormulations(G,T,Jd,roi,k,brainLabels)
%FINDS THE MATRICES NEEDED IN DIFFERENT OPTIMIZATION SCENARIOS
%Based on  the equations in the paper with the title "Optimized 
%Multi-electrode stimulation increases focality and intensity at the target 
%by Dmochowski et al, 2011. 

%%INPUTS;
%optimizationCase
% 0: least squares with no constraint. Sec 3.1 a)
% 1: weighted least squares with total current contraint. Sec 3.1 b)
% 2: weighted least squares with individual current constraint. Sec 3.2
% 3: Linearly constrained minimum variance. Sec 3.3
% 4: Intensity optimization. Sec 3.4
%parameters
% A: 
M = size(G,1)/3;
[N,L] = size(T);

if size(G,2) ~= N
    error('size of G matrix is wrong');
end
if size(roi,2) ~= M
    roi = roi';
end

%Variables

%%3.1 a)
%Q_brain = Q_ROI + Q_notROI is the square LxL matrix from which the current
%power in the brain is calculated. Used in 3.1 a)
%Q_brain = A_ROI' * A_ROI + A_nonROI' * A_nonROI (Q_ROI, Q_nonROI defined
%accordingly) = T' * G(ROI,:)' * G(ROI,:) * T + similar_for_nonROI
%v_brain = v_ROI = A_ROI' * Jd_ROI. = T' *  G(ROI,:)' * Jd_ROI.
roiIdx = find(roi==1);
nonroi = brainLabels;
nonroi(roi==1) = 0;
nonroiIdx = find(nonroi ==1);
expandedROI = sort([roiIdx*3 roiIdx*3-1 roiIdx*3-2]);
expandednonROI = sort([nonroiIdx*3 nonroiIdx*3-1 nonroiIdx*3-2]);

qROItemp = G(expandedROI,:)' * G(expandedROI,:);
qnonROItemp = G(expandednonROI,:)' * G(expandednonROI,:);
Q_ROI = T' * qROItemp * T;
Q_nonROI = T' * qnonROItemp * T;

Q_brain = Q_ROI + Q_nonROI;
if size(Jd,1) == 3*nnz(roi)
    Jd_ROI = Jd;
else
    Jd_ROI = Jd(expandedROI);
end
v_brain = T' * (G(expandedROI,:)' * Jd_ROI);
%%3.1 b)
%Q_brain_weighted = w * Q_ROI + wc * Q_nonROI
%v_brain_weighted = w * v_brain = w * v_ROI
%%
w = k*nnz(nonroi)/nnz(roi);
wC = 1/(w+1);
w = w/(w+1);
Q_brain_weighted = w * Q_ROI + wC * Q_nonROI;
v_brain_weighted = w * v_brain;
%%3.2 Similar to 3.1 b)

%%3.3
%Q_lcmv = Q_ROI + Q_nonROI
Q_lcmv = Q_brain;
variables = [];
%%3.4
%Nothing. 

% roiIndices = find(roi==1);
% diagVector = wC * ones(numel(roi)*3,1);
% diagVector([roiIndices*3 roiIndices*3-1 roiIndices*3-2]) = w;
% n = length(diagVector);
% W = spdiags(diagVector(:),0,n,n);



