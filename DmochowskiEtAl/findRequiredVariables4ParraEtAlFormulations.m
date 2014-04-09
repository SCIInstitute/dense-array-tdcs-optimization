function variables = findRequiredVariables4ParraEtAlFormulations(G,T,Jd,roi,brainLabels)
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
[N,~] = size(T);

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
variables.Q_ROI = T' * qROItemp * T;
variables.Q_nonROI = T' * qnonROItemp * T;

variables.Q_brain = variables.Q_ROI + variables.Q_nonROI;
if size(Jd,1) == 3*nnz(roi)
    variables.Jd_ROI = Jd;
else
    variables.Jd_ROI = Jd(expandedROI);
end
variables.v_brain = T' * (G(expandedROI,:)' * Jd_ROI);


