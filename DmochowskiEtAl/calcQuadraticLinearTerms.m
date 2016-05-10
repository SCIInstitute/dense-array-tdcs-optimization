function variables = calcQuadraticLinearTerms(G, T, Ed, ROI, brainLabels)
% Finds the qudratic and linear matrices that are needed for the
% optimizations in the paper by Dmochowski et al. 2011.
%
% Synopsis: [variables] = ( G, T, Ed, ROI )
%           [variables] = ( G, T, Ed, ROI, brainLabels )
%
% Input:    G      = transfer matrix from potential to electric field.
%           T      = transfer matrix from electrode currents to potential.
%           Ed     = desired electric field at the target.
%           ROI    = target region of interest.
%           brainLabels = brain labels.
%
% Output:   variables =   quadratic matrices and linear terms
%                .Q_ROI = ROI matrix
%                .Q_nonROI = nonROI matrix
%                .Q_brain = brain matrix (Q_ROI + Q_nonROI)
%                .v_brain = linear term for the brain
%                .Ed_ROI  = desired electric field at ROI.

% Notes:    1. Given a head model, an ROI, and a desired electric field at
%           the ROI, calculates certain matrices and vectors used in the
%           optimizations in Section 3 of the article:
%           Dmochowski et al., "Optimized multi-electrode stimulation 
%           increases focality and intensity at the target," 
%           Journal of neural engineering, 8.4 : 046011, 2011.
%
%           2. The relationship between the quadratic matrices and the original 
%           paper terminology:
%           Q_ROI    := A_ROI' * A_ROI = T' * G(ROI,:)' * G(ROI,:) * T
%           Q_nonROI := A_nonROI' * A_nonROI
%                     = T' * G(nonROI,:)' * G(nonROI,:) * T
%           Q_brain  := A' * A = A_ROI' * A_ROI + A_nonROI' * A_nonROI
%                     = Q_ROI + Q_nonROI
%
%           3. The reason we are calculating Q_ROI and Q_nonROI instead of
%           directly calculating Q_brain is because we would like to have
%           the flexibility of easily adjusting the weights for ROI and nonROI.

M = size(G,1)/3; % number of elements
[N,~] = size(T); % number of nodes 

if size(G,2) ~= N
    error('Size of G matrix is not correct.');
end

if size(ROI,2) ~= M
    ROI = ROI';
end


roiIdx = find(ROI==1);
eRoiIdx = sort([roiIdx*3 roiIdx*3-1 roiIdx*3-2]);

if (nargin == 5) %calculating quadratic terms require input brainLabels
    nonRoi = brainLabels;
    nonRoi(ROI==1) = 0;
    nonRoiIdx = find(nonRoi ==1);
    eNonRoiIdx = sort([nonRoiIdx*3 nonRoiIdx*3-1 nonRoiIdx*3-2]);
    
    %Quadratic matrices
    qROItemp = G(eRoiIdx,:)' * G(eRoiIdx,:); %multiply sparse matrices first
    qnonROItemp = G(eNonRoiIdx,:)' * G(eNonRoiIdx,:);
    variables.Q_ROI = T' * qROItemp * T;
    variables.Q_nonROI = T' * qnonROItemp * T;
    variables.Q_brain = variables.Q_ROI + variables.Q_nonROI;
end

%Linear terms
if size(Ed,1) == 3*nnz(ROI)
    variables.Ed_ROI = Ed;
else
    variables.Ed_ROI = Ed(eRoiIdx);
end
variables.v = T' * (G(eRoiIdx,:)' * variables.Ed_ROI);
