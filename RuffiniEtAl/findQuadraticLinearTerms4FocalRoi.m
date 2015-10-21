function [Q,b] = findQuadraticLinearTerms4FocalRoi(G, T, roi, cortex, surfNormal, E0)
%Finds quadratic and linear terms for the objective function for focal ROI.
%
%
% Synopsis: [Q,b] = findQuadraticLinearTerms4Objective(headMesh, mapV2E, ...
%                       roi, surfaceNormal)
%
%
% Input:
%           G               = mapping from potential to electric field
%           T               = transfer matrix
%           roi             = region of interest
%           cortex          = cortex definition
%           surfaceNormal   = surface normal for the head.
%           E0              = desired normal to the cortex electric field component.
%
%
% Output:   Q           = quadratic term in the objective
%           b           = linear coefficients in the objective

M = size(G,1)/3; % # elements
N = size(G,2); % # nodes
L = size(T,2); % # electrodes

if size(T,1) ~= N
    error('mismatch in transfer matrix size');
end
if size(roi,2) ~= M
    roi = roi';
end

if size(cortex,2) ~= M
    cortex = cortex';
    cortex(roi==1) = 0;
end

if (size(surfNormal,2) ~= M)
    error('size of surface normal should be of size #elements');
end

if size(E0,2) == M
    E0 = E0';
end
if numel(E0) == M
    E0 = E0(roi==1);
end

if numel(unique(E0)) == 1
    E0 = E0(1);
end

%% Calculate mapping matrices
roiSurfNormal = surfNormal(:,roi==1);
cortexSurfNormal = surfNormal(:,cortex==1);
clear surfNormal;

roiIdx = find(roi==1);
eRoiIdx = sort([3*roiIdx 3*roiIdx-1 3*roiIdx-2]);
cortexIdx = find(cortex==1);
eCortexIdx = sort([3*cortexIdx 3*cortexIdx-1 3*cortexIdx-2]);
nr = nnz(roi);
nc = nnz(cortex);
clear roiIdx cortexIdx roi cortex;
Gr = G(eRoiIdx,:);
Gc = G(eCortexIdx,:);
clear G;

surfNormalInCells4ROI = mat2cell(roiSurfNormal',ones(size(roiSurfNormal,2),1),3);
mapOntoSurfNormal4ROI = sparse(blkdiag(surfNormalInCells4ROI{:}));
clear surfNormalInCells4ROI;
surfNormalInCells4cortex = mat2cell(cortexSurfNormal',ones(size(cortexSurfNormal,2),1),3);
mapOntoSurfNormal4cortex = sparse(blkdiag(surfNormalInCells4cortex{:}));
clear surfNormalInCells4cortex;

%normal component of electric field (tmap weighted)
EwROI = 2 * mapOntoSurfNormal4ROI * Gr;
EwCortex = mapOntoSurfNormal4cortex * Gc;
clear Gr Gc;

%desired normal component of electric field 
Yw = 2 * E0; %Yw

k = N / (4 * nr + nc);

Qtemp = k * (4*EwROI' * EwROI + EwCortex' * EwCortex) ;
Q = T' * Qtemp * T;
btemp = -4 * k * E0 * sum(EwROI);
b = (btemp * T)';











