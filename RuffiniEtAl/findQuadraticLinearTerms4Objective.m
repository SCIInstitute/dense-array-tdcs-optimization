function [Q,b] = findQuadraticLinearTerms4Objective(mapV2E, T, roi, surfNormal, tmap, tMin, E0) 
%Finds quadratic and linear terms for the objective function 
%
%
% Synopsis: [Q,b] = findQuadraticLinearTerms4Objective(headMesh, mapV2E, ...
%                       roi, surfaceNormal)
%
%
% Input:    headMesh        = tetrahedral head mesh.
%               .cell       = tetrahedra definition matrix. Size: 4 x #cells
%               .node       = node coordinates.
%           mapV2E          = mapping from potential to electric field
%           T               = transfer matrix
%           roi             = region of interest
%           surfaceNormal   = surface normal for the ROI elements.
%           tmap            = tmap values for the ROI elements.
%           E0              = desired electric field.
%
%
% Output:   Q           = quadratic term in the objective
%           b           = linear coefficients in the objective

M = size(mapV2E,1)/3; % # elements
N = size(mapV2E,2); % # nodes
L = size(T,2); % # electrodes

if size(T,1) ~= N
    error('mismatch in transfer matrix size');
end
if size(roi,2) ~= M
    roi = roi';
end

if (size(surfNormal,2) == M)
    surfNormal = surfNormal(:,roi==1);
end

if (numel(tmap) == M)
    tmap = tmap(roi==1);
end

if isempty(tMin)
    tMin = 2;
end

if numel(E0) == M
    E0 = E0(roi==1);
end
    

%% Calculate mapping matrices

tmap(abs(tmap)<tMin) = 0; %thresholding tmap
W = abs(tmap);

roiIdx = find(roi==1);
eRoiIdx = sort([3*roiIdx 3*roiIdx-1 3*roiIdx-2]);

surfNormalInCells = mat2cell(surfNormal',ones(size(surfNormal,2),1),3);
mapOntoSurfNormal = sparse(blkdiag(surfNormalInCells{:}));

%normal component of electric field (tmap weighted) 
Ewtemp = (sparse(diag(W)) * mapOntoSurfNormal * mapV2E(eRoiIdx,:));
Ew = Ewtemp * T;

%desired normal component of electric field (tmap weighted)
Yw = diag(tmap) * E0; %Yw

k = numel(W) / sum(W.^2);

Q = k * (Ew' * Ew) ;
b = -2 * k * (Yw' * Ew)';











