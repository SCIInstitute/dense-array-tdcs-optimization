function [ w, Q ]  = linearQuadraticCoefficientCalculation(ROI, avoidRegions, desiredDirection, T, G, elemVolumes)
%% Calculates the linear and quadratic terms for the objective function and
% power constraint of hd-tDCS electrode current optimization. The objective
% function is linear and power constraint is quadratic functions of electrode
% currents.
%
% Written by: Seyhmus Guler 4/29/14.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M := # Elements in the mesh.
% N := # Nodes in the mesh.
% L := # Electrodes (reference excluded so L = #electrodes on the head - 1).
%INPUTS:
%ROI: Region of interest. Logical vector of size 1 x M.
%avoidRegions: critical regions for which we impose power constraints.
%Logical matrix of size (# AvoidRegions) x M.
%desiredDirection: desired modulation direction in the ROI.
% Size: either 3x1 or 3 x (# ROI elements)
%T: Transfer matrix from electrode current array to node potentials.
% This matrix is precalculated using forward problem solution.
% Size: N x L
%G: matrix linking the cell current density to node potentials
% Size: 3M x N
%elemVolumes: element volumes
% Size 1 x M
%OUTPUTS:
%w: The row vector of the linear coefficients of objective function.
% Size: 1 x L
%sqrtQ: Square root of the quadratic term matrices. Size: cell(1, #
% Avoid regions). For each avoid region we have LxL matrix in the
% corresponding cell.
%NOTES:
% Inputs G and elemVolumes can be calculated using
% anisomappingFromNodePotentialsToCurrentDensity() function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

%%check the sizes. If no match, give an error.
if size(G,2) ~= size(T,1) || ~isequal(numel(ROI),size(avoidRegions,2),numel(elemVolumes))
    error('Mismatch in matrix sizes.');
end
if size(desiredDirection,2) == numel(ROI)
  desiredDirection = desiredDirection(:,ROI==1);
end
% expanded ROI,expanded notRoi indices are needed since by convention we chose
% current density vector to be size of (3 #elements) x 1
% Current density = [x1; y1; z1; x2; y2; z2; ... ; xM; yM; zM]

roiIdx = find(ROI==1);
expandedROI = sort([3*roiIdx 3*roiIdx-1 3*roiIdx-2]);
numOfAvoid = size(avoidRegions,1);
expandedNotROI = cell(1,numOfAvoid);
for i = 1:numOfAvoid
    avoidRegionIdx = find(avoidRegions(i,:)==1);
    expandedNotROI{i} = sort([3*avoidRegionIdx 3*avoidRegionIdx-1 3*avoidRegionIdx-2]);
end


%Linear term: w * J = a * G * T * I where J = G * T * I is the current
% density and w is the linear coefficients. We know G, T (inputs). We need to
% calculate a, which is inner product of roi element volumes and desired
% modulation direction.
a = weightedInnerProductSumOverROI(desiredDirection,elemVolumes(roiIdx));

% Quadratic term: I' * Q * I = I' * (T' * G' * S * G * T) * I where S is
% the matrix that weights each element according to its volume. We
% calculate S here, for each avoidRegion.
S = s2NormIntegrationoverNotROI(avoidRegions,elemVolumes);


%% Combination of matrices
%We first multiply sparse matrices.
wTtemp = a * G(expandedROI,:); %intermadiate step
clear expandedROI a;
%We then multiply with full matrix T.
w = wTtemp * T;
clear wTtemp;
fprintf('%s%f%s\n','w is calculated in ',toc,' seconds.');

%
Q = cell(1,numOfAvoid);
for i = 1:numOfAvoid
    %First sparse terms.
    Qtemp = G(expandedNotROI{i},:)' * S{i} * G(expandedNotROI{i},:);
    %Then full terms.
    Q{i} = T' * Qtemp * T;
end

fprintf('%s%f%s\n','sqrtQ is calculated in ',toc,' seconds.');
end

function w = weightedInnerProductSumOverROI(desiredDirection,ROIelemVolumes)
%% FINDS THE LINEAR WEIGHTS FOR THE INTEGRAL OF THE DIRECTIONAL CURRENT
%DENSITY IN THE ROI
%
%Written by: Seyhmus Guler, Revised: Moritz Dannhauer
%Last edit: 4/28/13 by Guler,S
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
%desiredDirection: The predefined direction matrix.
% size: either 3x1 if the desired direction is fixed through ROI
%  or 3 x (# elements in ROI) if the desired direction changes
%  through ROI
%ROIelemVolumes: ROI element volumes.
% size: 1 x # (elements in the ROI).
%OUTPUTS:
%w: The row vector for weights. If dotted with the current density
%of the ROI, it will give weighted integral of directional
%current density in the ROI
%size: 1 x (3 #elements in the ROI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
nZ = numel(ROIelemVolumes); % number of ROI elements

if ~isequal(size(desiredDirection),[3 nZ]) && ~isequal(size(desiredDirection),[3 1])
    error('Mismatch between ROI size and desired modulation direction size.');
end


w = zeros(1,3 * nZ);
idx = 3:3:3*nZ;

w(idx-2) = desiredDirection(1,:) .* ROIelemVolumes;
w(idx-1) = desiredDirection(2,:) .* ROIelemVolumes;
w(idx) = desiredDirection(3,:) .* ROIelemVolumes;
end

function S = s2NormIntegrationoverNotROI(avoidRegions,elemVolumes)
%% THIS FUNCTION IS USED TO FIND THE INTEGRAL OF NORM 2 SQUARED OF CURRENT
% DENSITY OVER THE AVOIDANCE REGIONS
%
%Written by: Seyhmus Guler, Revised: Moritz Dannhauer
%Last edit: 10/7/13 by Guler,S
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
%avoidRegions: matrix, each row determining different avoidance region
%elemVolumes: volumes of the elements
%OUTPUTS:
%S:Volumes of the elements in avoidRegions, purmutated with replacement
%to be used in integral. size: cell array, each cell containing
%(3 #elements in avoid region x) x (3 # elements in avoid region x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADDITIONAL NOTES:
%By convention, the current density array is 3 #elements x 1 where
%(3i-2)th, (3i-1)th and (3i)th values give the current density of the i'th element.
%Assume that we know current density array.
%To find s2norm integration over avoid regions, we need to rescale current density of each element
%by the volume of that element since the intagration is just a sum when
%the function to be integrated is piece-wise constant. The resulting matrix of this function
%is used for rescaling purpose. And by multiplying this matrix with current
%density (of avoid region elements) transpose from left and current density (of avoid region elements)
%from right, we will get the current power in that avoid region.
numOfAvoid = size(avoidRegions,1);
S = cell(1,numOfAvoid);
for i = 1:numOfAvoid
    nZ = nnz(avoidRegions(i,:));
    notROIElementVolumes = elemVolumes(avoidRegions(i,:) ==1);
    S{i} = sparse(1:3*nZ,1:3*nZ,reshape(repmat(notROIElementVolumes,3,1),1,[]),3*nZ,3*nZ);
end
end
