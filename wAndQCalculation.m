function [w,Q] = wAndQCalculation(ROI,avoidRegions,desiredDirection,T,G,elemVolumes)
%CALCULATES THE MATRICES NEEDED IN THE OPTIMIZATION PROBLEM
%RUNS ONLY ONCE BEFORE THE OPTIMIZATION
%
%Written by: Seyhmus Guler, Revisited: Moritz Dannhauer
%Last edit: 3/26/13 by Guler,S
%
%Needs faster ways to do multiplications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
    %ROI: region of interest
    %avoidRegions: critical regions for which we impose power constraints
    %desiredDirection: desired modulation direction in the ROI. 
        % The component of the current density along this particular 
        % direction in the ROI will be the objective function
    %T: Transfer matrix from current array to potential. This matrix is 
        %precalculated using forward problem solution
    %G: matrix linking the current density to node potentials 
    %elemVolumes: element volumes
        %Note that JfromU and elemVolumes can be calculated using 
        %anisomappingFromNodePotentialsToCurrentDensity() function
%OUTPUTS:
    %w: The row vector from current array to directional current
        %density in the ROI. size: 1 x #electrodes
    %Q: the transfer matrices from current array to norm 2 integral
        % of current density over avoid regions. These matrices will be
        %used to set power constraints over different regions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
%%check if the sizes of the inputs match,if not give an error
if size(G,2) ~= size(T,1) || ~isequal(numel(ROI),size(avoidRegions,2),numel(elemVolumes))
    error('mismatch in matrix sizes');
end

% calculate matrices 
% expanded ROI,notRoi indices are needed since by convention we chose
% current density vector to be size of (3 #elements) x 1
expandedROI = sort([3*find(ROI==1) 3*find(ROI==1)-1 3*find(ROI==1)-2]);
nnRoi = size(avoidRegions,1);
expandedNotROI = cell(1,nnRoi);
for i = 1:nnRoi
expandedNotROI{i} = sort([3*find(avoidRegions(i,:)==1) 3*find(avoidRegions(i,:)==1)-1 3*find(avoidRegions(i,:)==1)-2]);
end


%integration vector for objective
a = weightedInnerProductSumOverROI(desiredDirection,ROI,elemVolumes);

%integration matrix for power constraint
S = s2NormIntegrationoverNotROI(avoidRegions,elemVolumes);

% combination of matrices
%whole transfer for directional integration is a * J(roiElementIndices,:) * T 
wTtemp = a * G(expandedROI,:); %intermadiate step
clear expandedROI a;
w = wTtemp * T;
clear wTtemp;
fprintf('%s%f%s\n','w is calculated in ',toc,' seconds.');

Q = cell(1,nnRoi);
for i = 1:nnRoi
Qtemp = G(expandedNotROI{i},:)' * S{i} * G(expandedNotROI{i},:); %intermadiate step
Qtemp2 = T' * Qtemp * T;
Q{i} = chol(Qtemp2);
end

fprintf('%s%f%s\n','sqrtQ is calculated in ',toc,' seconds.');


