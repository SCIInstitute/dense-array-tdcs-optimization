function w = weightedInnerProductSumOverROI(desiredDirection,ROIelemVolumes)
%FINDS THE LINEAR WEIGHTS FOR THE INTEGRAL OF THE DIRECTIONAL CURRENT 
%DENSITY IN THE ROI
%
%Written by: Seyhmus Guler, Revisited: Moritz Dannhauer
%Last edit: 3/24/13 by Guler,S
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


if size(desiredDirection,2) ~= nZ && size(desiredDirection) ~= 1
    error('Mismatch between ROI size and desired modulation direction size.');
end

w = zeros(1,3 * nZ);
idx = 3:3:3*nZ;

w(idx-2) = desiredDirection(1,:) .* ROIelemVolumes;
w(idx-1) = desiredDirection(2,:) .* ROIelemVolumes;
w(idx) = desiredDirection(3,:) .* ROIelemVolumes;

fprintf('%s%f%s\n','Weights for objective function is created in ',toc,' seconds.');
