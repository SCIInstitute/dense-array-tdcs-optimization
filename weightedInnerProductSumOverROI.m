function w = weightedInnerProductSumOverROI(desiredDirection,ROI,elemVolumes)
%FINDS THE MATRIX USED TO FIND THE INTEGRAL OF THE DIRECTIONAL CURRENT 
%DENSITY IN THE ROI
%
%Written by: Seyhmus Guler, Revisited: Moritz Dannhauer
%Last edit: 10/7/13 by Guler,S
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
    %desiredDirection: The predefined direction matrix. 
        % size: either 3x1 if the desired direction is fixed through ROI 
        %  or 3 x (# elements in ROI) if the desired direction changes
        %  through ROI
    %ROI: vector determining the area of interest.
        % size: 1 x #elements. the i'th value is 1 if i'th element is
            %  in roi or 0 if i'th element is not in roi
    %elemVolumes: the element volumes. size: 1 x #elements
%OUTPUTS:
    %w: The integration row vector. If multiplied with current density
        %array corresponding to ROI, it will give weighted integral of inner 
        %product of desired direction and current density over ROI
        %size: 1 x (3 #elements in ROI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
ROIElementVolumes = elemVolumes(ROI==1);
nZ = nnz(ROI);
if size(desiredDirection,2) == 1
    desiredDirection = reshape(repmat(desiredDirection,1,nZ),1,[]);
else
    if size(desiredDirection,2) ~= nZ
        error('mismatch between number of roi elements and desired directions');
    else
        desiredDirection = reshape(desiredDirection,1,[]);
    end
end

w = reshape(repmat(ROIElementVolumes,3,1),1,[]) .* desiredDirection;
fprintf('%s%f%s\n','integration vector for objective function is created in ',toc,' seconds.');
