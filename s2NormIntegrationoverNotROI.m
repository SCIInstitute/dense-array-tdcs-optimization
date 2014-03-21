function S = s2NormIntegrationoverNotROI(avoidRegions,elemVolumes)
%THIS FUNCTION IS USED TO FIND THE INTEGRAL OF NORM 3 SQUARED OF CURRENT 
% DENSITY OVER THE AVOIDANCE REGIONS
%
%Written by: Seyhmus Guler, Revisited: Moritz Dannhauer
%Last edit: 10/7/13 by Guler,S
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
    %avoidRegions: matrix, each row determining different avoidance region
    %elemVolumes: volumes of the elements
%OUTPUTS:
    %S:Volumes of the elements in avoidRegions, purmutated with replacement
        %to be used in integral. size: cell array, each cell containing 
        %(3 #elements in notROIr) x (3 #elements in notROI) matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%ADDITIONAL NOTES:
%By convention, the current density array is 3 #elements x 1 where
%(3i-2)'th to 3i'th row is the current density in the i'th element.
%Assume that we know current density array.
%To find s2norm integration over notRoi we need to rescale current density of each element
%by the volume of that element since the intagration is just a sum when 
%the function to be integrated is piece-wise linear. The resulting matrix of this function 
%is used for rescaling purpose. And by multiplying this matrix with current
%density (of notROI elements) transpose from left and current density (of notROI elements) 
%from right, we will get the power constraint integration
S = cell(1,size(avoidRegions,1));
for i = 1:size(avoidRegions,1)    
nZ = nnz(avoidRegions(i,:));
notROIElementVolumes = elemVolumes(avoidRegions(i,:) ==1);
S{i} = sparse(1:3*nZ,1:3*nZ,reshape(repmat(notROIElementVolumes,3,1),1,[]),3*nZ,3*nZ);
end