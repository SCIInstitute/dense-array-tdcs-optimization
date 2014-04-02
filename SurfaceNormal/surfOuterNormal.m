function [surfNormal, surfElements] = surfOuterNormal(elem,node,volSurface)
%Finds the outer surface normal on the surface of a mesh. 
%Also finds surface element centers and indices for the surface triangles to be used as 
%nodes for surface normal mesh.
%elem: 4xM, node: 3xN
if nargin < 3
volSurface = volface(elem');
end

nodeDifference21 = node(:,volSurface(:,2))-node(:,volSurface(:,1));
nodeDifference31 = node(:,volSurface(:,3))-node(:,volSurface(:,1));
volElementIdx = zeros(1,size(volSurface,1));
crossProduct = zeros(3, size(volSurface,1));
for i=1:size(volSurface,1)    
volElementIdx(i) = find(sum(ismember(elem,volSurface(i,:))) >= 3);
crossProduct(:,i) = cross(nodeDifference21(:,i),nodeDifference31(:,i));
crossProduct(:,i) = crossProduct(:,i)/norm(crossProduct(:,i));
end
surfElemC = elemCenter(node,elem(:,volElementIdx));
surfTriC = elemCenter(node,volSurface');
surfNormal = crossProduct .* repmat(sign(-sum(crossProduct .* (surfElemC-surfTriC))),3,1);
surfElements.elemCenters = surfTriC;
surfElements.elemIndices = volElementIdx; 
