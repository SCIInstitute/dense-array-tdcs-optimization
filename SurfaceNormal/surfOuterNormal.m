function [surfNormal, surfElements] = surfOuterNormal(elem,node,volSurface)
%Finds the outer surface normal on the surface of a mesh.
%Also finds surface triangle element centers and indices to be used as
%nodes for surface normal direction field.
%elem: 4xM, node: 3xN

%if surface elements are not defined, find them using iso2mesh toolbox
%volface func
if nargin < 3
    volSurface = volface(elem');
end

%We find the normal direction by cross product:
%Normal of a triangle is the cross of its two edges, we will use edges
%(vectors) 12 and 13.
nodeDifference21 = node(:,volSurface(:,2))-node(:,volSurface(:,1));
nodeDifference31 = node(:,volSurface(:,3))-node(:,volSurface(:,1));
volElementIdx = zeros(1, size(volSurface,1));
crossProduct = zeros(3, size(volSurface,1));

%For each surface triangle find the normal direction via for loop. There
%may be computationally more efficient ways.
for i=1:size(volSurface,1)
    %find which cell that triangle belongs to
    volElementIdx(i) = find(sum(ismember(elem,volSurface(i,:))) >= 3);
    %cross producto to find the normal direction and normalization.
    crossProduct(:,i) = cross(nodeDifference21(:,i),nodeDifference31(:,i));
    crossProduct(:,i) = crossProduct(:,i)/norm(crossProduct(:,i));
end

%We will need to make sure the direction is outwards by following script.
surfElemC = elemCenter(node,elem(:,volElementIdx)); %tetra centers.
surfTriC = elemCenter(node,volSurface'); %triangle centers.

%Convert the direction if it is inwards. 
surfNormal = crossProduct .* repmat(sign(sum(crossProduct .* (surfTriC-surfElemC))),3,1);
surfElements.elemCenters = surfTriC;
surfElements.elemIndices = volElementIdx;
