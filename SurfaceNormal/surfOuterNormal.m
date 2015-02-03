function [surfNormal, surfElements] = surfOuterNormal(elem,node,volSurface,mapping)
%Finds the outer surface normal on the surface of a mesh.
%Also finds surface triangle element centers and indices to be used as
%nodes for surface normal direction field.
%elem: 4xM, node: 3xN,

%if surface elements are not defined, find them using iso2mesh toolbox
%volface func
if nargin < 4
    if nargin < 3
    volSurface = volface(elem')';
    end
    permutations = perms(1:4);
    permutations(:,end) = [];
    mapping = zeros(size(volSurface,2),1);
    for i =1:size(permuations,1)
        [~,lcb] = ismember(volSurface',elem(permutations(i,:),:)','rows');
        mapping = mapping + lcb;
    end
end

if(numel(mapping) ~= size(volSurface,2))
    error('mapping and surface sizes dont match');
end

%We find the normal direction by cross product:
%Normal of a triangle is the cross of its two edges, we will use edges
%(vectors) 12 and 13.
u = node(:,volSurface(2,:))-node(:,volSurface(1,:));
v = node(:,volSurface(3,:))-node(:,volSurface(1,:));

%USE intersect function
%volElementIdx = ceil(find(ismember(faces',volSurface,'rows'))/4);
%calculate cross product analytically
crossProduct = zeros(3, size(volSurface,2));
crossProduct(1,:) = u(2,:).*v(3,:) - u(3,:).*v(2,:);
crossProduct(2,:) = u(3,:).*v(1,:) - u(1,:).*v(3,:);
crossProduct(3,:) = u(1,:).*v(2,:) - u(2,:).*v(1,:);
crossProduct = normc(crossProduct);

%mapping from triangles to elements, combined with sign multiplication to
%make sure the direction is towards outside of the gm.

%We will need to make sure the direction is outwards by following script.
surfElemC = elemCenter(node,elem(:,mapping)); %tetra centers.
surfTriC = elemCenter(node,volSurface); %triangle centers.

%Convert the direction if it is inwards.
surfNormal = crossProduct .* repmat(sign(sum(crossProduct .* (surfTriC-surfElemC))),3,1);
surfElements.triCenters = surfTriC;
surfElements.elemIndices = mapping;
surfElements.tetCenters = surfElemC;
