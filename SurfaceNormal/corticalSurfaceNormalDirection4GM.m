%% Finds the cortical surface normal direction field for gray matter elements.
%
%
% Synopsis: [surfNormD] = corticalSurfaceNormalDirection4GM(headMesh, opts )
%
%
% Input:    headMesh  = tetrahedral head mesh.
%           opts      = If 0, it interpolates from the CSF-GM interface. 
%                       If 1, it interpolates from both CSF-GM and GM-WM.
%
%
% Output:   surfNormD = cortical surface normal direction for GM

% Notes:    1. The assumption is that gray matter is labeled '4', CSF '3', 
%              and white matter '5' in the mesh field.
%           2. Needs iso2mesh package for finding the face neighbors and
%           surface of a volume mesh.

function surfNormD = corticalSurfaceNormalDirection4GM(headMesh,opts)
%% Reading inputs

grayMatter = headMesh.cell(:,headMesh.field==4);

a = finddisconnsurf(grayMatter');
if size(a,2) > 1
    for i = 1:size(a,1)
        ai(i) = size(a{i},1);
    end
    grayMatter = a{find(ai==max(ai))}';
end
whiteMatter = headMesh.cell(:,headMesh.field==5);
csf = headMesh.cell(:,headMesh.field==3);
node = headMesh.node;
clear headMesh;


%% Defining the surfaces from which the interpolation will be initiated

gmSurface = volface(grayMatter')';
gmCsfInterface = false(size(gmSurface,2),1);
gmWmInterface = false(size(gmSurface,2),1);
gmVolume2SurfaceMapping = zeros(size(gmSurface,2),1);

permutations = perms(1:4);
permutations(:,end) = [];

for i =1:size(permutations,1)
    [liaCsf] = ismember(gmSurface', csf(permutations(i,:),:)','rows');
    gmCsfInterface = gmCsfInterface | liaCsf;
    
    if opts == 1
    [liaWm] = ismember(gmSurface', whiteMatter(permutations(i,:),:)','rows');
    gmWmInterface = gmWmInterface | liaWm;
    end
    
    [~,locGm] = ismember(gmSurface', grayMatter(permutations(i,:),:)','rows');
    gmVolume2SurfaceMapping = gmVolume2SurfaceMapping + locGm;
end
gmWmInterface(gmCsfInterface) = 0;


%% Defining surface normal direction for surface elements

u = node(:,gmSurface(2,:)) - node(:,gmSurface(1,:));
v = node(:,gmSurface(3,:)) - node(:,gmSurface(1,:));
crossProduct = zeros(3, size(gmSurface,2));
crossProduct(1,:) = (u(2,:) .* v(3,:)) - (u(3,:) .* v(2,:));
crossProduct(2,:) = (u(3,:) .* v(1,:)) - (u(1,:) .* v(3,:));
crossProduct(3,:) = (u(1,:) .* v(2,:)) - (u(2,:) .* v(1,:));
crossProduct = normc(crossProduct);

%the direction is outwards
surfTetraCenter = elemCenter(node,grayMatter(:,gmVolume2SurfaceMapping)); 
surfTriangleCenter = elemCenter(node,gmSurface);

surfaceOuterNormal = crossProduct .* repmat(sign(sum(crossProduct .* (surfTriangleCenter-surfTetraCenter))),3,1);

%surfaceElements.triangleCenter = surfTriangleCenter;
%surfaceElements.tetrahedraCenter = surfTetraCenter;

%% interpolation towards inner elements

surfNormD = zeros(3, size(grayMatter,2)+1);
faceNeighbors4GM = faceneighbors(grayMatter');
faceNeighbors4GM(end+1,:) = [0 0 0 0];
faceNeighbors4GM(faceNeighbors4GM==0) = size(grayMatter,2)+1;

%SOLVE THE UNIQUENESS PROBLEM.
outerGmCsfSurfaceNormal = surfaceOuterNormal(:,gmCsfInterface);
innerGmWmSurfaceNormal = surfaceOuterNormal(:,gmWmInterface) * -1;
surfNormD(:,gmVolume2SurfaceMapping(gmCsfInterface)) = outerGmCsfSurfaceNormal;
surfNormD(:,gmVolume2SurfaceMapping(gmWmInterface)) = innerGmWmSurfaceNormal;

surfaceIndices = gmVolume2SurfaceMapping(gmCsfInterface | gmWmInterface);
doneIndices = surfaceIndices;
while (numel(doneIndices) < size(grayMatter,2)+1)
    newSetIndices = setdiff(unique(faceNeighbors4GM(doneIndices,:)),doneIndices);
    avIdx = faceNeighbors4GM(newSetIndices,:);
    surfNormD(:,newSetIndices) = normc(surfNormD(:,avIdx(:,1)) + ...
        surfNormD(:,avIdx(:,2))+surfNormD(:,avIdx(:,3))+surfNormD(:,avIdx(:,4)));
    surfNormD = normc(surfNormD);
    doneIndices = [doneIndices; newSetIndices];
    disp(num2str(numel(doneIndices)));
end

%Filling 0's in a matrix
innerElements = setdiff(1:size(grayMatter,2), surfaceIndices);
for i = 1:100
    disp(num2str(i));
    ordering = randperm(numel(innerElements));
    innerElements(ordering) = innerElements;
    blockSize = floor(numel(innerElements)/100);
    for k = 1:100
        if k ==100
            idx = innerElements((k-1)*blockSize+1:end);
        else
            idx = innerElements((k-1)*blockSize+1:k*blockSize);
        end
        avIdx = faceNeighbors4GM(idx,:);
        surfNormD(:,idx) = normc(surfNormD(:,avIdx(:,1))+surfNormD(:,avIdx(:,2))+surfNormD(:,avIdx(:,3))+...
            surfNormD(:,avIdx(:,4)));
    end
    surfNormD = normc(surfNormD);
end
surfNormD(:,end) = [];
end

function centers = elemCenter(node,elem)
%% finds the centers of the elements

[nD,~] = size(node);
[eD,M] = size(elem);

expandedCenters = zeros(1,nD*M);
for i=1:eD
    expandedCenters = expandedCenters + reshape(node(:,elem(i,:)),1,[]);
end

centers = reshape(expandedCenters/eD,nD,[]);
end
