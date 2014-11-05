function sfn = corticalSurfaceNormalDirection(headMesh,roi)
% Finds the cortical surface normal desired modulation direction for the 
% ROI elements. 
%
%
% Synopsis: [surfNormD] = corticalSurfaceNormalDirection(headMesh, roi )
%
%
% Input:    headMesh  = tetrahedral head mesh.
%           roi       = logical vector showing ROI elements.
%
%
% Output:   surfNormD =   cortical surface normal direction for ROI 

% Notes:    1. The function uses the elements in the ROI vicinity
%              in order to create a reasonable cortical surface normal.
%              It then interpolates the calculated cortical surface normal
%              towards the interior elements and returns the direction
%              found for the ROI elements.

%%Reading inputs
elem = headMesh.cell;
node = headMesh.node;
field = headMesh.field;
clear headMesh;
if size(roi,2) ~= size(elem,2)
    roi = roi';
end

%find the elements around ROI.
roiElements = elem(:,roi==1);
graymatterRow = elem(:,field==4)';
faceNeighbor = faceneighbors(graymatterRow);
roiInGrayMatter = ismember(graymatterRow,roiElements','rows');

surrElements = roiElements';
for i = 1:20
    surrElements = unique([surrElements; graymatterRow(nonzeros(faceNeighbor(ismember(graymatterRow,surrElements,'rows')==1,:)),:)],'rows');
end
surrElements = surrElements';


surfaceNormalD = zeros(3,size(surrElements,2));
surrSurface = volface(surrElements');

% Defining the surfaces from which the interpolation will be initiated 
surrCSFtouching = find(sum(ismember(surrSurface',elem(:,field==3))) >= 3);
surrWMtouching = find(sum(ismember(surrSurface',elem(:,field==5))) >= 3);

% csfSurface.face = roiSurface(roiCSFtouching,:)';
% csfSurface.node = node;
% save('csfSurfaceInterpolation','csfSurface');
% wmSurface.face = roiSurface(roiWMtouching,:)';
% wmSurface.node = node;
% save('wmSurfaceInterpolation','wmSurface');

%Defining the direction for the surface elements
[rSONfield,rSONelements] = surfOuterNormal(surrElements,node,surrSurface);
rSONfield = rSONfield(:,[surrCSFtouching surrWMtouching]);
rSONfield(:,1:numel(surrCSFtouching)) = rSONfield(:,1:numel(surrCSFtouching))*-1;
rSONindices = rSONelements.elemIndices;
rSONindices = rSONindices([surrCSFtouching surrWMtouching]);

% Defining neighborhood relationships
facenb = faceneighbors(surrElements');
doneIndices = rSONindices;

%interpolation
surfaceNormalD(:,doneIndices) = rSONfield;
while (numel(doneIndices) < size(surrElements,2))
    newSetIndices = unique(facenb(doneIndices,:));
    deleteSetIndices = [find(newSetIndices ==0); find(ismember(newSetIndices,doneIndices) == 1)];
    newSetIndices(deleteSetIndices) = [];
    for i=1:numel(newSetIndices)
        averagingIndices = facenb(newSetIndices(i),:);
        averagingIndices(averagingIndices==0) = [];
        surfaceNormalD(:,newSetIndices(i)) = sum(surfaceNormalD(:,averagingIndices),2);
        surfaceNormalD(:,newSetIndices(i)) = surfaceNormalD(:,newSetIndices(i)) / norm(surfaceNormalD(:,newSetIndices(i)));
    end
    doneIndices = [doneIndices newSetIndices'];
end

%10 random permutations to minimize the variance.
for i = 1:10
    ordering = randperm(size(surrElements,2));
    for k = 1:numel(ordering)
        if ~ismember(ordering(k), rSONindices)
            averagingIndices = facenb(ordering(k),:);
            averagingIndices(averagingIndices==0) = [];
            surfaceNormalD(:,ordering(k)) = sum(surfaceNormalD(:,averagingIndices),2);
            surfaceNormalD(:,ordering(k)) = surfaceNormalD(:,ordering(k))/norm(surfaceNormalD(:,ordering(k)));
        end
    end
end
[lia,locb] = ismember(surrElements',roiElements','rows');
sfn.surfaceNormalD = surfaceNormalD;
if(isequal(surrElements(:,lia==1),roiElements(:,locb(lia==1))))
    ROIinterSurfD = surfaceNormalD(:,lia==1);
    sfn.ROIsurfaceNormalD = ROIinterSurfD(:,locb(lia==1));
end
sfn.ROIcells = roiElements;
sfn.surrCells = surrElements;



 