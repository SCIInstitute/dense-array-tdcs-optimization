function surfaceNormalD = surfaceNormalInterpolation(headMesh,roi)
%Finds the surface normal direction for the elements in the ROI via
%interpolation. First finds the direction for the elements that are on the
%CSF-ROI or WM-ROI interface and then interpolates through the ROI. 

%Note: This function is terribly scripted, thus it is not guaranteed to be fast.

%%Reading inputs
elem = headMesh.cell;
node = headMesh.node;
field = headMesh.field;
clear headMesh;
if size(roi,2) ~= size(elem,2)
    roi = roi';
end

roiElements = elem(:,roi==1);
surfaceNormalD = zeros(3,size(roiElements,2));
roiSurface = volface(roiElements');

% Defining the surfaces from which the interpolation will be initiated 
roiCSFtouching = find(sum(ismember(roiSurface',elem(:,field==3))) >= 3);
roiWMtouching = find(sum(ismember(roiSurface',elem(:,field==5))) >= 3);

% csfSurface.face = roiSurface(roiCSFtouching,:)';
% csfSurface.node = node;
% save('csfSurfaceInterpolation','csfSurface');
% wmSurface.face = roiSurface(roiWMtouching,:)';
% wmSurface.node = node;
% save('wmSurfaceInterpolation','wmSurface');

%Defining the direction for the surface elements
[rSONfield,rSONelements] = surfOuterNormal(roiElements,node,roiSurface);
rSONfield = rSONfield(:,[roiCSFtouching roiWMtouching]);
rSONfield(:,1:numel(roiCSFtouching)) = rSONfield(:,1:numel(roiCSFtouching))*-1;
rSONindices = rSONelements.elemIndices;
rSONindices = rSONindices([roiCSFtouching roiWMtouching]);

% Defining neighborhood relationships
facenb = faceneighbors(roiElements');
doneIndices = rSONindices;

%interpolation
surfaceNormalD(:,doneIndices) = rSONfield;
while (numel(doneIndices) < size(roiElements,2))
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
    ordering = randperm(size(roiElements,2));
    for k = 1:numel(ordering)
        if ~ismember(ordering(k), rSONindices)
            averagingIndices = facenb(ordering(k),:);
            averagingIndices(averagingIndices==0) = [];
            surfaceNormalD(:,ordering(k)) = sum(surfaceNormalD(:,averagingIndices),2);
            surfaceNormalD(:,ordering(k)) = surfaceNormalD(:,ordering(k))/norm(surfaceNormalD(:,ordering(k)));
        end
    end
end

 