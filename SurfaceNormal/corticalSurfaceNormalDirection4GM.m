function sfn = corticalSurfaceNormalDirection4GM(headMesh,mapping,opts)
% Finds the cortical surface normal desired modulation direction for the
% gray matter elements.
%
%
% Synopsis: [surfNormD] = corticalSurfaceNormalDirection4GM(headMesh, opts )
%
%
% Input:    headMesh  = tetrahedral head mesh.
%           opts       = options. If 0, it interpolates from the CSF-GM
%           interface. If 1, it interpolates from both CSF-GM and GM-WM
%           directions.
%
%
% Output:   surfNormD =   cortical surface normal direction for ROI


%%Reading inputs

graymatterRow = headMesh.cell(:,headMesh.field==4)';
wm = headMesh.cell(:,headMesh.field==5);
csf = headMesh.cell(:,headMesh.field==3);
node = headMesh.node;
clear headMesh;
sfn = zeros(3,size(graymatterRow,1));


faceNeighborhood4GM = faceneighbors(graymatterRow);

gmSurface = volface(graymatterRow)';

% Defining the surfaces from which the interpolation will be initiated
GM_CSF_Interface = find(sum(ismember(gmSurface,csf)) >= 3);
if opts == 1
    GM_WM_Interface = find(sum(ismember(gmSurface,wm)) >= 3);
else
    GM_WM_Interface = [];
end

% csfSurface.face = roiSurface(roiCSFtouching,:)';
% csfSurface.node = node;
% save('csfSurfaceInterpolation','csfSurface');
% wmSurface.face = roiSurface(roiWMtouching,:)';
% wmSurface.node = node;
% save('wmSurfaceInterpolation','wmSurface');

% %Defining the direction for the surface elements
[rSONfield,~] = surfOuterNormal(graymatterRow',node,gmSurface,mapping);
rSONfield = rSONfield(:,[GM_CSF_Interface GM_WM_Interface]);
rSONfield(:,1:numel(GM_CSF_Interface)) = rSONfield(:,1:numel(GM_CSF_Interface))*-1;
%

doneIndices = mapping([GM_CSF_Interface GM_WM_Interface]);
%interpolation
sfn(:,doneIndices) = rSONfield;
while (numel(doneIndices) < size(graymatterRow,1))
    % for i=1:50
    %      sfn = normc(sfn(:,faceNeighborhood4GM(:,1))+...
    %                             sfn(:,faceNeighborhood4GM(:,2))+...
    %                             sfn(:,faceNeighborhood4GM(:,3))+...
    %                             sfn(:,faceNeighborhood4GM(:,4)));
    %      sfn(:,cidx) = surfaceNormal;
    newSetIndices = unique(faceNeighborhood4GM(doneIndices,:));
    deleteSetIndices = [find(newSetIndices ==0); find(ismember(newSetIndices,doneIndices) == 1)];
    newSetIndices(deleteSetIndices) = [];
    for i=1:numel(newSetIndices)
        averagingIndices = faceNeighborhood4GM(newSetIndices(i),:);
        averagingIndices(averagingIndices==0) = [];
        sfn(:,newSetIndices(i)) = sum(sfn(:,averagingIndices),2);
    end
    sfn = normc(sfn);
    doneIndices = [doneIndices newSetIndices'];
end
%Filling 0's in a matrix

for i = 1:5
    ordering = randperm(size(graymatterRow,1));
    for k = 1:numel(ordering)
        if ~ismember(ordering(k), rSONindices)
            averagingIndices = faceNeighborhood4GM(ordering(k),:);
            averagingIndices(averagingIndices==0) = [];
            sfn(:,ordering(k)) = normc(sum(sfn(:,averagingIndices),2));
        end
    end
end



