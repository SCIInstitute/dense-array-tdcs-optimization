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
sfn = zeros(3,size(graymatterRow,1)+1);
faceNeighborhood4GM = faceneighbors(graymatterRow);
faceNeighborhood4GM(end+1,:) = [0 0 0 0];
faceNeighborhood4GM(faceNeighborhood4GM==0) = size(graymatterRow,1)+1;

gmSurface = volface(graymatterRow)';

% Defining the surfaces from which the interpolation will be initiated
GM_CSF_Interface = find(sum(ismember(gmSurface,csf)) >= 3);
if opts == 1
    GM_WM_Interface = find(sum(ismember(gmSurface,wm)) >= 3);
    GM_WM_Interface = setdiff(GM_WM_Interface,GM_CSF_Interface);
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
[rSONfield] = surfOuterNormal(graymatterRow',node,gmSurface,mapping);
rSONfield = rSONfield(:,[GM_CSF_Interface GM_WM_Interface]);
rSONfield(:,1:numel(GM_CSF_Interface)) = rSONfield(:,1:numel(GM_CSF_Interface))*-1;
%

doneIndices = mapping([GM_CSF_Interface GM_WM_Interface]);
%interpolation
sfn(:,doneIndices) = rSONfield;
while (numel(doneIndices) < size(graymatterRow,1)+1)
    % for i=1:50
    %      sfn = normc(sfn(:,faceNeighborhood4GM(:,1))+...
    %                             sfn(:,faceNeighborhood4GM(:,2))+...
    %                             sfn(:,faceNeighborhood4GM(:,3))+...
    %                             sfn(:,faceNeighborhood4GM(:,4)));
    %      sfn(:,cidx) = surfaceNormal;
    newSetIndices = setdiff(unique(faceNeighborhood4GM(doneIndices,:)),doneIndices);
    avIdx = faceNeighborhood4GM(newSetIndices,:);
    sfn(:,newSetIndices) = normc(sfn(:,avIdx(:,1))+sfn(:,avIdx(:,2))+sfn(:,avIdx(:,3))+...
        sfn(:,avIdx(:,4)));
    sfn = normc(sfn);
    doneIndices = [doneIndices; newSetIndices];
    disp(num2str(numel(doneIndices)));
end
%Filling 0's in a matrix
innerElements = setdiff(1:size(graymatterRow,1),[GM_CSF_Interface GM_WM_Interface]);
for i = 1:10
    ordering = randperm(numel(innerElements));
    innerElements(ordering) = innerElements;
    blockSize = floor(numel(innerElements)/100);
    for k = 1:100
        if k ==100
            idx = innerElements((k-1)*blockSize+1:end);
        else
            idx = innerElements((k-1)*blockSize+1:k*blockSize);
        end
            avIdx = faceNeighborhood4GM(idx,:);
            sfn(:,idx) = normc(sfn(:,avIdx(:,1))+sfn(:,avIdx(:,2))+sfn(:,avIdx(:,3))+...
        sfn(:,avIdx(:,4)));
    end
    sfn = normc(sfn);
end
sfn(:,end) = [];



