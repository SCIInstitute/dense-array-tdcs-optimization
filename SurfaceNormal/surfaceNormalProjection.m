function surfaceNormalD = surfaceNormalProjection(headMesh,roi)
%%FINDS SURFACE NORMAL BASED ON THE PROJECTION TO THE SURFACE
%Written by: Seyhmus Guler, 3/8/14
% INPUTS
% OUTPUTS
%%surfaceNormalD: The surface normal desired modulation direction
% Size: 3 x (#roi Elements)

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


%%Defining the surfaces to which the elements will be projected
%We find CSF and WM touching of the ROI.
roiCSFtouching = find(sum(ismember(roiSurface',elem(:,field==3))) >= 3);
roiWMtouching = find(sum(ismember(roiSurface',elem(:,field==5))) >= 3);

% Save surface elements if needed. 
% csfSurface.face = roiSurface(roiCSFtouching,:)';
% csfSurface.node = node;
% save('csfSurfaceProjection','csfSurface');
% wmSurface.face = roiSurface(roiWMtouching,:)';
% wmSurface.node = node;
% save('wmSurfaceProjection','wmSurface');

%%Define surface normals for the surface elements
[rSONfield,rSONelements] = surfOuterNormal(roiElements,node,roiSurface);
rSONcenters = rSONelements.elemCenters;
%Keep only normals that are on the CSF or WM boundary.
rSONfield = rSONfield(:,[roiCSFtouching roiWMtouching]);
%Convert the direction for the CSF-ROI interface.
rSONfield(:,1:numel(roiCSFtouching)) = rSONfield(:,1:numel(roiCSFtouching))*-1;
rSONcenters = rSONcenters(:,[roiCSFtouching roiWMtouching]);

%%Finding projection to the surfaces and assigning surface normals
elementCenters = elemCenter(node,roiElements);
idx = knnsearch(rSONcenters',elementCenters');

%For each surface element, find the closest elements in the ROI and assign
%them with the same normal direction. 
for i=1:size(rSONcenters,2)
    surfaceNormalD(:,idx == i) = repmat(rSONfield(:,i),1,nnz(idx==i));
end