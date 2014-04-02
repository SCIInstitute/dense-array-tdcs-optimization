function dir = createNormalDirection(headMesh,roi)

M = size(headMesh.cell,2);
headMesh.cell = double(headMesh.cell); 
if(size(roi,2) ~= M)
	roi = roi';
end
dir = cell(1,size(roi,1));
for i = 1:size(roi,1)
  dir{i}{1} = surfaceNormalInterpolation(headMesh,roi(i,:),i);
  dir{i}{2} = surfaceNormalProjection(headMesh,roi(i,:),i);
end
