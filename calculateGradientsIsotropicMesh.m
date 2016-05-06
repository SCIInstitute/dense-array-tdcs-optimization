function [sigmaGradV, cellVolume] = calculateGradientsIsotropicMesh(headMesh, sigmaLookupTable)
%Finds the matrix that links the potential (at mesh nodes) to current
%density (at mesh cells). It also provides the volumes of mesh cells
%
%
% Synopsis: [G,V] = calculateGradientsIsotropicMesh(headMesh, sigmaLookUp )
%
%
% Input:    headMesh    = tetrahedral head mesh. Contains three variables:
%               .cell   = tetrahedra definition matrix. Size: 4 x #cells
%               .node   = node coordinates.
%               .field  = tissue labels
%           sigmaLookup = conductivity look up table
%
%
% Output:   G           = mapping matrix from potential to current density
%           V           = cell volumes
%
%

% Notes:    1.  This function is for the isotropic meshes. If you would
%               like to find the mapping for anisotropic mesh, please check
%               the function anisomappingFromNodePotentialsToCurrentDensity

tic;
elem = headMesh.cell;
node = headMesh.node;
field = headMesh.field;
if numel(sigmaLookupTable) == numel(unique(headMesh.field))
sigma = sigmaLookupTable(field);
elseif numel(sigmaLookupTable) == numel(headMesh.field);
    sigma = sigmaLookupTable';
end
clear headMesh;

if size(elem,1) ~= 4
    elem = elem';
end
if size(node,1) ~= 3
    node = node';
end
if size(field,1) ~= 1
    field = field';
end

if size(field,2) ~= size(elem,2)
    error('mismatch between element and field sizes');
end

fprintf('Calculating the transfer matrix from potentials to current densities...\n');
M = size(elem,2); %number of elements
N = size(node,2); %number of nodes


A = reshape(node(:,elem),12,[]);
clear node;

detA =   A(4,:) .* (A(8,:).*A(12,:) - A(9,:) .* A(11,:)) - A(7,:) .* (A(5,:) .* A(12,:) - A(6,:) .*A (11,:)) + A(10,:) .* (A(5,:) .* A(9,:) - A(6,:) .* A(8,:))...
       - A(1,:) .* (A(8,:).*A(12,:) - A(9,:) .* A(11,:)) + A(7,:) .* (A(2,:) .* A(12,:) - A(3,:) .*A (11,:)) - A(10,:) .* (A(2,:) .* A(9,:) - A(3,:) .* A(8,:))...
       + A(1,:) .* (A(5,:).*A(12,:) - A(6,:) .* A(11,:)) - A(4,:) .* (A(2,:) .* A(12,:) - A(3,:) .*A (11,:)) + A(10,:) .* (A(2,:) .* A(6,:) - A(3,:) .* A(5,:))...
       - A(1,:) .* (A(5,:).*A(9,:)  - A(6,:) .* A(8,:))  + A(4,:) .* (A(2,:) .* A(9,:)  - A(3,:) .*A (8,:))  - A(7,:)  .* (A(2,:) .* A(6,:) - A(3,:) .* A(5,:));

cellVolume = detA/6;

detA = -sigma ./ detA;
B = zeros(12,M);
B(1,:)  = ( - (A(8,:) .* A(12,:) - A(9,:) .* A(11,:)) + (A(5,:) .* A(12,:) - A(6,:) .* A(11,:)) - (A(5,:) .* A(9,:) - A(6,:) .* A(8,:)) ) .* detA;
B(2,:)  = (   (A(8,:) .* A(12,:) - A(9,:) .* A(11,:)) - (A(2,:) .* A(12,:) - A(3,:) .* A(11,:)) + (A(2,:) .* A(9,:) - A(3,:) .* A(8,:)) ) .* detA;
B(3,:)  = ( - (A(5,:) .* A(12,:) - A(6,:) .* A(11,:)) + (A(2,:) .* A(12,:) - A(3,:) .* A(11,:)) - (A(2,:) .* A(6,:) - A(3,:) .* A(5,:)) ) .* detA;
B(4,:)  = (   (A(5,:) .* A(9,:)  - A(6,:) .* A(8,:))  - (A(2,:) .* A(9,:)  - A(3,:) .* A(8,:))  + (A(2,:) .* A(6,:) - A(3,:) .* A(5,:)) ) .* detA;
B(5,:)  = (   (A(7,:) .* A(12,:) - A(9,:) .* A(10,:)) - (A(4,:) .* A(12,:) - A(6,:) .* A(10,:)) + (A(4,:) .* A(9,:) - A(6,:) .* A(7,:)) ) .* detA;
B(6,:)  = ( - (A(7,:) .* A(12,:) - A(9,:) .* A(10,:)) + (A(1,:) .* A(12,:) - A(3,:) .* A(10,:)) - (A(1,:) .* A(9,:) - A(3,:) .* A(7,:)) ) .* detA;
B(7,:)  = (   (A(4,:) .* A(12,:) - A(6,:) .* A(10,:)) - (A(1,:) .* A(12,:) - A(3,:) .* A(10,:)) + (A(1,:) .* A(6,:) - A(3,:) .* A(4,:)) ) .* detA;
B(8,:)  = ( - (A(4,:) .* A(9,:)  - A(6,:) .* A(7,:))  + (A(1,:) .* A(9,:)  - A(3,:) .* A(7,:))  - (A(1,:) .* A(6,:) - A(3,:) .* A(4,:)) ) .* detA;
B(9,:)  = ( - (A(7,:) .* A(11,:) - A(8,:) .* A(10,:)) + (A(4,:) .* A(11,:) - A(5,:) .* A(10,:)) - (A(4,:) .* A(8,:) - A(5,:) .* A(7,:)) ) .* detA;
B(10,:) = (   (A(7,:) .* A(11,:) - A(8,:) .* A(10,:)) - (A(1,:) .* A(11,:) - A(2,:) .* A(10,:)) + (A(1,:) .* A(8,:) - A(2,:) .* A(7,:)) ) .* detA;
B(11,:) = ( - (A(4,:) .* A(11,:) - A(5,:) .* A(10,:)) + (A(1,:) .* A(11,:) - A(2,:) .* A(10,:)) - (A(1,:) .* A(5,:) - A(2,:) .* A(4,:)) ) .* detA;
B(12,:) = (   (A(4,:) .* A(8,:)  - A(5,:) .* A(7,:))  - (A(1,:) .* A(8,:)  - A(2,:) .* A(7,:))  + (A(1,:) .* A(5,:) - A(2,:) .* A(4,:)) ) .* detA;

rowIdx = reshape(repmat(1:3*M,4,1),1,[]);
colIdx = double(reshape(repmat(elem,3,1),1,[]));
clear elem;
sigmaGradV = sparse(rowIdx,colIdx,B,3*M,N);

fprintf('%s%f%s\n','The mapping matrix from potential to current density and element volumes are calculated in ',toc,' seconds');

%%

%   //fidx : the index of x1 of node 1 for the corresponding element
%   // A =  1     enode[fidx]     enode[fidx+1]   enode[fidx+2]
%   //      1     enode[fidx+3]   enode[fidx+4]   enode[fidx+5]
%   //      1     enode[fidx+6]   enode[fidx+7]   enode[fidx+8]
%   //      1     enode[fidx+9]   enode[fidx+10]  enode[fidx+11];
%   //   =  a11   a12     a13     a14
%   //      a21   a22     a23     a24
%   //      a31   a32     a33     a34
%   //      a41   a42     a43     a44
%   //   =  1   x1     y1     z1
%   //      1   x2     y2     z2
%   //      1   x3     y3     z3
%   //      1   x4     y4     z4
%   //   enode(:,cidx) = [x1; y1; z1; x2; y2; z2; x3; y3; z3; x4; y4; z4];
%   //                   [1   2   3   4   5   6   7   8   9   10  11  12];
%   //   Then;
%   //  *   1   2   3
%   //  *   4   5   6
%   //  *   7   8   9
%   //  *   10  11  12
%
%   a12 <-> A(1,:)
%   a13 <-> A(2,:)
%   a14 <-> A(3,:)
%   a22 <-> A(4,:)
%   a23 <-> A(5,:)
%   a24 <-> A(6,:)
%   a32 <-> A(7,:)
%   a33 <-> A(8,:)
%   a34 <-> A(9,:)
%   a42 <-> A(10,:)
%   a43 <-> A(11,:)
%   a44 <-> A(12,:)
%   // Note that a11 = a21 = a31 = a41 = 1 so we take out these variables from the formula
%   //lets calculate determinant, we may do this on V variable directly if we don't want to define new variable
%   double detA = a22 * (a33*a44 - a34*a43) - a32 * (a23*a44 - a24*a43) + a42 * (a23*a34 - a24*a33) //b11
%     	      - a12 * (a33*a44 - a34*a43) + a32 * (a13*a44 - a14*a43) - a42 * (a13*a34 - a14*a33) //b12
%   	      + a12 * (a23*a44 - a24*a43) - a22 * (a13*a44 - a14*a43) + a42 * (a13*a24 - a14*a23) //b13
%     	      - a12 * (a23*a34 - a24*a33) + a22 * (a13*a34 - a14*a33) - a32 * (a13*a24 - a14*a23);//b14
%   //detA = b11 + b12 + b13 + b14 because a*1 = 1; *=1,2,3,4
%   detA = 1/detA;
%

%
%   b21 = ( - (a33*a44 - a34*a43) + (a23*a44 - a24*a43) - (a23*a34 - a24*a33) ) * detA;
%   b22 = (   (a33*a44 - a34*a43) - (a13*a44 - a14*a43) + (a13*a34 - a14*a33) ) * detA;
%   b23 = ( - (a23*a44 - a24*a43) + (a13*a44 - a14*a43) - (a13*a24 - a14*a23) ) * detA;
%   b24 = (   (a23*a34 - a24*a33) - (a13*a34 - a14*a33) + (a13*a24 - a14*a23) ) * detA;
%   b31 = (   (a32*a44 - a34*a42) - (a22*a44 - a24*a42) + (a22*a34 - a24*a32) ) * detA;
%   b32 = ( - (a32*a44 - a34*a42) + (a12*a44 - a14*a42) - (a12*a34 - a14*a32) ) * detA;
%   b33 = (   (a22*a44 - a24*a42) - (a12*a44 - a14*a42) + (a12*a24 - a14*a22) ) * detA;
%   b34 = ( - (a22*a34 - a24*a32) + (a12*a34 - a14*a32) - (a12*a24 - a14*a22) ) * detA;
%   b41 = ( - (a32*a43 - a33*a42) + (a22*a43 - a23*a42) - (a22*a33 - a23*a32) ) * detA;
%   b42 = (   (a32*a43 - a33*a42) - (a12*a43 - a13*a42) + (a12*a33 - a13*a32) ) * detA;
%   b43 = ( - (a22*a43 - a23*a42) + (a12*a43 - a13*a42) - (a12*a23 - a13*a22) ) * detA;
%   b44 = (   (a22*a33 - a23*a32) - (a12*a33 - a13*a32) + (a12*a23 - a13*a22) ) * detA;
%
%

