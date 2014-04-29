function [curDensFromPot,elemVolumes] = cellCurrDensFromNodePotentials(tetrahedralMesh,conductivity)
%FINDS THE MAPPING BETWEEN THE POTENTIAL AT THE MESH NODES TO CURRENT
%DENSITY OF EACH ELEMENT OF THE MESH.  
%
%Written by: Seyhmus Guler,
%Last edit: 2/14/14 by Guler,S.
%
%Works for linear tetrahedral mesh and isotropic conductivities.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
    %tetrahedralMesh: tetrahedral mesh. A structure with following variables:
        %'cell','node','field'. 
        %cell: (4 x #elements), each column containing the indices of
            %nodes of corresponding element.
        %node: (3 x #nodes), each column containing the position of
            %corresponding node. 
        %field: (1 x #elements), each column containing the material 
            %label of the corresponding element. 
    %conductivity: the volume conductor.
        % It is a vector of size 1 x #elements. (Isotropic conductivity)
%OUTPUTS:
    %curDensFromPot: the mapping from node potentials to current density 
        % size: (3 #elements) x (#nodes)
    %elemVolumes: the volumes of the elements. size: 1 x #elements    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;    
elem = double(tetrahedralMesh.cell); %tetrahedral element list
node = tetrahedralMesh.node; %node positions
%field = tetrahedralMesh.field; %Field variable not needed for now
clear tetrahedralMesh;

if size(conductivity,2) ~= size(elem,2) 
    error('mismatch between element and volume conductor sizes');
end

fprintf('Calculating the mapping matrix from node potentials to cell current densities...\n');
M = size(elem,2); %number of elements
N = size(node,2); %number of nodes
Jtemp = sparse(3*N,M);
elemVolumes = zeros(1,M);

for i =1:M %for each element in the mesh
    A = inv([ones(1,4); node(:,elem(:,i))])'; %inverse the position matrix
    elemVolumes(i) = abs(1/6/det(A)); %volume of element in terms of node potentials
    Jtemp(:,i) = sparse([elem(:,i); elem(:,i)+N; elem(:,i)+2*N],1,reshape((-conductivity(i)*A(2:4,:))',12,1),3*N,1);
end
curDensFromPot = reshape(Jtemp,N,[])';
fprintf('%s%f%s\n','The mapping matrix from potential to current density and element volumes are calculated in ',toc,' seconds.\n');



        