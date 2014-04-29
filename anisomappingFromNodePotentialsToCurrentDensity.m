function [curDensFromPot,elemVolumes] = anisomappingFromNodePotentialsToCurrentDensity(mesh,conductivity)
%FINDS THE MAPPING BETWEEN THE POTENTIAL AT THE MESH NODES TO CURRENT
%DENSITY OF EACH ELEMENT OF THE MESH.  
%
%Written by: Seyhmus Guler, Revisited: Moritz Dannhauer
%Last edit: 10/7/13 by Guler,S.
%
%Works only for linear tetrahedra mesh for now.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
    %mesh: tetrahedral mesh. A structure with following variables:
        %'cell','node','field'. 
        %mesh.cell: (4 x #elements), each column containing the indices of
            %nodes of corresponding element.
        %mesh.node: (3 x #nodes), each column containing the position of
            %corresponding node. 
        %mesh.field:(1 x #elements), each column containing the material 
            %label of the corresponding element. 
    %conductivity: the volume conductor.
        %It is a matrix of size  #elements x 9.
        % 3x3 conductivity tensor of each element is reshaped to 1x9 vector 
        % and put in corresponding row in the conductivity matrix.
%OUTPUTS:
    %curDensFromPot: the mapping from node potentials to current density 
        % size: (3 #elements) x (#nodes)
    %elemVolumes: the volumes of the elements. size: 1 x #elements    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;    
elem = mesh.cell;
node = mesh.node;
field = mesh.field;
clear mesh;

if size(elem,1) ~= 4
    elem = elem';
end
if size(node,1) ~= 3
    node = node';
end
if size(field,1) ~= 1
    field = field';
end
if size(conductivity,2) ~= 9
    conductivity = conductivity';
end

if size(field,2) ~= size(elem,2)
    error('mismatch between element and field sizes');
end
if size(conductivity,1) ~= size(elem,2) 
    error('mismatch between element and volume conductor sizes');
end

fprintf('Calculating the transfer matrix from potentials to current densities...\n');
M = size(elem,2); %number of elements
N = size(node,2); %number of nodes

nL = 1e5;
elemVolumes = zeros(1,M);
expandedNode = reshape(node(:,elem),12,[]);
expandedElem = double([elem; elem+N; elem+2*N]);
Jtemp = sparse(3*N,M);
K = floor(M/nL);
remain = M - K*nL;

temUtoJ = spalloc(3*N,nL,12*nL);
temUtoJremain = sparse(3*N,remain);
tempVol = zeros(1,nL);
tempVolRemain = zeros(1,remain);
matlabpool close force
matlabpool open;
for k=1: K+1
    if k <= K
        idx = 1+(k-1)*nL:k*nL;
        tempElem = expandedElem(:,idx);
        tempCond = conductivity(idx,:);
        tempNode = expandedNode(:,idx);
        parfor i = 1:nL
            A = inv([ones(1,4); reshape(tempNode(:,i),3,4)])';
            tempVol(i) = abs(1/6/det(A));
            temUtoJ(:,i) = sparse(tempElem(:,i),1,reshape((-reshape(tempCond(i,:),3,3)*A(2:4,:))',12,1),3*N,1);
        end
        Jtemp(:,idx) = temUtoJ;
        elemVolumes(idx) = tempVol;
    else
        idx = 1+(k-1)*nL:M;
        tempElem = expandedElem(:,idx);
        tempCond = conductivity(idx,:);
        tempNode = expandedNode(:,idx);
        parfor i = 1:remain
            A = inv([ones(1,4); reshape(tempNode(:,i),3,4)])';
            tempVolRemain(i) = abs(1/6/det(A));
            temUtoJremain(:,i) = sparse(tempElem(:,i),1,reshape((-reshape(tempCond(i,:),3,3)*A(2:4,:))',12,1),3*N,1);
        end
        Jtemp(:,1+(k-1)*nL:M) = temUtoJremain;
        elemVolumes(idx) = tempVolRemain;
    end
    disp([num2str(k) '/' num2str(K+1)])
end
matlabpool close;
clear temUtoJ temUtoJremain eelem enode K remain tempVol tempVolRemain;
curDensFromPot = reshape(Jtemp,N,[])';
clear Jtemp;
fprintf('%s%f%s\n','The mapping matrix from potential to current density and element volumes are calculated in ',toc,' seconds');



        