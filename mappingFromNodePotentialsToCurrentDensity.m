function [curDensFromPot,elemVolume] = mappingFromNodePotentialsToCurrentDensity(mesh,conductivity)
%finds the mapping from potentials at the nodes of each element to the current
%density of each element. 

%Works only for linear tetrahedra mesh for now.
%it would be faster if matlabs parallel toolbox is used here. Make sure
%that matlabpool is open before calling this function

%inputs:
    %tetMesh: tetrahedral mesh. Should have variables
        %'cell','node','field'
    %lut: the conductivities of each material.
        %It is a cell array of size 1 x #different materials of mesh.
        %Conductivity of each material is either 3x3 tensor or a scalar.
%outputs:
    %curDensFromPot: the mapping from node potentials to current density matrix 
        % size: (3 #elements) x (#nodes)
    %elemVolume: the volumes of the elements. size: 1 x #elements
tic;    
elem = mesh.cell;
node = mesh.node;
field = mesh.field;
clear tetMesh;

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
if numel(unique(field)) ~= numel(conductivity)
    error('mismatch between the number of materials in the mesh and conductivity table');
end

fprintf('Calculating the transfer matrix from potentials to current densities...\n');
M = size(elem,2); %number of elements
N = size(node,2); %number of nodes

nL = 2e5;
elemVolume = zeros(1,M);
expandedNode = reshape(node(:,elem),12,[]);
expandedElem = double([elem; elem+N; elem+2*N]);
Jtemp = sparse(3*N,M);
K = floor(M/nL);
remain = M - K*nL;

temUtoJ = spalloc(3*N,nL,12*nL);
temUtoJremain = sparse(3*N,remain);
tempVol = zeros(1,nL);
tempVolRemain = zeros(1,remain);
matlabpool;
for k=1:K+1
    if k <= K
        idx = 1+(k-1)*nL:k*nL;
        tempElem = expandedElem(:,idx);
        tempField = field(idx);
        tempNode = expandedNode(:,idx);
        parfor i = 1:nL
            A = inv([ones(1,4); reshape(tempNode(:,i),3,4)])';
            tempVol(i) = abs(1/6/det(A));
            temUtoJ(:,i) = sparse(tempElem(:,i),1,reshape((-conductivity{tempField(i)}*A(2:4,:))',12,1),3*N,1);
        end
        Jtemp(:,idx) = temUtoJ;
        elemVolume(idx) = tempVol;
    else
        idx = 1+(k-1)*nL:M;
        tempElem = expandedElem(:,idx);
        tempField = field(idx);
        tempNode = expandedNode(:,idx);
        parfor i = 1:remain
            A = inv([ones(1,4); reshape(tempNode(:,i),3,4)])';
            tempVolRemain(i) = abs(1/6/det(A));
            temUtoJremain(:,i) = sparse(tempElem(:,i),1,reshape((-conductivity{tempField(i)}*A(2:4,:))',12,1),3*N,1);
        end
        Jtemp(:,1+(k-1)*nL:M) = temUtoJremain;
        elemVolume(idx) = tempVolRemain;
    end
    disp([num2str(k) '/' num2str(K+1)])
end
matlabpool close;
clear temUtoJ temUtoJremain eelem enode K remain tempVol tempVolRemain;
curDensFromPot = reshape(Jtemp,N,[])';
clear Jtemp;
fprintf('%s%f%s\n','The mapping matrix from potential to current density and element volumes are calculated in ',toc,' seconds');



        