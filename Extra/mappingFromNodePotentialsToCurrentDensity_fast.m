function [curDensFromPot,elemVolume] = mappingFromNodePotentialsToCurrentDensity_fast(mesh,conductivity)
% THIS IS A REWRITTEN VERSION OF anisomappingFromNodePotentialsToCurrentDensity.m
% IT IS OPTIMIZED FOR COMPUTATIONAL SPEED AND IS KNOWN TO BE LESS ACCURATE
% AS ITS SLOW ORIGINAL SINCE MATRIX INVERSES ARE COMPUTED IN A NOT
% NUMERICALLY WAY
% WRITTEN BY MORITZ DANNHAUER, 04/17/13

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

if size(field,2) ~= size(elem,2)
    error('mismatch between element and field sizes');
end
if numel(field) ~= size(conductivity,3)
    error('mismatch between the number of materials in the mesh and conductivity table');
end

fprintf('Calculating the transfer matrix from potentials to current densities...\n');
M = size(elem,2); %number of elements
N = size(node,2); %number of nodes

expandedNode = reshape(node(:,elem),12,[]);

clear nodes;


nr_of_workers=12;
stepsize=floor(M/nr_of_workers);

for k=1:nr_of_workers
   if(k~=nr_of_workers)
    idx{k} = 1+(k-1)*stepsize:k*stepsize;
  else
    idx{k} = 1+(k-1)*stepsize:M;  
   end
  tempNode{k} = expandedNode(:,idx{k});
  A{k}=ones(4,4,length(idx{k}));
  B{k}=ones(4,3,length(idx{k}));
  inve{k}=zeros(4,4,length(idx{k}));
  elemVolume{k} = zeros(1,length(idx{k}));
  Jtemp{k} = [];
  tempCond{k} = -conductivity(:,:,idx{k});
  tmp= double([elem(:,idx{k}); elem(:,idx{k})+N; elem(:,idx{k})+2*N]);
  tmp=tmp(:);
  indexe{k} = tmp;
end
matlabpool close force
matlabpool open;
 
parfor k=1:nr_of_workers
  A{k}(2:4,1:4,:)=reshape(tempNode{k},3,4,length(idx{k}));
  inve{k}(1,1,:) = prod([ A{k}(2,2,:); A{k}(3,3,:) ; A{k}(4,4,:) ]) - prod([  A{k}(2,2,:)  ; A{k}(3,4,:) ; A{k}(4,3,:)]) - prod([  A{k}(3,2,:)  ; A{k}(2,3,:)  ; A{k}(4,4,:)]) + prod([  A{k}(3,2,:)  ; A{k}(2,4,:)  ; A{k}(4,3,:)]) + prod([  A{k}(4,2,:) ; A{k}(2,3,:)  ; A{k}(3,4,:) ])- prod([  A{k}(4,2,:) ; A{k}(2,4,:)  ; A{k}(3,3,:)]);
  inve{k}(1,2,:) = -prod([ A{k}(2,1,:); A{k}(3,3,:) ; A{k}(4,4,:)]) + prod([ A{k}(2,1,:)  ; A{k}(3,4,:) ; A{k}(4,3,:)]) + prod([ A{k}(3,1,:)  ; A{k}(2,3,:)  ; A{k}(4,4,:)]) - prod([  A{k}(3,1,:)  ; A{k}(2,4,:)  ; A{k}(4,3,:)]) - prod([  A{k}(4,1,:) ; A{k}(2,3,:)  ; A{k}(3,4,:)]) + prod([  A{k}(4,1,:) ; A{k}(2,4,:)  ; A{k}(3,3,:)]);
  inve{k}(1,3,:) = prod([ A{k}(2,1,:); A{k}(3,2,:) ; A{k}(4,4,:)]) - prod([ A{k}(2,1,:)  ; A{k}(3,4,:) ; A{k}(4,2,:)]) - prod([  A{k}(3,1,:)  ; A{k}(2,2,:) ; A{k}(4,4,:)]) + prod([  A{k}(3,1,:)  ; A{k}(2,4,:) ; A{k}(4,2,:)]) + prod([ A{k}(4,1,:) ; A{k}(2,2,:) ; A{k}(3,4,:)]) - prod([ A{k}(4,1,:) ; A{k}(2,4,:) ; A{k}(3,2,:)]);
  inve{k}(1,4,:) = -prod([ A{k}(2,1,:); A{k}(3,2,:) ; A{k}(4,3,:)]) + prod([ A{k}(2,1,:)  ; A{k}(3,3,:) ; A{k}(4,2,:)]) +prod([ A{k}(3,1,:)  ; A{k}(2,2,:) ; A{k}(4,3,:)]) - prod([  A{k}(3,1,:)  ; A{k}(2,3,:) ; A{k}(4,2,:)]) - prod([  A{k}(4,1,:) ; A{k}(2,2,:) ; A{k}(3,3,:)]) + prod([  A{k}(4,1,:) ; A{k}(2,3,:) ; A{k}(3,2,:)]);

  inve{k}(2,1,:) = -prod([  A{k}(3,3,:) ; A{k}(4,4,:)]) + prod([    A{k}(3,4,:) ; A{k}(4,3,:)]) + prod([  A{k}(3,2,:)  ;  A{k}(4,4,:)]) - prod([  A{k}(3,2,:)  ;   A{k}(4,3,:)]) - prod([  A{k}(4,2,:) ;   A{k}(3,4,:)]) + prod([  A{k}(4,2,:) ;   A{k}(3,3,:)]);
  inve{k}(2,2,:) = prod([  A{k}(3,3,:) ; A{k}(4,4,:)]) - prod([    A{k}(3,4,:) ; A{k}(4,3,:)]) - prod([ A{k}(3,1,:)  ;  A{k}(4,4,:)]) + prod([ A{k}(3,1,:)  ;   A{k}(4,3,:)]) + prod([ A{k}(4,1,:) ;   A{k}(3,4,:)]) - prod([ A{k}(4,1,:) ;   A{k}(3,3,:)]);
  inve{k}(2,3,:) = -prod([  A{k}(3,2,:) ; A{k}(4,4,:)]) + prod([   A{k}(3,4,:) ; A{k}(4,2,:)]) + prod([ A{k}(3,1,:)  ;   A{k}(4,4,:)]) - prod([ A{k}(3,1,:)  ;   A{k}(4,2,:)]) - prod([ A{k}(4,1,:) ;   A{k}(3,4,:)]) + prod([ A{k}(4,1,:) ;   A{k}(3,2,:)]);
  inve{k}(2,4,:) = prod([  A{k}(3,2,:) ; A{k}(4,3,:)]) - prod([    A{k}(3,3,:) ; A{k}(4,2,:)]) - prod([ A{k}(3,1,:)  ;   A{k}(4,3,:)]) + prod([ A{k}(3,1,:)  ;   A{k}(4,2,:)]) + prod([ A{k}(4,1,:) ;   A{k}(3,3,:)]) - prod([ A{k}(4,1,:) ;   A{k}(3,2,:)]);

  inve{k}(3,1,:) = prod([  A{k}(2,3,:) ; A{k}(4,4,:)]) - prod([    A{k}(2,4,:) ; A{k}(4,3,:)]) - prod([ A{k}(2,2,:)  ;  A{k}(4,4,:)]) + prod([ A{k}(2,2,:)  ;   A{k}(4,3,:)]) + prod([ A{k}(4,2,:) ;   A{k}(2,4,:)]) - prod([ A{k}(4,2,:) ;   A{k}(2,3,:)]);
  inve{k}(3,2,:) = -prod([ A{k}(2,3,:) ; A{k}(4,4,:)]) + prod([    A{k}(2,4,:) ; A{k}(4,3,:)]) + prod([ A{k}(2,1,:)  ;  A{k}(4,4,:)]) - prod([ A{k}(2,1,:)  ;   A{k}(4,3,:)]) - prod([ A{k}(4,1,:) ;   A{k}(2,4,:)]) + prod([ A{k}(4,1,:) ;   A{k}(2,3,:)]);
  inve{k}(3,3,:) = prod([  A{k}(2,2,:) ; A{k}(4,4,:)]) - prod([    A{k}(2,4,:) ; A{k}(4,2,:)]) - prod([ A{k}(2,1,:)  ;   A{k}(4,4,:)]) + prod([ A{k}(2,1,:)  ;   A{k}(4,2,:)]) + prod([ A{k}(4,1,:) ;   A{k}(2,4,:)]) - prod([ A{k}(4,1,:) ;   A{k}(2,2,:)]);
  inve{k}(4,3,:) = -prod([  A{k}(2,2,:) ; A{k}(3,4,:)]) + prod([   A{k}(2,4,:) ; A{k}(3,2,:)]) + prod([ A{k}(2,1,:) ;   A{k}(3,4,:)]) - prod([  A{k}(2,1,:) ;   A{k}(3,2,:)]) - prod([  A{k}(3,1,:) ;   A{k}(2,4,:)]) + prod([  A{k}(3,1,:) ;   A{k}(2,2,:)]);

  inve{k}(4,1,:) = -prod([ A{k}(2,3,:) ; A{k}(3,4,:)]) + prod([   A{k}(2,4,:) ; A{k}(3,3,:)]) + prod([ A{k}(2,2,:) ;   A{k}(3,4,:)]) - prod([ A{k}(2,2,:) ;   A{k}(3,3,:)]) - prod([ A{k}(3,2,:) ;   A{k}(2,4,:)]) + prod([ A{k}(3,2,:) ;   A{k}(2,3,:)]);
  inve{k}(4,2,:) = prod([  A{k}(2,3,:) ; A{k}(3,4,:)]) - prod([   A{k}(2,4,:) ; A{k}(3,3,:)]) - prod([ A{k}(2,1,:) ;   A{k}(3,4,:)]) + prod([  A{k}(2,1,:) ;   A{k}(3,3,:)]) + prod([  A{k}(3,1,:) ;   A{k}(2,4,:)]) - prod([ A{k}(3,1,:) ;   A{k}(2,3,:)]);
  inve{k}(3,4,:) = -prod([ A{k}(2,2,:) ; A{k}(4,3,:)]) + prod([    A{k}(2,3,:) ; A{k}(4,2,:)]) + prod([ A{k}(2,1,:)  ;   A{k}(4,3,:)]) - prod([ A{k}(2,1,:)  ;  A{k}(4,2,:)]) - prod([  A{k}(4,1,:) ;   A{k}(2,3,:)]) + prod([  A{k}(4,1,:) ;   A{k}(2,2,:)]);
  inve{k}(4,4,:) = prod([  A{k}(2,2,:) ; A{k}(3,3,:)]) - prod([   A{k}(2,3,:) ; A{k}(3,2,:)]) - prod([ A{k}(2,1,:) ;   A{k}(3,3,:)]) + prod([ A{k}(2,1,:) ;   A{k}(3,2,:)]) + prod([ A{k}(3,1,:) ;   A{k}(2,3,:)]) - prod([ A{k}(3,1,:) ;   A{k}(2,2,:)]);  
  det = inve{k}(1,1,:) +  inve{k}(1,2,:) +  inve{k}(1,3,:) +  inve{k}(1,4,:); 
  det = 1./det; 
  inve{k} = bsxfun(@times,inve{k},det);
  
  elemVolume{k} = abs(1/6/squeeze(det));
  B{k}(1,1,:)=prod([tempCond{k}(1,1,:) inve{k}(2,1,:)]) + prod([tempCond{k}(1,2,:) inve{k}(3,1,:)]) + prod([tempCond{k}(1,3,:) inve{k}(4,1,:)]); 
  B{k}(2,1,:)=prod([tempCond{k}(1,1,:) inve{k}(2,2,:)]) + prod([tempCond{k}(1,2,:) inve{k}(3,2,:)]) + prod([tempCond{k}(1,3,:) inve{k}(4,2,:)]);
  B{k}(3,1,:)=prod([tempCond{k}(1,1,:) inve{k}(2,3,:)]) + prod([tempCond{k}(1,2,:) inve{k}(3,3,:)]) + prod([tempCond{k}(1,3,:) inve{k}(4,3,:)]);
  B{k}(4,1,:)=prod([tempCond{k}(1,1,:) inve{k}(2,4,:)]) + prod([tempCond{k}(1,2,:) inve{k}(3,4,:)]) + prod([tempCond{k}(1,3,:) inve{k}(4,4,:)]);
  B{k}(1,2,:)=prod([tempCond{k}(2,1,:) inve{k}(2,1,:)]) + prod([tempCond{k}(2,2,:) inve{k}(3,1,:)]) + prod([tempCond{k}(2,3,:) inve{k}(4,1,:)]);
  B{k}(2,2,:)=prod([tempCond{k}(2,1,:) inve{k}(2,2,:)]) + prod([tempCond{k}(2,2,:) inve{k}(3,2,:)]) + prod([tempCond{k}(2,3,:) inve{k}(4,2,:)]);
  B{k}(3,2,:)=prod([tempCond{k}(2,1,:) inve{k}(2,3,:)]) + prod([tempCond{k}(2,2,:) inve{k}(3,3,:)]) + prod([tempCond{k}(2,3,:) inve{k}(4,3,:)]);
  B{k}(4,2,:)=prod([tempCond{k}(2,1,:) inve{k}(2,4,:)]) + prod([tempCond{k}(2,2,:) inve{k}(3,4,:)]) + prod([tempCond{k}(2,3,:) inve{k}(4,4,:)]);  
  B{k}(1,3,:)=prod([tempCond{k}(3,1,:) inve{k}(2,1,:)]) + prod([tempCond{k}(3,2,:) inve{k}(3,1,:)]) + prod([tempCond{k}(3,3,:) inve{k}(4,1,:)]);
  B{k}(2,3,:)=prod([tempCond{k}(3,1,:) inve{k}(2,2,:)]) + prod([tempCond{k}(3,2,:) inve{k}(3,2,:)]) + prod([tempCond{k}(3,3,:) inve{k}(4,2,:)]);
  B{k}(3,3,:)=prod([tempCond{k}(3,1,:) inve{k}(2,3,:)]) + prod([tempCond{k}(3,2,:) inve{k}(3,3,:)]) + prod([tempCond{k}(3,3,:) inve{k}(4,3,:)]);
  B{k}(4,3,:)=prod([tempCond{k}(3,1,:) inve{k}(2,4,:)]) + prod([tempCond{k}(3,2,:) inve{k}(3,4,:)]) + prod([tempCond{k}(3,3,:) inve{k}(4,4,:)]);
  B{k}=B{k}(:); B{k}=reshape(B{k},12,max(size(idx{k}))); 
  tmp=repmat(idx{k},12,1); 
  tmp=tmp(:); 
  tmp2=B{k}(:);
  Jtemp{k} = [indexe{k} tmp tmp2];
end
matlabpool close;
ind=0;
for k=1:nr_of_workers
   ind=ind+max(size(Jtemp{k}));
end
tmp1=zeros(ind,3);
from=1;
to=size(Jtemp{1},1);
for k=1:nr_of_workers
    tmp1(from:to,:)=Jtemp{k};
   if k<nr_of_workers
    from=to+1;
    to=to+size(Jtemp{k+1},1);
   else
     from=to+1;
     to=size(tmp1,1);
   end
end
curDensFromPot = spconvert([tmp1; 3*N M 0]);
curDensFromPot = reshape(curDensFromPot,N,[])';
clear Jtemp;

fprintf('%s%f%s\n','The mapping matrix from potential to current density and element volumes are calculated in ',toc,' seconds');
toc
