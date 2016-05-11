function [results] = optimizeWithLimitedSourcesAndElectrodes(w,Q,tot,ind,pmax,Te,nSources,nElectrodes,Ithresh,zhat,vOrder)
%% Finds the optimal electrode current stimulus pattern  with limited number of
%  current sources and electrodes using a branch and bound algorithm
%
%
% Synopsis: [results] = optimizeWithLimitedSourcesAndElectrodes(w,...
%             Q, tot, ind, pmax, Te, nSources, nElectrodes, Ithresh, zhat, vOrder)
%
%
% Input:    w           =   objective function coefficients
%           Q           =   matrix for the quadratic constraint on the current
%                           power in the brain outside the ROI
%           tot         =   bound on the total injected current
%           ind         =   bound on the individual electrode currents
%           pmax        =   bound on the current power in the brain outside
%                           the ROI
%           Te          =   matrix linking electrode currents to electrode
%                           potentials
%           nSources    =   number of available current sources
%           nElectrodes =   number of available electrodes
%           Ithresh     =   threshold current; smaller values are set to 0
%           vOrder      =   vertical ordering of the  electrodes in the
%                           branch and bound algorithm. Choices: 'high margin'
%                           (default), 'absolute', 'descend',
%                           'contribution'
%
%
% Output:   results     =   optimal electrode current stimulus pattern

% Notes:   1. Potential unit is V.

%% Reading inputs and checking sizes
tic;
L = numel(w);
pp = numel(Q);

if size(ind,1) == L %In case reference electrode bound is not defined
    ind(L+1,:) = ind(1,:); %reference electrode bound = bound on first electrode
end

if size(ind,2) == 1 %Lower bound = - Upper bound
    ind = [-ind ind];
end

if nSources > 20
    error('The number of available current sources should be less than 21.');
end

if isempty(nElectrodes)
    nElectrodes = L;
end

if isempty(vOrder)
    vOrder = 'high margin';
end

if isempty(zhat)
    zhat = -inf;
end

%% Find the unconstrained solution
[ca,fval,dv] = optimizationUsingCvxToolbox(w, Q, tot, ind, pmax);

results.unconstrainedSolution.currentArray = ca;
results.unconstrainedSolution.potential = Te * ca;
results.unconstrainedSolution.objectiveValue = fval;
results.unconstrainedSolution.dualVariables4Constraints = dv;

%% Reduce the problem size by setting small currents to 0.
if ~isempty(Ithresh)
    [newVar,percentLoss] = ridElectrodesWithSmallCurrents(w,Q,tot,ind,pmax,Ithresh);
    if percentLoss >= 5
        error('Percentage loss too high.');
    elseif percentLoss >= 1
        warning('Percentage loss by setting low currents to 0 is %f%s\n',...
            percentLoss,'.');
    else
        fprintf('%s%f%s\n','Percentage loss by setting low currents to 0 is ',...
            percentLoss,'.');
    end
    w = newVar.w;
    L = numel(w);
    Q = newVar.Q;
    Te = Te(newVar.idx,newVar.idx);
    ca = ca(newVar.idx);
    results.newOptimizationVariables = newVar;
end

sqrtQ = cell(1,pp);
for i=1:pp
    %Cholesky factor for reasons explained above
    [sqrtQ{i},p] = chol(Q{i});
    if p>0
        warning('quadratic matrix is not positive definite.');
        minEig = min(eig(Q{i}));
        sqrtQ{i} = chol(Q{i}+(-minEig+eps)*eye(size(Q{i})));
    end
end

nStates = nSources + 1; %'not connected' is a state

%% Initialize states for the electrodes.

%Use k-means clustering to start the initialization.
[idx0,c0] = kmeans(Te*ca,nSources);
configuration0{1} = [];
for i=2:nStates
    configuration0{i} = find(idx0 == i-1);
end

%calculate zhat.
[results.initialState.x0,zhat] = solveRelaxedProblem(w,sqrtQ,tot,ind,pmax,Te,configuration0);

if strcmp(vOrder,'high margin')
    [Dist,cIdx] = pdist2(c0,Te*ca,'euclidean','Smallest',1);
    %IDEA: use also the center indices in the node selection process.
    [~,idxOrder] = sort(Dist,'descend');
end

fprintf('%s\n','Initializating clustering of electrodes into states.');

initSet = cell(nStates,1);
for i = 1:nStates
    initSet{i} = []; % All states are initialized empty.
end
unknownSet = setdiff(1:L,cell2mat(initSet));

%Ordering of the unknown set for the branch and bound algorithm
%unknownSetOrder = f(unknownSet,optimalSolution,linearWeights)
switch vOrder
    case 'absolute'
        [~,idxOrder] = sort(abs(ca));
    case 'descend'
        [~,idxOrder] = sort(abs(w),'descend');
    case 'contribution'
        [~,idxOrder] = sort(abs(w .* ca'));
end
unknownSetOrder = unknownSet(idxOrder);

%% BRANCH AND BOUND ALGORITHM

%   zhat: objective function value of the best feasible solution found so far
%   activeSet: List of the unfathomed subsets
%   t: number of branches generated so far
%   F0: the set of all feasible solutions
%   r : the index of the subset selected for branching
%   pr : the number of branches generated from Fr
%   xi : the optimal solution of the relaxed subproblem defined on Ri
%   zi : the upper bound of the objective function on subset Fi
%   activeSet + Fi : operation of adding Fi to the list activeSet
%   activeSet - Fi : operation of deleting Fi from the list activeSet

%Initialization
clear x;
activeSet{1} = '';
%zhat = -inf;
z(1) = zhat;
t = 0;
totalActSetSize = 1;
labelAlphabet = '0123456789ABCDEFGHIJK';
percentDone = 0;
results.xhatBAB = [];
results.zhatBAB = [];
results.tBAB = [];
results.branchBAB = [];

zTotalDepth = zeros(L);
nCreatedBranchesDepth = zeros(L);
%Define electrode ordering here or in the loop (if we would like to have
%different ordering of the importance of the electrodes at each branch

%% Main loop of the branch and bound algorithm.
while ~isempty(activeSet)
    
    fprintf('%.2f%s\t%d\t%d\t',percentDone,'%',t,numel(activeSet));
    totalActSetSize = totalActSetSize + numel(activeSet);
    %100*(1-sum(1./(activeSet-mod(activeSet,nStates)))),'%') is another way
    %to calculate the ratio of finished branches vs total
    
    
    %Choose the next branch to check.
    %[~,idx] = max(z);
    idx = numel(activeSet);
    nextBranch = activeSet{idx};
    activeSet(idx) = [];
    z(idx) = [];
     
    %Determination of p(r) and branching: Fr = R1 U R2 U ... U Rpr.
    pr = max(2,min(numel(unique(nextBranch(nextBranch ~= '0')))+2,nStates));
    fprintf('%-24s\t', nextBranch);
    branchDepth = numel(nextBranch)+1;
    for i = 1:pr
        %% Setting the feasible set for the branch
        %  Ft+i = Fr n Ri
        
        %% Calculating the xt+i, zt+i for the branch
        
        % Assign electrode states according to branch number
        % Extract electrode states from the assignment vector and ordering
        elecAssign = [nextBranch labelAlphabet(i)];
        nConnectedElectrodes = numel(elecAssign(elecAssign ~= '0'));
        
        %continue if # connected electrodes exceed available number of electrodes
        if(nConnectedElectrodes == nElectrodes)
            elecAssign = [elecAssign repmat('0',1,L-numel(elecAssign))];
        end
        
        %Don't create a branch if the assignment is not sorted, to avoid
        %solving identical problems
        [~,ia,~] = unique(elecAssign(elecAssign ~= '0'));
        if issorted(ia)
            for j = 1:nStates
                jStateElec = elecAssign == labelAlphabet(j);
                conf{j} = [initSet{j} unknownSetOrder(jStateElec)];
            end
            [x,zr,pot] = solveRelaxedProblem(w,sqrtQ,tot,ind,pmax,Te,conf);
            zTotalDepth(branchDepth) = zTotalDepth(numel(elecAssign))+zr;
            nCreatedBranchesDepth(branchDepth) = nCreatedBranchesDepth(branchDepth)+1;
            if zhat < zr
                if numel(unique(round(pot(setdiff(1:L,conf{1}))))) < nStates
                    %update the best solution
                    zhat = zr;
                    fprintf('%s\t','Z');
                    results.xhatBAB(:,end+1) = x;
                    results.zhatBAB(end+1) = zhat;
                    results.tBAB(end+1) = t;
                    results.branchBAB{end+1} = elecAssign;
                elseif numel(elecAssign) < numel(unknownSetOrder)
                    %we add the children branch to the active set
                    activeSet{end+1} = elecAssign;
                    z(end+1) = zr;
                    fprintf('%s\t','B');
                else
                    fprintf('%s\t','L');
                end
            else
                %calculate percentage. The branch is pruned.
                percentDone = percentDone +100*(nStates^-numel(elecAssign));
                fprintf('%s\t','F');
            end
        else
            percentDone = percentDone +100*(nStates^-numel(elecAssign)-nSources);
            fprintf('%s\t','S');
        end
    end
    t = t + pr;
    fprintf('%.2f\n',zhat);
    if numel(activeSet) > 1e3
        warning('The number of branches exceeded 1000.')
    end
end
results.t = t;
results.convTime = toc;
results.avgActSetSize = totalActSetSize/t*pr;
results.percentLoss = 100*(results.origObj-zhat)/results.origObj;
fprintf('%s%f%s\n','Limited current sources optimization is solved in ',toc,' seconds.');
end

function [currentArray,fval,pot] = solveRelaxedProblem(w,sQ,tot,ind,pmax,Te,conf)
% Solves the relaxed problem on the electrode set ('not connected' set excluded')
nZeros = numel(conf{1});
L = numel(w);
normW = norm(w);
w = w/normW;
Te = Te/normW;

%% CVX, high precision, sedumi as solver
cvx_begin quiet
%cvx_precision high
cvx_solver sedumi

variable u(L-nZeros)

expression x(L)
x(conf{1}) = 0;
x(setdiff(1:L,conf{1})) = u;
expression y(L+1)
expression pot(L)
y = [x; -sum(x)];
pot = Te * x;

maximize w*x

subject to
norm(y,1) <= 2*tot;
ind(:,1) <= y;
y <= ind(:,2);
for pp = 1:numel(sQ)
    norm(sQ{pp}*x) <= sqrt(pmax(pp));
end

for k = 2:numel(conf)
    for ks = 2:numel(conf{k})
        pot(conf{k}(ks)) == pot(conf{k}(1));
    end
end
cvx_end

%% Not solved, medium precision is tried
if ~strcmp(cvx_status,'Solved')
    fprintf('%s','.');
    cvx_begin quiet
    cvx_precision low
    cvx_solver sedumi
    
    variable u(L-nZeros)
    
    expression x(L)
    x(conf{1}) = 0;
    x(setdiff(1:L,conf{1})) = u;
    expression y(L+1)
    expression pot(L)
    y = [x; -sum(x)]; %reference elec current = - (sum of the rest)
    pot = Te * x;
    
    maximize w*x
    
    subject to
    norm(y,1) <= 2*tot;
    ind(:,1) <= y;
    y <= ind(:,2);
    for pp=1:numel(sQ)
        norm(sQ{pp}*x) <= sqrt(pmax(pp));
    end
    
    for k = 2:numel(conf)
        for ks = 2:numel(conf{k})
            pot(conf{k}(ks)) == pot(conf{k}(1));
        end
    end
    cvx_end
end

%% Still not solved, SDPT3 is tried with low precision
if ~strcmp(cvx_status,'Solved')
    fprintf('%s','.');
    cvx_begin quiet
    cvx_precision low
    cvx_solver SDPT3
    variable u(L-nZeros)
    
    expression x(L)
    x(conf{1}) = 0;
    x(setdiff(1:L,conf{1})) = u;
    expression y(L+1)
    expression pot(L)
    y = [x; -sum(x)]; %reference elec current = - (sum of the rest)
    pot = Te * x;
    
    maximize w*x
    
    subject to
    norm(y,1) <= 2*tot;
    ind(:,1) <= y;
    y <= ind(:,2);
    for pp=1:numel(sQ)
        norm(sQ{pp}*x) <= sqrt(pmax(pp));
    end
    
    for k = 2:numel(conf)
        for ks = 2:numel(conf{k})
            pot(conf{k}(ks)) == pot(conf{k}(1));
        end
    end
    cvx_end
end
if ~strcmp(cvx_status,'Solved')
    warning('Not solved.');
end
fval = normW*cvx_optval;
currentArray = x;
end

% function [elecAssign,r,chunkSize] = dec2baseInfDigits(nodeIdx,nStates)
% elecAssign = '';
% while nodeIdx ~= 0
% elecAssign = [num2str(rem(nodeIdx,nStates)) elecAssign];
% nodeIdx = floor(nodeIdx/nStates);
% end
% end
