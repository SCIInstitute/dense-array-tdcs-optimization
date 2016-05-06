function [StimPattern] = branchBoundAlgorithm4Tdcs(w,Q,tot,ind,pmax,Te,nSources,Ith,zhat,vOrder)
%% Approximates the optimal solution with fewer current sources using branch
%  and bound algorithm. Stimulus pattern is denoted by x in the script. 
%
%
% Synopsis: [StimPattern] = ...
%           branchBoundAlgorithm4Tdcs(w, sqrtQ, tot, ind, pmax, Te, nSources, Ith)
%
%
% Input:    w           =   linear objective function coefficients
%                           ( objective: maximize w' * x )
%           Q           =   quadratic power constraint matrix
%                           ( power constraint: x' * Q * x <= pmax )
%           tot         =   total injected current bound
%           ind         =   individual electrode current bounds
%           pmax        =   power constraint bound
%           Te          =   matrix linking electrode currents to electrode 
%                           potential (Electrode potentials = Te * x)
%           nSources    =   number of available current sources
%           Ith         =   Threshold current value to be set to 0
%           vOrder      =   Vertical ordering of electrodes. Choose: 'high
%                           margin', 'absolute', 'descend', 'contribution'
%
%
% Output:   StimPattern =   optimal electrode current stimulus pattern 

% NOTES: 
%   1. Potential unit is V

%% Reading inputs and checking sizes
tic;
L = numel(w);
pp = numel(Q);

if size(ind,1) == L %In case reference electrode bound is not defined
    ind(L+1,:) = ind(end,:); %Make ref elec have the same const as last one
end

if size(ind,2) == 1 %Lower bound = - Upper bound
    ind = [-ind ind];
end

if nSources > 20
    error('The number of current sources is less than 21.');
end

sQ = cell(1,pp);
for i=1:pp
    [sQ{i},p] = chol(Q{i});
    if p>0
        warning('quadratic matrix is not positive definite.');
        minEig = min(eig(Q{i}));
        sQ{i} = chol(Q{i}+(-minEig+eps)*eye(size(Q{i})));
    end
end

%% First optimization to get the general solution
[ca,fval,dv] = optimizationUsingCvxToolbox(w, Q, tot, ind, pmax);

StimPattern.origCurrent = ca;
StimPattern.origPot = Te * ca;
StimPattern.origObj = fval;
StimPattern.origDV = dv;

%% Reduce the problem size by setting small currents to 0.
[newVar,percentLoss] = ...
    ridElectrodesWithSmallCurrents(w,Q,tot,ind,pmax,Ith);
if percentLoss >= 1e-3
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

sQ = cell(1,pp);
for i=1:pp
    [sQ{i},p] = chol(Q{i});
    if p>0
        warning('quadratic matrix is not positive definite.');
        minEig = min(eig(Q{i}));
        sQ{i} = chol(Q{i}+(-minEig+eps)*eye(size(Q{i})));
    end
end

StimPattern.newVar = newVar;

%% Find an initial set of states for the electrodes.
%
nStates = nSources+1;

if isempty(zhat)
    %Use kmeans clustering to start the initialization.
    [idx0,c0] = kmeans(Te*ca(newVar.idx),nSources);
    conf0{1} = [];
    for i=2:nStates
        conf0{i} = find(idx0 == i-1);
    end
    [StimPattern.x0,zhat] = solveRelaxedProblem(w,sQ,tot,ind,pmax,Te,conf0);
    pots2 = Te*ca(newVar.idx);
    bigCurrents = ca(newVar.idx) >= 2*Ith;
    bigCurrentIdx = find(bigCurrents == 1);
    [idx1,c1] = kmeans(pots2(bigCurrents),nSources);
    conf1{1} = find(bigCurrents==0);
    for i=2:nStates
        conf1{i} = bigCurrentIdx(find(idx1 == i-1));
    end
    [StimPattern.x1,zhat1] = solveRelaxedProblem(w,sQ,tot,ind,pmax,Te,conf1);
    if zhat1 > zhat
        fprintf('%s\n','Using not connected in the initialization is better.');
        fprintf('%s%f%s\n','Percentage increase for the initialization is ',100*(zhat1-zhat)/zhat,'.');
        zhat = zhat1;
        disp(zhat1);
    end
    
    %Use the indices in the node selection process.
    if strcmp(vOrder,'high margin')
    [Dist,cIdx] = pdist2(c0,Te*ca(newVar.idx),'euclidean','Smallest',1);
    [~,idxOrder] = sort(Dist,'descend');
    end
    if strcmp(vOrder,'high margin') && zhat1 == zhat
    [Dist,cIdx] = pdist2(c1,Te*ca(newVar.idx),'euclidean','Smallest',1);
    [~,idxOrder] = sort(Dist,'descend');
    end
end

% An initial state assignment could be performed here.
fprintf('%s\n','Initializating clustering of electrodes into states.');
initSet = cell(nStates,1);
for i = 1:nStates
    initSet{i} = []; % All states are initialized empty.
end
unknownSet = setdiff(1:L,cell2mat(initSet));

%Vertical ordering of electrodes in the branch and bound algorithm.
switch vOrder
    case 'absolute'
[~,idxOrder] = sort(abs(ca(newVar.idx)));
    case 'descend'
[~,idxOrder] = sort(abs(w),'descend');
    case 'contribution'
[~,idxOrder] = sort(abs(w .* ca(newVar.idx)'));
end
unknownSetOrder = unknownSet(idxOrder);

%% BRANCH AND BOUND ALGORITHM

%   zhat: objective function value of the best feasible solution so far
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
z(1) = zhat;
t = 0;
activeSetSizeCount = 1;
labelAlphabet = '0123456789ABCDEFGHIJK';
percentDone = 0;
StimPattern.xhatBAB = [];
StimPattern.zhatBAB = [];
StimPattern.tBAB = [];
StimPattern.branchBAB = [];

%% Branch and bound main loop.
while ~isempty(activeSet)
    
    fprintf('%.2f%s\t',percentDone,'%');
    activeSetSizeCount = activeSetSizeCount + numel(activeSet);
    %100*(1-sum(1./(activeSet-mod(activeSet,nStates)))),'%') is another way
    %to calculate the ratio of finished branches vs total
    
    
    %Choose the next branch to check. CHANGE!
    %[~,idx] = max(z);
    idx = numel(activeSet);
    nextBranch = activeSet{idx};
    activeSet(idx) = [];
    z(idx) = [];
    
    %Determination of p(r) and branching: Fr = R1 U R2 U ... U Rpr. CHANGE!
    pr = max(2,min(numel(unique(nextBranch(nextBranch ~= '0')))+2,nStates));
%     if nStates*nextBranch+pr-1 > 2^52
%         parentelecAssign = dec2baseInfDigits(nStates*nextBranch,nStates);
%     elseif numel(nextBranch) >= 2
%         %When the branch idx is too high,i.e. higher than 2^52 already.
%         for ir = 1:numel(nextBranch)
%             parentelecAssign = [dec2base(nextBranch(ir),nStates) parentelecAssign];
%         end
%     else
%         parentelecAssign = dec2base(nextBranch,nStates);
%     end
    fprintf('%-36s\t', nextBranch);
    for i = 1:pr
        %% Setting the feasible set for the branch
        %  Ft+i = Fr n Ri
        
        %% Calculating the xt+i, zt+i for the branch
        % Assign electrode states according to branch number
        % Extract electrode states from the assignment vector and ordering
        elecAssign = [nextBranch labelAlphabet(i)];
        %Eliminate the branch if it is not sorted.
        [~,ia,~] = unique(elecAssign(elecAssign ~= '0'));
        if issorted(ia)            
            for j = 1:nStates
                jStateElec = elecAssign == labelAlphabet(j);
                conf{j} = [initSet{j} unknownSetOrder(jStateElec)];
            end
            [x,zr,pot] = solveRelaxedProblem(w,Q,tot,ind,pmax,Te,conf);
            
            if zhat < zr
                if numel(unique(round(pot(setdiff(1:L,conf{1}))))) < nStates
                    %calculate percentage. The brach is pruned (best solution
                    % of that branch is found).
                    zhat = zr;
                    fprintf('%s\t','Z');
                    StimPattern.xhatBAB(:,end+1) = x;
                    StimPattern.zhatBAB(end+1) = zhat;
                    StimPattern.tBAB(end+1) = t;
                    StimPattern.branchBAB{end+1} = elecAssign;
                elseif numel(elecAssign) < numel(unknownSetOrder)
                    activeSet{end+1} = elecAssign;
                    z(end+1) = zr;
                    fprintf('%s\t','B');
                else
                    fprintf('%s\t','L');
                end
            else
                %calculate percentage. The branch is pruned.
                percentDone = percentDone + 100*(nStates^-numel(elecAssign));
                fprintf('%s\t','F');
            end
        else
            percentDone = percentDone + 100*(nStates^-numel(elecAssign));
            fprintf('%s\t','S');
        end
    end
    t = t + pr;
    fprintf('%.2f\n',zhat);
    if numel(activeSet) > 1e3
        warning('The number of branches exceeded 1000.')
    end
end
StimPattern.t = t;
StimPattern.convTime = toc;
StimPattern.avgActSetSize = activeSetSizeCount/t*pr;
StimPattern.percentLoss = 100*(StimPattern.origObj-zhat)/StimPattern.origObj;
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
