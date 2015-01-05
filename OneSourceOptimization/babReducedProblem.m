function [currentArray] = babReducedProblem(w,sQ,tot,ind,pmax,Te,nSources,Ith,zhat)
%% Finds the optimal solution with limited number of current sources
%  on the dominant electrode set in the original solution using branch
%  and bound algorithm
%
%
% Synopsis: [x,fval,dv] = ...
%           babReducedProblem(w, sqrtQ, tot, ind, pmax, Te, nSources, Ith)
%
%
% Input:    w           =   objective function coefficients
%           sqrtQ       =   square root matrix of quadratic const. matrix
%           tot         =   total current bound
%           ind         =   individual electrode current bounds
%           pmax        =   power constraint bound
%           Te          =   matrix linking elec currents to elec potential
%           nSources    =   number of current sources
%           Ith         =   Threshold current value to be set to 0
%
%
% Output:   x       =   optimal electrode current stimulus pattern
%           fval    =   optimized directional current density in the ROI
%           dv      =   dual variables for the constraints

%% Reading inputs and checking sizes
tic;
L = numel(w);
Ltemp = L;
pp = numel(sQ);

if size(ind,1) == L %In case reference electrode bound is not defined
    ind(L+1,:) = ind(end,:); %Make ref elec have the same const as last one
end

if size(ind,2) == 1 %Lower bound = - Upper bound
    ind = [-ind ind];
end

if nSources > 14
    error('The number of current sources is less than 15.');
end

%% First optimization to get the general unconstrained (i.e. there may be
%  as many current sources as number of electrodes) solution
[ca,fval,dv] = optimizationUsingCvxToolbox(w, sQ, tot, ind, pmax);

currentArray.origCurrent = ca;
currentArray.origPot = Te * ca;
currentArray.origObj = fval;
currentArray.origDV = dv;


if nargin <=8
    zhat = fval*9/10;
end

%% Reduce the problem size by setting small currents to 0.
[newVar,percentLoss] = ...
    ridElectrodesWithSmallCurrents(w,sQ,tot,ind,pmax,Ith);
if percentLoss >= 1e-3
    warning('Percentage loss by setting low currents to 0 is %f%s\n',...
        percentLoss,'.');
else
    fprintf('%s%f%s\n','Percentage loss by setting low currents to 0 is ',...
        percentLoss,'.');
end
w = newVar.w;
L = numel(w);
sQ = newVar.sQ;
Te = Te(newVar.idx,newVar.idx);

currentArray.newVar = newVar;
%% Find an initial set of states for the electrodes.
%
% IDEA: Use maybe k clustering to start the initialization.
% IDEA: Ordering of the electrodes for bAb algorithm may be initialized.
fprintf('%s\n','Initializating clustering of electrodes into states.');

nStates = nSources+1;

initSet = cell(nStates,1);
for i = 1:nStates
    initSet{i} = []; % All states are initialized empty.
end
unknownSet = setdiff(1:L,cell2mat(initSet));

%Ordering of the unknown set for the branch and bound algorithm
%unknownSetOrder = f(unknownSet,optimalSolution,linearWeights)
[~,idxOrder] = sort(abs(ca(newVar.idx)));
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
%[~,xhat,zhat,~] = babOne(w,sQ,tot,ind,pmax,Te);
activeSet(1) = 1;
%zhat = -inf;
z(1) = zhat;
t = 0;
totalactSetSize = 1;
labelfet = '0123456789ABCDE';
percentDone = 0;
currentArray.xhatBAB = [];
currentArray.zhatBAB = [];
currentArray.tBAB = [];
currentArray.branchBAB = [];
%Define electrode ordering here or in the loop (if we would like to have
%different ordering of the importance of the electrodes at each branch

%Main loop of branc and bound. Once finished, all the branches will have
%been fathomed.

while ~isempty(activeSet)
    
    fprintf('%.2f%s\t%d\t%d\t',percentDone,'%',t,numel(activeSet));
    totalactSetSize = totalactSetSize + numel(activeSet);
    %100*(1-sum(1./(activeSet-mod(activeSet,nStates)))),'%') is another way
    %to calculate the ratio of finished branches vs total
    
    
    %Choose the next branch. CHANGE!
    %[~,idx] = max(z);
    idx = numel(activeSet);
    r = activeSet(idx);
    activeSet(idx) = [];
    z(idx) = [];
    
    %Determination of p(r) and branching: Fr = R1 U R2 U ... U Rpr. CHANGE!
    pr = nStates;
    parentelecAssign = dec2base(nStates*r,nStates);
    parentelecAssign(1) = [];
    fprintf('%-24s\t',parentelecAssign(1:end-1));
    for i = 1:pr
        %% Setting the feasible set for the branch
        %  Ft+i = Fr n Ri
        
        %% Calculating the xt+i, zt+i for the branch
        % Assign electrode states according to branch number
        % Extract electrode states from the assignment vector and ordering
        elecAssign = dec2base(nStates*r+i-1,nStates);
        elecAssign(1) = [];
        %fprintf('%s\t\t\t',elecAssign);
        
        for j = 1:nStates
            jStateElec = elecAssign == labelfet(j);
            combinedStates{j} = [initSet{j} unknownSetOrder(jStateElec)];
        end
        nZeros = numel(combinedStates{1});
        
        % CVX, CHANGE! (higher the precision, better.)
        cvx_begin quiet
        %cvx_precision high
        cvx_solver sedumi
        
        variable u(L-nZeros)
        
        expression x(L)
        expression y(L+1)
        expression pot(L)
        x(combinedStates{1}) = 0; % 'Not connected' electrodes
        x(setdiff(1:L,combinedStates{1})) = u;
        y = [x; -sum(x)]; %reference elec current = - (sum of the rest)
        pot = 1e-3*Te * x;
        
        dual variable totConstV
        dual variable indConstLB
        dual variable indConstUB
        dual variable powConstV{pp}
        
        maximize w*x
        
        subject to
        totConstV : norm(y,1) <= 2*tot;
        indConstLB : ind(:,1) <= y;
        indConstUB : y <= ind(:,2);
        for il = 1:pp
            powConstV{il} : norm(sQ{il}*x) <= sqrt(pmax(il));
        end
        
        for k = 2:numel(combinedStates)
            for ks = 2:numel(combinedStates{k})
                pot(combinedStates{k}(ks)) == pot(combinedStates{k}(1));
            end
        end
        for kk = 3:numel(combinedStates)
            if ~isempty(combinedStates{kk}) && ~isempty(combinedStates{kk-1})
                x(combinedStates{kk}(1)) >= x(combinedStates{kk-1}(1))
            end
        end
        
        cvx_end
        
        %currentArray.electrodeCurrents(:,r) = x;
        if ~strcmp(cvx_status,'Solved')
            warning('Not solved. Trying different solver.');
            cvx_begin quiet
            %cvx_precision high
            cvx_precision low
            cvx_solver SDPT3
            
            variable u(L-nZeros)
            
            expression x(L)
            expression y(L+1)
            expression pot(L)
            x(combinedStates{1}) = 0; % 'Not connected' electrodes
            x(setdiff(1:L,combinedStates{1})) = u;
            y = [x; -sum(x)]; %reference elec current = - (sum of the rest)
            pot = 1e-3*Te * x;
            
            dual variable totConstV
            dual variable indConstLB
            dual variable indConstUB
            dual variable powConstV{pp}
            
            maximize w*x
            
            subject to
            totConstV : norm(y,1) <= 2*tot;
            indConstLB : ind(:,1) <= y;
            indConstUB : y <= ind(:,2);
            for il = 1:pp
                powConstV{il} : norm(sQ{il}*x) <= sqrt(pmax(il));
            end
            
            for k = 2:numel(combinedStates)
                for ks = 2:numel(combinedStates{k})
                    pot(combinedStates{k}(ks)) == pot(combinedStates{k}(1));
                end
            end
            for kk = 3:numel(combinedStates)
                if ~isempty(combinedStates{kk}) && ~isempty(combinedStates{kk-1})
                    x(combinedStates{kk}(1)) >= x(combinedStates{kk-1}(1))
                end
            end
            cvx_end
        end
        
        zr = cvx_optval;
        %fprintf('%f\t',zr);
        if zhat < zr
            if numel(unique(round(pot(setdiff(1:L,combinedStates{1}))))) < nStates
                %calculate percentage. The brach is pruned (best solution
                % of that branch is found).
                zhat = zr;
                fprintf('%s\t','Z');
                currentArray.xhatBAB(:,end+1) = x; 
                currentArray.zhatBAB(end+1) = zhat;
                currentArray.tBAB(end+1) = t;
                currentArray.branchBAB{end+1} = elecAssign;
            elseif numel(elecAssign) < numel(unknownSetOrder)
                %we add the children branch to the active set
                activeSet(end+1) = [nStates*r+i-1];
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
    end
    t = t + pr;
    fprintf('%.2f\n',zhat);
    if numel(activeSet) > 1e3
        warning('The number of branches exceeded 1000.')
    end
end
currentArray.t = t;
currentArray.convTime = toc;
currentArray.avgActSetSize = totalactSetSize/t*pr;
fprintf('%s%f%s\n','Limited current sources optimization is solved in ',toc,' seconds.');
