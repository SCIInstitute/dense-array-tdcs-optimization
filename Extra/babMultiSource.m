function [currentArray,xhat,zhat,t] = babMultiSource(w,sQ,tot,ind,pmax,Te,nSources)
%% Finds the optimal solution with limited number of current sources
%  using branch and bound algorithm
%
%
% Synopsis: [x,fval,dv] = babMultiSource(w, sqrtQ, tot, ind, pmax, Te)
%
%
% Input:    w           =   objective function coefficients
%           sqrtQ       =   square root matrix of quadratic const. matrix
%           tot         =   total current bound
%           ind         =   individual electrode current bounds
%           pmax        =   power constraint bound
%           nSources    =   number of current sources
%
%
% Output:   x       =   optimal electrode current stimulus pattern
%           fval    =   optimized directional current density in the ROI
%           dv      =   dual variables for the constraints

%% Reading inputs and checking sizes
tic;
if nSources > 7
    error('The number of current sources is less than 7.');
end

%% First optimization to get the general unconstrained (i.e. there may be
%  as many current sources as number of electrodes) solution
[ca,fval,dv] = optimizationUsingCvxToolbox(w, sQ, tot, ind, pmax);

currentArray.origCurrent = ca;
currentArray.origPot = Te * ca;
currentArray.origObj = fval;
currentArray.origDV = dv;

%% Reduce the problem size by setting small currents to 0.
[newVar,percentLoss] = ...
    ridElectrodesWithSmallCurrents(w,sQ,tot,ind,pmax,1e-6);

%% Find an initial set of states for the electrodes.
%  
% IDEA: Use maybe k clustering to start the initialization.
% IDEA: Ordering of the electrodes for bAb algorithm may be initialized.
fprintf('%s\n','Initializating clustering of electrodes into states.');

nStates = 2*nSources+1;

minX = min(x);
maxX = max(x);
nIdx = find(x < 0.99* minX);
pIdx = find(x > 0.99* maxX);

initSet = cell(nStates,1);
initSet{1} = find(abs(x) < 1e-6)'; % Not connected set: {i: I_i < 1e-9 A}
for i = 2:nStates
    initSet{i} = []; %Other states are initialized empty.
end
unknownSet = setdiff(1:L,cell2mat(initSet));

%Ordering of the unknown set for the branch and bound algorithm
%unknownSetOrder = f(unknownSet,optimalSolution,linearWeights)
unknownSetOrder = unknownSet;

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
zhat = -inf;
xhat = [];
activeSet(1) = 1;
z(1) = -inf; 
t = 0;
labelfet = '0123456789ABCDE';
%Define electrode ordering here or in the loop (if we would like to have
%different ordering of the importance of the electrodes at each branch

%Main loop of branc and bound. Once finished, all the branches will have
%been fathomed.
while ~isempty(activeSet)
    
    %Choose the next branch. CHANGE!
    [~,idx] = max(z);
    r = activeSet(idx);
    activeSet(idx) = [];
    z(idx) = [];
    
    %Determination of p(r) and branching: Fr = R1 U R2 U ... U Rpr. CHANGE!
    pr = nStates;
    
    for i = 1:pr
        %% Setting the feasible set for the branch
        %  Ft+i = Fr n Ri 
        
        %% Calculating the xt+i, zt+i for the branch 
        % Assign electrode states according to branch number 
        % Extract electrode states from the assignment vector and ordering
        elecAssign = dec2base(nStates*r+i-1,nStates);
        elecAssign(1) = [];
        disp(elecAssign);
        
        for j = 1:nStates
            jStateElec = elecAssign == labelfet(j);
            combinedStates{j} = [initSet{j} unknownSetOrder(jStateElec)];
        end
        nZeros = numel(combinedStates{1});
        
        % CVX, CHANGE! (higher the precision, better.) 
        cvx_begin quiet
        cvx_precision low
        cvx_solver SDPT3
        
        %optimization variable
        variable u(L-nZeros)
        
        %expression of electrode current array, reference added
        expression x(L)
        expression y(L+1)
        expression pot(L)
        x(combinedStates{1}) = 0; % 'Not connected' electrodes
        x(setdiff(1:L,combinedStates{1})) = u; 
        y = [x; -sum(x)]; %reference elec current = - (sum of other electrodes)
        pot = Te * x;
        
        %dual variables to check active/passive constraints
        dual variable totConstV
        dual variable indConstLB
        dual variable indConstUB
        dual variable powConstV{pp}
        
        %maximize current density in the ROI along the desired direction
        maximize w*x
        
        %subject to total and ind electrode currents, and the power constraints
        subject to
        totConstV : norm(y,1) <= 2*tot; %#ok<*VUNUS> %total current bound
        indConstLB : ind(:,1) <= y; %lower bound on individual currents
        indConstUB : y <= ind(:,2); %#ok<*BDSCI> %upper bound on individual currents
        %power in the avoid regions is bounded:
        for il = 1:pp
            powConstV{il} : norm(sQ{il}*x) <= sqrt(pmax(il));
        end
        
        for k = 1:numel(combinedStates)
            for ks = 2:numel(combinedStates{k})
            pot(combinedStates{k}(ks)) == pot(combinedStates{k}(1)); %#ok<EQEFF>
            end
        end
        cvx_end
        
        %currentArray.electrodeCurrents(:,r) = x;
        zr = cvx_optval;
        if zhat < zr
            if unique(round(pot(unknownSetOrder)*1e9) < nStates)
                zhat = zr;
                disp(zhat);
                xhat = x;
            elseif numel(elecAssign) < numel(unknownSetOrder)-1 
                activeSet(end+1) = [nStates*r+i-1];
                z(end+1) = zr;
            end
        end
    end
    t = t + pr;
    if numel(activeSet) > 1e3
        warning('The number of branches exceeded 1000.')
    end
end
fprintf('%s%f%s\n','Limited current sources optimization is solved in ',toc,' seconds.');
