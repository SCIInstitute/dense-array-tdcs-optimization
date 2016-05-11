function [currentArray,xhat,zhat,t] = babOne(w,sQ,tot,ind,pmax,Te)
%% Finds the optimal solution with one current source using branch and
% bound algorithm
%
%
% Synopsis: [x,fval,dv] = babOne(w, sqrtQ, tot, ind, pmax, Te)
%
%
% Input:    w       =   objective function coefficients
%           sQ      =   square root matrix of quadratic constraint matrix
%           tot     =   total current bound
%           ind     =   individual electrode current bounds
%           pmax    =   power constraint bound
%           Te      =   matrix linking elec currents to  elec potential
%
%
% Output:   x       =   optimal electrode current stimulus pattern
%           fval    =   optimized directional current density in the ROI
%           dv      =   dual variables for the constraints

%% Reading inputs and checking sizes
tic;
L = numel(w); %number of electrodes
pp = numel(sQ); % Number of power constraints

if size(ind,1) == L %In case reference electrode bound is not defined
    ind(L+1,:) = ind(end,:); %Make ref elec have the same const as last one
end

if size(ind,2) == 1 %Lower bound = - Upper bound
    ind = [-ind ind];
end

%% First optimization to get the general unconstrained (i.e. there may be
%  as many current sources as number of electrodes) solution
[ca, fval, dv] = optimizationUsingCvxToolbox(w, sQ, tot, ind, pmax);
fprintf('%s%f%s\n','Electrode current optimization is solved in ',toc,' seconds.');
currentArray.origCurrent = ca;
currentArray.origPot = Te * ca;
currentArray.origObj = fval;
currentArray.origDV = dv;

%% Find an initial set of 'not connected', 'source', 'sink' electrodes
fprintf('%s\n','checking the closest one current source solution...');

minX = min(ca);
maxX = max(ca);

nIdx = find(ca < 0.99* minX);
pIdx = find(ca > 0.99* maxX);
zeroIdx = find(abs(ca) < 0.001* min(maxX,-minX));

unknownNegatives = find(ca > 0.99*minX & ca < -0.001* min(maxX,-minX));
unknownPositives = find(ca < 0.99*maxX & ca > 0.001* min(maxX,-minX));
unknownElectrodes = [unknownPositives;unknownNegatives];
nN = numel(unknownNegatives);
nP = numel(unknownPositives);
if (numel(zeroIdx) + nN + nP + numel(nIdx) + numel(pIdx) ~= L)
    error('numbers mismatch');
end

%% BRANCH AND BOUND ALGORITHM

%   zhat: objective function value of the best feasible solution found so far
%   ListSet: List of the unfathomed subsets
%   t: number of branches generated so far
%   F0: the set of all feasible solutions
%   r : the index of the subset selected for branching
%   pr : the number of branches generated from Fr
%   xi : the optimal solution of the relaxed subproblem defined on Ri
%   zi : the upper bound of the objective function on subset Fi
%   ListSet + Fi : operation of adding Fi to the list ListSet
%   ListSet - Fi : operation of deleting Fi from the list ListSet

clear x;
zhat = -inf;
xhat = [];
activeSet(1) = 1;
t = 0;

while ~isempty(activeSet)
    r = activeSet(end);
    activeSet(end) = [];
    t = t + 1;
    
    %Calculate x(r) and z(r)
    elecAssign = dec2bin(r);
    elecAssign(1) = [];
    disp(elecAssign);
    onElectrodes = find(elecAssign == '1');
    offElectrodes = find(elecAssign == '0');
    extraSourcesIdx = unknownElectrodes(onElectrodes(onElectrodes <= nP));
    extraSinksIdx = unknownElectrodes(onElectrodes(onElectrodes > nP));
    extraOffIdx = unknownElectrodes(offElectrodes);
    sourceIndices = [pIdx; extraSourcesIdx];
    sinkIndices = [nIdx; extraSinksIdx];
    zeroIndices = [zeroIdx; extraOffIdx];
    nZeros = numel(zeroIndices);
    
    cvx_begin quiet
    %cvx_precision high
    cvx_solver sedumi
    
    %optimization variable
    variable u(L-nZeros)
    
    %expression of electrode current array, reference added
    expression x(L)
    expression y(L+1)
    expression pot(L)
    x(zeroIndices) = 0;
    ii = 1:L;
    ii(zeroIndices) = [];
    x(ii) = u;
    expression y(L+1)
    y = [x; -sum(x)]; %reference elec current = - (sum of other electrodes)
    pot = 1e-3*Te * x;
    
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
    for i = 1:pp
        powConstV{i} : norm(sQ{i}*x) <= sqrt(pmax(i));
    end
    
    for k = 2:numel(sourceIndices)
        pot(sourceIndices(k)) == pot(sourceIndices(1)); %#ok<EQEFF>
    end
    
    for ks = 2:numel(sinkIndices)
        pot(sinkIndices(ks)) == pot(sinkIndices(1)); %#ok<EQEFF>
    end
    cvx_end
    z(r) = cvx_optval;
    if zhat < z(r)
        if numel(sinkIndices)+numel(sourceIndices)+numel(zeroIndices) == L
            zhat = z(r);
            disp(zhat);
            xhat = x;
        else
            activeSet(end+1:end+2) = [2*r 2*r+1];
        end
    end
end
