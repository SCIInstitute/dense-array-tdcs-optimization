function ec =  testerCode4LimitedCurrentSources(w,sQ,tot,ind,pmax,Te,Ith)

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

ec.caorig = ca;
ec.potorig = Te * ca;
ec.objorig = fval;
ec.dvorig = dv;

%% One current source solution
%[~,xhat,zhat,~] = babOne(w,sQ,tot,ind,pmax,Te);
%ec.ca1 = xhat;
%ec.obj1 = zhat;
%ec.pot1 = Te*xhat;

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
w = newVar.w/norm(newVar.w);
normW = norm(newVar.w);
L = numel(w);
sQ = newVar.sQ;
Te = Te(newVar.idx,newVar.idx);

%% TESTER CODE
clear x;
labelfet = '0123456789ABCDE';
elecAssign = '1111111222222222';
disp(elecAssign);
nStates = 3;

initSet = cell(nStates,1);
for i = 1:nStates
    initSet{i} = []; % All states are initialized empty.
end
unknownSet = setdiff(1:L,cell2mat(initSet));

%Ordering of the unknown set for the branch and bound algorithm
%unknownSetOrder = f(unknownSet,optimalSolution,linearWeights)
[~,idxOrder] = sort(abs(ca(newVar.idx)));
unknownSetOrder = unknownSet(idxOrder);
for j = 1:nStates
    jStateElec = elecAssign == labelfet(j);
    combinedStates{j} = [initSet{j} unknownSetOrder(jStateElec)];
end
nZeros = numel(combinedStates{1});

% CVX, CHANGE! (higher the precision, better.)
cvx_begin
cvx_precision high
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
cvx_end

%currentArray.electrodeCurrents(:,r) = x;
zr = normW*cvx_optval;

fprintf('%s%f%s\n','Limited current sources optimization is solved in ',toc,' seconds.');
