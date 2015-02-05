function [newVar,percentLoss] = ...
    ridElectrodesWithSmallCurrents(w, sQ, tot, ind, pmax, Ith)
%% Reformulates the hd-tdcs optimization problem such that the electrodes 
% with small currents in the original solution are set to 0. 
%
% Written by: Seyhmus Guler, 11/30/14
%
%
% Synopsis: [wHat, sqHat, tot, ind, pmax] =
%           ridElectrodesWithSmallCurrents(w, sQ, tot, ind, pmax, Ith)
%
%
% Input:    w           =   objective function coefficients
%           sQ          =   square root matrix of quadratic const. matrix
%           tot         =   total current bound
%           ind         =   individual electrode current bounds
%           pmax        =   power constraint bound
%           Ith         =   Threshold current
%
%
% Output:   newVar          =   New variables.
%               .w          =   reduced problem objective coefficients
%               .sQ         =   reduced problem sqrt matrix of quad matrix
%               .tot        =   total current bound
%               .ind        =   individual electrode current bounds
%               .pmax       =   power constraint bound
%           percentLoss     =   Percentage loss in the obj func when the
%                               small electrode currents are set to 0.

% Notes: 
%   1.  This function is useful when we want to - as next step - find
%   limited current sources solution. For that, we are required to solve
%   a mixed integer program (integers: the states of electrodes, i.e. which
%   electrode is connected to which source, reals: the current value
%   applied from each source) and thus it is useful to reduce the number of
%   electrodes over which we optimize the input currents.

%% Reading the inputs

%% Solving the original problem using CVX toolkit
[ca, fval, dv] = optimizationUsingCvxToolbox(w,sQ,tot,ind,pmax);

%% Converting the problem to reduced form
idx = abs(ca) >= Ith;
fprintf('%s%d%s\n','The number of remaining electrode set is ',nnz(idx),'.');
newVar.w = w(idx);
for i = 1:numel(sQ)
newVar.sQ{i} = sQ{i}(:,idx);
end
newVar.tot = tot;
newVar.ind = ind;
newVar.pmax = pmax;
newVar.idx= idx;
percentLoss = (fval - newVar.w*ca(idx))/fval*100;

end
