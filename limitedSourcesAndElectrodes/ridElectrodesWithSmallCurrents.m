function [newVar,percentLoss] = ...
    ridElectrodesWithSmallCurrents(w, Q, tot, ind, pmax, Ith)
%% Reformulates dense array tDCS optimization problem such that the electrodes 
% with small currents in the original solution are set to 0. 
%
% Synopsis: [newVar, percentLoss] =
%           ridElectrodesWithSmallCurrents(w, Q, tot, ind, pmax, Ith)
%
%
% Input:    w           =   objective function coefficients
%           Q           =   quadratic constraint matrix
%           tot         =   total current bound
%           ind         =   individual electrode current bounds
%           pmax        =   power constraint bound
%           Ith         =   threshold current
%
%
% Output:   newVar          =   new parameters for the optimization
%               .w          =   reduced problem objective coefficients
%               .Q          =   reduced problem quadratic matrix
%               .tot        =   total current bound
%               .ind        =   individual electrode current bounds
%               .pmax       =   power constraint bound
%               .idx        =   indices of the electrodes kept in the
%                               reduced problem
%           percentLoss     =   Percentage loss in the obj func when the
%                               small electrode currents are set to 0.

% Notes:    1.  This function is useful when we want to find solutions with
%               fewer current sources. For that, we are required to solve
%               a mixed integer program (integers: the states, i.e. which
%               electrode is connected to which source, reals: the current value
%               applied from each source) and thus we reduce here the number of
%               electrodes over which we optimize the input currents.

%% Reading the inputs

%% Solving the original problem using CVX
[ca, fval, dv] = optimizationUsingCvxToolbox(w, Q, tot, ind, pmax);

%% Converting the problem to reduced form
idx = abs(ca) >= Ith;
fprintf('%s%d%s\n','The number of the remaining electrode set is ',nnz(idx),'.');
newVar.w = w(idx);
for i = 1:numel(Q)
newVar.Q{i} = Q{i}(idx,idx);
end
newVar.tot = tot;
newVar.ind = ind;
newVar.pmax = pmax;
newVar.idx= idx;
[~,fval2,~] = optimizationUsingCvxToolbox(newVar.w,newVar.Q,newVar.tot,newVar.ind,newVar.pmax);
percentLoss = (fval - fval2)/fval*100;
end
