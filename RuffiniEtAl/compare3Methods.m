function [getal, detal, retal, stats] = compare3Methods(dmoc, ruff, guler, G, T, field, roi, cortex, sfn, vole)
% [getal,detal,retal,stats] = compare3Methods(dmoc,ruff,guler,G,T,field,roi,sfn,vole)
if isempty(dmoc)
    Ed= 0.4;
    k = 0.01;
    variables = calcQuadraticLinearTerms(G, T, (reshape(Ed*sfn,1,[]))', roi, field==4 | field==5);
    dmoc.variables = variables;
    [dmoc.Q,dmoc.v] = findWeightedQuadratic(variables.Q_ROI,variables.Q_nonROI,variables.v,k,nnz(roi),nnz(field==4|field==5));
    dmoc.ind = 1;
    save(['dmoc' num2str(k) 'Ed' num2str(Ed*10)],'dmoc');
    
end
if isempty(ruff)
    E0 = 0.3;
    [ruff.Q, ruff.b] = findQuadraticLinearTerms4FocalRoi(G, T, roi, cortex, sfn , E0);
    save(['ruffE0' num2str(E0*10)],'ruff');
end
if isempty(guler)
    [ guler.w, guler.Q ]  = linearQuadraticCoefficientCalculation(roi, field==4 | field==5, sfn, T, G, vole);
    save('guler','guler');
end

%% optimization using Dmochowski et al.
[detal.currentArray, detal.fval, detal.dv] = weightedLeastSquaresL1ConstraintByParraEtAl(dmoc.Q, dmoc.v, dmoc.ind);

%% calculate bounds for the additional constraints on the other two methods.
rtot = (norm(detal.currentArray,1)+abs(sum(detal.currentArray)))/2;
ind = dmoc.ind;

%% optimization using Ruffini et al.
[retal.currentArray, retal.fval, retal.dv] = weightedLeastSquaresTotalIndividual(ruff.Q, ruff.b, rtot, ind);

%% optimization using our method.
gtot = min((norm(detal.currentArray,1) + abs(sum(detal.currentArray)))/2,(norm(retal.currentArray,1) + abs(sum(retal.currentArray)))/2);
for i =1:numel(guler.Q)
    pBrain(i) = max(detal.currentArray'*guler.Q{i}*detal.currentArray,retal.currentArray'* guler.Q{i}*retal.currentArray);
end
[getal.currentArray, getal.fval, getal.dv] = optimizationUsingCvxToolbox(guler.w,guler.Q,gtot,ind,pBrain);


stats.constraints.tot.detal = rtot;
stats.constraints.tot.retal = rtot;
stats.constraints.tot.getal = gtot;
stats.constraints.ind = ind;
stats.constraints.pow = pBrain;
stats.constraints.powDmoc = detal.currentArray'*guler.Q{i}*detal.currentArray;
stats.constraints.powRuff = retal.currentArray'*guler.Q{i}*retal.currentArray;

eca(:,1) = [detal.currentArray; -sum(detal.currentArray)];
eca(:,2) = [retal.currentArray; -sum(retal.currentArray)];
eca(:,3) = [getal.currentArray; -sum(getal.currentArray)];
stats.eca = eca;

if nargin >= 4
    eField = G * (T * getal.currentArray);
    eIntensity = sqrt(sum(reshape(eField.*eField,3,[])));
    clear eField;
    save('gSolution','eIntensity','getal');
    
    stats.gcmax = max(eIntensity(field~=8));
    stats.gcmaxb = max(eIntensity(field==4 | field ==5));
    stats.gcaveb = sum(eIntensity(field==4 | field== 5) .* vole(field ==4 | field==5))/sum(vole(field==4 | field ==5));
    stats.gcmedb = median(eIntensity(field==4 | field ==5));
    stats.gcaveROI = sum(eIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1));
    
    
    eField = G * (T* detal.currentArray);
    eIntensity = sqrt(sum(reshape(eField.*eField,3,[])));
    clear eField;
    save('dSolution','eIntensity','detal');
    
    stats.dcmax = max(eIntensity(field~=8));
    stats.dcmaxb = max(eIntensity(field==4 | field ==5));
    stats.dcaveb = sum(eIntensity(field==4 | field== 5) .* vole(field ==4 | field==5))/sum(vole(field==4 | field ==5));
    stats.dcmedb = median(eIntensity(field==4 | field ==5));
    stats.dcaveROI = sum(eIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1));
    
    
    eField = G * (T*retal.currentArray);
    eIntensity = sqrt(sum(reshape(eField.*eField,3,[])));
    clear eField;
    save('rSolution','eIntensity','retal');
    
    stats.rcmax = max(eIntensity(field~=8));
    stats.rcmaxb = max(eIntensity(field==4 | field ==5));
    stats.rcaveb = sum(eIntensity(field==4 | field== 5) .* vole(field ==4 | field==5))/sum(vole(field==4 | field ==5));
    stats.rcmedb = median(eIntensity(field==4 | field ==5));
    stats.rcaveROI = sum(eIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1));
    
    save('currentStatsFinal','retal','detal','getal','stats');
end

end

