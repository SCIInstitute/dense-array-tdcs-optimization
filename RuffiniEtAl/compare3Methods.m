function [ours, detal, retal, stats] = compare3Methods(sqrtQ_weighted, v_weighted, ind, w, sqrtQ, G, T, field, roi, vole)

[detal.Solution, detal.Fval, detal.tdV] = weightedLeastSquaresL1ConstraintByParraEtAl(sqrtQ_weighted,v_weighted,ind);

pBrain = norm(cell2mat(sqrtQ)*detal.Solution)^2;
tot = (norm(detal.Solution,1)+abs(sum(detal.Solution)))/2;

[retal.Solution, retal.Fval, retal.dV] = weightedLeastSquaresTotalIndividual(Q, b, tot, ind);

wScale = norm(w);
stats.pBrain = pBrain;
stats.tTot = tot;
[ours.Solution, ours.Fval, ours.dV] = optimizationUsingCvxToolbox(w/wScale,sqrtQ,tot,ind,pBrain);
stats.oTot = ( norm(ours.Solution,1) + abs(sum(ours.Solution))) /2;
ours.Fval= ours.Fval*wScale;
dVect = G * (T*ours.Solution);
currentIntensity = sqrt(sum(reshape(dVect.*dVect,3,[])));
clear dVect;
save('ourSolution','currentIntensity','ours');

stats.ocmax = max(currentIntensity(field~=8));
stats.ocmaxb = max(currentIntensity(field==4 | field ==5));
stats.ocaveb = sum(currentIntensity(field==4 | field== 5) .* vole(field ==4 | field==5))/sum(vole(field==4 | field ==5));
stats.ocmedb = median(currentIntensity(field==4 | field ==5));
stats.ocaveROI = sum(currentIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1));


dVect = G * (T*detal.Solution);
currentIntensity = sqrt(sum(reshape(dVect.*dVect,3,[])));
clear dVect;
save('dSolution','currentIntensity','theirs');

stats.dcmax = max(currentIntensity(field~=8));
stats.dcmaxb = max(currentIntensity(field==4 | field ==5));
stats.dcaveb = sum(currentIntensity(field==4 | field== 5) .* vole(field ==4 | field==5))/sum(vole(field==4 | field ==5));
stats.dcmedb = median(currentIntensity(field==4 | field ==5));
stats.dcaveROI = sum(currentIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1));


dVect = G * (T*retal.Solution);
currentIntensity = sqrt(sum(reshape(dVect.*dVect,3,[])));
clear dVect;
save('rSolution','currentIntensity','ours');

stats.rcmax = max(currentIntensity(field~=8));
stats.rcmaxb = max(currentIntensity(field==4 | field ==5));
stats.rcaveb = sum(currentIntensity(field==4 | field== 5) .* vole(field ==4 | field==5))/sum(vole(field==4 | field ==5));
stats.rcmedb = median(currentIntensity(field==4 | field ==5));
stats.rcaveROI = sum(currentIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1));

save('currentStats','theirs','ours','stats');