function [theirs,ours,stats] = compareMethods(sqrtQ_weighted, v_weighted, ind, w, sqrtQ, G, T, field, roi, vole)

[theirs.Solution, theirs.Fval, theirs.tdV] = weightedLeastSquaresL1ConstraintByParraEtAl(sqrtQ_weighted,v_weighted,ind);

pBrain = norm(sqrtQ*electrodeCurrent)^2;
tot = (norm(electrodeCurrent,1)+abs(sum(electrodeCurrent)));
wScale = norm(w);

[ours.Solution, ours.Fval, ours.dV] = optimizationUsingCvxToolbox(w/wScale,sqrtQ,tot,ind,pBrain);

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


dVect = G * (T*theirs.Solution);
currentIntensity = sqrt(sum(reshape(dVect.*dVect,3,[])));
clear dVect;
save('theirSolution','currentIntensity','theirs');

stats.tcmax = max(currentIntensity(field~=8));
stats.tcmaxb = max(currentIntensity(field==4 | field ==5));
stats.tcaveb = sum(currentIntensity(field==4 | field== 5) .* vole(field ==4 | field==5))/sum(vole(field==4 | field ==5));
stats.tcmedb = median(currentIntensity(field==4 | field ==5));
stats.tcaveROI = sum(currentIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1));

save('currentStats','t','o','stats');






