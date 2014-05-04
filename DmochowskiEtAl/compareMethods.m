function [t,o,stats] = compareMethods(sqrtQ_weighted,v_weighted,Smaxi,w,sqrtQ,G,T,field,roi,vole)

[t.Solution,t.Fval,t.tdV] = weightedLeastSquaresL1ConstraintByParraEtAl(sqrtQ_weighted,v_weighted,Smaxi);

pBrain = norm(sqrtQ*electrodeCurrent)^2;
tot = (norm(electrodeCurrent,1)+abs(sum(electrodeCurrent)));
wScale = norm(w);

[o.Solution,o.Fval,o.dV] = optimizationUsingCvxToolbox(w/wScale,sqrtQ,tot,Smaxi,pBrain);

o.Fval= o.Fval*wScale;
dVect = G * (T*o.Solution);
currentIntensity = sqrt(sum(reshape(dVect.*dVect,3,[])));
clear dVect;
save('ourSolution','currentIntensity','o');

stats.ocmax = max(currentIntensity(field~=8));
stats.ocmaxb = max(currentIntensity(field==4 | field ==5));
stats.ocaveb = sum(currentIntensity(field==4 | field== 5) .* vole(field ==4 | field==5))/sum(vole(field==4 | field ==5));
stats.ocmedb = median(currentIntensity(field==4 | field ==5));
stats.ocaveROI = sum(currentIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1));


dVect = G * (T*t.Solution);
currentIntensity = sqrt(sum(reshape(dVect.*dVect,3,[])));
clear dVect;
save('theirSolution','currentIntensity','t');

stats.tcmax = max(currentIntensity(field~=8));
stats.tcmaxb = max(currentIntensity(field==4 | field ==5));
stats.tcaveb = sum(currentIntensity(field==4 | field== 5) .* vole(field ==4 | field==5))/sum(vole(field==4 | field ==5));
stats.tcmedb = median(currentIntensity(field==4 | field ==5));
stats.tcaveROI = sum(currentIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1));

save('currentStats','t','o','stats');






