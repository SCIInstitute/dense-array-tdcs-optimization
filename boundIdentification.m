function [stats,fval,dualVariables] = boundIdentification(w,sqrtQ,G,T,field,roi,vole)
%ANALYSIS OF THE CONSTRAINT EFFECTS ON THE OBJECTIVE FUNCTION AND OTHER
%STATISTICS
%Written by: Seyhmus Guler 3/9/14

%INPUTS
    %w: Linear coefficients for the objective function
    %sqrtQ: Quadratic matrices representing avoid region power constraints.
    %G: Matrix representing -conductivity * Gradient() operator. If
    %multiplied with potential vector, results in current density in the 
    %domain
    %T: Lead field matrix. Matrix linking electrode currents to potential
    %in the domain. 
    %field: field vector showing the element labels.
%OUTPUTS
    %stats: Statistics about the current intensity found by optimization
    %fval: objective function values for optimized electrode currents
    %actSet: Active constraints at the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = numel(w);

%bound
smax = 1;
smaxi = 0.06:0.01:0.40;
pmax = 10.^(3.4:-0.1:1);

%initialization
fval = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
stats.cmax = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
stats.cmaxb = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
stats.caveb = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
stats.cmedb = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
stats.caveROI = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
dualVariables = zeros(size(smax,2),size(smaxi,2),size(pmax,2),3+2*L+numel(sqrtQ));
wScale = norm(w);

for ss = 1:size(smax,2)
    for si = 1: size(smaxi,2)
        for pi = 1:size(pmax,2)
            [x,f,dualVariables(ss,si,pi,:)] = optimizationUsingCvxToolbox(w/wScale,sqrtQ,smax(ss),smaxi(:,si),pmax(:,pi)); 
            fval(ss,si,pi) = f*wScale;
            dVect = G * (T*x);
            currentIntensity = sqrt(sum(reshape(dVect.*dVect,3,[])));
            clear dVect;
            stats.cmax(ss,si,pi) = max(currentIntensity(field~=8));
            stats.cmaxb(ss,si,pi) = max(currentIntensity(field==4 | field ==5));
            stats.caveb(ss,si,pi) = sum(currentIntensity(field==4 | field== 5) .* vole(field ==4 | field==5))/sum(vole(field==4 | field ==5));
            stats.cmedb(ss,si,pi) = median(currentIntensity(field==4 | field ==5));
            stats.caveROI(ss,si,pi) = sum(currentIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1)); 
        end
    end
end

save([pwd 'boundInformation.mat'],'w','sqrtQ','field','stats','fval','dualVariables','smax','smaxi','pmax'); 


