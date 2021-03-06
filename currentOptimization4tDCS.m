function currentOptimization4tDCS(mesh,T,conductivity,directions,ROIs,avoidRegions,tot,ind,pmax,G,V)
%CALCULATES THE BEST ELECTRODE CURRENT STIMULUS CONFIGURATION FOR 
%TARGETED AND DIRECTIONAL TDCS. 
%
%Maximizes the directional current density in the ROIs while imposing 
%constraints on the total injected current, each individual electrode 
%current and elecrical power in the avoidance areas.
%
%Written by: Seyhmus Guler, Revised: Moritz Dannhauer
%Last edit: 4/29/13 by Guler,S
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M := # Elements in the mesh.
% N := # Nodes in the mesh.
% L := # Electrodes (reference excluded so L = #electrodes on the head - 1).
%INPUTS: 
    %mesh: Tetrahedral realistic head mesh. A structure consisting of 
        %'cell','node','field' variables.
    %T: transfer matrix found by solving the forward problem. It is the
        %mapping from electrode currents to potential field at the nodes.
        %size: N x L
    %conductivity: volume conductor model. matrix of size M x 9
    %directions: desired directions. It is a cell array of cell arrays. 
        %Each cell consists of desired direction(s) for corresponding roi. 
    %ROIs: region(s) of interest. Each row represents a different ROI.
    %avoidRegions: avoidance region(s). Each row represents a avoidance
        %region
    %tot: the upper bound on the total current entering the head.
    %ind: uppers bound(s) on the individual electrode currents.
    %pmax: upper bound(s) on the electrical power in avoidRegions.
%OUTPUTS:
    %The best electrode current array, the corresponding potential field, 
    % current density field and current intensity field for each 
    % optimization scheme are saved in the path so there are no outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str2 = char(date); 
if(nargin <= 9)
    [G,V] = calculateGradientsIsotropicMesh(mesh,conductivity);
    save -v7.3 sigmaGradV.mat G
    save -v7.3 vole.mat V
end

M = size(mesh.cell,2);
if isempty(avoidRegions)
    avoidRegions = false(1,M);
    avoidRegions(mesh.field == 4 | mesh.field ==5) = 1;
end
if(size(avoidRegions,2) ~= M)
    avoidRegions = avoidRegions';
end
if(size(ROIs,2) ~= M)
    ROIs = ROIs';
end
if(size(pmax,1) ~= size(avoidRegions,1))
    error('mismatch in matrix sizes');
end

if (size(G,2) ~= size(T,1) || ~isequal(size(ROIs,2),size(avoidRegions,2),numel(V)))
    error('mismatch in matrix sizes');
end

for r = 1:size(ROIs,1)
    avoidRegionR = avoidRegions;
    mkdir([str2 '/roi' num2str(r)]);
    ROIr = ROIs(r,:);
    avoidRegionR(:,ROIr ==1) = 0;
    if isempty(directions)
        sfn = corticalSurfaceNormalDirection(mesh,ROIr);
        save('surfaceNormalInfo','sfn');
        desiredDir4ROIr{1} = sfn.ROIsurfaceNormalD; 
    else
    desiredDir4ROIr = directions{r};
    end
    for d = 1:numel(desiredDir4ROIr)
        desiredDirection = desiredDir4ROIr{d};
        [w,Q] = linearQuadraticCoefficientCalculation(ROIr,avoidRegionR,desiredDirection,T,G,V);
        save([pwd '/' str2 '/roi' num2str(r) '/variablesAndBounds' num2str(d) '.mat'],'w','Q','tot','ind','pmax','desiredDirection','ROIr');
        wScale = norm(w,2);
        w = w/wScale;
        for ss = 1: numel(tot)
            for si = 1:size(ind,2)
                for pi = 1:size(pmax,2)
                    [currentArray,fval,dualV] = optimizationUsingCvxToolbox(w,Q,tot(ss),ind(:,si),pmax(:,pi));
                    currentArrayReferenceAdded = [currentArray; -sum(currentArray)];
                    if ~isnan(fval)
                        fval = fval*wScale;
                        potential = T*currentArray;
                        dVect = G * potential;
                        potential = potential';
                        currentIntensity = sqrt(sum(reshape(dVect.*dVect,3,[])));
                        logIntensity = log2(currentIntensity);
                        currentDensity = reshape(dVect,3,[]);
                        save([pwd '/' str2  '/roi' num2str(r) '/fDual'  num2str(d) num2str(ss) num2str(si) num2str(pi) '.mat'],'fval','dualV');
                        save([pwd '/' str2  '/roi' num2str(r) '/elecCurrent'  num2str(d) num2str(ss) num2str(si) num2str(pi) '.mat'],'currentArrayReferenceAdded');
                        save([pwd '/' str2  '/roi' num2str(r) '/potential'  num2str(d) num2str(ss) num2str(si) num2str(pi) '.mat'],'potential');
                        save([pwd '/' str2  '/roi' num2str(r) '/intensity'  num2str(d) num2str(ss) num2str(si) num2str(pi) '.mat'],'currentIntensity');
                        save([pwd '/' str2  '/roi' num2str(r) '/logintensity'  num2str(d) num2str(ss) num2str(si) num2str(pi) '.mat'],'logIntensity');
                        save([pwd '/' str2  '/roi' num2str(r) '/density'  num2str(d) num2str(ss) num2str(si) num2str(pi) '.mat'],'currentDensity');
                    end
                end
            end
        end
    end
end