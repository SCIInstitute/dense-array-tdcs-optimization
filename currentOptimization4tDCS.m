function currentOptimization4tDCS(mesh,T,conductivity,directions,ROIs,avoidRegions,totCurBound,indCurBound,powerBound,G,V)
%CALCULATES THE BEST ELECTRODE CURRENT STIMULUS CONFIGURATION FOR 
%TARGETED AND DIRECTIONAL TDCS. 
%
%Maximizes the directional current density in the ROIs while imposing 
%constraints on the total injected current, each individual electrode 
%current and elecrical power in the avoidance areas.
%
%Written by: Seyhmus Guler, Revisited: Moritz Dannhauer
%Last edit: 10/7/13 by Guler,S
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS: 
    %mesh: Tetrahedral realistic head mesh. A structure consisting of 
        %'cell','node','field' variables.
    %T: transfer matrix found by solving the forward problem. It is the
        %mapping from electrode currents to potential field at the nodes.
        %size: (#mesh nodes) x (#electrodes)
    %conductivity: volume conductor model. matrix of size (9 x #elements)
    %directions: desired directions. It is a cell array of cell arrays. 
        %Each cell consists of desired direction(s) for corresponding roi. 
    %ROIs: region(s) of interest. Each row represents a different ROI.
    %avoidRegions: avoidance region(s). Each row represents a avoidance
        %region
    %totCurBound: the upper bound on the total current entering the head.
    %indCurBound: uppers bound(s) on the individual electrode currents.
    %powerBound: upper bound(s) on the electrical power in avoidRegions.
%OUTPUTS:
    %The best electrode current array, the corresponding potential field, 
    % current density field and current intensity field for each 
    % optimization scheme are saved in the path so there are no outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str2 = char(date); 
if(nargin <= 9)
    [G,V] = mappingFromNodePotentialsToCurrentDensity(mesh,conductivity);
    save -v7.3 G.mat G
    save -v7.3 V.mat V
end

M = size(mesh.cell,2);
if(size(avoidRegions,2) ~= M)
    avoidRegions = avoidRegions';
end
if(size(ROIs,2) ~= M)
    ROIs = ROIs';
end
if(size(powerBound,1) ~= size(avoidRegions,1))
    error('mismatch in matrix sizes');
end

if (size(G,2) ~= size(T,1) || ~isequal(size(ROIs,2),size(avoidRegions,2),numel(V)))
    error('mismatch in matrix sizes');
end
for r = 1:10
    mkdir([str2 '\roi' num2str(r)]);
    ROIr = ROIs(r,:);
    avoidRegions(:,ROIr ==1) = 0;
    desiredDir4ROIr = directions{r};
    for d = 1:numel(desiredDir4ROIr)
        desiredDirection = desiredDir4ROIr{d};
        [w,Q] = wAndQCalculation(ROIr,avoidRegions,desiredDirection,T,G,V);
        save([pwd '/' str2 '/roi' num2str(r) '/wMatrix' num2str(d) '.mat'],'w');
        wScale = norm(w,2);
        w = w/wScale;
        for ss = 1: numel(totCurBound)
            for si = 1:size(indCurBound,2)
                for pi = 1:size(powerBound,2)
                    [currentArray,fval] = optimizationUsingCvxToolbox(w,Q,totCurBound(ss),indCurBound(:,si),powerBound(:,pi));
                    %potential = T*currentArray;
                    %dVect = JfromU * potential;
                    %potential = potential';
                    %currentIntensity = sqrt(sum(reshape(dVect.*dVect,3,[])));
                    %currentDensity = reshape(dVect,3,[]);
                    save([pwd '/' str2  '/roi' num2str(r) '/elecCurrent'  num2str(d) num2str(ss) num2str(si) num2str(pi) '.mat'],'currentArray');
                    %save([pwd '/' str2  '/roi' num2str(r) '/potential'  num2str(d) num2str(ss) num2str(si) num2str(pi) '.mat'],'potential');
                    %save([pwd '/' str2  '/roi' num2str(r) '/intensity'  num2str(d) num2str(ss) num2str(si) num2str(pi) '.mat'],'currentIntensity');
                    %save([pwd '/' str2  '/roi' num2str(r) '/density'  num2str(d) num2str(ss) num2str(si) num2str(pi) '.mat'],'currentDensity');
                end
            end
        end
    end
end