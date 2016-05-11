function currentArray = oneSourceOptimization(w,Q,totCurBound,indCurBound,powerBound,Te)
%FINDS THE BEST ELECTRODE CURRENT STIMULUS CONFIGURATION USING CVX
%ASSUMES THERE IS ONLY ONE CURRENT SOURCE
%
%Written by: Seyhmus Guler
%Last edit: 12/29/13 by Guler,S
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
    %w: Linear weights for the objective function.
    %Q: Quadratic constraint matrices for each avoid region
    %totCurBound: the total current bound
    %indCurBound: individual electrode current bounds
    %powerBound: the bound on the electrical power in avoidance regions
    %Te: transfer matrix for the electrode potentials. Used in the 
        %equalities in the fewer source cases
%OUTPUTS:
    %x: best solution for electrode current array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
L = numel(w); %number of electrodes

cvx_begin quiet
variable x(L)
minimize w*x
subject to
-indCurBound <= x <= indCurBound; 
norm(x,1) <= 2*totCurBound;
sum(x) == 0; %#ok<EQEFF>
for i = 1:numel(Q)
quad_form(x,Q{i}) <= powerBound(i);
end
cvx_end
fprintf('%s%f%s\n','optimization is solved in ',toc,' seconds.');
currentArray.origCurrent = x;
currentArray.origPot = Te * x;
currentArray.origObj = cvx_optval;
fprintf('%s\n','checking the closest one current source solution...');

minX = min(x);
maxX = max(x);

nIdx = find(x < 0.95* minX);
pIdx = find(x > 0.95* maxX);
zeroIdx = find(abs(x) < 0.001* min(maxX,-minX));

unknownNegatives = find(x > 0.95*minX & x < -0.001* min(maxX,-minX));
unknownPositives = find(x < 0.95*maxX & x > 0.001* min(maxX,-minX));

nN = numel(unknownNegatives);
nP = numel(unknownPositives);
if (numel(zeroIdx) + nN + nP + numel(nIdx) + numel(pIdx) ~= L)
    error('numbers mismatch');
end
objFnc = zeros(1,2^(nN+nP));
nStats = cell(1,2^(nN+nP));

for i = 1:1
    elecAssign = dec2bin(i-1,nN+nP);
    elecAssign = fliplr(elecAssign);
    sourceAssign = false(1,nP);
    sinkAssign = false(1,nN);
    for k = 1:nN+nP
        if(k>nP)
            sinkAssign(k-nP) = strncmp(elecAssign(k),'1',1);
        else
            sourceAssign(k) = strncmp(elecAssign(k),'1',1);
        end
    end
    sourceIndices = [pIdx; unknownPositives(sourceAssign)];
    sinkIndices = [nIdx; unknownNegatives(sinkAssign)];
    zeroIndices = [zeroIdx; unknownPositives(~sourceAssign); unknownNegatives(~sinkAssign)];
    disp(x(sourceIndices));
    disp(x(sinkIndices));
    nStats{i}.sinks = sinkIndices;
    nStats{i}.sources = sourceIndices;
    nZeros = numel(zeroIndices);
    
    cvx_begin
    variable u(L-nZeros)
    expression x(L)
    x(zeroIndices) = 0;
    ii = 1:L;
    ii(zeroIndices) = [];
    x(ii) = u;
    cvx_precision low;
    cvx_solver Sedumi
    minimize w*x
    subject to
    -indCurBound <= x <= indCurBound; %#ok<VUNUS>
    norm(x,1) <= 2*totCurBound; %#ok<VUNUS>
     sum(x) == 0; %#ok<EQEFF>
    for ip = 1:numel(Q)
        quad_form(x,Q{ip}) <= powerBound(ip); %#ok<VUNUS>
    end
    for k = 2:numel(sourceIndices)
        Te(sourceIndices(k),:)*x == Te(sourceIndices(1),:)*x; %#ok<EQEFF>
    end
    for ks = 2:numel(sinkIndices)
        Te(sinkIndices(ks),:)*x == Te(sinkIndices(1),:)*x; %#ok<EQEFF>
    end
    cvx_end
    nStats{i}.current = x;
    objFnc(i) = cvx_optval;
    if(mod(i,100)==0)
         disp(num2str(i));
    end
end

bestConfig = find(objFnc == min(objFnc));

currentArray.bestConfig = bestConfig;
currentArray.limitedObj = objFnc(bestConfig);
currentArray.limitedPot = Te*nStats{bestConfig}.current;
currentArray.LimitedCurrent = nStats{bestConfig}.current;
currentArray.obj = objFnc;
currentArray.nStats = nStats;
end


