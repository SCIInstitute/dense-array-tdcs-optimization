function [stats,fval,dualVariables] = boundIdentification(w,sqrtQ,G,T,field,roi,vole)
%ANALYSIS OF THE CONSTRAINT EFFECTS ON THE OBJECTIVE FUNCTION AND OTHER
%STATISTICS
%Written by: Seyhmus Guler 3/9/14
%Last edit: 4/29/15

%INPUTS
    %w: Linear coefficients for the objective function
    %sqrtQ: Quadratic matrices representing avoid region power constraints.
    %G: Matrix representing -conductivity * Gradient() operator. If
    %multiplied with potential vector, results in current density in the 
    %domain
    %T: Lead field matrix. Matrix linking electrode currents to potential
    %in the domain. 
    %field: field vector showing the element labels.
    %roi: roi labels
    %vole: element volumes
%OUTPUTS
    %stats: Statistics about the current intensity found by optimization
    %fval: objective function values for optimized electrode currents
    %actSet: Active constraints at the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = numel(w);

%bound
%smax = 1;
%smaxi = 0.50:0.005:0.80;
%pmax = 10.^(3.4:-0.05:1);

%bounds
tot = 2;
ind = 0.10:0.02:1;
pmax = 10.^(5:-0.1:1);

%initialization
fval = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
stats.cmax = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
stats.cmaxb = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
stats.caveb = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
stats.cmedb = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
stats.caveROI = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
stats.cmedROI = zeros(size(smax,2),size(smaxi,2),size(pmax,2));
dualVariables = zeros(size(smax,2),size(smaxi,2),size(pmax,2),3+2*L+numel(sqrtQ));
wScale = norm(w);

for ss = 1:size(smax,2)
    for si = 1: size(smaxi,2)
        for pi = 1:size(pmax,2)
            [x,f,dualVariables(ss,si,pi,:)] = optimizeBound(w/wScale,sqrtQ,smax(ss),smaxi(:,si),pmax(:,pi)); 
            fval(ss,si,pi) = f*wScale;
            dVect = G * (T*x);
            currentIntensity = sqrt(sum(reshape(dVect.*dVect,3,[])));
            clear dVect;
            stats.cmax(ss,si,pi) = max(currentIntensity(field~=8));
            stats.cmaxb(ss,si,pi) = max(currentIntensity(field==4 | field ==5));
            stats.caveb(ss,si,pi) = sum(currentIntensity(field==4 | field== 5) .* vole(field ==4 | field==5))/sum(vole(field==4 | field ==5));
            stats.cmedb(ss,si,pi) = median(currentIntensity(field==4 | field ==5));
            stats.caveROI(ss,si,pi) = sum(currentIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1));
            stats.vmedROI(ss,si,pi) = sum(currentIntensity(roi==1) .* vole(roi==1))/sum(vole(roi==1));
        end
    end
end

save([pwd '/boundInformation2.mat'],'w','sqrtQ','field','stats','fval','dualVariables','smax','smaxi','pmax'); 
end

function [electrodeCurrents,fval,dualVariables] = optimizeBound(w,sqrtQ,tot,ind,pmax)
%% FINDS THE BEST ELECTRODE CURRENT STIMULUS CONFIGURATION USING CVX TOOLKIT
% IN ADDITION IT FINDS THE FUNCTION VALUE AT SOLUTION AND DUAL VARIABLES (FOR ACTIVE SET 
% INVESTIGATIN PURPOSES) IF REQUIRED.
%
%Written by: Seyhmus Guler, 4/29/14
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
%w: Linear weights for the objective function.
%sqrtQ: Square roots of matrices for quadratic constraint.
%tot: The total current bound
%ind: Individual electrode current bound.
%pmax: The bound on the electrical power in avoidance regions
%
%OUTPUTS:
%electrodeCurrents: best solution for electrode current array.
%fval: Objective function, the current in the ROI along the desired
    %direction.
%DualVariables: dual variables for the constraints at solution. The
    % variables corresponding to first total bound, second individual 
    % lower bound, third individual upper bound and lastly power bound 
    % is stored in an array of size #constraints x 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
L = numel(w); %Number of electrodes
pp = numel(sqrtQ); % Number of power constraints

if size(ind,1) == L %In case reference electrode bound is not defined
    ind(L+1,:) = ind(end,:); %Make ref elec have the same const as last one
end

if size(ind,2) == 1 %Lower bound = - Upper bound
    ind = [-ind ind];
end

cvx_begin quiet
cvx_precision high
cvx_solver sedumi

%optimization variable
variable x(L)

%expression of electrode current array, reference added
expression y(L+1)
y = [x; -sum(x)]; %reference elec current = - (sum of other electrodes)

%dual variables to check active/passive constraints
dual variable totConstV
dual variable indConstLB
dual variable indConstUB
dual variable powConstV{pp}

%maximize current density in the ROI along the desired direction
maximize w*x

%subject to total and ind electrode currents, and the power constraints
subject to
totConstV : norm(y,1) <= 2*tot; %#ok<*VUNUS> %total current bound
indConstLB : ind(:,1) <= y; %lower bound on individual currents
indConstUB : y <= ind(:,2); %#ok<*BDSCI> %upper bound on individual currents
%power in the avoid regions is bounded:
for i = 1:pp
    powConstV{i} : norm(sqrtQ{i}*x) <= sqrt(pmax(i)); 
end
cvx_end

if ~strcmp(cvx_status,'Solved')
    fprintf('%s\n','No solution with high precision, trying lower precision settings');
    cvx_begin quiet
    cvx_precision medium
    cvx_solver sedumi
    
    variable x(L)
   
    expression y(L+1)
    y = [x; -sum(x)];
    
    dual variable totConstV
    dual variable indConstLB
    dual variable indConstUB
    dual variable powConstV{pp}
    
    %maximize current density in the ROI along the desired direction
    maximize w*x
    
    %subject to total curr., ind elec currents and the power constraints
    subject to
    totConstV : norm(y,1) <= 2*tot; %total current bound
    indConstLB : ind(:,1) <= y; %lower bound on individual currents
    indConstUB : y <= ind(:,2); %upper bound on individual currents
    for i = 1:pp
        powConstV{i} : norm(sqrtQ{i} * x) <= sqrt(pmax(i)); %power in the avoid regions is bounded
    end
    cvx_end
end

if ~strcmp(cvx_status,'Solved')
    fprintf('%s\n','Still no solution, trying SDPT3 as solver');
    
    cvx_begin quiet
    cvx_solver SDPT3
    variable x(L)
   
    expression y(L+1)
    y = [x; -sum(x)];
    
    dual variable totConstV
    dual variable indConstLB
    dual variable indConstUB
    dual variable powConstV{pp}
    
    %maximize current density in the ROI along the desired direction
    maximize w*x
    
    %subject to total curr., ind elec currents and the power constraints
    subject to
    totConstV : norm(y,1) <= 2*tot; %total current bound
    indConstLB : ind(:,1) <= y; %lower bound on individual currents
    indConstUB : y <= ind(:,2); %upper bound on individual currents
    for i = 1:pp
        powConstV{i} : norm(sqrtQ{i}*x) <= sqrt(pmax(i)); %power in the avoid regions is bounded
    end
    cvx_end
end

electrodeCurrents = x;
fval = cvx_optval;
dualVariables = [totConstV; indConstLB; indConstUB; cell2mat(powConstV)]; 

fprintf('%s%f%s\n','Electrode current optimization is solved in ',toc,' seconds.');
end



