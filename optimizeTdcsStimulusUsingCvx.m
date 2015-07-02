% ==============================================================================
% Finds the optimal electrode current stimulus pattern for HD-tDCS using
% cvx toolkit
%
%
% Synopsis: StimulusPattern = optimizeTdcsStimulusUsingCvx(w, Q, SafetyBound)
%           
%
% Input:    w           =   linear objective function coefficients
%           Q           =   cell of matrices for power constraints on
%                           avoid regions.  
%           SafetyBound =   safety constraint bounds
%           	.tot    =   total applied current bound
%               .ind    =   individual electrode current bounds
%               .pow    =   power in the avoid regions constraint bound 
%
% Output:   StimulusPattern     =   optimal stimulus pattern
%               .currentArray   =   electrode currents, reference at the end of
%                                   the array by default.
%               .fval           =   directional current density in the ROI
%               .dualVariables  =   dual variables for the constraints

% Written by: Seyhmus Guler, 5/1/15
% Notes:
% ==============================================================================
function StimulusPattern = optimizeTdcsStimulusUsingCvx(w,Q,SafetyBound)

tic;
L = numel(w); %Number of (free) electrodes
pp = numel(Q); % Number of power constraints
totalBound = SafetyBound.tot;
ind = SafetyBound.ind;
pmax = SafetyBound.pow;

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
totConstV : norm(y,1) <= 2*totalBound; %#ok<*VUNUS> %total current bound
indConstLB : ind(:,1) <= y; %lower bound on individual currents
indConstUB : y <= ind(:,2); %#ok<*BDSCI> %upper bound on individual currents
%power in the avoid regions is bounded:
for i = 1:pp
    powConstV{i} : norm(Q{i}*x) <= sqrt(pmax(i)); 
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
    totConstV : norm(y,1) <= 2*totalBound; %total current bound
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
    totConstV : norm(y,1) <= 2*totalBound; %total current bound
    indConstLB : ind(:,1) <= y; %lower bound on individual currents
    indConstUB : y <= ind(:,2); %upper bound on individual currents
    for i = 1:pp
        powConstV{i} : norm(Q{i}*x) <= sqrt(pmax(i)); %power in the avoid regions is bounded
    end
    cvx_end
end

StimulusPattern.currentArray = x;
StimulusPattern.fval = cvx_optval;
StimulusPattern.dualVariables = [totConstV; indConstLB; indConstUB; cell2mat(powConstV)]; 

fprintf('%s%f%s\n','Electrode current optimization is solved in ',toc,' seconds.');
end
