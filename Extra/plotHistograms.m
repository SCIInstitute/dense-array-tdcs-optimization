function plotHistograms(legendOpts)
load('C:\Users\Seyhmus\Desktop\PROJECTS\HD_TDCS_Optimization\ProcessedData\refinedHeadMesh\field.mat');
load('C:\Users\Seyhmus\Desktop\PROJECTS\HD_TDCS_Optimization\ProcessedData\refinedHeadMesh\vole.mat');
load('C:\Users\Seyhmus\Desktop\PROJECTS\HD_TDCS_Optimization\ProcessedData\refinedHeadMesh\gmroi.mat');

for i = 1:4
    load(['C:\Users\Seyhmus\Desktop\PROJECTS\HD_TDCS_Optimization\ProcessedData\refinedHeadMesh\06-Jul-2015-NoEyeConstraint\roi' num2str(i) '\intensity1111.mat']);
    switch i
        case 1
            roiName = 'MFC';
        case 2
            roiName = 'ACC';
        case 3
            roiName = 'PHCG';
        case 4
            roiName = 'PC';
    end
    visualizeIntensityHistogram(currentIntensity,field,gmroi(:,i),vole);
    title(roiName,'FontSize',14);
    if nargin >= 1 && legendOpts == 1
        aa = legend('Scalp','Skull','CSF','GM','ROI');
        set(aa,'FontSize',14);
    end
    print(gcf, '-dpdf', '-r150', [roiName '.pdf']);
end