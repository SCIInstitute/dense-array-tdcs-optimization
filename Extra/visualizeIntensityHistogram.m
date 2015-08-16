function visualizeIntensityHistogram(intensity, field, roi, vole)
if size(roi,2) ~= numel(field)
    roi = roi';
end
if size(vole,2) ~= numel(field)
    vole = vole';
end

NBARS = 40;
logintensity = log10(intensity);

f = cell(9,1);
x = cell(9,1);
v = zeros(9,1);
for i=1:9
    if i==9
        idx = roi;
        nbars = NBARS;
    else
        if i==4
            idx = field == i & roi ~= 1;
            nbars = NBARS;
        else
            idx = field ==i;
            nbars = NBARS;
        end
    end
    [~, x{i}] = hist(logintensity(idx), nbars);
    xdist(i) = x{i}(2) - x{i}(1);
    [~, bins] = histc(logintensity(idx), [-inf x{i}(2:end-1) inf]);
    v(i) = sum(vole(idx));
    pv = vole(idx);
    for j = 1:nbars
        f{i}(j) = sum(pv(bins == j-1))/v(i);
    end
end

% range = max(logintensity(roi)) - min(logintensity(roi));

a = figure;
stairs(x{1},f{1},':k','LineWidth',1.5,'Marker','none');
hold on;
stairs(x{2},f{2},'-.m','LineWidth',1.5,'Marker','.');
stairs(x{3},f{3},'--g','LineWidth',1.5);
stairs(x{4},f{4},'b','LineWidth',2,'Marker','d','MarkerSize',3);
stairs(x{9},f{9},'r','LineWidth',2);
xlim([-4 1])
ylim([0 0.2])

%set(gca,'XTick',x)
%set(gca,'XTickLabel',sprintf('%3.4f|',x))
set(gca,'YTick',[0 0.05 0.10 0.15 0.20])
%set(gca,'YTickLabel',sprintf('%1.2f|',y))

title('PC','FontSize',14);
xlabel('log_{10}(Current Intensity) [A/m^2]','FontSize',14);
ylabel('Normalized volume','FontSize',14);

aa = legend('Skin','Skull','CSF','GM','ROI');
set(aa,'FontSize',14);
set(gca,'FontSize',14);

pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)
%%



%xlabel, ylabel and title
set(get(gca,'xlabel'),'FontSize', 18, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', 18, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', 18, 'FontWeight', 'Bold');

% box and ticks
set(gca,'box','off');
set(gca,'TickDir','out')

%Axis fontsize, weight, line properties
set(gca,'FontSize',16);
set(gca,'FontWeight','Bold');
set(gca,'LineWidth',2);

%set the paper size, figure position 
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [8 8]);
set(gcf,'PaperPosition',[0.5 0.5 7 5])
set(gcf,'PaperPositionMode','Manual');
%set(gcf,'color','w');

%print the figure as pdf
print(gcf, '-dpdf', '-r150', 'MFC.pdf');


%saveas(a,['roi' num2str(r) '\Hist' num2str(d) num2str(s) num2str(ss) num2str(p) '.png'],'png');
%close(a);


