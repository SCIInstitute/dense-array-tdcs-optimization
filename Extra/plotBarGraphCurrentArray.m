function plotBarGraphCurrentArray(eca,threshold,electrodeLabels)
L = size(eca,1);
if nargin <= 2
    electrodeLabels = 1:L;
end

idx = find(eca > threshold | eca < -threshold);
idx = unique([idx(idx<=L);idx(idx>L & idx <= 2*L)-L; idx(idx>2*L)-2*L]); 

b = bar(eca(idx,:));
b(1).FaceColor = 'b';
b(1).EdgeColor = 'k';
b(2).FaceColor = 'r';
b(2).EdgeColor = 'k';
%b(3).FaceColor = 'c';
%b(3).EdgeColor = 'k';
b(1).LineWidth = 0.05;
b(2).LineWidth = 0.05;
%b(3).LineWidth = 0.05;

title('Comparison of two optimal current patterns');
b_legend = legend('Ruffini et al.','Guler et al.');
set(b_legend,'FontSize',10);
set(b_legend,'Location','best');
ylabel('Electrode current [mA]');
te = text(numel(idx)-1,0.08,'Electrodes'); 
%te = text([numel(idx)+1 numel(idx)+1],[0.04 -0.04],{'Electrode' 'Incex'});
set(te,'FontSize',12);

set(gca,'box','off');
set(gca,'FontSize',12);
set(gca,'color','none');
set(gca,'Xtick',[]);
%set(gca,'Ycolor',get(gca,'color'));
set(gca,'Xcolor',get(gca,'color'));
xlim([0 numel(idx)+1]); 

%print(gcf,'-dpdf','-r150','comparison3Methods.pdf');

%ii = find(eca > threshold | eca < -threshold);
%pidx = find( ca > 0.01 | cagul >0.01);
%[~,locp] = ismember(pidx,ii);
%tp = text(locp,-0.01*ones(1,numel(locp)),electrodeAxis(ii(locp))'); 
%set(tp,'HorizontalAlignment','center');
%set(tp,'FontSize',14);