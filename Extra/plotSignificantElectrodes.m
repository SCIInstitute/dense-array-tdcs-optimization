function plotSignificantElectrodes(eca)

sorted = sort(abs(eca),'descend');

h = plot(sorted(1:20,:));

h(1).LineWidth = 1;
h(2).LineWidth = 1;
h(3).LineWidth = 1;
set(h(1),'color','b');
set(h(2),'color','c');
set(h(3),'color','r');

h_legend = legend('Dmochowski et al.','Ruffini et al.','Guler et al.');
set(h_legend,'FontSize',10);
set(h_legend,'Location','best');

set(gca,'box','off');
set(gca,'color','none');
set(gca,'FontSize',12);
title('Highest current magnitudes in each optimal pattern'); 
ylabel('Electrode current magnitude [mA]');
xlabel('Electrodes');

print(gcf,'-dpdf','-r150','highest20magnitudes.pdf');