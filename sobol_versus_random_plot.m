
figure
D = 1000;
sobolseq2 = net(sobolset(1),(D));
plot(sobolseq2(:, 1),  '*', 'linewidth',3)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 
set(findall(gcf,'-property','FontSize'),'FontSize',28, 'FontName', 'consolas')


figure
r = rand(D,1);
plot(r(:, 1), 'r*', 'linewidth',3)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 
set(findall(gcf,'-property','FontSize'),'FontSize',28, 'FontName', 'consolas')