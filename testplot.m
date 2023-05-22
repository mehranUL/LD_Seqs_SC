%heatmap(sim_faure,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Faure Sequence','XLabel','Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f','');
% heatmap(sim_faure_shift,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Shifted Faure Sequence','XLabel','# of shifts at each Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_gold,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Gold Sequence','XLabel','Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_gold_shift,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Shifted Gold Sequence','XLabel','# of shifts at each Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_hadamard,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Hadamard Sequence','XLabel','Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_hadamard_shift,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Shifted Hadamard Sequence','XLabel','# of shifts at each Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_halton,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Halton Sequence','XLabel','Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_halton_shift,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Shifted Halton Sequence','XLabel','# of shifts at each Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_hammersley,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Hammersley Sequence','XLabel','Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_hammersley_shift,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Shifted Hammersley Sequence','XLabel','# of shifts at each Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_nd,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Niederreiter Sequence','XLabel','Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_nd_shift,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Shifted Niederreiter Sequence','XLabel','# of shifts at each Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_R2,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of R2 Sequence','XLabel','Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_R2_shift,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Shifted R2 Sequence','XLabel','# of shifts at each Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_sobol,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Sobol Sequence','XLabel','Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_sobol_shift,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Shifted Sobol Sequence','XLabel','# of shifts at each Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_vd,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of VanDerCorput Sequence','XLabel','Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_vd_shift,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Shifted VanDerCorput Sequence','XLabel','# of shifts at each Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_weyl,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Weyl Sequence','XLabel','Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');
% heatmap(sim_weyl_shift,'FontSize',25,'FontName', 'consolas','Title','Cosine Similarity of Shifted Weyl Sequence','XLabel','# of shifts at each Sequence Number','YLabel','Sequence Number','Colormap',colorcube,'CellLabelFormat','%.3f');

  hp = heatmap(MAE_ps,'FontSize',50,'FontName', 'consolas','Colormap',parula);
  % Removing the Data labels
  Ax = gca; 
  Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
  Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
%  hp.Title = 'Cosine Similarity of Faure';
 %title = addTitle(hp,'Sim_Faure','Color','red');

%  histogram(MAE_faure_add);
%  title('OR-based Stochastic Addition Absolute Error of Faure sequence')
%  fontname('consolas')




figure
%bir = bar(x,[g1_NO;g1_YES]);
h = heatmap(MAE_R2_add);
set(findall(gcf,'-property','FontSize'),'FontSize',18, 'FontName', 'consolas')
%set(gca,'')
% set(gca,'xtick',1:8,'xticklabel',longbar) 
%xlabel("Bit-stream Size");

%ylabel("MAE");
legend("With Recon.", "Without Recon.")
%bir(1).FaceColor = 'r';
%bir(2).FaceColor = 'g';
box off
% bir(1).EdgeColor = 'none';
% bir(2).EdgeColor = 'none';