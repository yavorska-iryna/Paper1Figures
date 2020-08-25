% make_figure_4bc

cd ('/Users/mikewehr/Resilio Sync/Paper1Figures/figure4/matlab figs/')
open RS_FS_LaserEffect_by_Layer.fig
subplot(211)
set(gca, 'fontsize', 18)
set(gca, 'ytick', [-.4:.2:.2])
subplot(212)
set(gca, 'fontsize', 18)
set(gca, 'ytick', [[-.4:.2:.4]])


f=findobj('type', 'line')
set(f, 'linewidth', 2)
f=findobj('marker', '*')
set(f, 'marker', 'o', 'color', 'g', 'markersize', 12, 'LineWidth', 2)
f=findobj('marker', 'o')
set(f, 'marker', 'o', 'markersize', 12, 'LineWidth', 2)

print -dpsc2 fig4bc -bestfit