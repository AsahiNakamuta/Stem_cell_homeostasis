duration_heatmap = importdata("duration_heatmap.mat");
shannon_heatmap = importdata("shannon_heatmap.mat");


x_heatmap = ([0 1000]);
y_heatmap = ([0 1]);
x_contour = 10:10:1000;
y_contour = 0.01:0.01:1;
figure
hold on
imagesc(x_heatmap,y_heatmap,shannon_heatmap)
[M,c] = contour(x_contour,y_contour,shannon_heatmap, [3.64,3.64], 'r');
c.LineWidth = 3;
colorbar
hold off
saveas(gcf,'Fig6B.png')

figure
imagesc(x_heatmap,y_heatmap,duration_heatmap)
colorbar
saveas(gcf,'Fig6C.png')