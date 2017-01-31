function [h]=plotingsolution_1(G,W,name,x)
% [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
% clf;
% Plot the grid
%h1=subplot(2,2,n);
plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
h=plotCellData(G,x);

% Plot the wells
% plotWell(G, W);
view(0,90)
axis equal tight; colormap(jet(128));
title(name);
