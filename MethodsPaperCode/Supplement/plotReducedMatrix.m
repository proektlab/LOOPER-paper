

finalStream = saveData.FinalStream;

[pcaBasis, ~] = pca(finalStream, 'NumComponents', 3);

figureHandle = figure(10);
figureHandle.Renderer='Painters';
clf;
plot3(finalStream*pcaBasis(:,1), finalStream*pcaBasis(:,2), finalStream*pcaBasis(:,3));
scatter3(finalStream*pcaBasis(:,1), finalStream*pcaBasis(:,2), finalStream*pcaBasis(:,3));

