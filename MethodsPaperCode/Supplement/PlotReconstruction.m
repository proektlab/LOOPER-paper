

app.SavedData = saveData;

% finalStream = app.SavedData.FinalStream;
finalStream = reshape(finalTrials, size(finalTrials,1), size(finalTrials,2)*size(finalTrials,3));

stateMap = [];
newEmissions = [];
for i = 1:size(app.SavedData.BestLoopAssignments,1)
    thisIndicies = find(ismember(app.SavedData.BestStateMap,app.SavedData.BestLoopAssignments(i,:),'rows'));
    
    newEmissions(i,:) = mean(finalStream(thisIndicies,:),1);
    
    stateMap(thisIndicies) = i;
end

plotIndicies = 1:size(app.SavedData.RawData,1);

reconstruction = newEmissions(stateMap,plotIndicies);
originalData = finalStream(:,plotIndicies);

correlations = [];
labels = {};
for i = plotIndicies
    correlations(i) = corr(reconstruction(:,i), originalData(:,i));
    labels{i} = ['Dimension ' num2str(i)];
end

zOriginalData = zscore(originalData);
zReconstruction = zscore(reconstruction);


numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
times = 1:trialLength;
trialIndicies = repmat(times, [1, 6]);
trialIndicies = trialIndicies + kron((0:10:numTrial-1)*trialLength, ones(size(times)));

if length(plotIndicies) > 3
    figureHandle = figure(1);
    figureHandle.Renderer='Painters';
    clf;
    colormap(parula);
    h(1) = subplot(3,1,1);
    imagesc(zOriginalData(trialIndicies,:)');
    coloraxis = caxis;
    h(1).XTick = [];
    h(1).YTick = [];
    title('Original data');
    h(2) = subplot(3,1,2);
    imagesc(zReconstruction(trialIndicies,:)');
    caxis(coloraxis)
    h(2).XTick = [];
    h(2).YTick = [];
    linkaxes(h, 'x');
    title('Reconstructed data');
    h = subplot(3,1,3);
    bar(correlations)
    xticklabels(labels);
else
    figureHandle = figure(1);
    figureHandle.Renderer='Painters';
    clf;
    colormap(parula);
    numPlots = length(plotIndicies);
    clear h;
    for i = plotIndicies
        h(i) = subplot(numPlots+1,1,i); 
        hold on;
        plot(originalData(:,i));
        plot(reconstruction(:,i));
        h(i).XTick = [];
        h(i).YTick = [];
        title(labels{i});
    end
    linkaxes(h, 'x');
    h = subplot(numPlots+1,1,numPlots+1);
    bar(correlations)
    xticklabels(labels);
end



