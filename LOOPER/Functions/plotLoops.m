%% Display results

PLOT_DELAY_EMBEDDED = 1;

numTrial = max(saveData.TrialData);
trialLength = size(saveData.FinalStream,1) / numTrial;
times = 1:trialLength;

if max(times) > trialLength
    times = times(1):trialLength;
end

trialIndicies = repmat(times, [1, numTrial]);
trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(times)));

finalStream = saveData.FinalStream;

[trialPCABasis, ~] = pca(saveData.ClusterMeans, 'NumComponents',3);
% trialPCABasis = eye(size(finalStream,2), 3);

loopStarts = (0:numTrial-1)*trialLength+1;

finalStream(loopStarts,:) = nan;

figureHandle = figure(10000);
figureHandle.Renderer='Painters';
clf;
hold on;

colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;
% colors(4,:) = colors(3,:);

lineColors = 1:saveData.BestLoopCount;
lineColors = colors(lineColors,:);

allClusterMeans = [];
for i = 1:saveData.BestLoopCount
    thisLoopIDs = find(saveData.BestLoopAssignments(:,1) == i);
    thisLoopClusters = saveData.BestLoopAssignments(thisLoopIDs,2);
    
    allIDs = find(saveData.BestStateMap(:,1) == i);
    trialIDs = floor((allIDs-1)/trialLength) + 1;
    conditionID = i;
    
    badIDs = find(saveData.BestStateMap(:,1) == i);
    badIDs = setdiff(1:size(finalStream,1), badIDs);
    
    plotStream = finalStream;
    plotStream(badIDs,:) = nan;
    
    if PLOT_DELAY_EMBEDDED
        h = plot3(plotStream(trialIndicies,:)*trialPCABasis(:,1), plotStream(trialIndicies,:)*trialPCABasis(:,2), plotStream(trialIndicies,:)*trialPCABasis(:,3), 'LineWidth', 1, 'Color', lineColors(conditionID,:));
    else
        h = plot3(plotStream(trialIndicies,1), plotStream(trialIndicies,2), plotStream(trialIndicies,3), 'LineWidth', 1, 'Color', lineColors(conditionID,:));
    end
    h.Color(4) = 0.9;
    
    meanTimes = [];
    clusterLengths = [];
    clusterMeans = [];
    clusterSTDs = [];
    for j = 1:length(thisLoopClusters)
        thisIndices = find(saveData.BestStateMap(:,1) == i & saveData.BestStateMap(:,2) == thisLoopClusters(j));

        meanTimes(j) = mode(mod(thisIndices, trialLength)+1);
        clusterLengths(j) = length(thisIndices);
        
        clusterMeans(j,:) = nanmean(plotStream(thisIndices, :), 1);
        clusterSTDs(j,:) = nanstd(plotStream(thisIndices, :), [], 1);
        
        allClusterMeans(i, thisLoopClusters(j),:) = clusterMeans(j,:);
    end
    
    thisLoopIDs = 1:length(thisLoopClusters);
    clusterOrder = 1:length(thisLoopClusters);
    
    startCluster = 1;
    
    sortedOrder = [startCluster:length(clusterOrder) 1:startCluster-1];
    thisLoopIDs = thisLoopIDs(sortedOrder);
    
    meanTimes = meanTimes(sortedOrder);
    
    badTimes = find(meanTimes <= min(times) | meanTimes >= max(times));
    goodTimes = setdiff(sortedOrder, badTimes);
    
    clusterSTDs = clusterSTDs(thisLoopIDs(goodTimes),:);
    thisTrace = clusterMeans(thisLoopIDs(goodTimes),:);
    
    if length(goodTimes) <= 1
        continue;
    end
    
    thisTrace = filterData(thisTrace', 0.5, [], 1, 0)';
    
    if PLOT_DELAY_EMBEDDED
        tubeMean = thisTrace*trialPCABasis(:,1:3);
        projectedSTDs = clusterSTDs*trialPCABasis(:,1:3);
    else
        tubeMean = [thisTrace(:,1:3) zeros(size(thisTrace,1), 1)];
        projectedSTDs = [clusterSTDs(:,1:3) zeros(size(thisTrace,1), 1)];
    end
    
    projectedDerivative = diff(tubeMean);
    projectedDerivative = [projectedDerivative(1,:); projectedDerivative];
    
    orthogonalSTDs = zeros(size(projectedSTDs,1),1);
    for j = 1:size(projectedDerivative,1)
        orthogonalSTDs(j) = norm(projectedSTDs(j,:));
    end
    
    orthogonalSTDs = filterData(orthogonalSTDs, 0.5, [], 1, 2);
    
    [X,Y,Z,V] = tubeplot(tubeMean(:,1), tubeMean(:,2), tubeMean(:,3), orthogonalSTDs, ones(1,size(tubeMean,1)),10,[0 0 1]);

    t = surf(X,Y,Z,V, 'FaceColor', lineColors(conditionID,:));
    t.EdgeColor = lineColors(conditionID,:);
    t.EdgeAlpha = 0.5;
    t.FaceAlpha = 0.2;
end