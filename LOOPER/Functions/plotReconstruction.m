%% Reconstruction

numTrial = max(saveData.TrialData);
trialLength = size(saveData.FinalStream,1) / numTrial;
plotTimes = 1:trialLength;

rawUnits = size(saveData.RawData,1);

clusterMeans = [];
clusterSTD = [];
clusterLengths = [];
meanTimes = [];
stateLookup = [];
for i = 1:size(saveData.BestLoopAssignments,1)
    loopID = saveData.BestLoopAssignments(i,1);
    phaseID = saveData.BestLoopAssignments(i,2);
    
    thisIndices = find(saveData.BestStateMap(:,1) == loopID & saveData.BestStateMap(:,2) == phaseID);
    
    trialTimes = mod(thisIndices, trialLength)+1;
    meanTimes(i,:) = mode(trialTimes);
    clusterLengths(i,:) = length(thisIndices);

    clusterMeans(i,:) = nanmean(saveData.FinalStream(thisIndices, :), 1);
    clusterSTD(i,:) = nanstd(saveData.FinalStream(thisIndices, :), [], 1);
    
    stateLookup(loopID, phaseID) = i;
end

thisLoopIDs = saveData.BestStateMap(:,1);
thisPhaseIDs = saveData.BestStateMap(:,2);

allIDs = [];
reconstructedStream = [];
for j = 1:length(thisLoopIDs)
    thisID = stateLookup(thisLoopIDs(j), thisPhaseIDs(j));
    
    if thisID < 1 || thisID > size(clusterMeans,1) || isnan(thisID)
        reconstructedStream(j,:) = nan(1, rawUnits);
        allIDs(j) = nan;
    else
        reconstructedStream(j,:) = clusterMeans(thisID,1:rawUnits)';
        allIDs(j) = thisID;
    end
end
reconstructedStream = reconstructedStream';

processedStream = saveData.FinalStream(:,1:rawUnits)';

figure(1001);
clf;
hold on;
if size(processedStream,1) <= 3
    for i = 1:size(processedStream,1)
        subplot(size(processedStream,1),1,i)
        
        plot(processedStream, 'r')
        plot(reconstructedStream, 'b')
    end
else
    subplot(2,1,1);
    imagesc(processedStream);
    subplot(2,1,2);
    imagesc(reconstructedStream);
end

validIDs = ~isnan(allIDs);

testStream = processedStream(:,validIDs);
testReconstruction = reconstructedStream(:,validIDs);

Rsquared = corr(testStream(:), testReconstruction(:));

suptitle(['Reconstruction R^2 = ' num2str(Rsquared)]);



