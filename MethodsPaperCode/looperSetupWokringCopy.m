
SMOTHING_TIME = 0.1; %in s
FRAME_RATE = 20;

classes = [result{2:end,1}];
hits = [result{2:end,3}];
maxClass = max(classes);
allSpikes = {};
longestTrial = 0;
numTrials = max([result{2:end,2}]);
for i = 2:size(result,1)
    spikes = [result{i,6}];
    
    for j = 1:length(spikes)-2
        allSpikes{i,j} = spikes{j};
        
        longestTrial = max([longestTrial, spikes{j}]);
    end
end

TOTAL_TRIAL_TIME = ceil(FRAME_RATE * (longestTrial / 1000 + SMOTHING_TIME*3));

timeVector = (1:TOTAL_TRIAL_TIME) / FRAME_RATE;

smoothedAll = zeros(size(allSpikes, 1), TOTAL_TRIAL_TIME, size(allSpikes, 2));
for i = 1:size(allSpikes, 1)
    for neuronID = 1:size(allSpikes, 2)
        for spikeID = 1:length(allSpikes{i, neuronID})            
            gaussianVector = exp(-(timeVector - allSpikes{i, neuronID}(spikeID)/1000).^2/(2*SMOTHING_TIME^2));

            smoothedAll(i, :, neuronID) = smoothedAll(i, :, neuronID) + reshape(gaussianVector, [1, length(gaussianVector)]);
        end
    end
end

%%
hits = [result{2:end,3}];
conditions = [result{2:end,1}];
frequencies = [result{2:end,4}];
frequencies2 = [result{2:end,5}];
uniqueFrequencies = unique(frequencies);

plotNeuron = 7;
colors = parula(length(uniqueFrequencies));

figure(3);
clf;
hold on
colorCounter = 1;
labels = {};
for frequencyID = uniqueFrequencies
    thisIDs = find(frequencies == frequencyID);
    thisIDs(hits(thisIDs) == 0) = [];
    thisIDs(frequencies(thisIDs) > frequencies2(thisIDs)) = [];
    
    firstCondition = min(conditions(thisIDs));
    
    thisIDs(conditions(thisIDs) > firstCondition) = [];
    
    thisIDs = thisIDs + 1;
    thisTrace = mean(smoothedAll(thisIDs, :, plotNeuron), 1);
    
    labels{colorCounter} = num2str(frequencyID);
    
    plot(timeVector, thisTrace, 'Color', colors(colorCounter,:), 'LineWidth', 2);
    colorCounter = colorCounter + 1;
end

legend(labels);

%%
frequencies = [result{2:end,4}];
frequencies2 = [result{2:end,5}];

trials = [result{2:end,2}];
maxTrial = max(trials);
conditions = [result{2:end,1}];

% plotNeuron = 3;
condition = 3;


figure(2);
clf;
hold on;
colorCounter = 1;

conditionIDs = find(conditions == condition);

colors = parula(length(conditionIDs));

for trialID = 1:length(conditionIDs)
%     thisIDs = find(conditions == condition);
%     thisIDs(hits(thisIDs) == 0) = [];
%     thisIDs(frequencies(thisIDs) < frequencies2(thisIDs)) = [];
    
%     firstCondition = min(conditions(thisIDs));
%     
%     thisIDs(conditions(thisIDs) > firstCondition) = [];
%     
    thisID = conditionIDs(trialID) + 1;
    thisTrace = smoothedAll(thisID , :, plotNeuron);
    
%     labels{colorCounter} = num2str(frequencyID);
    
    plot(timeVector, thisTrace, 'Color', colors(colorCounter,:), 'LineWidth', 2);
    colorCounter = colorCounter + 1;
end
title(['Stim 1: ' num2str(frequencies(conditionIDs(1))) ' Stim 2: ' num2str(num2str(frequencies2(conditionIDs(1))))]);

%% From PCA

% X = zscore(firingRatesAverage(:,:)');
X = (firingRatesAverage(:,:)');
times = 1:900;
doZScore = 0;
usePCA = 1;
numDelays = 0;
delayTime = 0;
figurePlot = 1;

[pcaBasis, ~, ~, ~, explained] = pca(X, 'NumComponents', 20);
Z = X * pcaBasis;
dataDim = size(firingRatesAverage);

Zfull = reshape(Z', [size(pcaBasis,2) dataDim(2:end)]);

colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;

plotComponents = [1,2,3];
figure(1);
clf;
hold on;
for f=1:size(Zfull,2)
    for d = 1:size(Zfull,3)
        if d == 1
            lineType = '.';
        else
            lineType = '-';
        end
        
        thisTrace = squeeze(Zfull(plotComponents, f, d, :));

        plot3(thisTrace(1,:), thisTrace(2,:), thisTrace(3,:), lineType, 'color', colors(f,:), 'LineWidth', 2)
    end
end

%% From dPCA

%Run romoAnalysisPipeline and first bit of pfcRomo_dpcaPipeline for dPCA

X = firingRatesAverage(:,:)';
Xcen = bsxfun(@minus, X, mean(X));
dataDim = size(firingRatesAverage);
Z = Xcen * W;

componentsToPlot = 1:50;

Zfull = reshape(Z(:,componentsToPlot)', [length(componentsToPlot) dataDim(2:end)]);

colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;

stimMarg = find(strcmp(margNames, 'Stimulus'));
decisionMarg = find(strcmp(margNames, 'Decision'));
interactionMarg = find(strcmp(margNames, 'Interaction'));

topStim = find(whichMarg == stimMarg, 1);
topDecision = find(whichMarg == decisionMarg, 1);
topInteraction = find(whichMarg == interactionMarg, 1);

% plotComponents = [topStim, topDecision, topInteraction];
plotComponents = [topStim, topDecision];%, topInteraction];
% plotComponents = [topStim, 12, 18];
times = 1:900;
numDelays = 5;
delayTime = 5;
doZScore = 0;
usePCA = 0;
figurePlot = 2;

allTraces = [];
for f=1:size(Zfull,2)
    for d = 1:size(Zfull,3)
        thisTrace = squeeze(Zfull(plotComponents, f, d, times));        
        if doZScore
            thisTrace = zscore(thisTrace, [], 2);
        end
        thisTrace = delayEmbed(thisTrace, numDelays, delayTime);
        
        allTraces = [allTraces thisTrace];
    end
end

if size(allTraces,1) < 3
    pcaBasis = eye(size(allTraces,1),3);
else
    [pcaBasis, ~] = pca(allTraces', 'NumComponents', 3);
end

figure(2);
clf;
hold on;
for f=1:size(Zfull,2)
    for d = 1:size(Zfull,3)
        if d == 1
            lineType = '.';
        else
            lineType = '-';
        end
        
        thisTrace = squeeze(Zfull(plotComponents, f, d, times));
        if doZScore
            thisTrace = zscore(thisTrace, [], 2);
        end
        thisTrace = delayEmbed(thisTrace, numDelays, delayTime);
        
%         [pcaBasis, ~] = pca(thisTrace', 'NumComponents', 3);

        plot3(thisTrace'*pcaBasis(:,1), thisTrace'*pcaBasis(:,2), thisTrace'*pcaBasis(:,3), lineType, 'color', colors(f,:), 'LineWidth', 2)
%         scatter3(thisTrace'*pcaBasis(:,1), thisTrace'*pcaBasis(:,2), thisTrace'*pcaBasis(:,3), 32, 1:size(thisTrace,2))
    end
end


%% Bootstrap trials

% descision 1 -> f1 < f2
% descision 2 -> f1 > f2

BOOTSTRAP_AMOUNT = 1/2;

f1stimulationset = [10 14 18 24 30 34; 10 14 18 24 30 34];
f2stimulationset = [ 6  8 10 16 22 26; 18 22 26 32 38 44];

useF1s = [1,3,6];

% usePCA = 0;

if usePCA
    plotComponents = 1:size(pcaBasis,2);
end

testCounter = 1;

dPCACount = 30;
numBootStraps = 10;
firingRatesAverageSize = size(firingRatesAverage);
allBootstrappedFiringRates = zeros([length(plotComponents), length(useF1s), firingRatesAverageSize(3), length(times), numBootStraps]);
allInputs = zeros([1, length(useF1s), firingRatesAverageSize(3), length(times), numBootStraps]);
figure(figurePlot);
clf;
for bootStrapID = 1:numBootStraps
    bootStrappedFiringRates = zeros([firingRatesAverageSize(1), length(useF1s), firingRatesAverageSize(3), length(times)]);
    for i = 1:size(bootStrappedFiringRates,1)
        for j = 1:length(useF1s)%1:size(bootStrappedFiringRates,2)
            for k = 1:size(bootStrappedFiringRates,3)
                trialCount = trialNum(i,j,k);
                trialCount = ceil(trialCount*BOOTSTRAP_AMOUNT);
trialSizes(testCounter) = size(thisData,2);
testCounter = testCounter + 1;
                thisData = squeeze(firingRates(i,useF1s(j),k,times,1:trialCount));
%                 thisData = squeeze (firingRates(i,useF1s(j),k,times,1:2:trialNum(i,j,k)));
%                 thisData = squeeze (firingRates(i,useF1s(j),k,times,trialCount+1:trialNum(i,j,k)));

                bootNum = trialCount;
%                 bootNum = 1;
                bootstrappedIndices = randsample(size(thisData,2), bootNum, true);

                bootStrappedFiringRates(i,j,k,:) = mean(thisData(:,bootstrappedIndices),2);
                
                if i == 1
                    inputs = zeros(size(time));
                    inputs(time >= 0 & time <= 0.5) = f1stimulationset(k,useF1s(j));
                    inputs(time >= 3.5 & time <= 4) = f2stimulationset(k,useF1s(j));
                    
                    allInputs(1, j, k, :, bootStrapID) = inputs(times);
                end
            end
        end
    end

    if ~usePCA
        X = bootStrappedFiringRates(:,:)';
        Xcen = bsxfun(@minus, X, mean(X));
        dataDim = size(bootStrappedFiringRates);
        Z = Xcen * W;

        Zfull = reshape(Z(:,componentsToPlot)', [length(componentsToPlot) dataDim(2:end)]);

        allBootstrappedFiringRates(:,:,:,:,bootStrapID) = Zfull(plotComponents,:,:,:);
    else
        X = bootStrappedFiringRates(:,:)';

%         [pcaBasis, ~] = pca(X, 'NumComponents', 20);
        Z = X * pcaBasis;
        dataDim = size(bootStrappedFiringRates);

        Zfull = reshape(Z', [size(pcaBasis,2) dataDim(2:end)]);

        allBootstrappedFiringRates(:,:,:,:,bootStrapID) = Zfull(:,:,:,:);
    end
    
    hold on;
    for f=1:size(Zfull,2)
        for d = 1:size(Zfull,3)
            if d == 1
                lineType = '.';
            else
                lineType = '-';
            end

            thisTrace = squeeze(Zfull(plotComponents, f, d, times));
            if doZScore
                thisTrace = zscore(thisTrace, [], 2);
            end
            thisTrace = delayEmbed(thisTrace, numDelays, delayTime);

    %         [pcaBasis, ~] = pca(thisTrace', 'NumComponents', 3);
    
            if ~usePCA
                plotData = thisTrace'*pcaBasis;
            else
                plotData = thisTrace';
            end

            plot3(plotData(:,1), plotData(:,2), plotData(:,3), lineType, 'color', colors(useF1s(f),:), 'LineWidth', 2)
    %         scatter3(thisTrace'*pcaBasis(:,1), thisTrace'*pcaBasis(:,2), thisTrace'*pcaBasis(:,3), 32, 1:size(thisTrace,2))
        end
    end
end

%% Setup for LOOPER

decimateAmount = 5;

matSize = size(allBootstrappedFiringRates);
matData = permute(allBootstrappedFiringRates, [1, 4, 2, 3, 5]);
allTrials = reshape(matData, [matSize(1), matSize(4), matSize(2) * matSize(3) * matSize(5)]);

inputData = permute(allInputs, [1, 4, 2, 3, 5]);
allInputs = reshape(inputData, [1, matSize(4), matSize(2) * matSize(3) * matSize(5)]);


finalTrials = [];
finalInputs = [];
for i = 1:size(allTrials,1)
    for j = 1:size(allTrials,3)
        finalTrials(i,:,j) = decimate(squeeze(allTrials(i,:,j)),decimateAmount);
        
        if i == 1
        	finalInputs(1,:,j) = decimate(squeeze(allInputs(1,:,j)),decimateAmount);
        end
    end
end


figure(2);
clf;
hold on;
pcaData = reshape(finalTrials, size(finalTrials, 1),[]);
if size(pcaData,1) < 3
    trialPCABasis = eye(size(pcaData,1),3);
else
    [trialPCABasis, ~] = pca(pcaData', 'NumComponents',3);
end
for i = 1:size(finalTrials,3)
    thisTrace = finalTrials(:,:,i)'*trialPCABasis;
    
    plot3(thisTrace(:,1),thisTrace(:,2),thisTrace(:,3));
end

%% From RNN

NUM_TRIALS = 10;

useClasses = [1,2,3,6,7,8];

trialData = permute(dynamics, [3, 1, 2]);

allData = reshape(trialData, size(trialData,1),[]);

[pcaBasis,~] = pca(allData', 'NumComponents', 20);

% pcaBasis = eye(size(size(trialData,1)));

figure(1);
clf;
colors = lines(length(useClasses));
hold on;

finalTrials = [];
finalInputs = [];
trialCounter = 1;
for i = 1:length(useClasses)
    thisIndices = find(classes == useClasses(i));
    
    for j = 1:NUM_TRIALS    
        thisData = (squeeze(trialData(:,:,thisIndices(j)))' * pcaBasis)';
        for k = 1:size(thisData,1)
            finalTrials(k,:,trialCounter) = decimate(thisData(k,:),2);
        end
        finalInputs(:,:,trialCounter) = decimate(inputs(:,:,thisIndices(j)),2);
        
        trialCounter = trialCounter + 1;
    end
    
    plot(mean(finalInputs(1,:,end-9:end),3), 'Color', colors(i,:));
end



%% Display results

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
times = 1:400;

if max(times) > trialLength
    times = times(1):trialLength;
end

trialIndicies = repmat(times, [1, numTrial]);
trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(times)));

finalStream = app.SavedData.FinalStream;

[trialPCABasis, ~] = pca(finalStream(trialIndicies,:), 'NumComponents',3);
% trialPCABAsis = pcaBasis;

loopStarts = (0:numTrial-1)*trialLength+1;

finalStream(loopStarts,:) = nan;

figure(1);
clf;
hold on;
colors = jet(app.SavedData.BestLoopCount);
% h = plot3(finalStream(loopStarts+1,:)*trialPCABasis(:,1), finalStream(loopStarts+1,:)*trialPCABasis(:,2), finalStream(loopStarts+1,:)*trialPCABasis(:,3), 'o', 'LineWidth', 0.5);

% for i = [6, 7]%app.SavedData.BestLoopCount    
% for i = [4, 5]%app.SavedData.BestLoopCount    
% for i = [1,2,5,6]%app.SavedData.BestLoopCount    
% for i = [3,7,8,10]%app.SavedData.BestLoopCount    
for i = 1:app.SavedData.BestLoopCount    
    thisLoopIDs = find(app.SavedData.BestLoopAssignment(:,1) == i);
    thisLoopClusters = app.SavedData.BestLoopAssignment(thisLoopIDs,2);
    
    badIDs = find(app.SavedData.BestStateMap(:,1) == i);
    badIDs = setdiff(1:size(finalStream,1), badIDs);
    
    plotStream = finalStream;
    plotStream(badIDs,:) = nan;
    
    h = plot3(plotStream(trialIndicies,:)*trialPCABasis(:,1), plotStream(trialIndicies,:)*trialPCABasis(:,2), plotStream(trialIndicies,:)*trialPCABasis(:,3), 'LineWidth', 0.5, 'Color', colors(i,:));
    h.Color(4) = 0.5;
    
    meanTimes = [];
    clusterLengths = [];
    for j = 1:length(thisLoopClusters)
        thisIndices = find(app.SavedData.BestStateMap(:,1) == i & app.SavedData.BestStateMap(:,2) == thisLoopClusters(j));

        meanTimes(j) = mean(mod(thisIndices, trialLength)+1);
        clusterLengths(j) = length(thisIndices);
    end
    
%     startPoints = find(app.SavedData.BestStateMap(loopStarts,1) == i);
%     bestStartPoint = mode(app.SavedData.BestStateMap(loopStarts(startPoints),2));
%     
%     
    clusterOrder = app.SavedData.BestLoopAssignment(thisLoopIDs,2);
    
%     testTimes = meanTimes;
%     testTimes(clusterLengths < 5) = 1000;
    [~, startCluster] = min(meanTimes);
    
    sortedOrder = [startCluster:length(clusterOrder) 1:startCluster-1];
    thisLoopIDs = thisLoopIDs(sortedOrder);
    
    for j = 1:length(clusterOrder)
        thisIndices = find(app.SavedData.BestLoopAssignment(:,1) == i & app.SavedData.BestLoopAssignment(:,2) == clusterOrder(sortedOrder(j)));
    end
    
    meanTimes = meanTimes(sortedOrder);
    
    badTimes = find(meanTimes < min(times) | meanTimes > max(times));
    goodTimes = setdiff(sortedOrder, badTimes);
    
    thisTrace = app.SavedData.BestEmission(thisLoopIDs(goodTimes), :);
    
    plot3(thisTrace*trialPCABasis(:,1), thisTrace*trialPCABasis(:,2), thisTrace*trialPCABasis(:,3), 'LineWidth', 2, 'Color', colors(i,:))
end

%%

STATE_SMOOTH = 2;
FLUX_CUTOFF = 3;
MINIMUM_STATE_TIME = 0 *(2+1);
IS_RNN = 0;

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
times = 1:900;

if max(times) > trialLength
    times = times(1):trialLength;
end

trialIndicies = repmat(times, [1, numTrial]);
trialIndicies = trialIndicies + kron((0:numTrial-1)*trialLength, ones(size(times)));

finalStream = app.SavedData.FinalStream;

[trialPCABasis, ~] = pca(finalStream(trialIndicies,:), 'NumComponents',3);
% trialPCABAsis = pcaBasis;

loopStarts = (0:numTrial-1)*trialLength+1;

finalStream(loopStarts,:) = nan;

loopIDs = app.SavedData.BestStateMap(:,1);
if MINIMUM_STATE_TIME > 0
    loopIDs = colfilt(loopIDs, [MINIMUM_STATE_TIME 1], 'sliding', @mode);
end

if IS_RNN
    loopOrder = [1, 6, 4, 3, 2, 5, 7, nan];

    lineStyles = {'r', 'g', 'b', 'r.', 'g.', 'b.'};
    lineColors = {'r', 'g', 'b', 'r', 'g', 'b'};
else
    loopOrder = [3, 6, 4, 5, 1, 2, nan];

    lineStyles = {'r', 'g', 'b', 'r.', 'g.', 'b.'};
    lineColors = {'r', 'g', 'b', 'r', 'g', 'b'};
end

figureHandle = figure(2);
figureHandle.Renderer='Painters';
clf;
hold on;
for i = 1:numTrial
    trialIDs = (0:trialLength-1)+(i-1)*trialLength+1;
    
    if IS_RNN
        conditionID = floor((i-1)/10)+1;
    else
        conditionID = mod((i-1),6)+1;
    end
    
    IDs = loopIDs(trialIDs);
    IDs(IDs == 0) = length(loopOrder);
    
    h = plot(0:trialLength-1, loopOrder(IDs) + normrnd(0,0.03,size(trialIDs)), lineStyles{conditionID});
    h.Color(4) = 0.2;
    
    transitionCheck = colfilt(loopOrder(IDs), [1 2], 'sliding', @mean);
    transitionCheck(transitionCheck == loopOrder(IDs)) = 0;
    transitionPoints = find(transitionCheck);
    transitionPoints = [transitionPoints, transitionPoints+1];
    transitionPoints(transitionPoints > trialLength) = [];
    
    transitionData = nan(size(loopOrder(IDs)));
    transitionData(transitionPoints) = loopOrder(IDs(transitionPoints));
    
    h = plot(0:trialLength-1, transitionData + normrnd(0,0.03,size(trialIDs)), lineColors{conditionID});
    h.Color(4) = 0.2;
end

%%

% transitionTimes = {};
% for i = 1:app.SavedData.BestLoopCount
%     for j = 1:app.SavedData.BestLoopCount
%         transitionTimes{i,j} = 
%     end
% end


figure(2);
clf;
hold on;
colors = jet(app.SavedData.BestLoopCount);

clusterLengths = [];
for i = 1:size(app.SavedData.BestModel,1)
    loopID = app.SavedData.BestLoopAssignments(i,1);
    clusterID = app.SavedData.BestLoopAssignments(i,2);
    
    thisIndices = find(app.SavedData.BestStateMap(:,1) == loopID & app.SavedData.BestStateMap(:,2) == clusterID);

    clusterTimes = zeros(size(app.SavedData.BestStateMap(:,1)));
    clusterTimes(thisIndices) = 1;
    
    [~,trajectoryCounts] = findpeaks(clusterTimes);
    
    clusterLengths(i) = length(trajectoryCounts);
end

fluxMatrix = app.SavedData.BestModel .* clusterLengths';

loopTransitions = [];
loopTransitionStrength = [];
denosiedFlux = zeros(size(fluxMatrix));
for i = 1:app.SavedData.BestLoopCount    
    thisLoopIDs = find(app.SavedData.BestLoopAssignments(:,1) == i);
    
    for j = 1:app.SavedData.BestLoopCount
        if i ~= j        
            otherLoopIDs = find(app.SavedData.BestLoopAssignments(:,1) == j);

            transitionsSmoothed = imgaussfilt(fluxMatrix(thisLoopIDs, otherLoopIDs),STATE_SMOOTH); 
            transitionsSmoothed = transitionsSmoothed;
            transitionPeaks = imregionalmax(transitionsSmoothed);
    %         transitionRegions = transitionsSmoothed;
    %         transitionRegions = transitionRegions > 0.015; % just under smallest pixel for smoothing gaussian

            transitionWatershed = watershed(-transitionsSmoothed);

            transitionRegions = regionprops(transitionPeaks);

            finalTransitions = zeros(size(transitionsSmoothed));
            for k = 1:length(transitionRegions)
                centriod = round(transitionRegions(k).Centroid);

                finalTransitions(centriod(2),centriod(1)) = sum(sum(transitionsSmoothed(transitionWatershed == transitionWatershed(centriod(2), centriod(1)))));
            end
            
            finalTransitions(finalTransitions < FLUX_CUTOFF) = 0;

            loopTransitions{i,j} = [];
            loopTransitionStrength{i,j} = [];
            if ~isempty(finalTransitions)
                peakPositions = find(finalTransitions > 0);
                for k = 1:length(peakPositions)
                    loopTransitions{i,j}(k,:) = ind2sub(peakPositions(k), size(finalTransitions));
                    loopTransitionStrength{i,j}(k) = finalTransitions(peakPositions(k));
                end
            end
            denosiedFlux(thisLoopIDs, otherLoopIDs) = finalTransitions;
        end
    end
end

loopTimes = [];
clusterTimes = [];
for i = 1:app.SavedData.BestLoopCount
    thisLoopIDs = find(app.SavedData.BestLoopAssignments(:,1) == i);
    clusterID = app.SavedData.BestLoopAssignments(thisLoopIDs,2);
    
    meanTimes = [];
    for j = 1:length(clusterID)
        thisIndices = find(app.SavedData.BestStateMap(:,1) == i & app.SavedData.BestStateMap(:,2) == clusterID(j));

        meanTimes(j) = mode(mod(thisIndices, trialLength)+1);
    end
    
    loopTimes{i} = meanTimes;
    clusterTimes = [clusterTimes, meanTimes];
end


figure(3);
clf;
hold on;
colors = jet(app.SavedData.BestLoopCount);

for i = 5:7%app.SavedData.BestLoopCount    
    thisLoopIDs = find(app.SavedData.BestLoopAssignments(:,1) == i);
    thisLoopClusters = app.SavedData.BestLoopAssignments(thisLoopIDs,2);
    
    badIDs = find(app.SavedData.BestStateMap(:,1) == i);
    badIDs = setdiff(1:size(finalStream,1), badIDs);
    
    plotStream = finalStream;
    plotStream(badIDs,:) = nan;
    
%     h = plot3(plotStream(trialIndicies,:)*trialPCABasis(:,1), plotStream(trialIndicies,:)*trialPCABasis(:,2), plotStream(trialIndicies,:)*trialPCABasis(:,3), 'LineWidth', 0.5, 'Color', colors(i,:));
%     h.Color(4) = 0.5;
    
    clusterOrder = app.SavedData.BestLoopAssignments(thisLoopIDs,2);
    
    [~, startCluster] = min(loopTimes{i});
    
    sortedOrder = [startCluster:length(clusterOrder) 1:startCluster-1];
    thisLoopIDs = thisLoopIDs(sortedOrder);
    
    for j = 1:length(clusterOrder)
        thisIndices = find(app.SavedData.BestLoopAssignments(:,1) == i & app.SavedData.BestLoopAssignments(:,2) == clusterOrder(sortedOrder(j)));
    end
    
    sortedMeanTimes = loopTimes{i}(sortedOrder);
    
    badTimes = find(sortedMeanTimes < min(times) | sortedMeanTimes > max(times));
    goodTimes = setdiff(sortedOrder, badTimes);
    
%     thisTrace = app.SavedData.BestEmission(thisLoopIDs(goodTimes), :);
    
    plot(sortedMeanTimes, ones(size(sortedMeanTimes,2),1)*i, 'LineWidth', 2, 'Color', colors(i,:))

    for j = 1:size(loopTransitions,2)
        for k = 1:size(loopTransitions{i,j},1)
            clusterTransition = loopTransitions{i,j}(k,:);
            transitionTimes = [loopTimes{i}(clusterTransition(1)), loopTimes{j}(clusterTransition(2))];
            transitionLoops = [i, j];
            transitionTimes
            plot(transitionTimes, transitionLoops, 'LineWidth', 2, 'Color', colors(i,:));
        end
    end
end
