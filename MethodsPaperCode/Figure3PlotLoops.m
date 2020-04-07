%% Run once with SHOULD_VALDIATE = 0 and once with SHOULD_VALIDATE = 1

STATE_SMOOTH = 2;
FLUX_CUTOFF = 3;
MINIMUM_STATE_TIME = 2 *(2+1);
IS_RNN = 0;
decimateTime = 10;
% SHOULD_VALIDATE = 0;

app.SavedData = saveData;



numTrial = max(app.SavedData.TrialData);
originalLength = size(app.SavedData.RawData,2)/numTrial;
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


if SHOULD_VALIDATE
    matSize = size(allBootstrappedFiringRates);
    matData = permute(allBootstrappedFiringRates, [1, 4, 2, 3, 5]);
    allTrials = reshape(matData, [matSize(1), matSize(4), matSize(2) * matSize(3) * matSize(5)]);

    inputData = permute(allInputs, [1, 4, 2, 3, 5]);
    allInputs = reshape(inputData, [1, matSize(4), matSize(2) * matSize(3) * matSize(5)]);


    finalTrials = [];
    finalInputs = [];
    for i = 1:size(allTrials,1)
        for j = 1:size(allTrials,3)
            finalTrials(i,:,j) = decimate(squeeze(allTrials(i,:,j)),decimateTime);

            if i == 1
                finalInputs(1,:,j) = decimate(squeeze(allInputs(1,:,j)),decimateTime);
            end
        end
    end
    
    rawData = convertToCell(finalTrials);

    lastEnd = 0;
    rawTrialData = [];
    for i = 1:length(rawData)
        rawTrialData(lastEnd + (1:size(rawData{i}, 2))) = i;

        lastEnd = lastEnd + size(rawData{i}, 2);
    end

    rawData = mergeData(rawData);

    [tempData, trialData, procesedTrialSwitches] = preprocessData(rawData, [], 0, [], 0, rawTrialData, 0, true, app.SavedData.PreprocessData.Smoothing, app.SavedData.PreprocessData.ZScore, app.SavedData.PreprocessData.DelayTime, app.SavedData.PreprocessData.DelayCount, app.SavedData.DataMean, app.SavedData.DataSTD);
    badIndicies = [procesedTrialSwitches procesedTrialSwitches+1];
    badIndicies(badIndicies > size(tempData,2)) = [];
    tempData(:,badIndicies) = [];
    
%     tempData = app.SavedData.FinalStream';

%     figure(4);
%     clf;
%     plot(tempData(1,:));
%     hold on;
%     plot(app.SavedData.FinalStream(:,1)');

    
    clusterDistances = [];
    finalStream = app.SavedData.FinalStream;
    for i = 1:size(app.SavedData.BestEmission,1)
        thisLoopPosition = app.SavedData.BestLoopAssignments(i,:);
        thisIndices = find(ismember(app.SavedData.BestStateMap, thisLoopPosition, 'rows'));
        
        stds = std(finalStream(thisIndices,:), [], 1);
        
        thisEmission = mean(finalStream(thisIndices,:), 1);
        
        clusterStream = (thisEmission - tempData') ./ stds;           
        
        clusterDistances(i,:) = sum(clusterStream.^2,2);% ./ sqrt(length(thisIndices));
    end
    
%     bestClusters = [];
%     clusterDistances = pdist2(app.SavedData.BestEmission, tempData');
    [~, bestClusters] = min(clusterDistances, [], 1);

    loopIDs = app.SavedData.BestLoopAssignments(bestClusters,1);
    
%     MINIMUM_STATE_TIME = 5;
    
    figureID = 3;
else
    figureID = 2;
    
    loopIDs = app.SavedData.BestStateMap(:,1);
end
if MINIMUM_STATE_TIME > 0
    loopIDs = colfilt(loopIDs, [MINIMUM_STATE_TIME 1], 'sliding', @mode);
end

if exist('tempIDs')
    loopOrder = tempIDs;

    lineStyles = {'r', 'g', 'b', 'r.', 'g.', 'b.'};
    lineColors = {'r', 'g', 'b', 'r', 'g', 'b'};
elseif IS_RNN
    loopOrder = [1, 6, 4, 3, 2, 5, 7, nan];

    lineStyles = {'r', 'g', 'b', 'r.', 'g.', 'b.'};
    lineColors = {'r', 'g', 'b', 'r', 'g', 'b'};
else
    loopOrder = [6, 5, 1, 2, 3, 4, 7, nan];
%     loopOrder = [4, 6, 5, 3, 1, 2, nan];

    lineStyles = {'r', 'g', 'b', 'r.', 'g.', 'b.'};
    lineColors = {'r', 'g', 'b', 'r', 'g', 'b'};
    
    decimateAmount = decimateTime;
    originalTime = [-2500 10000];
    originalTime = originalTime(1)/1000 : 0.01 : originalTime(2)/1000;
end

plotTime = originalTime((decimateAmount:decimateAmount:trialLength*decimateAmount) + (originalLength - trialLength)*decimateAmount);
% plotTime = 1:trialLength;


trialLoopIDs = [];
trialConditionIDs = [];

figureHandle = figure(figureID);
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
    
    trialLoopIDs(:,i) = loopOrder(IDs);
    trialConditionIDs(i) = conditionID;
    
    h = plot(plotTime, loopOrder(IDs) + normrnd(0,0.03,size(trialIDs)), lineStyles{conditionID});
    h.Color(4) = 0.2;
    
    transitionCheck = colfilt(loopOrder(IDs), [1 2], 'sliding', @mean);
    transitionCheck(transitionCheck == loopOrder(IDs)) = 0;
    transitionPoints = find(transitionCheck);
    transitionPoints = [transitionPoints, transitionPoints+1];
    transitionPoints(transitionPoints > trialLength) = [];
    
    transitionData = nan(size(loopOrder(IDs)));
    transitionData(transitionPoints) = loopOrder(IDs(transitionPoints));
    
    h = plot(plotTime, transitionData + normrnd(0,0.03,size(trialIDs)), lineColors{conditionID});
    h.Color(4) = 0.2;
end

yLimit = ylim;

F1Time = 0;
F2Time = 3.5;

plot([F1Time F1Time], yLimit);
plot([F2Time F2Time], yLimit);


