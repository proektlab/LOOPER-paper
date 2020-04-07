
load('F:\Dropbox\MethodsPaperCode\Data\monkeyDataForFigure3.mat');
originalData = saveData;

saveDir = 'F:\Dropbox\MethodsPaperCode\Data\parameterTests\';

chanceRate = 0.25;

%% "true" data

load('monkeyDataForFigure3')

SHOULD_VALIDATE = 0;

clear tempIDs;
tempIDs = 1:8;

Figure3PlotLoops

originalLength = size(saveData.RawData,2)/numTrial;
trialLength = size(saveData.FinalStream,1) / numTrial;
originalTime = [-2500 10000];
originalTime = originalTime(1)/1000 : 0.01 : originalTime(2)/1000;
truePlotTime = originalTime((decimateAmount:decimateAmount:trialLength*decimateAmount) + (originalLength - trialLength)*decimateAmount);

bestLoopCounts = [];
for i = 1:size(trialLoopIDs,1)    
    bestLoopCounts(i) = length(unique(trialLoopIDs(i,:)));
end

F1StartTime = 2.5;
F1StartTime = truePlotTime(find(bestLoopCounts == 3, 1));
F1EndTime = 3;
F2StartTime = truePlotTime(find(bestLoopCounts == 6, 1));
F2EndTime = truePlotTime(find(bestLoopCounts(1:end-1) == 6 & bestLoopCounts(2:end) ~= 6, 1));

F1Indices = find(truePlotTime >= F1StartTime & truePlotTime <= F1EndTime);
F2Indices = find(truePlotTime >= F2StartTime & truePlotTime <= F2EndTime);

F1Indices(F1Indices > trialLength) = [];
F2Indices(F2Indices > trialLength) = [];

targetF1Counts = mode(bestLoopCounts(F1Indices)); %3
targetF2Counts = mode(bestLoopCounts(F2Indices)); %6


%% Nearest neighbors
totalNeighborRatesF1 = [];
totalNeighborRatesF2 = [];
allNeighbors = 2:14;
counter = 1;
for nn = allNeighbors
    saveName = [saveDir 'testNN' num2str(nn) '.mat'];
    
    params = [];
%     params.PreprocessData.DelayCount = 7;
    params.NearestNeighbors = nn;

    if ~exist(saveName, 'file')
        LOOPER(originalData, true, [], [], [], params);
        save(saveName, 'saveData');
    else
        load(saveName);
    end    
    
    %% Do stuff.
    
    originalLength = size(saveData.RawData,2)/numTrial;
    trialLength = size(saveData.FinalStream,1) / numTrial;
    originalTime = [-2500 10000];
    originalTime = originalTime(1)/1000 : 0.01 : originalTime(2)/1000;
    plotTime = originalTime((decimateAmount:decimateAmount:trialLength*decimateAmount) + (originalLength - trialLength)*decimateAmount);
    
    loopIDs = saveData.BestStateMap(:,1);
    trialLength = length(loopIDs) / length(saveData.TrialSwitches);
    conditionIDs = [];
    for i = 1:6
        thisConditionIndices = repmat(1:trialLength, [10, 1]) + ((0:9)'*6*trialLength + (i-1)*trialLength);
        conditionIDs(:,:,i) = loopIDs(thisConditionIndices);
    end
    
    bestLoopCounts = [];
    for i = 1:size(trialLoops,1)    
        bestLoopCounts(i) = length(unique(trialLoops(i,:)));
    end
    
    F1Indices = find(plotTime >= F1StartTime & plotTime <= F1EndTime);
    F2Indices = find(plotTime >= F2StartTime & plotTime <= F2EndTime);
    
    F1Indices(F1Indices > trialLength) = [];
    F2Indices(F2Indices > trialLength) = [];
    
%     F1Rate = mean(bestLoopCounts(F1Indices) == targetF1Counts);
%     F2Rate = mean(bestLoopCounts(F2Indices) == targetF2Counts);

    F1Sets = [];
    for i = 1:3
        F1conditions = [i i+3];
        allF1IDs = conditionIDs(:,:,F1conditions);
        allF1IDs = permute(allF1IDs, [2,1,3]);
        
        F1Sets(:,:,i) = reshape(allF1IDs, [size(allF1IDs,1), size(allF1IDs,2)*size(allF1IDs,3)]);
    end
    
    F1overlaps = [];
    for i = 1:3
        thisSet = F1Sets(:,:,i);
        
        otherIDs = 1:3;
        otherIDs(i) = [];
        
        otherSets = F1Sets(:,:,otherIDs);
        otherSets = reshape(otherSets, [size(otherSets,1), size(otherSets,2)*size(otherSets,3)]);
        
        for t = 1:size(otherSets,1)
            F1overlaps(i,t) = sum(ismember(otherSets(t,:), thisSet(t,:))) / size(otherSets,2);
        end
    end
    
    F1overlaps = (1 - F1overlaps)/3;
    
    F2Sets = [];
    for i = 1:2
        F2conditions = [1:3] + (i-1)*3;
        allF1IDs = conditionIDs(:,:,F2conditions);
        allF1IDs = permute(allF1IDs, [2,1,3]);
        
        F2Sets(:,:,i) = reshape(allF1IDs, [size(allF1IDs,1), size(allF1IDs,2)*size(allF1IDs,3)]);
    end
    
    F2overlaps = [];
    for i = 1:2
        thisSet = F2Sets(:,:,i);
        
        otherIDs = 1:2;
        otherIDs(i) = [];
        
        otherSets = F2Sets(:,:,otherIDs);
        otherSets = reshape(otherSets, [size(otherSets,1), size(otherSets,2)*size(otherSets,3)]);
        
        for t = 1:size(otherSets,1)
            F2overlaps(i,t) = sum(ismember(otherSets(t,:), thisSet(t,:))) / size(otherSets,2);
        end
    end
    
    F2overlaps = (1 - F2overlaps)/2;
    
    F1Rate = mean(sum(F1overlaps(:,F1Indices)) - sum(F2overlaps(:,F1Indices)));
    F2Rate = mean(sum(F2overlaps(:,F2Indices)));
    
    if isnan(F1Rate)
        F1Rate = 0;
    end
    if isnan(F2Rate)
        F2Rate = 0;
    end
    
    totalNeighborRatesF1(counter) = F1Rate*F2Rate;
    totalNeighborRatesF2(counter) = F2Rate;
    counter = counter + 1;
end

figure(1);
clf;
hold on;
plot(allNeighbors,totalNeighborRatesF1, 'k');
% plot(allNeighbors,totalNeighborRatesF2, 'g');
plot(7,totalNeighborRates(find(allNeighbors==7)), 'rx');
xlabel('Neighbor counts');
ylabel('Bullshit made up decoding rate');

%% Delay embeds
totalDelayRatesF1 = [];
totalDelayRatesF2 = [];
allDelays = 0:30;
counter = 1;
for delays = allDelays
    saveName = [saveDir 'testDelays' num2str(delays) '.mat'];
    
    params = [];
    params.PreprocessData.DelayCount = delays;
%     params.NearestNeighbors = 4;

    if ~exist(saveName, 'file')
        LOOPER(originalData, true, [], [], [], params);
        save(saveName, 'saveData');
    else
        load(saveName);
    end
    
    %% Do stuff.
    
    originalLength = size(saveData.RawData,2)/numTrial;
    trialLength = size(saveData.FinalStream,1) / numTrial;
    originalTime = [-2500 10000];
    originalTime = originalTime(1)/1000 : 0.01 : originalTime(2)/1000;
    plotTime = originalTime((decimateAmount:decimateAmount:trialLength*decimateAmount) + (originalLength - trialLength)*decimateAmount);
    
    loopIDs = saveData.BestStateMap(:,1);
    trialLength = length(loopIDs) / length(saveData.TrialSwitches);
    conditionIDs = [];
    for i = 1:6
        thisConditionIndices = repmat(1:trialLength, [10, 1]) + ((0:9)'*6*trialLength + (i-1)*trialLength);
        conditionIDs(:,:,i) = loopIDs(thisConditionIndices);
    end
    
    bestLoopCounts = [];
    for i = 1:size(trialLoops,1)    
        bestLoopCounts(i) = length(unique(trialLoops(i,:)));
    end
    
    F1Indices = find(plotTime >= F1StartTime & plotTime <= F1EndTime);
    F2Indices = find(plotTime >= F2StartTime & plotTime <= F2EndTime);
    
    F1Indices(F1Indices > trialLength) = [];
    F2Indices(F2Indices > trialLength) = [];
    
%     F1Rate = mean(bestLoopCounts(F1Indices) == targetF1Counts);
%     F2Rate = mean(bestLoopCounts(F2Indices) == targetF2Counts);

    F1Sets = [];
    for i = 1:3
        F1conditions = [i i+3];
        allF1IDs = conditionIDs(:,:,F1conditions);
        allF1IDs = permute(allF1IDs, [2,1,3]);
        
        F1Sets(:,:,i) = reshape(allF1IDs, [size(allF1IDs,1), size(allF1IDs,2)*size(allF1IDs,3)]);
    end
    
    F1overlaps = [];
    for i = 1:3
        thisSet = F1Sets(:,:,i);
        
        otherIDs = 1:3;
        otherIDs(i) = [];
        
        otherSets = F1Sets(:,:,otherIDs);
        otherSets = reshape(otherSets, [size(otherSets,1), size(otherSets,2)*size(otherSets,3)]);
        
        for t = 1:size(otherSets,1)
            F1overlaps(i,t) = sum(ismember(otherSets(t,:), thisSet(t,:))) / size(otherSets,2);
        end
    end
    
    F1overlaps = (1 - F1overlaps)/3;
    
    F2Sets = [];
    for i = 1:2
        F2conditions = [1:3] + (i-1)*3;
        allF1IDs = conditionIDs(:,:,F2conditions);
        allF1IDs = permute(allF1IDs, [2,1,3]);
        
        F2Sets(:,:,i) = reshape(allF1IDs, [size(allF1IDs,1), size(allF1IDs,2)*size(allF1IDs,3)]);
    end
    
    F2overlaps = [];
    for i = 1:2
        thisSet = F2Sets(:,:,i);
        
        otherIDs = 1:2;
        otherIDs(i) = [];
        
        otherSets = F2Sets(:,:,otherIDs);
        otherSets = reshape(otherSets, [size(otherSets,1), size(otherSets,2)*size(otherSets,3)]);
        
        for t = 1:size(otherSets,1)
            F2overlaps(i,t) = sum(ismember(otherSets(t,:), thisSet(t,:))) / size(otherSets,2);
        end
    end
    
    F2overlaps = (1 - F2overlaps)/2;
    
    F1Rate = mean(sum(F1overlaps(:,F1Indices)) - sum(F2overlaps(:,F1Indices)));
    F2Rate = mean(sum(F2overlaps(:,F2Indices)));
    
    if isnan(F1Rate)
        F1Rate = 0;
    end
    if isnan(F2Rate)
        F2Rate = 0;
    end
    
    totalDelayRatesF1(counter) = F1Rate*F2Rate;
    totalDelayRatesF2(counter) = F2Rate;
    counter = counter + 1;
end

figure(2);
clf;
hold on;
plot(allDelays,totalDelayRatesF1, 'k');
% plot(allDelays,totalDelayRatesF2, 'g');
plot(11,totalDelayRates(find(allDelays==11)), 'rx');
xlabel('Delay counts');
ylabel('Bullshit made up decoding rate');
