
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

bestLoopIDs = [];
for i = 1:size(trialLoopIDs,1)
    for j = 1:6
        thisTrials = find(trialConditionIDs == j);
        
        bestLoopIDs(i,j) = mode(trialLoopIDs(i,thisTrials));
    end
end

%% Nearest neighbors
totalNeighborRates = [];
allNeighbors = 2:20;
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
    
    SHOULD_VALIDATE = 0;
    Figure3PlotLoops
    
    originalLength = size(saveData.RawData,2)/numTrial;
    trialLength = size(saveData.FinalStream,1) / numTrial;
    originalTime = [-2500 10000];
    originalTime = originalTime(1)/1000 : 0.01 : originalTime(2)/1000;
    plotTime = originalTime((decimateAmount:decimateAmount:trialLength*decimateAmount) + (originalLength - trialLength)*decimateAmount);
    
    F1Indices = find(plotTime >= 0.5 & plotTime <= 3);
    F2Indices = find(plotTime >= 3.5 & plotTime <= 6);
    
    F1Indices(F1Indices > trialLength) = [];
    F2Indices(F2Indices > trialLength) = [];
    
    validationPercentCorrect = [];
    for i = 1:size(trialLoopIDs,1)
        validationPercentCorrect(i) = sum(trialLoopIDs(i,:) == bestLoopIDs(i, trialConditionIDs)) / length(trialConditionIDs);
    end

    overall = mean(validationPercentCorrect);
    F1 = mean(validationPercentCorrect(F1Indices));
    F1remember = min(validationPercentCorrect(F1Indices));
    F2 = mean(validationPercentCorrect(F2Indices));
    
    totalNeighborRates(counter) = nanmean([F1 F2]);
    counter = counter + 1;
end



%% Delay embeds
totalDelayRates = [];
allDelays = 0:40;
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
    
    SHOULD_VALIDATE = 0;
    Figure3PlotLoops
    
    originalLength = size(saveData.RawData,2)/numTrial;
    trialLength = size(saveData.FinalStream,1) / numTrial;
    originalTime = [-2500 10000];
    originalTime = originalTime(1)/1000 : 0.01 : originalTime(2)/1000;
    plotTime = originalTime((decimateAmount:decimateAmount:trialLength*decimateAmount) + (originalLength - trialLength)*decimateAmount);
    
    F1Indices = find(plotTime >= 0.5 & plotTime <= 3);
    F2Indices = find(plotTime >= 3.5 & plotTime <= 6);
    
    F1Indices(F1Indices > trialLength) = [];
    F2Indices(F2Indices > trialLength) = [];
    
    validationPercentCorrect = [];
    for i = 1:size(trialLoopIDs,1)
        validationPercentCorrect(i) = sum(trialLoopIDs(i,:) == bestLoopIDs(i, trialConditionIDs)) / length(trialConditionIDs);
    end

    overall = mean(validationPercentCorrect);
    F1 = mean(validationPercentCorrect(F1Indices));
    F1remember = min(validationPercentCorrect(F1Indices));
    F2 = mean(validationPercentCorrect(F2Indices));
    
    totalDelayRates(counter) = nanmean([F1 F2]);
    counter = counter + 1;
end


%% Display
figure(1);
clf;
hold on;
plot(allNeighbors,totalNeighborRates, 'k');
plot(7,totalNeighborRates(find(allNeighbors==7)), 'rx');
plot(allNeighbors, ones(size(allNeighbors))*chanceRate, 'r');
xlabel('Neighbor counts');
ylabel('Bullshit made up decoding rate');

figure(2);
clf;
hold on;
plot(allDelays,totalDelayRates, 'k');
plot(11,totalDelayRates(find(allDelays==11)), 'rx');
plot(allDelays, ones(size(allDelays))*chanceRate, 'r');
xlabel('Delay counts');
ylabel('Bullshit made up decoding rate');
