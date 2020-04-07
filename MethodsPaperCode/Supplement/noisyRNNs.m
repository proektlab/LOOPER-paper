%% Use this to preprocess RNN data for LOOPER
% load('F:\Dropbox\ConservationOfAgentDynamics\WorkingMemoryTask\Experiments\workingMemoryCustomBroken\lstm\1\networkTrace.mat');

try
    TAIL_LENGTH = 1;

    percentCorrect = [];
    for i = 1:size(outputs,1)
        thisTargets = squeeze(targets(1,:,i));
        thisOutputs = double(squeeze(outputs(i,:)));

        DEBUG = 0;
        if DEBUG
            figure(1);
            clf;
            hold on;
            plot(thisTargets)
            plot(thisOutputs)
        end

        outputDiff = (thisTargets - thisOutputs);

        targetIndices = find(thisTargets ~= 0);
        checkIndicies = min(targetIndices) - TAIL_LENGTH:max(targetIndices) + TAIL_LENGTH;

        percentCorrect(i) = sum(outputDiff(checkIndicies) == 0) / length(checkIndicies);
    end

    missTrials = percentCorrect < 0.3;
    sum(missTrials) / length(missTrials)

    % From RNN

    NUM_TRIALS = 10;

    times = 1:180;

    useClasses = [0,1,2,3,4,5];

    trialData = permute(dynamics, [3, 1, 2]);

    allData = reshape(trialData, size(trialData,1),[]);

    [pcaBasis,~,~,~,explained] = pca(allData', 'NumComponents', 20);
    percentExplained = cumsum(explained);
    percentExplained = percentExplained(20)

%     pcaBasis = eye(size(trialData,1));

%     figure(1);
%     clf;
%     colors = lines(length(useClasses));
%     hold on;

    trialCounts = [];
    finalTrials = [];
    trialID = 1;
    for i = 1:length(useClasses)
        thisIndices = find(classes == useClasses(i) & missTrials == 0);

%         F1ID = mod(useClasses(i), 3) + 1;
%         F2ID = floor(useClasses(i) / 3) + 1;


        for j = 1:NUM_TRIALS
            thisData = (squeeze(trialData(:,times,thisIndices(j)))' * pcaBasis)';
            for k = 1:size(thisData,1)
                finalTrials(k,:, trialID) = decimate(thisData(k,:),2);
            end
            
            trialID = trialID + 1;
        end

%         plot(mean(finalInputs(1,:,end-9:end),3), 'Color', colors(i,:));
    end
catch
end

%% Good RNN
% load('F:\Dropbox\ConservationOfAgentDynamics\WorkingMemoryTask\Experiments\workingMemoryCustomFrequenciesLong\lstm\1\networkTrace.mat');
% load('goodRNNResult');

figureNumber = 1;

PlotNoisyLoops

%% Find miss trials
PLOT_BROKEN = 1;

RNNAccuracies = {};

% for PLOT_BROKEN = 0:1

%     if PLOT_BROKEN == 0
%         % Good RNN
%         load('F:\Dropbox\ConservationOfAgentDynamics\WorkingMemoryTask\Experiments\workingMemoryCustomGood\lstm\1\networkTraceTests.mat')
%         figureNumber = 3;
%     elseif PLOT_BROKEN == 1
%         % Broken RNN
%         load('F:\Dropbox\ConservationOfAgentDynamics\WorkingMemoryTask\Experiments\workingMemoryCustomGood\lstm\1\networkTraceTests0.2.mat')
%         figureNumber = 4;
%     end
    figureNumber = 3;

    TAIL_LENGTH = 0;

    allOutputs = [];
    allTargets = [];
    percentCorrect = [];
    for i = 1:size(outputs,1)
        thisTargets = squeeze(targets(1,:,i));
        thisOutputs = double(squeeze(outputs(i,:)));

        DEBUG = 0;
        if DEBUG

            figure(1);
            clf;
            hold on;
            plot(thisTargets)
            plot(thisOutputs)
        end


        outputDiff = (thisTargets - thisOutputs);

        targetIndices = find(thisTargets ~= 0);
        checkIndicies = min(targetIndices) - TAIL_LENGTH:max(targetIndices) + TAIL_LENGTH;


        allOutputs(i,:) = thisOutputs(checkIndicies);
        allTargets(i,:) = thisTargets(checkIndicies);

        percentCorrect(i) = sum(outputDiff(checkIndicies) == 0) / length(checkIndicies);
    end

    missTrials = percentCorrect < 0.3;
    hitTrials = percentCorrect > 0.5;

    sum(missTrials) / length(missTrials)

    maxClass = max(classes);

    accuracies = [];
    hitAccuracy = [];
    targetOutputs = [];
    for i = 1:maxClass+1
        thisIDs = find(classes == i-1);

    %     hitAccuracy(i) = mean(percentCorrect(thisIDs));
        thisOutputs = allOutputs(thisIDs,:);

        thisTargets = allTargets(thisIDs,:);
        targetOutputs(i) = mode(thisTargets(:));

        hitAccuracy(i,:) = histcounts(thisOutputs(:), [-0.5, 0.5, 1.5, 2.5], 'Normalization', 'pdf');

        for j = 1:2        
            if j == 1
                accuracies(i,j) = length(find(classes == i-1 & hitTrials));
            else
                accuracies(i,j) = length(find(classes == i-1 & missTrials));
            end
        end
    end

    missConditions = find(sum(accuracies >= 10,2) >= 2);

    figureHandle = figure(figureNumber);
    figureHandle.Renderer='Painters';
    clf
    hold on;
    imagesc(hitAccuracy);
    yticks(1:length(hitAccuracy));

    if length(hitAccuracy) == 6
         yticklabels({   '10Hz -> 5Hz', ...
                        '20Hz -> 15Hz', ...
                        '40Hz -> 30Hz', ...
                        '10Hz -> 150Hz', ...
                        '20Hz -> 30Hz', ...
                        '40Hz -> 50Hz'});
    elseif length(hitAccuracy) > 15
        yticklabels({   '20Hz -> 5Hz', ...
                        '20Hz -> 15Hz', ...
                        '20Hz -> 25Hz', ...
                        '20Hz -> 30Hz', ...
                        '20Hz -> 50Hz', ...
                        '40Hz -> 5Hz', ...
                        '40Hz -> 15Hz', ...
                        '40Hz -> 25Hz', ...
                        '40Hz -> 30Hz', ...
                        '40Hz -> 50Hz', ...
                        '10Hz -> 5Hz', ...
                        '10Hz -> 15Hz', ...
                        '10Hz -> 25Hz', ...
                        '10Hz -> 30Hz', ...
                        '10Hz -> 50Hz', ...
                        '0Hz -> 5Hz', ...
                        '0Hz -> 15Hz', ...
                        '0Hz -> 25Hz', ...
                        '0Hz -> 30Hz', ...
                        '0Hz -> 50Hz'});
    elseif length(hitAccuracy) > 15
        yticklabels({   '30Hz -> 20Hz', ...
                        '30Hz -> 30Hz', ...
                        '30Hz -> 40Hz', ...
                        '30Hz -> 50Hz', ...
                        '30Hz -> 60Hz', ...
                        '40Hz -> 20Hz', ...
                        '40Hz -> 30Hz', ...
                        '40Hz -> 40Hz', ...
                        '40Hz -> 50Hz', ...
                        '40Hz -> 60Hz', ...
                        '50Hz -> 20Hz', ...
                        '50Hz -> 30Hz', ...
                        '50Hz -> 40Hz', ...
                        '50Hz -> 50Hz', ...
                        '50Hz -> 60Hz'});
    end
    xticks(1:3);
    xticklabels({'No output', 'Less than', 'Greater than'});
    scatter(targetOutputs + 1, 1:length(targetOutputs), 'rx');
    ylim([0.5 length(targetOutputs) + 0.5]);
    
    RNNAccuracies{end+1} = hitAccuracy;
% end