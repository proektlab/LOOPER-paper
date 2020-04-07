
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

if length(hitAccuracy) > 15
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
else
    yticklabels({   '10Hz -> 5Hz', ...
                    '40Hz -> 15Hz', ...
                    '20Hz -> 30Hz', ...
                    '10Hz -> 15Hz', ...
                    '40Hz -> 25Hz', ...
                    '20Hz -> 50Hz'});
end
xticks(1:3);
xticklabels({'No output', 'Less than', 'Greater than'});
scatter(targetOutputs + 1, 1:length(targetOutputs), 'rx');
ylim([0.5 length(targetOutputs) + 0.5]);
