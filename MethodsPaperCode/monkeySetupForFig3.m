forValidation = 1;

%% From PCA

% X = zscore(firingRatesAverage(:,:)');
X = (firingRatesAverage(:,:)');
times = 1:900;
doZScore = 0;
usePCA = 1;
numDelays = 0;
delayTime = 0;
figurePlot = 1;

[pcaBasis, ~, ~, ~, explained] = pca(X, 'NumComponents', 10);
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



%% Bootstrap trials

% descision 1 -> f1 < f2
% descision 2 -> f1 > f2

BOOTSTRAP_AMOUNT = 1/2;


f1stimulationset = [10 24 34; 10 24 34];
f2stimulationset = [ 6  16 26; 18 32 44];

useF1s = [1,2,3];

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
allRawFiringRates = zeros([size(firingRatesAverage,1), length(useF1s), firingRatesAverageSize(3), length(times), numBootStraps]);
figure(figurePlot);
clf;
for bootStrapID = 1:numBootStraps
    bootStrappedFiringRates = zeros([firingRatesAverageSize(1), length(useF1s), firingRatesAverageSize(3), length(times)]);
    for i = 1:size(bootStrappedFiringRates,1)
        for j = 1:length(useF1s)%1:size(bootStrappedFiringRates,2)
            for k = 1:size(bootStrappedFiringRates,3)
                trialCount = trialNum(i,j,k);
                trialCount = ceil(trialCount*BOOTSTRAP_AMOUNT);
                
%                 thisData = squeeze(firingRates(i,useF1s(j),k,times,1:trialCount));
                if forValidation
                    thisData = squeeze (firingRates(i,useF1s(j),k,times,2:2:trialNum(i,j,k)));
                else
                    thisData = squeeze (firingRates(i,useF1s(j),k,times,1:2:trialNum(i,j,k)));
                end
% %                 thisData = squeeze (firingRates(i,useF1s(j),k,times,trialCount+1:trialNum(i,j,k)));

                trialSizes(testCounter) = size(thisData,2);
                testCounter = testCounter + 1;

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
    allRawFiringRates(:,:,:,:,bootStrapID) = bootStrappedFiringRates;
    
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

decimateAmount = 10;

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

