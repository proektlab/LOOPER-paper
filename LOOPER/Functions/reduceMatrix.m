
 
if ~exist('clusterCounts', 'var') || isempty(clusterCounts)
    clusterCounts = [100 80 60 50 40 30 20 17 15 12 10 8];
end

if ~exist('distanceType', 'var') || isempty(distanceType)
    distanceType = 'correlation';
end

if ~exist('maxCheckTime', 'var') || isempty(maxCheckTime)
    maxCheckTime = 10;
end

hadError = 0;

try
    pdist2(1,1,distanceType);
catch
    msgbox('Distance measure must be one of: euclidean squaredeuclidean seuclidean mahalanobis cityblock minkowski chebychev cosine correlation hamming jaccard or spearman. See pdist2 documentation for details.');
    
    hadError = 1;
    
    return;
end

%%

% [eigenModes, eigenvalues] = eigs(asymmetricProbabilities, 20);
% eigenvalues = (diag(eigenvalues));
% eigenModes = eigenModes;
% 
% ANGLE_EPLISON = 0.01;
% ANGLE_EIGENMODE = find(imag(eigenvalues) > 0, 1);
% ANGLE_SIGMA = 2*pi/100;
% 
% TAIL_LENGTH = 1;
% 
% allPhases = angle(eigenModes(:,ANGLE_EIGENMODE))';
% 
% if isempty(allPhases)
%     allPhases = zeros(size(eigenModes(:,1)))';
% end

%%

showIndices = 1:size(finalDynamicsStream,1);

% figure(2);
% clf;
% colors = lines(8);
% scatter3(finalDynamicsStream(showIndices,:)*pcaBasis(:,1),finalDynamicsStream(showIndices,:)*pcaBasis(:,2),finalDynamicsStream(showIndices,:)*pcaBasis(:,3), 32, allPhases)
% xlabel('DPC1');
% ylabel('DPC2');
% zlabel('DPC3');
% colormap(hsv);

% %% Get detailed balance decomp
% [steadyState, ~] = eigs(asymmetricProbabilities', 1);
% steadyState = steadyState ./ nansum(steadyState);
% 
% if var(steadyState) == 0
%     steadyState = sum(asymmetricProbabilities^100,1) / size(asymmetricProbabilities,1);
% end
% 
% fluxMatrix = diag(steadyState) * asymmetricProbabilities;
% symmetricMatrix = (fluxMatrix+fluxMatrix.')/2;
% antisymmetricMatrix = (fluxMatrix-fluxMatrix.')/2;
% symmetricMatrix = symmetricMatrix ./ steadyState;
% antisymmetricMatrix = antisymmetricMatrix ./ steadyState;
% 
% %% Rescale diffusion map
% %         diffusedProbabilities = max(0, expm(symmetricMatrix) - eye(size(symmetricMatrix)) + antisymmetricMatrix);
% diffusedProbabilities = max(0, symmetricMatrix^2 + antisymmetricMatrix);
% % diffusedProbabilities = symmetricMatrix + antisymmetricMatrix;
% diffusedProbabilities = diffusedProbabilities ./ sum(diffusedProbabilities,2);

diffusedProbabilities = asymmetricProbabilities;

%% Get time evolution of true matrix
% asymmetricProbabilities = full(denseProbabilityMatrix);
% finalDynamicsStream = allPoints(1:size(asymmetricProbabilities,1),:);

MIN_MODE = 0.1;
MAX_TIMES = 50;
PHASE_SIGMA = pi/10;
PHASE_STEPS = 100;
% longestTime = round(log(MIN_MODE) / log(abs(eigenvalues(ANGLE_EIGENMODE))));
% timeStep = ceil(longestTime / MAX_TIMES);
%
%
% checkTimes = 1:timeStep:longestTime;
checkTimes = 1:maxCheckTime;
% deltaT = checkTimes(2) - checkTimes(1);
stepMatrix = diffusedProbabilities;

stateDistributions = [];
currentMatrix = diffusedProbabilities^checkTimes(1);

waitHandle = parfor_progressbar(length(checkTimes), 'Calculating true distributions');
for i = 1:length(checkTimes)
    stateDistributions(:,:,i) = currentMatrix*eye(size(diffusedProbabilities,2));
    
    currentMatrix = stepMatrix*currentMatrix;
    
    waitHandle.iterate(1);
end
close(waitHandle);




%% Build similarity of state matrix
% distanceType = 'correlation';

%         diffusionDistance = symmetricMatrix^(1/sigmaScales(bestSigmaValue));

allTimeMatrix = (asymmetricProbabilities);% * diffusionDistance;%expm(asymmetricProbabilities);
%         allTimeMatrix = allTimeMatrix .* (1 - eye(size(allTimeMatrix)));
similarities1 = pdist((allTimeMatrix), distanceType);
% matrixExponential = 1 - expm(symmetricMatrix);
% matrixExponential = matrixExponential .* (1 - eye(size(matrixExponential)));
% similarities2 = squareform(matrixExponential);
% angleDifferences = angleDiff(allPhases, allPhases');

%         indexChunk{1} = 1:size(asymmetricProbabilities,1);
%         indexChunk{2} = size(asymmetricProbabilities,1) + firstChuckIndices;


stateDistances = zeros(1, size(diffusedProbabilities,1) * (size(diffusedProbabilities,1)-1) / 2);
statePairIndex = 1;
waitHandle = parfor_progressbar(size(diffusedProbabilities,1), 'Calculating distances');
for i = 1:size(diffusedProbabilities,1)
    %             testMatrix = asymmetricProbabilities;
    %             KL1 = -log(testMatrix ./ testMatrix(i,:)) * testMatrix(i,:)';
    
    for j = i+1:size(diffusedProbabilities,1)
        %                 stateDistances(statePairIndex) = (KL1(j)) * exp(-angleDifferences(i,j)^2/(2*PHASE_SIGMA^2));
        %                 stateDistances(statePairIndex) = similarities1(statePairIndex) * similarities2(statePairIndex) * exp(-angleDifferences(i,j)^2/(2*PHASE_SIGMA^2));
        stateDistances(statePairIndex) = similarities1(statePairIndex);% + 0.0 * (similarities2(statePairIndex));% * similarities2(statePairIndex) * exp(-angleDifferences(i,j)^2/(2*PHASE_SIGMA^2));
        
        statePairIndex = statePairIndex + 1;
    end
    waitHandle.iterate(1);
end
close(waitHandle);

clustering = linkage(stateDistances, 'average');


%%
%     clusterCounts = [5000, 4000, 3000, 2000, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100, 90, 80, 70, 60, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 40, 35, 30, 25, 23, 21, 19, 17, 15, 13, 12, 11, 10, 9 ,8, 7, 6, 5, 4, 3, 2, 1];
%         clusterCounts = [200, 150, 120, 110, 105, 100, 95, 90, 80, 70, 60, 50, 40, 30, 20];

KLDivergences = [];
reducedMatrices = {};
cosineDistances = [];
waitHandle = parfor_progressbar(length(clusterCounts), 'Calculating reduced distributions');
for clusterIndex = 1:length(clusterCounts)
    
    maxClusters = clusterCounts(clusterIndex);
    clusterIDs = cluster(clustering,'maxclust',maxClusters);
    %     cutoff = median([clustering(end-maxClusters+1,3) clustering(end-maxClusters+2,3)]);
    %     dendrogram(clustering,'ColorThreshold',cutoff);
    
    disp(['Calculating response of ' num2str(maxClusters) ' clusters']);
    
    %% Build reduced matrix
    %     thisIndices = find(clusterIDs == i);
    %
    %     clusterSimilarities = similarities(thisIndices,thisIndices);
    %     clusterSimilarities = sum(clusterSimilarities, 2);
    %
    %     weights = [];
    %     for j = 1:length(clusterSimilarities)
    %         weights(j,1) = sum(clusterSimilarities > clusterSimilarities(j));
    %     end
    %     weights = weights ./ sum(weights);
    
    reducedMatrixTemp = zeros(max(clusterIDs), size(asymmetricProbabilities,2));
    for i = 1:size(reducedMatrixTemp,1)
        reducedMatrixTemp(i,:) = sum(asymmetricProbabilities(clusterIDs == i, :), 1) ./ sum(clusterIDs == i);
    end
    
    reducedMatrix = zeros(max(clusterIDs));
    for i = 1:size(reducedMatrix,1)
        reducedMatrix(:,i) = sum(reducedMatrixTemp(:, clusterIDs == i), 2);
    end
    
    reducedStateCounts = [];
    for j = 1:size(reducedMatrix,2)
        reducedStateCounts(j) = sum(clusterIDs == j);
    end
    
    testStateDistributions = [];
    currentMatrix = reducedMatrix^checkTimes(1);
    stepMatrix = reducedMatrix;
    
    for i = 1:length(checkTimes)
        reducedDistribution = currentMatrix*eye(size(reducedMatrix,2));
        %         reducedDistribution = currentMatrix*diag(size(asymmetricProbabilities,1) .* reducedStateCounts);
        
        expandedMatrixTemp = zeros(length(clusterIDs), size(reducedDistribution,2));
        for j = 1:size(expandedMatrixTemp,2)
            numClusters = sum(clusterIDs == j);
            expandedMatrixTemp(clusterIDs == j,:) = repmat(reducedDistribution(j, :), [numClusters,1]);
        end
        
        expandedMatrix = zeros(length(clusterIDs), length(clusterIDs));
        for j = 1:size(expandedMatrixTemp,2)
            numClusters = sum(clusterIDs == j);
            expandedMatrix(:,clusterIDs == j) = repmat(expandedMatrixTemp(:, j), [1,numClusters]) / numClusters;
        end
        
        testStateDistributions(:,:,i) = expandedMatrix;
        
        %         testStateDistributions(:,:,i) = zeros(size(asymmetricProbabilities));
        %         for j = 1:size(reducedDistribution,1)
        %             for k = 1:size(reducedDistribution,1)
        %                 indicesX = find(clusterIDs == j);
        %                 indicesY = find(clusterIDs == k);
        %
        %                 testStateDistributions(indicesX,indicesY,i) = reducedDistribution(j,k) ./ sum(clusterIDs == k);
        %             end
        %         end
        %         testStateDistributions(:,:,i) = testStateDistributions(:,:,i) ./ sum(testStateDistributions(:,:,i),2);
        
        
        currentMatrix = stepMatrix*currentMatrix;
    end
    
%     timeIndex = 3;
%     figure(1);
%     clf;
%     h(1) = subplot(1,2,1);
%     imagesc(sqrt(stateDistributions(:,:,timeIndex)));
%     colorAxis = caxis;
%     h(2) = subplot(1,2,2);
%     imagesc(sqrt(testStateDistributions(:,:,timeIndex)));
%     caxis(colorAxis);
%     linkaxes(h);
    
    %% Calculate KLs
    
    
    for i = 1:length(checkTimes)
        realDistribution = stateDistributions(:,:,i);
        reducedDistribution = testStateDistributions(:,:,i);
        
        elementwiseKL = reducedDistribution.*log(reducedDistribution ./ realDistribution);
        elementwiseKL(isnan(elementwiseKL)) = 0;
        elementwiseKL(isinf(elementwiseKL)) = 0;
        
        
        KLDivergences(clusterIndex,i) = quantile(sum(elementwiseKL,2), 0.95);
%         KLDivergences(clusterIndex,i) = max(sum(elementwiseKL,2));
        %         KLDivergences(clusterIndex,i) = mean(sum(elementwiseKL,2));
        
        %         similarities = pdist2(realDistribution, reducedDistribution, distanceMeasure);
        %         cosineDistances(clusterIndex,i) = max(diag(similarities));
    end
    
    figure(1);
    clf;
    plot(KLDivergences(clusterIndex,:));
    
    
    
    reducedMatrices{clusterIndex} = reducedMatrix;
    %     figure(2);
    %     clf;
    %     plot(cosineDistances(clusterIndex,:));
    waitHandle.iterate(1);
end
close(waitHandle);

%%

worstDivervenges = [];
KLPeaks = [];
for i = 1:size(KLDivergences,1)
    worstDivervenges(i) = max(KLDivergences(i,:));
%     [~,peaks] = findpeaks(KLDivergences(i,:));
%     if ~isempty(peaks)
%         KLPeaks(i) = peaks(1);
%     else
%         KLPeaks(i) = 1;
%     end
end

BICs = worstDivervenges*size(diffusedProbabilities,1) + (clusterCounts).^2/2*log(size(diffusedProbabilities,1)/(2*pi));

figure(1);
clf;
plot(KLDivergences');
strings = {};
for i = 1:length(clusterCounts)
    strings{i} = [num2str(clusterCounts(i)) ' clusters'];
end
legend(strings);
title('Time series of model fits');
xlabel('Model time steps');
ylabel('95% qunatile of KL divergence');

figure(2);
clf
plot(clusterCounts, BICs);
title('Model fitting');
xlabel('Number of clusters');
ylabel('BIC');

%% Build best model
[~, minIndex] = min(BICs);
maxClusters = clusterCounts(minIndex);
clusterIndex = find(clusterCounts == maxClusters);
clusterIDs = cluster(clustering,'maxclust',maxClusters);



validClusters = [];
countClusters = [];
clusterMeans = [];
clusterMeansPCA = [];
clusterSTDs = [];
clusterSTDsPCA = [];
% clusterPhaseMeans = [];
% clusterPhaseSTDs = [];
for i = 1:max(clusterIDs)
    thisIndices = find(clusterIDs == i);
    
    clusterMeans(i,:) = mean(finalDynamicsStream(thisIndices,:),1);
    clusterSTDs(i,:) = std(finalDynamicsStream(thisIndices,:), [], 1) / 3;
    
    clusterMeansPCA(i,:) = mean(finalDynamicsStream(thisIndices,:)*pcaBasis(:,1:3),1);
    clusterSTDsPCA(i,:) = std(finalDynamicsStream(thisIndices,:)*pcaBasis(:,1:3),[],1) / 3;
    
    countClusters(i) = length(thisIndices);
    validClusters(i) = countClusters(i) > size(diffusedProbabilities,1) / size(reducedMatrices{clusterIndex},1) / 10;
    
    
    figure(i);
    clf;
    hold on;
    if exist('DEBUG', 'var') && DEBUG
        thisTrace = finalDynamicsStream * pcaBasis;
        plot3(thisTrace(:,1), thisTrace(:,2), thisTrace(:,3));
        scatter3(thisTrace(thisIndices,1), thisTrace(thisIndices,2), thisTrace(thisIndices,3));
        scatter3(clusterMeansPCA(i,1), clusterMeansPCA(i,2), clusterMeansPCA(i,3), 300, 'rx', 'LineWidth', 10);
    end
    
%     meanValues = exp(1i*allPhases(thisIndices));
%     clusterPhaseMeans(i) = angle(mean(meanValues));
%     clusterPhaseSTDs(i) = std(allPhases(thisIndices));
end

finalReducedMatrix = reducedMatrices{clusterIndex};

% if shouldUseTerminalState        
%     clusterMeans(end+1,:) = nan(size(clusterMeans(1,:)));
%     clusterSTDs(end+1,:) = nan(size(clusterSTDs(1,:)));
%     clusterMeansPCA(end+1,:) = nan(size(clusterMeansPCA(1,:)));
%     clusterSTDsPCA(end+1,:) = nan(size(clusterSTDsPCA(1,:)));
%     countClusters(end+1) = 0;
%     validClusters(end+1) = 0;
%     clusterPhaseMeans(end+1) = nan;
%     clusterPhaseSTDs(end+1) = nan;
% end

% asymmetricReducedMatrix = circshift(finalReducedMatrix, 1, 2);
% 
% [clusterEigenVectors, clusterEigenValues] = eigs(asymmetricReducedMatrix, min(30, size(reducedMatrices{clusterIndex},1)));
% clusterEigenValues = diag(clusterEigenValues);


%%

% reducedMatrixTemp = zeros(max(clusterIDs), size(symmetricDistances,2));
% for i = 1:size(reducedMatrixTemp,1)
%     reducedMatrixTemp(i,:) = sum(symmetricDistances(clusterIDs == i, :), 1) ./ sum(clusterIDs == i);
% end
% 
% reducedMatrix = zeros(max(clusterIDs));
% for i = 1:size(reducedMatrix,1)
%     reducedMatrix(:,i) = sum(reducedMatrixTemp(:, clusterIDs == i), 2);
% end
% 
% reducedSymmetricDistances = reducedMatrix;
% 
% figure(6);
% clf
% imagesc(reducedSymmetricDistances);
% title('Reduced distance matrix');
% ylabel('From state');
% xlabel('To state');
% colormap(parula(256));

%%

% figure(7);
% clf;
% hold on;
% for i = 1:max(clusterIDs)
%     thisIndices = find(clusterIDs == i);
%     
%     scatter(allPhases(thisIndices), ones(size(allPhases(thisIndices))) * clusterPhaseMeans(i));
% end
% ylabel('Cluster mean');
% xlabel('Cluster distribution');

