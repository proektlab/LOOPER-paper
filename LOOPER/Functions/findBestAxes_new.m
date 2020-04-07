function [bestLocalProjections, bestSigmaValues, bestDiffSigmaValues, localDistances, bestNeighborCounts, distanceRatios] = findBestAxes(dynamicsStream, calcIndex, meanDistanceMatrix, verbose)
    %meanDistanceMatrix = pdist2(dynamicsStream, dynamicsStream);

    bestSigmaValues = [];
    bestNeighborCounts = [];
    bestLocalProjections = [];
    distanceRatios = [];
    bestDiffSigmaValues = [];
    localDistances = [];
    
%     dynamicsStream = [dynamicsStream(1:end-4,:), dynamicsStream(2:end-3,:), dynamicsStream(3:end-2,:), dynamicsStream(4:end-1,:), dynamicsStream(5:end-0,:) ];

    %% Near trajectories

    MIN_RETURN_TIME = 10;
    
    if ~exist('meanDistanceMatrix', 'var') || isempty(meanDistanceMatrix)
        meanDistanceMatrix = pdist2(dynamicsStream, dynamicsStream);
    end
    if ~exist('verbose', 'var') || isempty(verbose)
        verbose = 0;
    end
%     verbose=1
    
    test = meanDistanceMatrix(1:end,1:end);
    totalValues = test(calcIndex,:);

    [peaks, peakIDs] = findpeaks(-totalValues);
    peaks = -peaks;

    [sortedPeaks, sortIDs] = sort(peaks);

    sortedPeakTimes = peakIDs(sortIDs);

    i = 1;
    while i < length(sortedPeakTimes)
        tempPeakTimes = sortedPeakTimes(i+1:end);
        repeatIndices = find(tempPeakTimes > sortedPeakTimes(i) - MIN_RETURN_TIME & tempPeakTimes < sortedPeakTimes(i) + MIN_RETURN_TIME);

        sortedPeakTimes(repeatIndices + i) = [];


        i = i + 1;
    end

    %%

    MIN_NEIGHBOR_COUNT = 0;
    BEST_TRAJECTORY_COUNT = 1;

    currentPeakCount = 2;

    maxDistances = [];
    bestSigmas = [];
    bestProjections = [];
    bestDistances = [];
    peakIndices = [];
    peakValues = [];
    bestValues = [];
    worstDistances = [];
    minDistances = [];
    neighborhoodCounts = [];
    diffSigmas = [];
    
    for i = 1:10
        if i == 1                
            currentPeaks = [peakIDs(sortIDs(1:currentPeakCount)), peakIDs(sortIDs(1:currentPeakCount)) + 1];
        else
            currentPeaks = peakIDs(sortIDs(1:currentPeakCount));
        end

        currentPCAWeight = 1;
        globalMean = mean(dynamicsStream);
        globalSTD = std(dynamicsStream - globalMean);
        currentMean = mean(dynamicsStream(currentPeaks,:));
%         if calcIndex < size(dynamicsStream,1)
%             derivative = dynamicsStream(calcIndex,:) - dynamicsStream(calcIndex + 1,:);
%             currentSTD = sum([derivative; BEST_TRAJECTORY_COUNT*globalSTD; currentPCAWeight*std(dynamicsStream(currentPeaks,:) - currentMean)]) / (1 + currentPCAWeight + BEST_TRAJECTORY_COUNT);
%         else
            currentSTD = sum([BEST_TRAJECTORY_COUNT*globalSTD; currentPCAWeight*std(dynamicsStream(currentPeaks,:) - currentMean)]) / (currentPCAWeight + BEST_TRAJECTORY_COUNT);
%             currentSTD = std(dynamicsStream(currentPeaks,:) - currentMean);
%         end
        %         currentSTD = std(dynamicsStream(currentPeaks,:));
        currentSTD = currentSTD / sqrt(sum(currentSTD.^2));
%         currentSTD = ones(size(currentSTD));
        currentStream = dynamicsStream ./ currentSTD;
        currentStream(isinf(currentStream)) = 0.0001;

%         thisValues = pdist2(currentStream(calcIndex,:), currentStream);
        thisValues = pdist2(currentStream(currentPeaks,:), currentStream);
%         thisValues = thisValues(1,:);
%         nextValues = pdist2(currentStream(calcIndex+1,:), currentStream(2:end,:));
        
%         sortedValues = sort(thisValues);
%         thisProbabilities = cumsum(thisValues) / sum(thisValues);

        smoothSigma = 1;
        x = -ceil(3*smoothSigma):ceil(3*smoothSigma);
        kernel = -x.*exp(-x.^2/(2*smoothSigma^2))/(smoothSigma^3*sqrt(2*pi));

        currentDiff = filterData(currentStream, 0, kernel, 1, 0);
%         thisDiffDistances = pdist2(currentDiff(calcIndex,:), currentDiff, 'cosine');
        thisDiffDistances = pdist2(currentDiff(currentPeaks,:), currentDiff, 'cosine');
%         thisDiffDistances = thisDiffDistances(1,:);
%         sortedDiffDistances = sort(thisDiffDistances);
%         
%         diffProbabilities = cumsum(sortedDiffDistances) / sum(sortedDiffDistances);
        
        
%         diffCurvature = filterData(filterData(sortedDiffDistances, 0, kernel, 1, 0), 50);
%         [~, bestDiffSigmaIndex] = min(diffCurvature);
%         bestDiffSigma = sortedDiffDistances(ceil(bestDiffSigmaIndex/2));
%         guassianDiff = exp(-thisDiffDistances.^2/(2*bestDiffSigma^2));
%         guassianDiff(guassianDiff < 0.1) = 0.1;

%         clf
%         plot(sort(thisValues)/max(thisValues))
%         hold on
%         plot(sort(thisDiffDistances)/max(thisDiffDistances))

        totalValues = 1 - (1 - thisValues ./ max(thisValues,[],2)) .* (1 - thisDiffDistances ./ max(thisDiffDistances,[],2));
    
        allValues = totalValues(1,:);
        totalValues = totalValues(1,:);
        
%         totalValues = totalValues(1,:)
        
        allPeaks = [];
        allPeakIDs = [];
        for j = 1%:size(thisValues,1)
            [peaks, peakIDs] = findpeaks(-totalValues(j,:));            
            peakIDs = peakIDs + size(totalValues,2)*(j-1);
            
            allPeaks = [allPeaks, peaks];
            allPeakIDs = [allPeakIDs, peakIDs];
            
            allPeaks = -allPeaks;
        end
        
        [sortedPeaks, sortIDs] = sort(allPeaks);

        sortedPeakTimes = allPeakIDs(sortIDs);

%         j = 1;
%         while j < length(sortedPeakTimes)
%             tempPeakTimes = sortedPeakTimes(j+1:end);
%             repeatIndices = find(tempPeakTimes > sortedPeakTimes(j) - MIN_RETURN_TIME & tempPeakTimes < sortedPeakTimes(j) + MIN_RETURN_TIME);
% 
%             sortedPeakTimes(repeatIndices + j) = [];
% 
% 
%             j = j + 1;
%         end

        if currentPeakCount >= length(sortedPeakTimes)
            break;
        end
        
        sortedDistances = totalValues(1,sortedPeakTimes);

        totalErrors = [];
        sigmas = [];
        for j = 1:length(sortedPeakTimes)%min(length(sortedPeakTimes), currentPeakCount*2)            
            sigma = sortedDistances(j);
            
            gaussianDistances = exp(-allValues.^2/(2*sigma^2));
                
            error = 1 - gaussianDistances;
            inCluster = gaussianDistances(1,:) >= exp(-1);
            outCluster = gaussianDistances(1,:) < exp(-1);

            smallestValue = min(min(gaussianDistances));
            if sum(outCluster(1,:)) > 0
                totalError = median(gaussianDistances(:, inCluster)) / (mean(gaussianDistances(:, outCluster)) + exp(-1));
            else
                totalError = median(gaussianDistances(:, inCluster)) / (smallestValue + exp(-1));
            end
            
            totalErrors(j) = totalError;
            sigmas(j) = sortedDistances(j);
        end
        
        totalErrors = filterData(totalErrors, 1);
        [totalError, bestIndex] = max(totalErrors);
        
        sigma = sigmas(bestIndex);
        
%         bestTraceValues = -findpeaks(-allValues);
        
        gaussianDistances = exp(-allValues.^2/(2*sigma^2));
                
        error = 1 - gaussianDistances;
        inCluster = gaussianDistances(1,:) >= exp(-1);
        outCluster = gaussianDistances(1,:) < exp(-1);

        smallestValue = min(min(gaussianDistances));
        if sum(outCluster(1,:)) > 0
            totalError = median(gaussianDistances(:, inCluster)) / (mean(gaussianDistances(:, outCluster)) + exp(-1));
        else
            totalError = median(gaussianDistances(:, inCluster)) / (smallestValue + exp(-1));
        end
        
        if i == 7
            test = 1;
        end
        
%         slope = filterData(filterData(sortedDistances, 0, kernel, 1, 0), 50);
%         [~, bestDiffSigmaIndex] = min(slope);
%         sigma = sortedDiffDistances(ceil(bestDiffSigmaIndex/2));
        
%         if i == 1
%             sigma = thisValues(sortedPeakTimes(2));
%         else
%             sigma = thisValues(sortedPeakTimes(ceil(currentPeakCount/2+1)));
%         end
%         gaussianDistances = exp(-thisValues.^2/(2*sigma^2));


        diffDistances = diff(sortedDistances);

%         if i == 1
%             calcTimes = sortedPeakTimes(1);
%             calcTimes = [calcTimes calcTimes + 1];
%             worstDistance = pdist2(currentStream(calcTimes,:),currentStream(calcTimes,:));
% %         elseif sum(inCluster) == 1
% %             calcTimes = sortedPeakTimes(inCluster);
% %             calcTimes = [calcTimes calcTimes + 1];
% %             worstDistance = pdist2(currentStream(calcTimes,:),currentStream(calcTimes,:));
%         else
%             calcTimes = inCluster;
%             worstDistance = pdist2(currentStream(calcTimes,:),currentStream(calcTimes,:));
%         end
%         worstDistance(logical(eye(size(worstDistance)))) = max(max(worstDistance));
%         minDistance = min(min(worstDistance));
%         worstDistance = max(min(worstDistance));

        totalError = totalError;

        peakSortedDistances = totalValues(1,sortedPeakTimes);
        peakDistances = exp(-peakSortedDistances.^2/(2*sigma^2));
        newCount = sum(peakDistances >= exp(-1));

        bestValues{end+1} = gaussianDistances(1,:);
        peakIndices{end+1} = allPeakIDs;
        peakValues{end+1} = totalValues(allPeakIDs);
        maxDistances(end+1) = totalError;
%         worstDistances(end+1) = worstDistance;
%         minDistances(end+1) = minDistance;
        bestSigmas(end+1) = sigma;
        bestProjections(:,end+1) = currentSTD;
        bestDistances(:,end+1) = exp(-totalValues(1,:).^2/(2*sigma^2));
        neighborhoodCounts(end+1) = currentPeakCount;
%         diffSigmas(end+1) = bestDiffSigma;

        if verbose
            figure(3+i)
            clf;
            plot(sort(gaussianDistances, 'descend'));
            title(['k = ' num2str(currentPeakCount)]);
        end

        currentPeakCount = ceil(max([currentPeakCount+1, currentPeakCount*1.2, newCount]));
        
        if currentPeakCount >= length(sortedPeakTimes)
            break;
        end
    end
    
    finalScore = maxDistances;% .* minDistances ./ worstDistances;% .* -log(minDistances ./ worstDistances);

%     distances = filterData(finalScore(MIN_NEIGHBOR_COUNT+1:end), 1, 'gaussian', false, 0);
    distances = finalScore(MIN_NEIGHBOR_COUNT+1:end);
    if 0%length(distances) > 2
        [~,peakIndices] = findpeaks(distances);
        if isempty(peakIndices)
            peakIndices = length(finalScore) - 1;
        end
        bestIndex = peakIndices(1);
    else
        [~, bestIndex] = max(distances);
    end
    [~, overrideIndex] = max(distances);
    if overrideIndex == 1
        bestIndex = 1;
    end
    
    if verbose
        figure(3+bestIndex)
%         clf;
%         plot(sort(gaussianDistances, 'descend'));
        title(['Best k = ' num2str(neighborhoodCounts(bestIndex))]);
    end

    if verbose
        figure(3);
        clf;
        hold on;
        plot(finalScore(MIN_NEIGHBOR_COUNT+1:end) / max(finalScore(MIN_NEIGHBOR_COUNT+1:end)));
%         plot(distances(MIN_NEIGHBOR_COUNT+1:end) / max(distances(MIN_NEIGHBOR_COUNT+1:end)));
%         plot(bestSigmas(MIN_NEIGHBOR_COUNT+1:end) / max(bestSigmas(MIN_NEIGHBOR_COUNT+1:end)));
%         plot(worstDistances(MIN_NEIGHBOR_COUNT+1:end) / max(worstDistances(MIN_NEIGHBOR_COUNT+1:end)));
    end
    
    distanceValues = bestDistances(:,bestIndex+MIN_NEIGHBOR_COUNT);
    plotIndices = find(distanceValues > exp(-1));

    [plotBasis, ~] = pca(currentStream, 'NumComponents', min(4, size(currentStream,2)));

    if verbose
        figure(2);
        clf;
        hold on;
        colors = lines(8);
        h = plot3(currentStream*plotBasis(:,1),currentStream*plotBasis(:,2),currentStream*plotBasis(:,3));
        h.Color(4) = 0.2;
        scatter3(currentStream(plotIndices,:)*plotBasis(:,1),currentStream(plotIndices,:)*plotBasis(:,2),currentStream(plotIndices,:)*plotBasis(:,3), 32, distanceValues(plotIndices))
        xlabel('DPC1');
        ylabel('DPC2');
        zlabel('DPC3');
        colormap(jet(256));
    end
    
    test = 1;

    bestLocalProjections = bestProjections(:,bestIndex+MIN_NEIGHBOR_COUNT);
    bestSigmaValues = bestSigmas(bestIndex+MIN_NEIGHBOR_COUNT);
    localDistances = bestValues{bestIndex+MIN_NEIGHBOR_COUNT};
    bestDiffSigmaValues = 0;%minDistances(bestIndex+MIN_NEIGHBOR_COUNT) / worstDistances(bestIndex+MIN_NEIGHBOR_COUNT);
    bestNeighborCounts = neighborhoodCounts(bestIndex+MIN_NEIGHBOR_COUNT);
    distanceRatios = 0;%minDistances(bestIndex+MIN_NEIGHBOR_COUNT) / worstDistances(bestIndex+MIN_NEIGHBOR_COUNT);
    
end