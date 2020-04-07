
finalDynamicsStream = saveData.FinalStream;
reducedMatrix = saveData.ReducedMatrix;
clusterIDs = saveData.ClusterIDs;
clusterMeansPCA = saveData.ClusterMeansPCA;

trialLength = size(finalDynamicsStream,1)/length(saveData.TrialSwitches);
meanClusterTimes = [];
for i = 1:size(reducedMatrix,1)
    thisIndices = find(clusterIDs == i);
    thisTimes = mod((thisIndices-1),trialLength);
    
    meanClusterTimes(i) = mean(thisTimes);
end

% firstComplexMode = find(imag(clusterEigenValues) > 0, 1, 'first');
% for plotCount = 1:2
    plotCount = 2;
    if plotCount == 0
        phase = 0;
    else
%         phase = angle(clusterEigenVectors(:,firstComplexMode)) + pi;
        phase = meanClusterTimes / trialLength * 2*pi;
        
        phase(phase > pi) = phase(phase > pi) - 2*pi;
    end
    
    SHOW_PHASE = 1;
    PLOT_PHASE = pi+2*pi/100*48;
    PHASE_SIGMA = 1000;%2*pi/50;
%     phaseIndices = exp(-angleDiff(phase(clusterIDs), PLOT_PHASE).^2/(2*PHASE_SIGMA^2)) > exp(-1);
    clusterWeights = exp(-angleDiff(phase, PLOT_PHASE).^2/(2*PHASE_SIGMA^2));
    phaseClusters = clusterWeights > exp(-1);
    
    figure(3+plotCount-1);
    clf;
    if SHOW_PHASE
%         scatter3(clusterMeans(phaseClusters,1), clusterMeans(phaseClusters,2), clusterMeans(phaseClusters,3), 32, phase(phaseClusters));
        colormap(hsv(256));
        colors = hsv(256);
        caxis([-pi pi]);
        
        phaseColors = colors(floor((phase + pi) / (2*pi) * 255)+1, :);
    else
%         scatter3(clusterMeans(:,1), clusterMeans(:,2), clusterMeans(:,3), 70, 1:length(phase), 'x', 'LineWidth', 3);
        colormap(lines(256));
    end
%     colormap(lines(256));
    hold on;
%     h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
%     h.Color(4) = 0.2;
    for i = 1:size(reducedMatrix,1)
        if ~validClusters(i) || ~phaseClusters(i)
            continue
        end
        
        transitionProbabilities = reducedMatrix(i,:);
        transitionProbabilities(i) = 0;
        transitionProbabilities(~validClusters) = 0;
        [bestTransferProbability, bestNextState] = max(transitionProbabilities);
        
        thisIndices = find(clusterIDs == i);
        
        for j = 1:length(transitionProbabilities)
            if transitionProbabilities(j) > 0.1
                thisTransition = [clusterMeansPCA(i,:); clusterMeansPCA(j,:)];
                
                colors = jet(256);
                thisColor = ceil(transitionProbabilities(j) * 256);
                
                plot3(thisTransition(:,1), thisTransition(:,2), thisTransition(:,3), 'LineWidth', transitionProbabilities(j) * 10, 'Color', colors(thisColor,:));
            end
        end
        
        [x, y, z] = sphere;
        if SHOW_PHASE
            surf(x.*clusterSTDsPCA(i,1) + clusterMeansPCA(i,1),y.*clusterSTDsPCA(i,2) + clusterMeansPCA(i,2),z.*clusterSTDsPCA(i,3) + clusterMeansPCA(i,3), ones(size(x)) * phase(i))
            
            thisPlot = finalDynamicsStream;
            badIndices = 1:size(finalDynamicsStream,1);
            badIndices(thisIndices) = [];
            thisPlot(badIndices,:) = nan;
            h = plot3(thisPlot*pcaBasis(:,1), thisPlot*pcaBasis(:,2), thisPlot*pcaBasis(:,3), 'Color', phaseColors(i,:));
        else
            surf(x.*clusterSTDsPCA(i,1) + clusterMeansPCA(i,1),y.*clusterSTDsPCA(i,2) + clusterMeansPCA(i,2),z.*clusterSTDsPCA(i,3) + clusterMeansPCA(i,3), ones(size(x)) * i)
            h = scatter3(finalDynamicsStream(thisIndices,:)*pcaBasis(:,1), finalDynamicsStream(thisIndices,:)*pcaBasis(:,2), finalDynamicsStream(thisIndices,:)*pcaBasis(:,3), 32, i*ones(size(finalDynamicsStream(thisIndices,3))));
        end
        h.Color(4) = 0.8;
    end
    % scatter3(finalDynamicsStream(:,1), finalDynamicsStream(:,2), finalDynamicsStream(:,3), 32, allAngles);
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    if plotCount == 1
        title('Reduced dynamics - mean phase');
    else
        title('Reduced dynamics - reduced phase');
    end
% end

figure(6);
clf
imagesc(reducedMatrix);
title('Reduced markov matrix');
ylabel('From state');
xlabel('To state');
colormap(parula(256));
