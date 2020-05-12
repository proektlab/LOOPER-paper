
% finalDynamicsStream
% validationModel
% validationEmission

% returns scoreMean, scoreSTD

maxTime = 10;

scoreMean = 0;
scoreSTD = 0;

if size(validationEmission, 2) ~= size(finalDynamicsStream,2)
    disp(['Emission matrix must be STATES by CHANNELS (' num2str(size(validationModel,2)) ', ' num2str(size(finalDynamicsStream,2)) ')']);
    return;
end

if size(validationEmission, 1) ~= size(validationModel,2)
    disp(['Emission matrix must be STATES by CHANNELS (' num2str(size(validationModel,2)) ', ' num2str(size(finalDynamicsStream,2)) ')']);
    return;
end

disp('Validating...');

USE_PCA = 0;

noInputDynamics = finalDynamicsStream;%(:,1:100);
noInputEmissions = validationEmission;%(:,1:100);

dynamicsMean = mean(noInputDynamics,1);
dynamicsSTD = std(noInputDynamics - dynamicsMean,1);

noInputDynamics = (noInputDynamics - dynamicsMean) ./ dynamicsSTD;
noInputEmissions = (noInputEmissions - dynamicsMean) ./ dynamicsSTD;

noInputDynamics(isnan(noInputDynamics)) = 0;
noInputEmissions(isnan(noInputEmissions)) = 0;

if USE_PCA
    [finalPCABasis, ~, ~, ~, explained] = pca(noInputDynamics);
    requiredDimensions = round(find(cumsum(explained) > 95, 1, 'first'));

    finalPCABasis = finalPCABasis(:, 1:requiredDimensions);
else
    finalPCABasis = eye(size(noInputDynamics, 2));
end

emissionDistances = pdist2(noInputEmissions, noInputDynamics).^2;
emissionVelocities = repmat(permute(noInputEmissions, [3,1,2]), [size(noInputEmissions,1),1,1]) - repmat(permute(noInputEmissions, [1,3,2]), [1,size(noInputEmissions,1),1]);

predictTime = 1;

possibleTransitions = validationModel;
possibleTransitions = possibleTransitions^predictTime;
% currentTransitions = possibleTransitions;
% possibleTransitionsList = {};
% for i = 1:maxTime
%     possibleTransitionsList{i} = currentTransitions;
%     possibleTransitionsList{i}(possibleTransitionsList{i} > 1) = 1;
%     possibleTransitionsList{i}(possibleTransitionsList{i} == 0) = 100000;
%     
%     currentTransitions = currentTransitions * possibleTransitions;
% end

fromStates = [];
toStates = [];
scores = [];
for i = 1:size(noInputDynamics,1)-predictTime
    distances = emissionDistances(:,i);
    nextDistances = emissionDistances(:,i+predictTime);
    
%     min(1 - possibleTransitions(:))
    
    transitionScores = log(((distances .* repmat(nextDistances, [1, size(nextDistances,1)])) ./ possibleTransitions));
    transitionScores(isinf(transitionScores)) = 100000;
    
    scores(i) = nanmin(transitionScores(:));
    
    index = find(transitionScores == scores(i), 1);
    [currentState, nextState] = ind2sub(size(transitionScores), index);
    
    fromStates(i) = currentState;
    toStates(i) = nextState;
%     disp([num2str(currentState) ' -> ' num2str(nextState) ' = ' num2str(scores(i))]);
%     position = noInputDynamics(i,:);
%     velocity = noInputDynamics(i+1,:) - position;
%     
%     distances = emissionDistances(:,i);
%     velocityDistances = emissionVelocities - permute(velocity, [3,1,2]);
%     velocityDistances = sqrt(sum(velocityDistances.^2,3));
%     velocityDistances = min(velocityDistances,[],1);
%     
%     emissionsTotals = distances + velocityDistances;
%     
%     velocityDistance = pdist2(velocity);
end


scoreMean = mean(scores);
scoreSTD = std(scores);

%%

[pcaBasis, pcaOutputs] = pca(finalDynamicsStream, 'NumComponents', 3);

figure(1);
clf;
hold on;
h = plot3(finalDynamicsStream*pcaBasis(:,1), finalDynamicsStream*pcaBasis(:,2), finalDynamicsStream*pcaBasis(:,3), 'LineWidth', 0.5);
h.Color(4) = 0.2;

colors = jet(256);
scoreIndices = scores - min(scores);
scoreIndices = scoreIndices / max(scoreIndices);
scoreIndices = round(1 + scoreIndices * 255);

plotDynamics = finalDynamicsStream(1:end-predictTime,:);

scatter3(plotDynamics*pcaBasis(:,1), plotDynamics*pcaBasis(:,2), plotDynamics*pcaBasis(:,3), 32, scores);
colormap(jet)

scatter3(validationEmission*pcaBasis(:,1), validationEmission*pcaBasis(:,2), validationEmission*pcaBasis(:,3), 100, 'kx', 'LineWidth', 3);


for i = 1:5:length(toStates)
    lines = validationEmission([fromStates(i),toStates(i)],:)*pcaBasis(:,1:3);
    
    plot3(lines(:,1), lines(:,2), lines(:,3), 'k', 'LineWidth', 1);
end


