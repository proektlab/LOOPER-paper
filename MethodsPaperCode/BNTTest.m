
% app.SavedData = saveData;
% 
% numTrial = max(app.SavedData.TrialData);
% trialLength = size(app.SavedData.FinalStream,1) / numTrial;
% trialStarts = 1:trialLength:size(app.SavedData.FinalStream,1);
% 
% rawStream = app.SavedData.FinalStream;
% 
% data = [];
% for i = 1:numTrial
%     finalindices = trialStarts(i):trialStarts(i)+trialLength-1;
%     
%     data(:,:,i) = rawStream(finalindices,1:20)';
% end

data = finalTrials;


O = size(data,1);          %Number of coefficients in a vector 
T = size(data,2);         %Number of vectors in a sequence 
nex = size(data,3);        %Number of sequences 
M = 1;          %Number of mixtures 
Q = 70;          %Number of states 
cov_type = 'diag';

matrixPrior = zeros(Q,Q);
% matrixPrior = diag(ones(Q-1,1),1);
% matrixPrior(Q,1) = 1;
% matrixPrior = kron(eye(3),matrixPrior) * 5000/Q/3;
% matrixPrior = matrixPrior + eye(size(matrixPrior)) * 5000/Q/6;


intra = zeros(2);
intra(1,2) = 1; % node 1 in slice t connects to node 2 in slice t

inter = zeros(2);
inter(1,1) = 1; % node 1 in slice t-1 connects to node 1 in slice t

ns = [Q O];
dnodes = [1];
onodes = [2];

eclass1 = [1 2];
eclass2 = [3 2];
eclass = [eclass1 eclass2];

bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'observed', onodes, 'eclass1', eclass1, 'eclass2', eclass2);
prior0 = normalise(rand(Q,1));
% transmat0 = eye(Q,Q);
transmat0 = mk_stochastic(rand(Q,Q));
obsmat0 = mk_stochastic(rand(Q,O));
bnet.CPD{1} = tabular_CPD(bnet, 1, prior0);
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', obsmat0, 'cov_type', 'diag');
bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', transmat0);
% bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', transmat0, 'prior_type', 'dirichlet', 'dirichlet_type', ...
%    'prior', 'dirichlet_weight', matrixPrior);

engine = bk_inf_engine(bnet, 'clusters', 'exact');
% engine = hmm_inf_engine(bnet);

ss = 2;%slice size(ss)
max_iter=500;%iterations for EM
cases = cell(1, nex);
for i=1:nex
  cases{i} = cell(ss,T);
  for t = 1:T
    cases{i}{2,t} = squeeze(data(:, t, i));
  end
end

% fakeDataCount = floor(T/3/3);
% for i = 1:3
%     cases{end+1} = cell(ss,fakeDataCount);
%     for t = 1:fakeDataCount
%         cases{end}{1,t} = mod(t, 100) + 1 + 100*(i-1);
%     end
% end
[bnet2, LLtrace] = learn_params_dbn_em(engine, cases, 'max_iter', max_iter);


%% Build loop model

CPT1 = get_field(bnet2.CPD{1}, 'cpt');
means = get_field(bnet2.CPD{2}, 'mean');
covs = get_field(bnet2.CPD{2}, 'cov');
CPT3 = get_field(bnet2.CPD{3}, 'cpt');

% loglikHMM = mhmm_logprob(validationData, CPT1, CPT3, means, covs, ones(Q,1))

stateIDs = [];
for i = 1:size(data,3)
    B = mixgauss_prob(data(:,:,i), means, covs, ones(size(CPT3,1),1));
    
    stateIDs(:,i) = viterbi_path(CPT1, CPT3, B);
end

% figure(1);
% clf;
% hold on;
% for i = 1:size(data,3)
%     h = plot(data(1,:,i), data(2,:,i), 'k');
%     h.Color(4) = 0.3;
%     scatter(data(1,:,i), data(2,:,i), 32, stateIDs(:,i));
% end

stateSimilarities = [];
for i = 1:size(CPT3,1)
    varsi = diag(covs(:,:,i));
    
    for j = 1:size(CPT3,1)
        if i == j
            stateSimilarities(i,j) = 1;
        else
            varsj = diag(covs(:,:,j));

            stateDiff = abs(means(:,i) - means(:,j));
            normDiff = stateDiff ./ norm(stateDiff);

            totalVars = sum(abs((varsi + varsj) .* stateDiff));

            stateSimilarities(i,j) = exp(-(norm(stateDiff)).^2/(2*totalVars));
        end
    end
end


finalDynamicsStream = reshape(data, size(data,1), [])';
reducedMatrix = CPT3;
clusterIDs = reshape(stateIDs, 1, [])';
clusterMeans = means;

countClusters = [];
for i = 1:size(CPT3,1)
    countClusters(i) = sum(clusterIDs == i);
end

if size(finalDynamicsStream,2) > 2
    [~, pcaBasis] = pca(finalDynamicsStream, 'NumComponents', 3);
else
    pcaBasis = eye(size(clusterMeans,1),3);
end

clusterMeansPCA = (clusterMeans' * pcaBasis);

putativeLoopCounts = 6;  
shouldUseTerminalState = true;
totalClusters = 200;

selectingLoops = 0;

trialSwitchTimes = [];
endTime = 1;
for i = 1:size(data,3)
    trialSwitchTimes(i) = endTime;
    
    endTime = endTime + size(data,2);
end

buildMinimalModelFromMarkovMatrix;

% app.SavedData.BestStateCount = bestStateCount;
% app.SavedData.LogLikelihoods = allLikelihoods;
% app.SavedData.BestLoopCount = bestLoopCount;
% app.SavedData.BestModel = bestModel;
% app.SavedData.BestEmission = bestEmission;
% app.SavedData.BestLoopAssignments = bestLoopAssignments;
% app.SavedData.BestStateMap = stateMap;



%% Load LOOPER markov model

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
trialStarts = 1:trialLength:size(app.SavedData.FinalStream,1);

stateModel = app.SavedData.BestModel;

stateMap = [];
for i = 1:size(stateModel,1)
    loopID = app.SavedData.BestLoopAssignments(i,1);
    phaseID = app.SavedData.BestLoopAssignments(i,2);
    
    stateIndices = find(app.SavedData.BestStateMap(:,1) == loopID & app.SavedData.BestStateMap(:,2) == phaseID);

    stateMap(stateIndices) = i;
end

rawData = app.SavedData.RawData;

trialStates = [];
trialActivity = [];
for i = 1:numTrial
    thisIndices = find(app.SavedData.TrialData == i);
    finalindices = trialStarts(i):trialStarts(i)+trialLength-1;
    
    thisData = rawData(:, thisIndices);
    
    finalTimes = size(thisData,2)-trialLength+1-2:size(thisData,2)-2;
    
    trialActivity = [trialActivity thisData(:,finalTimes)];
    trialStates = [trialStates stateMap(finalindices)];
end

stateMeans = [];
stateVars = [];
for i = 1:size(stateModel,1)
    loopID = app.SavedData.BestLoopAssignments(i,1);
    phaseID = app.SavedData.BestLoopAssignments(i,2);
    
    stateIndices = find(app.SavedData.BestStateMap(:,1) == loopID & app.SavedData.BestStateMap(:,2) == phaseID);

    stateMeans(:,i) = mean(trialActivity(:,stateIndices),2);
    
    networkNoise = 2;
    stateVars(:,:,i) = diag(max(networkNoise, var(trialActivity(:,stateIndices),[],2)));
end

startDistribution = zeros(size(stateModel,1),1);
startState = mode(trialStates(trialStarts));
startDistribution(startState) = 1;


%% Load LOOPER diffusion map model

app.SavedData = saveData;

numTrial = max(app.SavedData.TrialData);
trialLength = size(app.SavedData.FinalStream,1) / numTrial;
trialStarts = 1:trialLength:size(app.SavedData.FinalStream,1);

stateModel = app.SavedData.ReducedMatrix;

stateMap = app.SavedData.ClusterIDs;

rawData = app.SavedData.RawData;

trialStates = [];
trialActivity = [];
for i = 1:numTrial
    thisIndices = find(app.SavedData.TrialData == i);
    finalindices = trialStarts(i):trialStarts(i)+trialLength-1;
    
    thisData = rawData(:, thisIndices);
    
    finalTimes = size(thisData,2)-trialLength+1-2:size(thisData,2)-2;
    
    trialActivity = [trialActivity thisData(:,finalTimes)];
    trialStates = [trialStates stateMap(finalindices)];
end

stateMeans = [];
stateVars = [];
for i = 1:size(stateModel,1)
    stateIndices = find(app.SavedData.ClusterIDs == i);

    stateMeans(:,i) = mean(trialActivity(:,stateIndices),2);
    
    networkNoise = 2;
    stateVars(:,:,i) = diag(max(networkNoise, var(trialActivity(:,stateIndices),[],2)));
end

startDistribution = zeros(size(stateModel,1),1);
startState = mode(trialStates(trialStarts));
startDistribution(startState) = 1;

%% Find loops in Markov model

app.SavedData = saveData;

CPT1 = get_field(bnet2.CPD{1}, 'cpt');
means = get_field(bnet2.CPD{2}, 'mean');
covs = get_field(bnet2.CPD{2}, 'cov');
CPT3 = get_field(bnet2.CPD{3}, 'cpt');

asymmetricProbabilities = app.SavedData.AsymmetricMap;
finalDynamicsStream = app.SavedData.FinalStream;
reducedMatrix = CPT3;
clusterIDs = app.SavedData.ClusterIDs;
stateHasNext = app.SavedData.StateHasNext;
stateValidities= app.SavedData.ValidStates;
clusterMeansPCA = app.SavedData.ClusterMeansPCA;
clusterMeans = app.SavedData.ClusterMeans;
countClusters = app.SavedData.CountClusters;

app.SavedData.PutativeLoopCounts = cell2mat(cellfun(@str2num,strsplit(app.PutativeloopcountsEditField.Value, ','),'uniform',0));    
putativeLoopCounts = app.SavedData.PutativeLoopCounts;

app.SavedData.UseTerminalState = app.TerminalStateCheckbox.Value;    
shouldUseTerminalState = app.SavedData.UseTerminalState;

app.SavedData.TotalStates = app.TargetstatecountEditField.Value;
totalClusters = app.SavedData.TotalStates;

selectingLoops = 0;

[pcaBasis, pcaOutputs] = pca(app.SavedData.FinalStream, 'NumComponents', 3);


trialSwitchTimes = app.SavedData.TrialSwitches;

buildMinimalModelFromMatrix;

app.SavedData.BestStateCount = bestStateCount;
app.SavedData.LogLikelihoods = allLikelihoods;
app.SavedData.BestLoopCount = bestLoopCount;
app.SavedData.BestModel = bestModel;
app.SavedData.BestEmission = bestEmission;
app.SavedData.BestLoopAssignments = bestLoopAssignments;
app.SavedData.BestStateMap = stateMap;

%% Compare models
CPT1 = get_field(bnet2.CPD{1}, 'cpt');
means = get_field(bnet2.CPD{2}, 'mean');
covs = get_field(bnet2.CPD{2}, 'cov');
CPT3 = get_field(bnet2.CPD{3}, 'cpt');

% figure(1)
% clf;
% imagesc(CPT3);

validationData = finalTrials;

loglikHMM = mhmm_logprob(validationData, CPT1, CPT3, means, covs, ones(Q,1))
loglikLOOPER = mhmm_logprob(validationData, startDistribution, stateModel, stateMeans, stateVars, ones(size(stateModel,1),1))


%% Sim markov model
CPT1 = get_field(bnet2.CPD{1}, 'cpt');
means = get_field(bnet2.CPD{2}, 'mean');
covs = get_field(bnet2.CPD{2}, 'cov');
CPT3 = get_field(bnet2.CPD{3}, 'cpt');

validationModel = CPT3;
validationEmissions = means';
validationCovs = covs;

% validationModel = stateModel;
% validationEmissions = stateMeans';
% validationCovs = stateVars;

numStates = length(CPT1);



Y = [];

for j = 1:100
    simulation(1) = randsample(numStates, 1, true, CPT1);

    for i = 2:size(data,2)
        simulation(i) = randsample(numStates, 1, true, CPT3(simulation(i-1),:));
    end

    for i = 1:length(simulation)
        Y(i,:,j) = normrnd(means(:,simulation(i)), sqrt(diag(validationCovs(:,:,simulation(i)))))';
    %     Y(i,:) = means(:,simulation(i))';
    end
end

plotData = data(:,:,1)';
% pcaData = pca(Y', 'NumComponents', 3);
% [pcaBasis, ~] = pca(plotData, 'NumComponents', 3);
pcaBasis = eye(2,3);

figure(1);
clf;
hold on;
for j = 1:100
    h = plot3(Y(:,:,j) * pcaBasis(:,1), Y(:,:,j) * pcaBasis(:,2), Y(:,:,j) * pcaBasis(:,3), 'LineWidth', 0.5, 'Color', 'k');
    h.Color(4) = 0.3;
end
% plot3(plotData * pcaBasis(:,1), plotData * pcaBasis(:,2), plotData * pcaBasis(:,3));

%%

CPT1 = get_field(bnet2.CPD{1}, 'cpt');
means = get_field(bnet2.CPD{2}, 'mean');
covs = get_field(bnet2.CPD{2}, 'cov');
CPT3 = get_field(bnet2.CPD{3}, 'cpt');

validationModel = CPT3;
validationEmissions = means';

plotData = data';
[pcaBasis, ~] = pca(plotData, 'NumComponents', 3);

figure(1);
clf;
h = plot3(plotData * pcaBasis(:,1), plotData * pcaBasis(:,2), plotData * pcaBasis(:,3));
h.Color(4) = 0.2;
hold on;
for i = 1:3
    states = (i-1)*100 + (1:100);
    plot3(validationEmissions(states,:) * pcaBasis(:,1), validationEmissions(states,:) * pcaBasis(:,2), validationEmissions(states,:) * pcaBasis(:,3), 'LineWidth', 2);
end

