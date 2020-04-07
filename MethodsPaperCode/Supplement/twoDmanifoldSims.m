%%

numFluxes = 5;

figureHandle = figure(1);
figureHandle.Renderer='Painters';
clf
hold on
for i = 1:numFluxes
    x = sin(0:2*pi/100:pi);
    y = 0:2*pi/100:pi;
    z = ones(size(x))*(i-0.5)*pi/numFluxes;
    
    plot3(x, y, z, 'LineWidth', 4, 'Color', colors(i,:));
end

[X,Y] = meshgrid(0:2*pi/100:pi, 0:2*pi/100:pi);
Z = sin(Y);
surf(Z,Y,X, ones(size(Z)))
light('Position',[2 3 3],'Style','local')
lighting gouraud
shading interp

zlim([min(X(:)) - 2, max(X(:)) + 2]);
xlim([min(Z(:)) - 2, max(Z(:)) + 2]);
ylim([min(Y(:)) - 2, max(Y(:)) + 2]);

%%
noise = 0.3;

starts = [-2 -1 0 1 2];
% starts = [-0 -0 0 0 0];
simTrials = 50;

traces = repmat(starts, [1 1 simTrials]);
for t = 2:100
    for i = 1:simTrials
        traces(t,:,i) = traces(t-1,:,i) + noise*normrnd(0,1,size(starts));
        
%         thisTraces = traces(t,:,i);
%         thisTraces(thisTraces > 4) = 4;
%         thisTraces(thisTraces < -4) = -4;
%         traces(t,:,i) = thisTraces;
    end
end

figureHandle = figure(3);
figureHandle.Renderer='Painters';
clf
subplot(4,1,1:3);
hold on
for i = 1:length(starts)
    meanValues = mean(squeeze(traces(:,i,:)),2);
    stdValues = std(squeeze(traces(:,i,:)),[],2);
    
    plotShadedCI(meanValues, stdValues, [], colors(i,:));
    
%     for j = 1:size(traces,3)
%         plot(traces(:,i,j), 'Color', colors(i,:));
%     end
end

data = squeeze(traces(1,:,:))';
startVector = data(:);

if ~exist('nullModel')
    nullModel = [];
    for i = 1:10000
        dataVector = rand(1, length(startVector));

        decodeRate = [];
        for j = 1:length(dataVector)
            greaterThanTarget = find(startVector > startVector(j));
            lessThanTarget = find(startVector < startVector(j));

            maxDecodes = length(greaterThanTarget) + length(lessThanTarget);

            greaterThan = find(dataVector > dataVector(j));
            lessThan = find(dataVector < dataVector(j));

            correctDecodes = sum(ismember(greaterThan, greaterThanTarget)) + sum(ismember(lessThan, lessThanTarget));

            decodeRate(j) = correctDecodes / maxDecodes;
        end

        nullModel(i) = mean(decodeRate);
    end
end

decodeRateTotal = [];
for i = 1:size(traces,1)
    data = squeeze(traces(i,:,:))';
    dataVector = data(:);
    
    decodeRate = [];
    for j = 1:length(dataVector)
        greaterThanTarget = find(startVector > startVector(j));
        lessThanTarget = find(startVector < startVector(j));

        maxDecodes = length(greaterThanTarget) + length(lessThanTarget);

        greaterThan = find(dataVector > dataVector(j));
        lessThan = find(dataVector < dataVector(j));
        
        correctDecodes = sum(ismember(greaterThan, greaterThanTarget)) + sum(ismember(lessThan, lessThanTarget));
        
        decodeRate(j) = correctDecodes / maxDecodes;
    end
    
    decodeRateTotal(i) = mean(decodeRate);
%     decodeRateTotal(i) = 1 - sum(mean(decodeRate) > nullModel) / length(nullModel);
end

subplot(4,1,4);
plot(decodeRateTotal)

%%
noise = 0.07;
restoreForce = 0.03;

starts = [-2 -1 0 1 2];

colors = parula(length(starts));


traces = starts;
for t = 2:100
    traces(t,:) = traces(t-1,:) - restoreForce*sign(traces(t-1,:)).*(abs(traces(t-1,:)))  + noise*normrnd(0,1,size(starts));
end


figureHandle = figure(3);
figureHandle.Renderer='Painters';
clf
hold on
for i = 1:length(starts)
    plot(traces(:,i), 'Color', colors(i,:));
end


traces = starts;
for t = 2:100
    traces(t,:) = traces(t-1,:) - restoreForce*sign(traces(t-1,:)).*(abs(traces(t-1,:)))  + 0.00*normrnd(0,1,size(starts));
end


figureHandle = figure(4);
figureHandle.Renderer='Painters';
clf
hold on
for i = 1:length(starts)
    plot(traces(:,i), 'Color', colors(i,:));
end


targets = starts/2;
traces = starts;
for t = 2:100
    traces(t,:) = traces(t-1,:) - restoreForce*sign(traces(t-1,:) - targets).*(abs(traces(t-1,:) - targets))  + noise*normrnd(0,1,size(starts));
end


figureHandle = figure(5);
figureHandle.Renderer='Painters';
clf
hold on
plot(traces);
for i = 1:length(starts)
    plot(traces(:,i), 'Color', colors(i,:));
end