
%% Simulate the two neuron system illustarted in the appendix
% Simulate without noise to find limit cycles
t1 = 1;
t2 = 6;
w11 = 8;
w12 = -6;
w21 = 16;
w22 = -2;
b1 = .34;
b2 = 2.5;

limX = [0 1];
limY = [0.2 1];

% limX(1) = limX(1) - 2;
% limX(2) = limX(2) + 2;
% limY(1) = limY(1) - 2;
% limY(2) = limY(2) + 2;

xProbe = linspace(limX(1),limX(2),2000);
yProbe = linspace(limY(1),limY(2),2000);
[xGrid,yGrid] = meshgrid(xProbe, yProbe);


fx = @(x)(-x(:,1) + sigmf(w11*x(:,1) + w12*x(:,2) - b1, [1 0]))/t1;
fy = @(x)(-x(:,2) + sigmf(w21*x(:,1) + w22*x(:,2) - b2, [1 0]))/t2;

dt = 1/100;

limit1x = 0.1;
limit1y = 0.5;
for t = 2:100000
    limit1x(t) = limit1x(t-1) + ((-limit1x(t-1) + sigmf(w11*limit1x(t-1) + w12*limit1y(t-1) - b1, [1 0]))/t1)*dt;
    limit1y(t) = limit1y(t-1) + ((-limit1y(t-1) + sigmf(w21*limit1x(t-1) + w22*limit1y(t-1) - b2, [1 0]))/t2)*dt;
end

limit2x = 0.5;
limit2y = 0.1;
for t = 2:100000
    limit2x(t) = limit2x(t-1) + ((-limit2x(t-1) + sigmf(w11*limit2x(t-1) + w12*limit2y(t-1) - b1, [1 0]))/t1)*dt;
    limit2y(t) = limit2y(t-1) + ((-limit2y(t-1) + sigmf(w21*limit2x(t-1) + w22*limit2y(t-1) - b2, [1 0]))/t2)*dt;
end

%% Simulate network with noise

MAX_TIME = 1000000;
NOISE_SIGMA = 0.1;
TRACE_LENGTH = 3000;

x = zeros(1,MAX_TIME);
y = zeros(1,MAX_TIME);
x(1) = 0.5;
y(1) = 0.1;
for t = 2:MAX_TIME
    x(t) = x(t-1) + ((-x(t-1) + sigmf(w11*x(t-1) + w12*y(t-1) - b1, [1 0]))/t1 + normrnd(0, NOISE_SIGMA))*dt;
    y(t) = y(t-1) + ((-y(t-1) + sigmf(w21*x(t-1) + w22*y(t-1) - b2, [1 0]))/t2 + normrnd(0, NOISE_SIGMA))*dt;
end

%% Plot data as a time series

finalTraceIndices = MAX_TIME-TRACE_LENGTH*10+1:10:MAX_TIME;

figure(1);
clf;
subplot(2,1,1)
plot(zscore(x(finalTraceIndices)), 'r','LineWidth',1);
ylabel('Activity of neuron A (au)');
xlabel('Time (timesteps)');
subplot(2,1,2)
plot(zscore(y(finalTraceIndices)), 'b','LineWidth',1);
ylabel('Activity of neuron B (au)');
xlabel('Time (timesteps)');
suptitle('Activity trace of toy system');

%% Plot data as a a phase plot

traceColor = lines(2);

figure(2);
clf;
xlim(limX);
ylim(limY);
hold on;
quiver(xGrid(:),yGrid(:),fx([xGrid(:),yGrid(:)]),fy([xGrid(:),yGrid(:)]),'k','linewidth',0.5);
plot(x(finalTraceIndices), y(finalTraceIndices), 'Color', traceColor(2,:),'LineWidth',2);
plot(limit1x(end-10000:end),limit1y(end-10000:end),'color','k','linewidth',4);
plot(limit2x(end-10000:end),limit2y(end-10000:end),'color','k','linewidth',4);
xlabel('Activity of neuron A (au)');
ylabel('Activity of neuron B (au)');
title('Phase protrait of toy system');

%%

options.bound = 'per';
mesh_size = 0.1;
% xp = -2.05:mesh_size:2.05;
% yp = -2.05:mesh_size:2.05;
% [x,y]=meshgrid(xp,yp);
% f = x.^2 - y.^2;
% [fx,fy]=gradient(f);

fxGrid = zeros(size(xGrid));
fyGrid = zeros(size(xGrid));
for i = 1:size(xGrid,1)
    for j = 1:size(xGrid,2)
        fxGrid(i,j) = fx([xGrid(i,j), yGrid(i,j)]);
        fyGrid(i,j) = fy([xGrid(i,j), yGrid(i,j)]);
    end
end
fn = sqrt(fxGrid.^2 +fyGrid.^2);
[a,b] = compute_hodge_decompositon(cat(3,fxGrid./fn,fyGrid./fn),options);
figure(1);
clf;
% a = cat(3, fxGrid, fyGrid);
[curlz,cav] = curl(xGrid,yGrid,a(:,:,2),a(:,:,1));
imagesc(curlz)
figure(2);
clf
% b = cat(3, fxGrid, fyGrid);
d = divergence(xGrid,yGrid,b(:,:,2),b(:,:,1));
imagesc(d)

%%

subSampling = 50;

figure(2);
clf;
xlim(limX);
ylim(limY);
hold on;
subSampledX = xGrid(1:subSampling:end,1:subSampling:end);
subSampledY = yGrid(1:subSampling:end,1:subSampling:end);
ax = a(1:subSampling:end,1:subSampling:end,1);
ay = a(1:subSampling:end,1:subSampling:end,1);
P = potential(V,X);
quiver(subSampledX(:),subSampledY(:),ax(:),ay(:),'k','linewidth',0.5);

