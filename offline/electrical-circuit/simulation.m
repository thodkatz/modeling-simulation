%%% Theodoros Katzalis AEM:9282

clear;
clc;

format long;

% samples
timeEnd = 5;
timeRange = 0:0.000001:timeEnd;
N = length(timeRange);

VR = zeros(1,N);
VC = zeros(1,N);

% collect data
index = 0;
for t = timeRange
    index = index + 1;
    Vout = v(t);
    VC(index) = Vout(1);
    VR(index) = Vout(2);
end

% inputs
input1 = @(t)2*sin(t);
input2 = ones(1,N);

% filter design
polesRange = [100];
%polesRange = 1:4:200;
simulSize = length(polesRange);
thetaLinear = zeros(6,simulSize);
filterCoeff = zeros(2,simulSize);
count = 0;

% parametric analysis
for i = polesRange
    count = count + 1;
    fprintf("iteation id %d/%d\n",count, simulSize);
    eig1 = -150;
    eig2 = -150;
    filterCoeff(:,count) = [-(eig1+eig2) ; eig1*eig2];
    thetaLinear(:,count) = estimator(eig1, eig2, N, VC, input1, input2, timeRange);
end

invRC1 = thetaLinear(1,:) + filterCoeff(1,:);
invLC1 = thetaLinear(2,:) + filterCoeff(2,:);
invRC2 = thetaLinear(3,:);
invRC3 = thetaLinear(5,:);
invLC2 = thetaLinear(6,:);

fprintf('1/RC1=%f\n1/RC2=%f\n1/RC3=%f\n1/LC1=%f\n1/LC2=%f\n', invRC1(end),... 
        invRC2(end), invRC3(end), invLC1(end), invLC2(end));

% calculate init conditions
Vdiff = v(1e-10);
Vinit = v(0);
x0(1) = Vinit(1);
x0(2) = (Vdiff(1) - Vinit(1))/1e-10; % approximate derivative

% solve ode and compare the results with data
[time,x] = ode45(@(time,x)firstOrder(time,x, invRC1(end), invLC1(end)), timeRange, x0);

%PLOTS

% error VR and VC estimated
figure();
plot(timeRange, VC(:) - x(:,1));
title("Estimation Error VC");
print('VCerror', '-dpng', '-r300');

VRhat = ones(1,length(time));
VRhat = VRhat + input1(time)' - x(:,1)';
figure();
plot(timeRange, abs(VR(:)) - abs(VRhat(:))); 
title("Estimation Error VR");

% poles parametric analysis
figure();
subplot(311);
plot(polesRange, invRC1);
title("1/RC offset");
subplot(312);
plot(polesRange, invRC2, 'color', 'magenta');
subplot(313);
plot(polesRange, invRC3, 'color', 'red');
xlabel('poles $\Lambda(s)$', 'interpreter', 'latex');
print('RCoffset', '-dpng', '-r300');

figure();
subplot(211);
plot(polesRange, invLC1); 
title("1/LC offset");
subplot(212);
plot(polesRange, invLC2, 'color', 'magenta');
xlabel('poles $\Lambda(s)$', 'interpreter', 'latex');
print('LCoffset', '-dpng', '-r300');

%plot output
figure();
plot(timeRange,VC, 'color','red');
title("VC");
xlabel("Time(s)");
ylabel("Voltage");
print('vc', '-dpng', '-r300');

figure();
plot(timeRange,VR, 'color', 'magenta');
title("VR");
xlabel("Time(s)");
ylabel("Voltage");
print('vr', '-dpng', '-r300');


%%%%%%% QUESTION 2

fprintf("\n...OUTLIERS...\n");
% add three outliers to three random positions in the data
rng(0,'twister');
randomIndex = randi([1 length(VC)], 1, 3);
%outlier  = randi([1e2 1e3], 1, 3);
outlier = [100 110 130];
VCnoised = VC;
VCnoised(randomIndex) = VCnoised(randomIndex) + outlier;

thetaLinear = estimator(-100, -100, N, VCnoised, input1, input2, timeRange);

invRC1 = thetaLinear(1) + 200;
invLC1 = thetaLinear(2) + 1e4;
invRC2 = thetaLinear(3);
invRC3 = thetaLinear(5);
invLC2 = thetaLinear(6);

meanRC = (invRC1 + invRC2 + invRC3) / 3;
meanLC = (invLC1 + invLC2) / 2;

fprintf('1/RC1=%f\n1/RC2=%f\n1/RC3=%f\n1/LC1=%f\n1/LC2=%f\n', invRC1,... 
        invRC2, invRC3, invLC1, invLC2);
fprintf("mean 1/RC = %f\nmean 1/LC = %f\n", meanRC, meanLC);
    
% calculate init conditions
Vdiff = v(1e-10);
Vinit = v(0);
x0(1) = Vinit(1);
x0(2) = (Vdiff(1) - Vinit(1))/1e-10; % approximate derivative

% solve ode and compare the results with data
[time,x] = ode45(@(time,x)firstOrder(time, x, meanRC, meanLC), timeRange, x0);

% error VR and VC estimated
figure();
plot(timeRange, VC(:) - x(:,1));
title("Estimation Error VC with outliers");
print('vcestnoised', '-dpng', '-r300');
figure();
plot(timeRange, x(:,1));
title("VC estimation with outliers");

figure();
plot(timeRange,VCnoised);
title("VC noised");
xlabel("Time(s)");
ylabel("Voltage");
axis([0 5 -2 10])
print('vcnoised', '-dpng', '-r300');

function dx = firstOrder(t,x,invRC,invLC)

% diff(Vc,t,2) + diff(Vc,t,1)*1/RC + Vc/LC = diff(u1,t,1)/RC + diff(u2,t,1)+ u2/LC
% u1(t) = 2sin(t)
% u2(t) = 1

dx(1) = x(2);
dx(2) = -x(2)*invRC - x(1)*invLC + 2*cos(t)*invRC + 0 + invLC;

dx = dx';

end

function theta = estimator(eig1, eig2, N, VC, input1, input2, timeRange)

filter = [1 -(eig1+eig2) eig1*eig2];

phi = zeros(N, 6);
phi(:,1) = lsim(tf(-[1 0], filter), VC, timeRange);
phi(:,2) = lsim(tf(-[0 1], filter), VC, timeRange);
phi(:,3) = lsim(tf([1 0],filter), input1(timeRange), timeRange);
phi(:,4) = lsim(tf([0 1],filter), input1(timeRange), timeRange);
phi(:,5) = lsim(tf([1 0],filter), input2, timeRange);
phi(:,6) = lsim(tf([0 1],filter), input2, timeRange); 

theta = VC * phi / (phi' * phi);

end