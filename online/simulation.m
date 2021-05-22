%% SISO

clc;
clear;
close all;

% the parameters I want to estimate
a = 2; b = 1;

% input
u = @(t) 5*sin(3*t);

% output
start = 0; step = 0.001; finish = 10; 
t = start:step:finish;
dx = @(t,x) -a*x + b*u(t);
[t, x] = ode45(dx,t,0);

% convert vector to function of time
x = @(t) x(round((t - start)/step + 1));

% PLOT
figure
hold on
tPlot = 0:0.01:min(10,finish);
plot(tPlot,x(tPlot))
plot(tPlot,u(tPlot))
legend('output','input')
print('in_out','-dpng','-r300')

% GRADIENT DESCENT

% the root of the stable filter
eig = -2;

% learning rate
gamma = 1;

[phi,aHat,bHat] = gradientDescent(x,u,t,eig,gamma);
phi1 = phi(:,1); phi2 = phi(:,2);

% convert vectors to function of time
theta1 = @(t) -eig -aHat(round((t - start)/step + 1));
theta2 = @(t) bHat(round((t - start)/step + 1));
phi1 = @(t) phi1(round((t - start)/step + 1));
phi2 = @(t) phi2(round((t - start)/step + 1));

% estimate output
xHat = @(t) theta1(t).*phi1(t) + theta2(t).*phi2(t);

figure
plot(t,xHat(t),t,x(t)),legend('estimation','real')

% plot
figure
subplot(211),plot(t,aHat)
subplot(212),plot(t,bHat)
sgtitle('Gradient descent')

% LYAPUNOV PARALLEL
[xHat,aHat,bHat] = lyapunovParallel(x,u,t);

% plot
figure
subplot(211)
hold on
plot(t,aHat),plot(t,bHat),legend('a','b'),title('Lyapunov Parallel')

% LYAPUNOV MIX
[xHat,aHat,bHat] = lyapunovMix(x,u,t);

% plot
subplot(212)
hold on
plot(t,aHat),plot(t,bHat),legend('a','b'),title('Lyapunov Mix')

% NOISE
amplitude = 1;
noise = @(t) 0.15*sin(2*pi*20*t);
x = @(t) x(t) + noise(t);

% LYAPUNOV PARALLEL
[xHat,aHat,bHat] = lyapunovParallel(x,u,t);

% plot
figure
subplot(211)
hold on
plot(t,aHat),plot(t,bHat),legend('a','b'),title('Lyapunov Parallel with noise')

% LYAPUNOV MIX
[xHat,aHat,bHat] = lyapunovMix(x,u,t);

% plot
subplot(212)
hold on
plot(t,aHat),plot(t,bHat),legend('a','b'),title('Lyapunov Mix with noise')

%% MIMO

clc;
clear;
close all;

% the parameters I want to estimate
a = [-0.25 3; -5 -1];
b = [1;2];

% input
u = @(t) 10*sin(2*t) + 5*sin(7.5*t);

% time range
t = 0:0.001:100;

% learning rate
gamma = [1 1];

% compare to the previous implementations, now we compute everything 
% (output,estimations) within one ode function in lyapunovParallelMIMO()
[xHat,aHat,bHat] = lyapunovParallelMIMO(u,t,a,b,gamma);

% plot
figure
subplot(211)
hold on
plot(t,aHat),legend('a11','a12','a21','a22'),title('Lyapunov Parallel')

subplot(212)
hold on
plot(t,bHat),legend('b1','b2'),title('Lyapunov Parallel')