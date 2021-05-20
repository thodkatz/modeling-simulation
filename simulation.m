clc;
clear;
close all;

% the parameters I want to estimate
a = 2; b = 1;

% input
u = @(t) 5*sin(3*t);

% output
start = 0; step = 0.1; 
t = start:step:1000;
dx = @(t,x) -a*x + b*u(t);
[t, x] = ode45(dx,t,0);

% convert vector to function of time
x = @(t) x(round((t - start)/step + 1));

% plot
figure,hold on
plot(t,x(t)),plot(t,u(t)),legend('output','input')

% learning rate
gamma = 0.1; 

% the root of the stable filter
eig = -1; 

% find input of linearized model
dPhi1 = @(t,phi) eig*phi + x(t);
dPhi2 = @(t,phi) eig*phi + u(t);
[t, phi1] = ode45(dPhi1,t,0);
[t, phi2] = ode45(dPhi2,t,0);

% convert vectors to function of time
phi1 = @(t) phi1(round((t - start)/step + 1));
phi2 = @(t) phi2(round((t - start)/step + 1));

[t, thetaHat] = ode45(@(t,theta)diffSystem(t,theta,phi1,phi2,x,gamma),t,[0 0]);
aHat = -eig -thetaHat(:,1);
bHat = thetaHat(:,2);

figure
hold on
plot(t,aHat),plot(t,bHat),legend('a','b')

function dtheta = diffSystem(t,theta,phi1,phi2,x,gamma)

% thetaDot = gamma*error*phi
error = x(t) - theta(1)*phi1(t) - theta(2)*phi2(t);
dtheta(1) = gamma*error*phi1(t);
dtheta(2) = gamma*error*phi2(t);

dtheta = dtheta';

end








