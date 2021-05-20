function [aHat,bHat] = gradientDescent(x,u,t)
% the root of the stable filter
eig = -1; 

% find input of linearized model
dPhi1 = @(t,phi) eig*phi + x(t);
dPhi2 = @(t,phi) eig*phi + u(t);
[t, phi1] = ode45(dPhi1,t,0);
[t, phi2] = ode45(dPhi2,t,0);

% convert vectors to function of time
start = t(1); step = t(2) - t(1);
phi1 = @(t) phi1(round((t - start)/step + 1));
phi2 = @(t) phi2(round((t - start)/step + 1));

% learning rate
gamma = 1;

[t, thetaHat] = ode45(@(t,theta)diffSystem(t,theta,phi1,phi2,x,gamma),t,[0 0]);
aHat = -eig -thetaHat(:,1);
bHat = thetaHat(:,2);
end

function dtheta = diffSystem(t,theta,phi1,phi2,x,gamma)
% thetaDot = gamma*error*phi
error = x(t) - theta(1)*phi1(t) - theta(2)*phi2(t);
dtheta(1) = gamma*error*phi1(t);
dtheta(2) = gamma*error*phi2(t);
dtheta = dtheta';
end