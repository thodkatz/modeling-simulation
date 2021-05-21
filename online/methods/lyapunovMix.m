function [xHat,aHat,bHat] = lyapunovMix(x,u,t)
gamma = [1 1];
thetaM = 1;
[t,y] = ode45(@(t,y)diffSystem(x,u,t,y,gamma,thetaM), t, [0 0 0]);
xHat = y(:,1);
aHat = y(:,2);
bHat = y(:,3);
end

function dy = diffSystem(x,u,t,y,gamma,thetaM)
xHat=y(1); theta(1) = y(2); theta(2) = y(3);

error = x(t) - xHat;
dxHat = -theta(1)*x(t) + theta(2)*u(t) + thetaM*error;
dtheta(1) = -gamma(1)*error*x(t);
dtheta(2) = gamma(2)*error*u(t);

dy(1) = dxHat; dy(2) = dtheta(1); dy(3) = dtheta(2);
dy = dy';
end