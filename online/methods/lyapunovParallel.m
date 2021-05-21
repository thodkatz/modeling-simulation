function [xHat,aHat,bHat] = lyapunovParallel(x,u,t)
gamma = [1 1];
[t,y] = ode45(@(t,y)diffSystem(x,u,t,y,gamma), t, [0 0 0]);
xHat = y(:,1);
aHat = y(:,2);
bHat = y(:,3);
end

function dy = diffSystem(x,u,t,y,gamma)
xHat=y(1); theta(1) = y(2); theta(2) = y(3);

dxHat = -theta(1)*xHat + theta(2)*u(t);
error = x(t) - xHat;
dtheta(1) = -gamma(1)*error*xHat;
dtheta(2) = gamma(2)*error*u(t);

dy(1) = dxHat; dy(2) = dtheta(1); dy(3) = dtheta(2);
dy = dy';
end