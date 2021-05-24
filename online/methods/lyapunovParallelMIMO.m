function [x,xHat,aHat,bHat] = lyapunovParallelMIMO(u,t,a,b,gamma)
[t,y] = ode45(@(t,y)diffSystem(t,y,u,a,b,gamma), t, zeros(10,1));
x    = [y(:,1) y(:,4)];
xHat = [y(:,3) y(:,4)];
aHat = [y(:,5) y(:,6) y(:,7) y(:,8)];
bHat = [y(:,9) y(:,10)];
end

function dy = diffSystem(t,y,u,a,b,gamma)
x    = [y(1);y(2)];
xHat = [y(3);y(4)];
aHat = [y(5) y(6);y(7) y(8)];
bHat = [y(9);y(10)];

dx = a*x + b*u(t);
e = [x(1) - xHat(1); x(2) - xHat(2)];
dxHat = aHat*xHat + bHat*u(t);
daHat = gamma(1) * (xHat*e');
dbHat = gamma(2)*u(t)*e';

dy = [dx(1); dx(2);dxHat(1); dxHat(2);...
    daHat(1);daHat(2);daHat(3);daHat(4);...
    dbHat(1);dbHat(2)];
end