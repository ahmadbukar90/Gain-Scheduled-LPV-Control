% Non linear system simulation
tspan=[0 5];
x0=[0 0];
U=[0.3,1.8,-1.8,20];
[t,x] = ode45(@(t,x) odefcn(t,x,U), tspan, x0);
plot(t,x(:,1),'-.',t,x(:,2),'-.');