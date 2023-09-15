% x(1) - prey population 1, x(2) - prey population 2, 
%x(3) - predator population 1, x(4) - predator population 2 
global a b c d %declaring non-dimensional parameters k,m,c global

%a - 2x1 array of a_i; a_i - prey growth rate
%b - 2x2 matrix of the effect of predators on prey
%c - 2x2 matrix of the effect of prey on predators
%d - 2x1 array of c_i; c_i - predator death rate
%parameters
a=2;
b=0.1;
c=0.4;
d=2;

options = odeset('RelTol',1e-4,'AbsTol',1e-5);
tspan=[0 50];
x0=rand(2,1);
[t1,x] = ode45(@LV2and2fun, tspan, x0);

y0 = [5, 20];
[t2, y] = ode45(@LV2and2fun, tspan, y0);

z0 = [4, 19];
[t3, z] = ode45(@LV2and2fun, tspan, z0);

clf;
figure(1)
set(gcf,'color','w');
plot(t1, x(:,1), t1,x(:,2))
grid on
title("Periodic solution for 2-species LV system")
xlabel("time")
ylabel("Population size")
legend("prey population", "predator population", location="southoutside", Orientation="horizontal")

figure(2)
set(gcf,'color','w');
subplot(2,1,1)
plot(t2, y(:,1), t2,y(:,2))
grid on
title("(a) x_0 = (5,20)")
xlabel("time")
ylabel("Population size")
ylim([0 25])

subplot(2,1,2)
plot(t3, z(:,1), t3,z(:,2))
grid on
title("(b) x_0 = (4,19)")
xlabel("time")
ylabel("Population size")
ylim([0 25])

figure(3)
set(gcf,'color','w');
plot(x(:,1), x(:,2))
grid on
title("Periodic solution for 2-species LV system")
xlabel("prey population")
ylabel("predator population")