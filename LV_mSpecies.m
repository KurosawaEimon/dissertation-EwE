
%M - the number of prey species and predator species
% x(1:M) - prey populations, x(M+1:2*M) - predator populations
global n a b c d
%a -  M x 1 array of a_i; a_i - prey growth rates
%b - M x M matrix of the effect of predators on prey
%c - M x M matrix of the effect of prey on predators
%d - M x 1 array of c_i; c_i - predator death rates

%declaring parameters
n=4;                   %number of species
a=0.1+0.2*rand(n,1);   %vector of random numbers uniformly distributed in [0.1,0.3]
b=0.15*rand(n,n);  %M x M matrix of rand. numbers from [0,0.15]
c=0.25*rand(n,n);   %M x M matrix of rand. numbers from [0,0.25
d=0.1+0.2*rand(n,1);   %vector of random numbers uniformly distributed in [0.1,0.3]

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
tspan=[0 500]; %time interval for solving the system
x0=0.5+0.5*rand(2*n,1); %ICs - M x 1 array of random numbers from [0.5,1.0]
[t,x] = ode45(@LV_2n_fun, tspan, x0);

%%%%%%%%%%%%%%%% Plotting the solution %%%%%%%%%%%%%%
figure(1); clf;
set(gcf,'color','w');
subplot(2,1,1)
for k=1:n
    plot(t,x(:,k))
    hold on;
end
grid;
title('(a) Dynamics of prey species');
xlabel('time');
ylabel("population size")
hold off;

subplot(2,1,2)
set(gcf,'color','w');
for k=1:n
    plot(t,x(:,n+k))
    hold on;
end
grid;
title('(b) Dynamics of predator species');
xlabel('time');
ylabel('population size')
hold off;

function dxdt = LV_2n_fun(t,x)
%M - the number of prey species and predator species 
% x(1:M) - prey populations; x(M+1:2*M) - predator populations

global n a b c d
%alpha - M x M matrix of the effect of predators on prey
%beta - M x M matrix of the effect of prey on predators
%a -  M x 1 array of a_i; a_i - prey growth rates
%b - M x 1 array of c_i; c_i - predator death rates

dxdt = zeros(2*n,1); 
dxdt(1:n) = a.*x(1:n)-diag(x(1:n))*b*x(n+1:2*n);
dxdt(n+1:2*n) = -d.*x(n+1:2*n)+diag(x(n+1:2*n))*c*x(1:n);

end