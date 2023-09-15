%EwE_MATsim
%integrate ECOPATH with ECOSIM (EwE) biomass differential equations through time for P primary
%producers and C consumers.
%plots biomass levels of producers in dotted lines, consumers in solid
%lines.
%Implement increased fishing mortality on one species by setting 'fishing'
%to 1


global P C Mo g_i v a r_i h_i fishing fished

P = 1;                  %num of primary producers
C = 15;                  %num of consumers
Mo = 0.2*ones(P+C, 1);  %Mortality other than consumption
g_i = 0.8*ones(C, 1);       %growth coeff of consumers
v = 0.3*ones(P+C, 1);   %vulnerability coeff (0.3 is default EwE setting)
a = 0.2+0.8*rand(P+C, P+C);     %catchability coeff
r_i = 1.5*ones(P, 1);
h_i = ones(P, 1);
fished = P+1;

%to increase fishing mortality by 5x on group P+1 set to 1
fishing = 1;

%declare ode params
tspan = [0 400];        %time span
x0 = rand(P+C,1);       %initial biomass levels

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

%solve system
[t, x] = ode45(@EwE_biomass, tspan, x0);

figure(1); clf; set(gcf,'color','w');
hold on

%plot producers
for p=1:P
    plot(t, x(:,p), '--')
end

%plot consumers
for c=P+1:P+C
    plot(t,x(:,c))
end
hold off
grid on;
title('Dynamics of species');
xlabel('Time');
ylabel('Biomass levels');
legend('Primary Producers');
hold off;