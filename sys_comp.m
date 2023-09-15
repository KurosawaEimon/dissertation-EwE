%sys_comp solves EwE biomass equations for 8 consumers and 1 producer
%Generates 100 runs with number of interactions removed in range 0:28:2
%Finds and plots average num and proportion of extinct species
%Finds and plots average prop of recovered/unrecovered systems
%Finds and plots average time taken to recover

global P C Mo g_i v a r_i h_i fishing fished

%PARAMS
C = 8;                 %number of consumers
P = 1;                 %num of primary producers
Mo = 0.2*ones(P+C, 1); %Mortality other than consumption
g_i = 0.8*ones(C, 1);  %growth coeff of consumers
v = 0.3*ones(P+C, 1);  %vulnerability coeff
r_i = 1.5*ones(P, 1);  %
h_i = ones(P, 1);      %
tspan = [0 500];       %time span
num_sims = 250;        %number of simulations to run

%store number of extinct species
extinct_data = zeros(15,num_sims);

%store proportion of recovered/unrecovered systems
avg_recovered = zeros(15,1);
avg_unrecovered = zeros(15,1);

%store average time for recovery
avg_rectime = zeros(15,1);

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

%set increased fishing
fishing = 1;
for n = 1:15
    
    %store number of recovered/unrecovered per iteration
    num_recovered = 0;
    num_unrecovered = 0;

    %store recovery times
    rectimes = zeros(num_sims,1);

    %number of a_ij to remove
    num_rem = 2*(n-1);

    for m = 1:num_sims %100 samples, take average
        %reset fished species each iteration
        fished = P+1;
        a = 0.2+0.8*rand(P+C, P+C);     %catchability coeff
        %for n=2:15 randomly remove 2*(n-1) a_ij coeffs
        if n > 1
            int_rem = randi([1,(P+C)^2], num_rem,1);
            a(int_rem) = 0;
        end

        %generate random initial conditions
        x0 = rand(P+C,1);

        %Solve system
        [t, x] = ode45(@EwE_biomass, tspan, x0);

        %Store final biomass levels
        biomass_fin2 = x(end,:);

        %Find number of extinct species
        num_extinct2 = find(biomass_fin2 < 0.001);

        %Store number of extinct species
        extinct_data(n,m) = length(num_extinct2);

        %Declare time of system's steady state
        t_steady = 148;
        %Store index of time nearest to 148
        [~,steady_ind] = min( abs( t-t_steady));

        %Biomass levels at steady state
        steady_state = x(steady_ind,:);

        %Change in biomass from steady state to end
        bmass_dif = abs(steady_state - biomass_fin2);

        %error margin for biomass levels
        b_error = 0.001;

        %Declare time fishing ends
        t_fStop = 200;
        %Find value and index closest to 200
        [fStop_val, fStop_ind] = min( abs(t-t_fStop));

        %Biomass matrix after fishing ends
        rec_period = x(fStop_ind:end,:);

        %If final biomass levels return to steady state (within error)
        if all( bmass_dif <= b_error)
            %Count number of recovered systems
            num_recovered = num_recovered + 1;
            
            %Create repeated matrix of steady state vector
            sState_matrix = repmat(steady_state, length(rec_period), 1);
            
            %Find difference in biomass levels from steady state
            rPeriod_Bmassdif = sum( abs( rec_period - sState_matrix), 2);

            %Find index of when system recovers
            rec_ind = find(rPeriod_Bmassdif < 0.05, 1);

            %Compute time taken from fishing end to recover
            rectimes(m) = t(rec_ind+fStop_ind) - t(fStop_ind);

        %if final biomass levels did not return to steady state
        elseif any( bmass_dif(:) > b_error)
            %Count number of unrecovered systems
            num_unrecovered = num_unrecovered + 1;
        end
    end

    %proportion of systems that did not return to steady state
    avg_unrecovered(n) = num_unrecovered/num_sims;
    %proportion of systems that returned to steady state
    avg_recovered(n) = num_recovered/num_sims;
    
    %average time for systems to recover
    avg_rectime(n) = sum(rectimes)./num_recovered;
end

%generate vector for num of removed a_ij coeffs to plot
a_rem = (0:2:28)';

%Compute average num of extinct species
fish_extinct = sum(extinct_data, 2)./num_sims;
%Compute average proportion of extinct species
fPropExtinct = fish_extinct./C;

clf;

%plot number of extinct species
figure(1)
set(gcf,'color','w');
plot(a_rem, fish_extinct, "b-*")
grid on
title("Average number of extinct species in system")
xlabel("Number of species interactions removed")
ylabel("Average number extinct")
xticks(0:2:42)

%plot proportion of extinct species
figure(2)
set(gcf,'color','w');
plot(a_rem, fPropExtinct, "b-*")
grid on
title("Average proportion of extinct species in system")
xlabel("Number of species interactions removed")
ylabel("Average proportion extinct")
xticks(0:2:42)

%plot proportion of unrecovered and recovered systems under fishing
figure(3)
set(gcf,'color','w');
plot(a_rem, avg_unrecovered, "b-*", a_rem, avg_recovered, "r--*")
grid on
title("Average proportion of systems recovered and unrecovered")
xlabel("Number of species interactions removed")
ylabel("Average unrecovered/recovered")
xticks(0:2:42)
legend("Unrecovered", "Recovered")

%plot average time to recover under fishing
figure(4)
set(gcf,'color','w');
plot(a_rem, avg_rectime, "b-*")
grid on
title("Average time for system to return to steady state")
xlabel("Number of species interactions removed")
ylabel("Average time taken for recovery")
xticks(0:2:42)