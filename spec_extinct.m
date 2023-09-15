%spec_extinct checks how many species go extinct according to number of species
%Generates 100 runs with number of speces in range 1:15
%Finds and plots average num and proportion of extinct species
%Finds and plots average prop of recovered/unrecovered systems
%Finds and plots average time taken to recover

global P C Mo g_i v a r_i h_i fishing fished

%PARAMS
%C   - number of consumers
%Mo  - Mortality other than consumption
%g_i - growth coeff of consumers
%v   - vulnerability coeff
%a   - consumer interaction coeffs
%fishing - set to 1 to increase fishing mort

P = 1;                %num of primary producers
r_i = 1.5*ones(P, 1); %
h_i = ones(P, 1);     %
tspan = [0 500];      %time span

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

%declare matrices to store number of extinct species
fished_data = zeros(15,100);    %with fishing increased
nonfished_data = zeros(15,100); %with base level of fishing

%declare vectors to store proportion of recovered/urecovered systems
avg_recovered = zeros(15,1);
avg_unrecovered = zeros(15,1);

%declare vector to store time taken for system to recover
avg_rectime = zeros(15,1);

%Initialise 100 simulations without fishing for 1 to 15 consumers
fishing = 0;
for n = 1:15 %1-15 consumers
    C = n;                  %number of consumers
    Mo = 0.2*ones(P+C, 1);  %Mortality other than consumption
    g_i = 0.8*ones(C, 1);       %growth coeff of consumers
    v = 0.3*ones(P+C, 1);   %vulnerability coeff

    for m = 1:100
        %reset fished species 
        fished = P+1;

        %declare random parameters
        a = 0.2+0.8*rand(P+C, P+C);     %catchability coeff matrix
        x0 = rand(P+C,1);               %initial conditions

        %solve system
        [t, x] = ode45(@EwE_biomass, tspan, x0);
        
        %store final biomass levels
        biomass_fin = x(end,:);

        %find number of extinct species
        num_extinct = find(biomass_fin < 0.001);
        nonfished_data(n,m) = length(num_extinct);
    end
end

%Initialise 100 simulations with fishing for 1 to 15 consumers
fishing = 1;
for n = 1:15 %1-15 consumers
    C = n;   %num of consumers
    Mo = 0.2*ones(P+C, 1);  %Mortality other than consumption
    g_i = 0.8*ones(C, 1);       %growth coeff of consumers
    v = 0.3*ones(P+C, 1);   %vulnerability coeff (0.3 is default EwE setting)

    %To store num of recovered/unrecovered systems
    num_recovered = 0;
    num_unrecovered = 0;

    %To store recovery times
    rectimes = zeros(100,1);

    %run 100 simulations 
    for m = 1:100
        %reset fished species
        fished = P+1;

        %declare random params
        a = 0.2+0.8*rand(P+C, P+C);  %catchability coeff
        x0 = rand(P+C,1);            %initial conditions

        %solve system
        [t, x] = ode45(@EwE_biomass, tspan, x0);

        %store final biomass levels
        biomass_fin2 = x(end,:);

        %find num of extinct species
        num_extinct2 = find(biomass_fin2 < 0.001);
        fished_data(n,m) = length(num_extinct2);

        %declare time of steady state
        t_steady = 148;
        %find closest t value and index to 148
        [steady_val,steady_ind] = min( abs( t-t_steady));

        %Find biomass levels at steady state
        steady_state = x(steady_ind,:);
        %compute change in biomass from steady state to end
        bmass_dif = abs(steady_state - biomass_fin2);

        %margin of error for biomass diffs
        b_error = 0.001;

        %declare time fishing stops
        t_fStop = 200;
        %find closest t value and ind to 200
        [fStop_val, fStop_ind] = min( abs(t-t_fStop));

        %Matrix of biomass levels after fishing stops
        rec_period = x(fStop_ind:end,:);

        %if final biomass levels return to steady state (within error)
        if all( bmass_dif <= b_error)
            %count num of recovered systems
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
            %count num of unrecovered systems
            num_unrecovered = num_unrecovered + 1;
        end
    end
    %proportion of systems that did not return to steady state
    avg_unrecovered(n) = num_unrecovered/100;
    %proportion of systems that returned to steady state
    avg_recovered(n) = num_recovered/100;
    
    %average time for systems to recover
    avg_rectime(n) = sum(rectimes)./num_recovered;
end

%declare vector for number of species to plot
species_num = (1:15)';

%Find average num and proportion of extinct species; base-level fishing
nonfished_extinct = sum(nonfished_data, 2)./100;
nfPropExtinct = nonfished_extinct./species_num;

%Find average num and prop of extinct species in fished systems
fished_extinct = sum(fished_data, 2)./100;
fPropExtinct = fished_extinct./species_num;

clf;

%plot number of extinct species
figure(1)
set(gcf,'color','w');
plot(species_num, nonfished_extinct, "b-*", species_num, fished_extinct, "r--*")
grid on
title("Average number of extinct species per total species in system")
xlabel("Number of species")
ylabel("Average number extinct")
legend("Non-fishing extinction","Fishing extinction", location = "northwest")

%plot proportion of extinct species
figure(2)
set(gcf,'color','w');
plot(species_num, nfPropExtinct, "b-*", species_num, fPropExtinct, "r--*")
grid on
title("Average proportion of extinct species per total species in system")
xlabel("Number of species")
ylabel("Average proportion extinct")
legend("Non-fishing extinction","Fishing extinction", location = "northwest")

%plot proportion of unrecovered and recovered systems under fishing
figure(3)
set(gcf,'color','w');
plot(species_num, avg_unrecovered, "b-*", species_num, avg_recovered, "r--*")
grid on
title("Average proportion of fished systems recovered and unrecovered")
xlabel("Number of species")
ylabel("Average unrecovered/recovered")
legend("Unrecovered", "Recovered", location="north")

%plot average time to recover under fishing
figure(4)
set(gcf,'color','w');
plot(species_num, avg_rectime, "b-*")
grid on
title("Average time for fished system to return to steady state")
xlabel("Number of species")
ylabel("Average time taken for recovery")