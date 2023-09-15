%input is biomass levels of each group, return P+C x 1 vector of biomass levels

function dBdt = EwE_biomass(t,x)

global P C Mo g_i v a r_i h_i fishing fished
%P is primary producer
%C is consumers
%x is biomass levels

%create vector to store biomass levels 
dBdt = zeros(P+C, 1);

%bug fix; check for negative biomass
x(x < 0) = 0;

%assign biomass levels variable
Bmass = x;

%check species is not extinct before fishing;
%considered extinct if biomass is < 0.005
if (t > 148 && t < 150)
    while Bmass(fished) < 0.005
        %Assign new species to be fished until biomass > 0.01
        fished = fished + 1;
        %bug fix; stop index exceeding max when all species are extinct
        if fished > P+C
            fished = fished - 1;
            break
        end
    end
end

%Declare fishing mortality
F_i = 0.1*ones(P+C, 1);

%increase fishing mortality by 5x for 150<t<200
if fishing == 1
    if (t < 200 && t > 150)
        F_i(fished) = F_i(fished)*5;
    end
end


%consumption matrix; element i,j is amount species j consumes species i
consumption = zeros(P+C, P+C);
for i = 1:P+C      %index rows
    for j = 1:P+C  %index columns
        if j <= P  %producers don't consume fish
            consumption(i,j) = 0;
        else
            consumption(i,j) = (a(i,j)*v(i)*Bmass(i)*Bmass(j))/(2*v(i)+a(i,j)*Bmass(j));
        end
    end
end

%production of primary producers
production = zeros(P,1);
for k = 1:P
    production(k) = (r_i(k)*Bmass(k))/(1+Bmass(k)*h_i(k));
end

%Biomass lost by each group
Bio_loss = Mo.*Bmass + F_i.*Bmass + sum(consumption, 2); %sum func sums each row of matrix

dBdt(1:P) = production - Bio_loss(1:P);
dBdt(P+1:P+C) = g_i.*sum(consumption(:,P+1:P+C))' - Bio_loss(P+1:P+C);

end