%% Problem 1 (done)

clear
clc
format shortg

%constants list
k = 400;                     %w/m*c_
t_plate = 90;                %celcius

h = 10;                      %w/m^2*c 
t_air = 20;                  %celcius

l = 0.06;                   %6cm in meters
element_length = l/3;
dia_fin = 0.002;            %0.2cm in meters
area = pi*dia_fin^2/4;
perimeter = pi * dia_fin;

%fea setup
elements = 3;
nodes = elements + 1;
k_matrix = zeros(nodes,nodes);
heat_matrix = zeros(nodes,1);

%2d conduction convection
c1 = (k * area)/element_length; %1d conduction coefficient
c2 = (h*perimeter*element_length)/6; %element convenction coefficient
c3 = (h*t_air*perimeter*element_length)/2; %heat coefficient

cond_stiff = c1 *[1,-1; -1,1 ];
conv_stiff = c2 *[2,1; 1,2];
element_heat_matrix = c3 *[1;1];

%create element relations
for n = 1:elements
        elations (n, 1) = n;
        elations (n, 2) = n+1;
end

%global matrix creation
for element_counter = 1:elements
    n1 = elations(element_counter,1);
    n2 = elations(element_counter,2);
    k_local = cond_stiff + conv_stiff;
    
    k_matrix(n1,n1) = k_matrix(n1,n1) + k_local(1,1);
    k_matrix(n1,n2) = k_matrix(n1,n2) + k_local(1,2);
    k_matrix(n2,n1) = k_matrix(n2,n1) + k_local(2,1);
    k_matrix(n2,n2) = k_matrix(n2,n2) + k_local(2,2);
    
    heat_matrix (n1,1) = heat_matrix(n1,1) + element_heat_matrix(1,1);
    heat_matrix (n2,1) = heat_matrix(n2,1) + element_heat_matrix(2,1);
end

%apply bc's and solve for remaining temps
k_mod = k_matrix;

%constant temp boundary condition at node 1
k_mod(1,:)=0;
k_mod(1,1)= 1;
heat_matrix(1,1) = t_plate;

%end open to convection boundary condition

free_end = nodes; %free end is last node hence this variable
k_mod(free_end,free_end) = k_mod(free_end,free_end) + h*area;
heat_matrix(free_end,1) = heat_matrix(free_end,1) + (h*area*t_air);

%find steady state nodal temps

nodal_temp = k_mod\heat_matrix;

%plot temp along the fin
x = linspace(0,l,nodes);
plot(x,nodal_temp)
%% Problem 2 (done)

clear
clc
format shortg

%constants
k = 60;         %w/m*c_
t_wall = 200;    %celcius

h = 15;          %w/m^2*c
t_air = 25;      %celcius

l_total = 0.2 + 0.3;
l_ins = 0.2;
l_ele = 0.3/2;
dia_fin = 0.025;
area = pi*dia_fin^2/4;
perimeter = pi * dia_fin;

%fea setup
elements = 3;
nodes = elements + 1;
k_matrix = zeros(nodes,nodes);
heat_matrix = zeros(nodes,1);

%create element relations
for n = 1:elements
        elations (n, 1) = n;
        elations (n, 2) = n+1;
end

for element_counter = 1:elements
    n1 = elations(element_counter,1);
    n2 = elations(element_counter,2);
    if element_counter == 1
        c1 = (k * area)/l_ins; %1d conduction coefficient
        cond_stiff = c1 *[1,-1; -1,1 ];
        k_local = cond_stiff;
        
        k_matrix(n1,n1) = k_matrix(n1,n1) + k_local(1,1);
        k_matrix(n1,n2) = k_matrix(n1,n2) + k_local(1,2);
        k_matrix(n2,n1) = k_matrix(n2,n1) + k_local(2,1);
        k_matrix(n2,n2) = k_matrix(n2,n2) + k_local(2,2);
    elseif element_counter >1
        c1 = (k * area)/l_ele; %1d conduction coefficient
        c2 = (h*perimeter*l_ele)/6; %element convenction coefficient
        c3 = (h*t_air*perimeter*l_ele)/2; %heat coefficient
        
        cond_stiff = c1 *[1,-1; -1,1 ];
        conv_stiff = c2 *[2,1; 1,2];
        element_heat_matrix = c3 *[1;1];

        k_local = cond_stiff + conv_stiff;
        
        k_matrix(n1,n1) = k_matrix(n1,n1) + k_local(1,1);
        k_matrix(n1,n2) = k_matrix(n1,n2) + k_local(1,2);
        k_matrix(n2,n1) = k_matrix(n2,n1) + k_local(2,1);
        k_matrix(n2,n2) = k_matrix(n2,n2) + k_local(2,2);
        
        heat_matrix (n1,1) = heat_matrix(n1,1) + element_heat_matrix(1,1);
        heat_matrix (n2,1) = heat_matrix(n2,1) + element_heat_matrix(2,1);
    end
end

%apply bc's and solve for remaining temps
k_mod = k_matrix;

%fixed temp at node 1
k_mod(1,:)=0;
k_mod(1,1)= 1;
heat_matrix(1,1) = t_wall;
%insulated at node 2 and node 4
heat_matrix(2,1) = 0;
heat_matrix(4,1) = 0;

%plot temp along the fin
nodal_temp = k_mod\heat_matrix;
x = [0, 0.2, 0.35, 0.5];

plot(x,nodal_temp)
%% Problem 3

clear
clc
format shortg

k = 43;
t_wall = 20;

h = 250;         %w/m^2*c
t_air = 1000;    %celcius

l = 3;
l_ele = l/3;
dia = 0.03; %3cm in meters
area = pi*dia^2/4;
perimeter = pi * dia;


%fea setup
elements = 3;
nodes = elements + 1;
k_matrix = zeros(nodes,nodes);
heat_matrix = zeros(nodes,1);
x = linspace(0,3,nodes);

%create element relations
for n = 1:elements
        elations (n, 1) = n;
        elations (n, 2) = n+1;
end

for element_counter = 1:elements
    n1 = elations(element_counter,1);
    n2 = elations(element_counter,2);
    
    c1 = (k * area)/l_ele; %1d conduction coefficient
    cond_stiff = c1 *[1,-1; -1,1 ];
    k_local = cond_stiff;

    k_matrix(n1,n1) = k_matrix(n1,n1) + k_local(1,1);
    k_matrix(n1,n2) = k_matrix(n1,n2) + k_local(1,2);
    k_matrix(n2,n1) = k_matrix(n2,n1) + k_local(2,1);
    k_matrix(n2,n2) = k_matrix(n2,n2) + k_local(2,2);
end

%apply bc's and solve for remaining temps
k_mod = k_matrix;

%constant temp boundary condition at node 1
k_mod(1,:)= 0;
k_mod(1,1)= 1;
heat_matrix(1,1) = t_wall;

%end open to convection boundary condition

free_end = nodes; %free end is last node hence this variable
k_mod(nodes,nodes) = k_mod(nodes,nodes) + h*area;
heat_matrix(nodes,1) = heat_matrix(nodes,1) + (h*area*t_air);

%find steady state nodal temps
nodal_temp = k_mod\heat_matrix;

%transient parts

a = 1.17*10^-5; %alpha
rho_cp = k/a;

t_i = 20 * ones(nodes,1);
t_i_1 = t_i * 0;
b = 0.5; %beta
delta_t = 1; %minutes
t_max = 1800; %30 hours in minutes

t = 0;
counter = 1;
while t<t_max
    t = t + delta_t;
    heat_matrix_RHS = zeros(nodes,1);
    k_matrix = zeros(nodes,nodes);

    for element_counter = 1:elements
        n1 = elations(element_counter,1);
        n2 = elations(element_counter,2);
        
        c1 = (k * area)/l_ele; %1d conduction coefficient
        cond_stiff = c1 *[1,-1; -1,1 ];
        
        m = (rho_cp * area * l_ele)/6 * [2 1; 1 2]; 
        k_local = (1/delta_t) * m + ((b) * cond_stiff);
        T_i_local = [t_i(n1);t_i(n2)];
        QQ = (1/delta_t*m - ((1-b)*(cond_stiff)))*T_i_local;

        k_matrix(n1,n1) = k_matrix(n1,n1) + k_local(1,1);
        k_matrix(n1,n2) = k_matrix(n1,n2) + k_local(1,2);
        k_matrix(n2,n1) = k_matrix(n2,n1) + k_local(2,1);
        k_matrix(n2,n2) = k_matrix(n2,n2) + k_local(2,2);
        
        heat_matrix_RHS(n1,1) = heat_matrix_RHS(n1,1) + QQ(1,1);
        heat_matrix_RHS(n2,1) = heat_matrix_RHS(n2,1) + QQ(2,1);


    end
    
    %apply boundary conditions
	k_mod = k_matrix;

    %constant temp boundary condition at node 1
    k_mod(1,:)= 0;
    k_mod(1,1)= 1;
    heat_matrix_RHS(1,1) = 20; %fixed temp

    %end open to convection boundary condition
    free_end = nodes; %free end is last node hence this variable
    
    k_mod(free_end,free_end) = k_mod(free_end,free_end) + h*area;
    heat_matrix_RHS(free_end,1) = heat_matrix_RHS(free_end,1) + (h*area*t_air);
    
    t_i_1 = k_mod\heat_matrix_RHS;

    plot(x,t_i_1,'*-');
    hold on
    t_i = t_i_1;
    Time_plot(counter)=t;
    Tem_plot(counter) = t_i(nodes-1);
    counter= counter+1;
end

%figure 
%plot(Time_plot,Tem_plot,'o-')
%% Problem 4

clear
clc
format shortg

k = 425;
h = 200;         %w/m^2*c
t_air = 25;    %celcius

l = 1;
l_ele = l/2;

%assuming a linear distribution, node 2 has an area of 0.2
%we take averages of the 2 values to obtain appropriate element areas
area(1,1) = (0.1 + 0.2) /2;
area(2,1) = (0.3 + 0.2) /2;

%fea setup
elements = 2;
nodes = elements + 1;
k_matrix = zeros(nodes,nodes);
heat_matrix = zeros(nodes,1);
x = linspace(0,3,nodes);

%create element relations
for n = 1:elements
        elations (n, 1) = n;
        elations (n, 2) = n+1;
end

for element_counter = 1:elements
    n1 = elations(element_counter,1);
    n2 = elations(element_counter,2);
    
    c1 = (k * area(element_counter,1))/l_ele; %1d conduction coefficient
    cond_stiff = c1 *[1,-1; -1,1 ];
    k_local = cond_stiff;

    k_matrix(n1,n1) = k_matrix(n1,n1) + k_local(1,1);
    k_matrix(n1,n2) = k_matrix(n1,n2) + k_local(1,2);
    k_matrix(n2,n1) = k_matrix(n2,n1) + k_local(2,1);
    k_matrix(n2,n2) = k_matrix(n2,n2) + k_local(2,2);
end

%apply bc's and solve for remaining temps
k_mod = k_matrix;

%constant heat boundary condition at node 1
heat_matrix(1,1) = 1000;

%end open to convection boundary condition

free_end = nodes; %free end is last node hence this variable
k_mod(free_end,free_end) = k_mod(free_end,free_end) + h*area(2,1);
heat_matrix(free_end,1) = heat_matrix(free_end,1) + (h*area(2,1)*t_air);

%find steady state nodal temps
nodal_temp = k_mod\heat_matrix;



%transient parts
t_i = 100 * ones(nodes,1);
t_i_1 = t_i * 0;
density = 10000;
specific_heat = 240;
b = 0.5; %beta
delta_t = 20; %seconds
t_max = 1200; %20 minutes in seconds

t = 0;
counter = 1;
while t<t_max
    t = t+delta_t;
    heat_matrix_RHS = zeros(nodes,1);
    k_matrix = zeros(nodes,nodes);

    for element_counter = 1:elements
        n1 = elations(element_counter,1);
        n2 = elations(element_counter,2);
        
        c1 = (k * area(element_counter,1))/l_ele; %1d conduction coefficient
        cond_stiff = c1 *[1,-1; -1,1 ];
        %consistent_capacitance_matrix
        m = (density * specific_heat * area(element_counter,1) * l_ele)/6 * [2 1; 1 2]; 
        k_local = (1/delta_t) * m + (b * cond_stiff);
        T_i_local = [t_i(n1);t_i(n2)];
        QQ = (1/delta_t*m- (1-b)*(cond_stiff))*T_i_local;

        k_matrix(n1,n1) = k_matrix(n1,n1) + k_local(1,1);
        k_matrix(n1,n2) = k_matrix(n1,n2) + k_local(1,2);
        k_matrix(n2,n1) = k_matrix(n2,n1) + k_local(2,1);
        k_matrix(n2,n2) = k_matrix(n2,n2) + k_local(2,2);
        
        heat_matrix_RHS(n1,1) = heat_matrix_RHS(n1,1) + QQ(1,1);
        heat_matrix_RHS(n2,1) = heat_matrix_RHS(n2,1) + QQ(2,1);


    end
    
    %apply boundary conditions
	k_mod = k_matrix;

    %constant heat boundary condition at node 1
    heat_matrix_RHS(1,1) = -1000;
    %insulated at node 2
    
    %end open to convection boundary condition
    free_end = nodes; %free end is last node hence this variable
    
    k_mod(free_end,free_end) = k_mod(free_end,free_end) + h*area(2,1);
    heat_matrix_RHS(free_end,1) = heat_matrix_RHS(free_end,1) + (h*area(2,1)*t_air);
    
    t_i_1 = k_mod\heat_matrix_RHS;

    plot(t_i_1,'*-');
    hold on
    t_i = t_i_1;
    Time_plot(counter)=t;
    Tem_plot(counter) = t_i(nodes-1);
    counter= counter+1;
end

figure 
plot(Time_plot,Tem_plot,'o-')