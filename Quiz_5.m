%% Quiz 5
%organiziation stuff
clear
clc
format shortg

%constants
%thermal constants
k = 15
h = 35
t_air = 20
t_fixed = 350
q = 100 %W/m
%geometrical constants
area = 500e-6; %500 mm^2 in m^2
peri = 120e-3; %120 mm in m
l_element = 100e-3 %100mm in m

%fea setup
elements = 2;
nodes = elements + 1;
k_matrix = zeros(nodes,nodes);
heat_matrix = zeros(nodes,1);
%create element relations
for n = 1:elements
        elations (n, 1) = n;
        elations (n, 2) = n+1;
end

%global matrix creation
for element_counter = 1:elements
    n1 = elations(element_counter,1);
    n2 = elations(element_counter,2);
    c1 = (k * area)/l_element; %1d conduction coefficient
    c2 = (h*peri*l_element)/6; %element convenction coefficient
    c3 = (h*t_air*peri*l_element)/2; %heat coefficient
    cond_stiff = c1 *[1,-1; -1,1 ];
    conv_stiff = c2 *[2,1; 1,2];
    element_heat_matrix = c3 *[1;1];
    
    if element_counter == 1
        k_local = cond_stiff;
        
        k_matrix(n1,n1) = k_matrix(n1,n1) + k_local(1,1);
        k_matrix(n1,n2) = k_matrix(n1,n2) + k_local(1,2);
        k_matrix(n2,n1) = k_matrix(n2,n1) + k_local(2,1);
        k_matrix(n2,n2) = k_matrix(n2,n2) + k_local(2,2);

    elseif element_counter >1 
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

%constant temp boundary condition at node 1
bc = 1
k_mod(bc,:)=0;
k_mod(bc,bc)= 1;
heat_matrix(bc,bc) = t_fixed;

%heat loss at nodes 1 and 2

net_heat = q * l_element
heat_matrix(1,1) = heat_matrix(1,1) - (2/3 * net_heat)
heat_matrix(2,1) = heat_matrix(2,1) - (1/3 * net_heat)


%convection at node 3
bc = 3
k_mod(bc,bc) = k_mod(bc,bc) + h*area;
heat_matrix(bc,1) = heat_matrix(bc,1) + (h*area*t_air);

%find steady state nodal temps
nodal_temp = k_mod\heat_matrix;

%plot temp along the fin
x = linspace(0,200,nodes);
plot(x,nodal_temp)