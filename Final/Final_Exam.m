%% Problem 2 (transient)
clear
clc
format shortg
% Steady State Part
%constants
k = 15;                     %w/m*c_
t_water = 100;                %celcius

h = 17;                      %w/m^2*c 
t_air = 25;                  %celcius

l = 18/100;                   %18cm in meters
element_length = l/3;
w = 0.2/100; % 0.2 cm in meters
h = 1.25/100; %1.25 cm in meters
area = w*h;
perimeter = 2 *(w+h);

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
heat_matrix(1,1) = t_water;

%end open to convection boundary condition
free_end = nodes; %free end is last node hence this variable
k_mod(free_end,free_end) = k_mod(free_end,free_end) + h*area;
heat_matrix(free_end,1) = heat_matrix(free_end,1) + (h*area*t_air);

%find steady state nodal temps
nodal_temp = k_mod\heat_matrix;

%plot temp along the fin
disp(nodal_temp)
x = linspace(0,l,nodes);
plot(x,nodal_temp)

% Transient Part
t1 = 25 * ones(nodes,1); %t1 = t_i , t2 = t_i+1
t2 = t1 * 0;
density = 8055;
specific_heat = 480;
b = 2/3; %beta
delta_t = 60; %seconds
t_max = 15 * 60 ; %15 minutes in seconds
rho_cp = density * specific_heat;

t = 0;
counter = 1;
while t<t_max
    t = t + delta_t;
    k_matrix = zeros(nodes,nodes);
    heat_matrix_RHS = zeros(nodes,1);
    
    for element_counter = 1:elements
        n1 = elations(element_counter,1);
        n2 = elations(element_counter,2);
        cond_stiff = ((k * area)/element_length) *[1,-1; -1,1 ]; 
        conv_stiff = ((h*perimeter*element_length)/6) *[2,1; 1,2];
        
        m = (rho_cp * area * element_length)/6 * [2 1; 1 2]; 
        k_local = (1/delta_t) * m + (b * (cond_stiff+conv_stiff));
        
        q_conv = (h*t_air*perimeter*element_length)/2 * [1;1];
        T_i_local = [t1(n1,1);t1(n2,1)];
        QQ = ((m* 1/delta_t) - ((1-b)*(cond_stiff+conv_stiff)))* T_i_local;
        
        k_matrix(n1,n1) = k_matrix(n1,n1) + k_local(1,1);
        k_matrix(n1,n2) = k_matrix(n1,n2) + k_local(1,2);
        k_matrix(n2,n1) = k_matrix(n2,n1) + k_local(2,1);
        k_matrix(n2,n2) = k_matrix(n2,n2) + k_local(2,2);
        
        heat_matrix_RHS(n1,1) = heat_matrix_RHS(n1,1) + QQ(1,1);
        heat_matrix_RHS(n2,1) = heat_matrix_RHS(n2,1) + QQ(2,1);
    end
    
    %apply boundary conditions
	k_mod = k_matrix;
    heat_mod = heat_matrix_RHS;
    
    %constant temp boundary condition at node 1
    bc = 1;
    k_mod(bc,:)= 0;
    k_mod(bc,bc)= 1;
    heat_mod(bc,1) = t_water; %fixed temp
      
    %convection boundary condition at node 4
    bc = 4;
    k_mod(bc,bc) = k_mod(bc,bc) + h*area;
    heat_mod(bc,1) = heat_mod(bc,1) + (h*area*t_air);
    t2 = k_mod\heat_mod;

    plot(x,t2,'*-');
    hold on
    t1 = t2;
    
    
    Time_plot(counter)=t;
    Tem_plot(counter) = t1(nodes-1);
    counter= counter+1;
end

disp(t1)

figure 
plot(Time_plot,Tem_plot,'o-')


%% Problem 2 (transient to steady state)
clear
clc
format shortg
% Steady State Part
%constants
k = 15;                     %w/m*c_
t_water = 100;                %celcius

h = 17;                      %w/m^2*c 
t_air = 25;                  %celcius

l = 18/100;                   %18cm in meters
element_length = l/3;
w = 0.2/100; % 0.2 cm in meters
h = 1.25/100; %1.25 cm in meters
area = w*h;
perimeter = 2 *(w+h);

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
heat_matrix(1,1) = t_water;

%end open to convection boundary condition
free_end = nodes; %free end is last node hence this variable
k_mod(free_end,free_end) = k_mod(free_end,free_end) + h*area;
heat_matrix(free_end,1) = heat_matrix(free_end,1) + (h*area*t_air);

%find steady state nodal temps
nodal_temp = k_mod\heat_matrix;

%plot temp along the fin
disp('steady state temps used for tolerance equation')
disp(nodal_temp)
x = linspace(0,l,nodes);
%plot(x,nodal_temp)

% Transient Part
t1 = 25 * ones(nodes,1); %t1 = t_i , t2 = t_i+1
t2 = t1 * 0;
density = 8055;
specific_heat = 480;
b = 2/3; %beta
delta_t = 60; %seconds
rho_cp = density * specific_heat;

error = [10;10;10;10];
tolerance = 1*10^-2;
t = 0;
counter = 1;
while any(error(:) >= tolerance) == 1 %things change here
    t = t + delta_t;
    k_matrix = zeros(nodes,nodes);
    heat_matrix_RHS = zeros(nodes,1);
    
    for element_counter = 1:elements
        n1 = elations(element_counter,1);
        n2 = elations(element_counter,2);
        cond_stiff = ((k * area)/element_length) *[1,-1; -1,1 ]; 
        conv_stiff = ((h*perimeter*element_length)/6) *[2,1; 1,2];
        
        m = (rho_cp * area * element_length)/6 * [2 1; 1 2]; 
        k_local = (1/delta_t) * m + (b * (cond_stiff+conv_stiff));
        
        q_conv = (h*t_air*perimeter*element_length)/2 * [1;1];
        T_i_local = [t1(n1,1);t1(n2,1)];
        QQ = ((m* 1/delta_t) - ((1-b)*(cond_stiff+conv_stiff)))* T_i_local;
        
        k_matrix(n1,n1) = k_matrix(n1,n1) + k_local(1,1);
        k_matrix(n1,n2) = k_matrix(n1,n2) + k_local(1,2);
        k_matrix(n2,n1) = k_matrix(n2,n1) + k_local(2,1);
        k_matrix(n2,n2) = k_matrix(n2,n2) + k_local(2,2);
        
        heat_matrix_RHS(n1,1) = heat_matrix_RHS(n1,1) + QQ(1,1);
        heat_matrix_RHS(n2,1) = heat_matrix_RHS(n2,1) + QQ(2,1);
    end
    
    %apply boundary conditions
	k_mod = k_matrix;
    heat_mod = heat_matrix_RHS;
    
    %constant temp boundary condition at node 1
    bc = 1;
    k_mod(bc,:)= 0;
    k_mod(bc,bc)= 1;
    heat_mod(bc,1) = t_water; %fixed temp
      
    %convection boundary condition at node 4
    bc = 4;
    k_mod(bc,bc) = k_mod(bc,bc) + h*area;
    heat_mod(bc,1) = heat_mod(bc,1) + (h*area*t_air);
    t2 = k_mod\heat_mod;
    error = abs((nodal_temp-t1)./nodal_temp);
    
    if rem(counter,5) == 0;
        plot(x,t2,'*-');
        hold on
    end
    
    %test = error > tolerance
    %test2 = any(error(:) > tolerance)
    
    t1 = t2;
    Time_plot(counter)=t;
    Tem_plot(1,counter) = t1(1);
    Tem_plot(2,counter) = t1(2);
    Tem_plot(3,counter) = t1(3);
    Tem_plot(4,counter) = t1(4);
    counter= counter+1;

end

disp('nodal temp distributions at earliest steady state temp')
disp(t1)
disp('Seconds taken to reach steady state within tolerance')
disp(t)
disp('Minutes taken to reach steady state within tolerance')
disp(t/60)

figure 
subplot(2,2,1)
plot(Time_plot,Tem_plot(1,:),'o-')
title('Temperature Field History for Node 1')
subplot(2,2,2)
plot(Time_plot,Tem_plot(2,:),'o-')
title('Temperature Field History for Node 2')
subplot(2,2,3)
plot(Time_plot,Tem_plot(3,:),'o-')
title('Temperature Field History for Node 3')
subplot(2,2,4)
plot(Time_plot,Tem_plot(4,:),'o-')
title('Temperature Field History for Node 4')



%% Problem 3-1 (need to transfer to PDF)
clear
clc
format shortg

syms f1(x,y,z)
f1(x,y,z) = x^3 - 10*x +3 + y -z
syms f2(x,y,z)
f2(x,y,z) = y^3 + 10*y - 2*x - 2*z - 5
syms f3(x,y,z)
f3(x,y,z) = x + y +5 - 10*z + 2*sin(z)

J = [diff(f1(x,y,z),x), diff(f1(x,y,z),y),diff(f1(x,y,z),z);
     diff(f2(x,y,z),x), diff(f2(x,y,z),y),diff(f2(x,y,z),z);
     diff(f3(x,y,z),x), diff(f3(x,y,z),y),diff(f3(x,y,z),z)];
F = [f1(x,y,z);f2(x,y,z);f3(x,y,z)];

syms r(x,y,z)
r(x,y,z) = J;

syms q(x,y)
q(x,y,z) = F;

x1 = [1;1;1]; %x1 = i, x2 = i+1
%j = r(x1(1),x1(2))
%s = q(x1(1),x1(2))

tolerance = 1*10^-4;
error = [10;10;10];

while any(error(:) >= tolerance) == 1
   j = double(r(x1(1),x1(2),x1(3)));
   s = double(q(x1(1),x1(2),x1(3)));
   x2 = x1 - (inv(j)*s);
   error = abs((x2-x1)./x2);
   disp('initial values')
   disp(x1)
   disp('function values using initial values')
   disp(s)
   disp('Jacobian Matrix')
   disp(j)
   x1 = x2;
   disp('new initial values')
   disp(x1)
   disp(' ')
   disp('below is the next iteration')
   disp(' ')
end

%% Problem 3-2 (done)

clear
clc
format shortg

syms f(x)
f(x) = exp(-x) * sin(x)^2
a = 0;
b = pi/4;

%evaluating directly
direct = int(f(x),a,b);
disp('evaluating directly:')
disp(double(direct))


%Newton-Cotes (NC) - Simpson's 1/3 rule
h = (b-a)/2;
f0 = double(f(a + 0*h)) ;
f1 = double(f(a + 1*h)) ;
f2 = double(f(a + 2*h)) ;
NC = (b-a) *(f0 + (4*f1) + f2)/6;
disp('evaluating using 3 point Newton-Cotes')
disp(NC)

%Gausian Quadrature (GQ) - 2 point (REFER TO TABLE FOR E and W)
syms r(g)
h = (b-a)/2 ;
r(g) = vpa(h * g + (b+a)/2);
dx = diff(r(g),g);
gp = 1/sqrt(3); %gauss point refer to table and reference to number of points needed
gp2 = - gp;
w = 1; %weights refer to table and refernece to number of points requested

e1 = r(gp);
e2 = r(gp2);
NC = vpa(h * ((w*f(e1)) + (w*f(e2))));

disp('evaluating using 2 point Gaussian Quadrature')
disp(double(NC))

