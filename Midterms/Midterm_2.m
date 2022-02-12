%% Problem 2 (thin = plane stress) (done)

clear
clc
format longg

nodes = 4;
elements = 2;
elations = [1,2,4; 2,3,4];

k_matrix = zeros(2*nodes,2*nodes);

e = 30e6; %in psi)
v = 0.3; %poisson ratio
t = 1; %thickness (in)
x_coor = [0,10,10,0] ; %coordinates in inches
y_coor = [0,0,5,5];

%assign external forces to matrix
s = -200;
f_const = (s* 10 * 1)/2;
f = -500;
force_matrix = [0; 0; 0; f; 0; f_const; 0; f_const];

for n = 1:elements
    %extract node assignments from element relations and their
    %corresponding coordinates
    i= elations(n,1);
    j= elations(n,2);
    m= elations(n,3);
    x_i = x_coor(i);
    x_j = x_coor(j);
    x_m = x_coor(m);
    y_i = y_coor(i);
    y_j = y_coor(j);
    y_m = y_coor(m);
    
    %correspond each node direction to a corresponding spot on the global
    %stiffness matrix
    u_i = 2*i-1;
    u_j = 2*j-1;
    u_m = 2*m-1;
    v_i = 2*i;
    v_j = 2*j;
    v_m = 2*m;

    %determine the area of the element
    area = 0.5* (x_i *(y_j -y_m)+ x_j *(y_m -y_i)+ x_m *(y_i -y_j));
    
    %figure out the greek letters
    beta_i  = y_j - y_m;
    beta_j  = y_m - y_i;
    beta_m  = y_i - y_j;
    gamma_i = x_m - x_j;
    gamma_j = x_i - x_m;
    gamma_m = x_j - x_i;
    %generate b matrix
    b_matrix(1,1) = beta_i; 
    b_matrix(1,3) = beta_j;
    b_matrix(1,5) = beta_m;
    b_matrix(2,2) = gamma_i;
    b_matrix(2,4) = gamma_j;
    b_matrix(2,6) = gamma_m;
    b_matrix(3,1) = gamma_i ;
    b_matrix(3,2) = beta_i  ;
    b_matrix(3,3) = gamma_j ;
    b_matrix(3,4) = beta_j  ;
    b_matrix(3,5) = gamma_m ;
    b_matrix(3,6) = beta_m  ;
    b_const = 1/(2*area);
    b_matrix = b_const .* b_matrix;
    %generate d matrix
    d_matrix = eye(3);
    d_matrix(1,2)= v;
    d_matrix(2,1) = v;
    d_matrix(3,3) = (1-v)/2;
    d_const = e/(1-(v^2));
    d_matrix = d_const .* d_matrix;
    %determine local stiffness matrix using previous ones
    k_local = t * area * transpose(b_matrix) * d_matrix * b_matrix;
    % add them to their corresponding spot on the global stiffness matrix
    k_matrix(u_i,u_i) = k_matrix(u_i,u_i)+ k_local (1,1);
    k_matrix(u_i,v_i) = k_matrix(u_i,v_i)+ k_local (1,2);
    k_matrix(u_i,u_j) = k_matrix(u_i,u_j)+ k_local (1,3);
    k_matrix(u_i,v_j) = k_matrix(u_i,v_j)+ k_local (1,4);
    k_matrix(u_i,u_m) = k_matrix(u_i,u_m)+ k_local (1,5);
    k_matrix(u_i,v_m) = k_matrix(u_i,v_m)+ k_local (1,6);

    k_matrix(v_i,u_i) = k_matrix(v_i,u_i)+ k_local (2,1);
    k_matrix(v_i,v_i) = k_matrix(v_i,v_i)+ k_local (2,2);
    k_matrix(v_i,u_j) = k_matrix(v_i,u_j)+ k_local (2,3);
    k_matrix(v_i,v_j) = k_matrix(v_i,v_j)+ k_local (2,4);
    k_matrix(v_i,u_m) = k_matrix(v_i,u_m)+ k_local (2,5);
    k_matrix(v_i,v_m) = k_matrix(v_i,v_m)+ k_local (2,6);

    k_matrix(u_j,u_i) = k_matrix(u_j,u_i)+ k_local (3,1);
    k_matrix(u_j,v_i) = k_matrix(u_j,v_i)+ k_local (3,2);
    k_matrix(u_j,u_j) = k_matrix(u_j,u_j)+ k_local (3,3);
    k_matrix(u_j,v_j) = k_matrix(u_j,v_j)+ k_local (3,4);
    k_matrix(u_j,u_m) = k_matrix(u_j,u_m)+ k_local (3,5);
    k_matrix(u_j,v_m) = k_matrix(u_j,v_m)+ k_local (3,6);

    k_matrix(v_j,u_i) = k_matrix(v_j,u_i)+ k_local (4,1);
    k_matrix(v_j,v_i) = k_matrix(v_j,v_i)+ k_local (4,2);
    k_matrix(v_j,u_j) = k_matrix(v_j,u_j)+ k_local (4,3);
    k_matrix(v_j,v_j) = k_matrix(v_j,v_j)+ k_local (4,4);
    k_matrix(v_j,u_m) = k_matrix(v_j,u_m)+ k_local (4,5);
    k_matrix(v_j,v_m) = k_matrix(v_j,v_m)+ k_local (4,6);

    k_matrix(u_m,u_i) = k_matrix(u_m,u_i)+ k_local (5,1);
    k_matrix(u_m,v_i) = k_matrix(u_m,v_i)+ k_local (5,2);
    k_matrix(u_m,u_j) = k_matrix(u_m,u_j)+ k_local (5,3);
    k_matrix(u_m,v_j) = k_matrix(u_m,v_j)+ k_local (5,4);
    k_matrix(u_m,u_m) = k_matrix(u_m,u_m)+ k_local (5,5);
    k_matrix(u_m,v_m) = k_matrix(u_m,v_m)+ k_local (5,6);

    k_matrix(v_m,u_i) = k_matrix(v_m,u_i)+ k_local (6,1);
    k_matrix(v_m,v_i) = k_matrix(v_m,v_i)+ k_local (6,2);
    k_matrix(v_m,u_j) = k_matrix(v_m,u_j)+ k_local (6,3);
    k_matrix(v_m,v_j) = k_matrix(v_m,v_j)+ k_local (6,4);
    k_matrix(v_m,u_m) = k_matrix(v_m,u_m)+ k_local (6,5);
    k_matrix(v_m,v_m) = k_matrix(v_m,v_m)+ k_local (6,6);
end

%modify matrix based on boundary conditions
k_modified = k_matrix;
fixed_nodes = [1,2,7]; %fixed x and y at node 1, fixed x at node 4
for fixed = fixed_nodes
    k_modified (fixed,:) = 0;
    k_modified(fixed,fixed) = 1;
end

%find nodal displacements
disp_matrix = k_modified\force_matrix;
%for each element
%find the stress components, planar principal stresses, associted principal
%angles, von mises and tresca effective stresses

for a = 1:elements
    i= elations(a,1);
    j= elations(a,2);
    m= elations(a,3);

    u_i = 2*i-1;
    u_j = 2*j-1;
    u_m = 2*m-1;
    v_i = 2*i;
    v_j = 2*j;
    v_m = 2*m;

    local_disp = zeros(6,1);
    local_disp(1,1) = disp_matrix(u_i,1);
    local_disp(2,1) = disp_matrix(v_i,1);
    local_disp(3,1) = disp_matrix(u_j,1);
    local_disp(4,1) = disp_matrix(v_j,1);
    local_disp(5,1) = disp_matrix(u_m,1);
    local_disp(6,1) = disp_matrix(v_m,1);

    x_i = x_coor(i);
    x_j = x_coor(j);
    x_m = x_coor(m);
    y_i = y_coor(i);
    y_j = y_coor(j);
    y_m = y_coor(m);
    area = 0.5* (x_i *(y_j -y_m)+ x_j *(y_m -y_i)+ x_m *(y_i -y_j));

    beta_i  = y_j - y_m;
    beta_j  = y_m - y_i;
    beta_m  = y_i - y_j;
    gamma_i = x_m - x_j;
    gamma_j = x_i - x_m;
    gamma_m = x_j - x_i;

    b_const = 1/(2*area);

    b_matrix(1,1) = beta_i; 
    b_matrix(1,3) = beta_j;
    b_matrix(1,5) = beta_m;
    b_matrix(2,2) = gamma_i;
    b_matrix(2,4) = gamma_j;
    b_matrix(2,6) = gamma_m;
    b_matrix(3,1) = gamma_i ;
    b_matrix(3,2) = beta_i  ;
    b_matrix(3,3) = gamma_j ;
    b_matrix(3,4) = beta_j  ;
    b_matrix(3,5) = gamma_m ;
    b_matrix(3,6) = beta_m  ;

    b_matrix = b_const .* b_matrix;

    d_matrix = eye(3);
    d_matrix(1,2)= v;
    d_matrix(2,1) = v;
    d_matrix(3,3) = (1-v)/2;
    d_const = e/(1-(v^2));
    d_matrix = d_const .* d_matrix;


    %2. find element stress components
    stress = d_matrix * b_matrix * local_disp;

    %3. find element planar principal stresses
    stress_avg = (stress(1,1) + stress(2,1))/2;
    useful = (stress(1,1) - stress(2,1))/2;
    r = sqrt((useful^2)+(stress(3,1)^2));
    principal = [(stress_avg+r) ; (stress_avg-r)];
    
    %4 find the respective principal angles
    placeholder = 0.5*atan2d(stress(3,1),useful);
    principal_angles = [placeholder; placeholder+90];
    
    %5. find von mises stress and tresca stress for the element
    s1 = stress(1,1)^2;
    s2 = stress(1,1) * stress(2,1);
    s3 = stress(2,1)^2;
    s4 = 3*stress(3,1)^2;
    vm = sqrt((s1-s2)+s3+s4);
    
    sigma_1 = principal(1,1);
    sigma_2 = principal(2,1);
    sigma_3 = 0; %plane stress condition states sigma 3 is 0

    t1 = abs(sigma_1 - sigma_2);
    t2 = abs(sigma_2 - sigma_3);
    t3 = abs(sigma_1 - sigma_3);
    tresca_criteria = [t1;t2;t3];
    
    disp('For element:')
    disp(a)
    disp('Element Stress Components')
    disp(stress) %first row = sigma_x, second row = sigma_y, third = tau_xy
    disp('Planar Principal Stresses')
    disp(principal)
    disp('Associated Principal Angles')
    disp(principal_angles)
    disp('Von Mises Effective Stress')
    disp(vm)
    disp('Tresca Effective Stress')
    disp(max(tresca_criteria))
    disp('-----------')
    
end
%% Problem 3 (thin = plane stress) (done)

clear
clc
format longg
%constants
lx = 50; %a
ly = 10; %b
x_nodes = 21; %number of nodes request from figure 3
y_nodes = 11;
v = 0.3; %poisson ratio
e = 30e6; %youngs in psi
t = 1; %thickness in inches

element_y = y_nodes-1; %number of elements in between the y min and y max
l_between = ly/element_y; %element length in the y axis

%create triangular elements and element relations
nodes = x_nodes * y_nodes;
k_matrix = zeros(2*nodes,2*nodes);
f_matrix = zeros(2*nodes,1);

q = 1000;
angle = 0;
f_mag = q * l_between * t/2; %magnitude of distributed load per element
f_mag_x = f_mag *cosd(angle);
f_mag_y = f_mag *sind(angle);

for N_y = 1:y_nodes;
    for N_x = 1:x_nodes;
        Node_number = N_x + (N_y -1)*x_nodes;
        x_coor (Node_number) = (N_x-1)*lx/(x_nodes- 1);
        y_coor (Node_number) = (N_y-1)*ly/(y_nodes- 1);
    end
end
elations = delaunay(x_coor,y_coor);
%triplot(elations,x_coor,y_coor)
[elements,verification_step] = size (elations);

%using element information with corresponding coordinates, do the local
%stiffness matrix stuff and add them to the global
for n = 1:elements
    %extract node assignments from element relations and their
    %corresponding coordinates
    i= elations(n,1);
    j= elations(n,2);
    m= elations(n,3);
    x_i = x_coor(i);
    x_j = x_coor(j);
    x_m = x_coor(m);
    y_i = y_coor(i);
    y_j = y_coor(j);
    y_m = y_coor(m);
    
    %correspond each node direction to a corresponding spot on the global
    %stiffness matrix
    u_i = 2*i-1;
    u_j = 2*j-1;
    u_m = 2*m-1;
    v_i = 2*i;
    v_j = 2*j;
    v_m = 2*m;
    
    %plot([x_i,x_j, x_m, x_i],[y_i,y_j,y_m, y_i],'bo-')
        %hold on 
        %x_c = 1/3*(x_i+x_j+x_m);
        %y_c = 1/3*(y_i+y_j+y_m);
        %text (x_c,y_c, num2str(n),'Color','r')
        %axis equal
    %axis off
    %axis equal 
    
    %determine the area of the element
    area = 0.5* (x_i *(y_j -y_m)+ x_j *(y_m -y_i)+ x_m *(y_i -y_j));
    
    %figure out the greek letters
    beta_i  = y_j - y_m;
    beta_j  = y_m - y_i;
    beta_m  = y_i - y_j;
    gamma_i = x_m - x_j;
    gamma_j = x_i - x_m;
    gamma_m = x_j - x_i;
    %generate b matrix
    b_matrix(1,1) = beta_i; 
    b_matrix(1,3) = beta_j;
    b_matrix(1,5) = beta_m;
    b_matrix(2,2) = gamma_i;
    b_matrix(2,4) = gamma_j;
    b_matrix(2,6) = gamma_m;
    b_matrix(3,1) = gamma_i ;
    b_matrix(3,2) = beta_i  ;
    b_matrix(3,3) = gamma_j ;
    b_matrix(3,4) = beta_j  ;
    b_matrix(3,5) = gamma_m ;
    b_matrix(3,6) = beta_m  ;
    b_const = 1/(2*area);
    b_matrix = b_const .* b_matrix;
    %generate d matrix
    d_matrix = [1, v, 0;
                v, 1, 0;
                0, 0, (1-v)/2];
    d_const = e/(1-(v^2));
    d_matrix = d_const .* d_matrix;
    %determine local stiffness matrix using previous ones
    k_local = t * area * transpose(b_matrix) * d_matrix * b_matrix;
    % add them to their corresponding spot on the global stiffness matrix
    k_matrix(u_i,u_i) = k_matrix(u_i,u_i)+ k_local (1,1);
    k_matrix(u_i,v_i) = k_matrix(u_i,v_i)+ k_local (1,2);
    k_matrix(u_i,u_j) = k_matrix(u_i,u_j)+ k_local (1,3);
    k_matrix(u_i,v_j) = k_matrix(u_i,v_j)+ k_local (1,4);
    k_matrix(u_i,u_m) = k_matrix(u_i,u_m)+ k_local (1,5);
    k_matrix(u_i,v_m) = k_matrix(u_i,v_m)+ k_local (1,6);

    k_matrix(v_i,u_i) = k_matrix(v_i,u_i)+ k_local (2,1);
    k_matrix(v_i,v_i) = k_matrix(v_i,v_i)+ k_local (2,2);
    k_matrix(v_i,u_j) = k_matrix(v_i,u_j)+ k_local (2,3);
    k_matrix(v_i,v_j) = k_matrix(v_i,v_j)+ k_local (2,4);
    k_matrix(v_i,u_m) = k_matrix(v_i,u_m)+ k_local (2,5);
    k_matrix(v_i,v_m) = k_matrix(v_i,v_m)+ k_local (2,6);

    k_matrix(u_j,u_i) = k_matrix(u_j,u_i)+ k_local (3,1);
    k_matrix(u_j,v_i) = k_matrix(u_j,v_i)+ k_local (3,2);
    k_matrix(u_j,u_j) = k_matrix(u_j,u_j)+ k_local (3,3);
    k_matrix(u_j,v_j) = k_matrix(u_j,v_j)+ k_local (3,4);
    k_matrix(u_j,u_m) = k_matrix(u_j,u_m)+ k_local (3,5);
    k_matrix(u_j,v_m) = k_matrix(u_j,v_m)+ k_local (3,6);

    k_matrix(v_j,u_i) = k_matrix(v_j,u_i)+ k_local (4,1);
    k_matrix(v_j,v_i) = k_matrix(v_j,v_i)+ k_local (4,2);
    k_matrix(v_j,u_j) = k_matrix(v_j,u_j)+ k_local (4,3);
    k_matrix(v_j,v_j) = k_matrix(v_j,v_j)+ k_local (4,4);
    k_matrix(v_j,u_m) = k_matrix(v_j,u_m)+ k_local (4,5);
    k_matrix(v_j,v_m) = k_matrix(v_j,v_m)+ k_local (4,6);

    k_matrix(u_m,u_i) = k_matrix(u_m,u_i)+ k_local (5,1);
    k_matrix(u_m,v_i) = k_matrix(u_m,v_i)+ k_local (5,2);
    k_matrix(u_m,u_j) = k_matrix(u_m,u_j)+ k_local (5,3);
    k_matrix(u_m,v_j) = k_matrix(u_m,v_j)+ k_local (5,4);
    k_matrix(u_m,u_m) = k_matrix(u_m,u_m)+ k_local (5,5);
    k_matrix(u_m,v_m) = k_matrix(u_m,v_m)+ k_local (5,6);

    k_matrix(v_m,u_i) = k_matrix(v_m,u_i)+ k_local (6,1);
    k_matrix(v_m,v_i) = k_matrix(v_m,v_i)+ k_local (6,2);
    k_matrix(v_m,u_j) = k_matrix(v_m,u_j)+ k_local (6,3);
    k_matrix(v_m,v_j) = k_matrix(v_m,v_j)+ k_local (6,4);
    k_matrix(v_m,u_m) = k_matrix(v_m,u_m)+ k_local (6,5);
    k_matrix(v_m,v_m) = k_matrix(v_m,v_m)+ k_local (6,6);
end

%apply boundary conditions and solve for local displacements

k_modified = k_matrix;
for disp_checker =1:elements
    for node_check =1:3 %for every element, check every node
        node = elations(disp_checker,node_check);
        x_L = x_coor(node);
        y_L = y_coor(node);
        
        if abs(x_L-0)<1e-6 %if one node is at the left side, apply the fixed boundary condition
            u_LNN = node*2-1;
            v_LNN = node*2;
            k_modified(u_LNN,:)=0;
            k_modified(u_LNN,u_LNN)=1;
            k_modified(v_LNN,:)=0;
            k_modified(v_LNN,v_LNN)=1;       
        end
    end
end

for LNN = 1:nodes %check every node for boundary condition
    x_L = x_coor(LNN);
    y_L = y_coor(LNN);
        if abs(x_L-lx)<1e-6 % if the node lies on the end of the bar apply the distributed force boundary condition
            u_LNN = LNN*2-1;
            v_LNN = LNN*2;
            f_matrix(u_LNN) = f_mag_x;
            f_matrix(v_LNN) = f_mag_y;
        end
end
disp_matrix = k_modified\f_matrix;

sigma_x_all = zeros(elements,1);
for a = 1:elements
    i= elations(a,1);
    j= elations(a,2);
    m= elations(a,3);

    u_i = 2*i-1;
    u_j = 2*j-1;
    u_m = 2*m-1;
    v_i = 2*i;
    v_j = 2*j;
    v_m = 2*m;

    local_disp = zeros(6,1);
    local_disp(1,1) = disp_matrix(u_i,1);
    local_disp(2,1) = disp_matrix(v_i,1);
    local_disp(3,1) = disp_matrix(u_j,1);
    local_disp(4,1) = disp_matrix(v_j,1);
    local_disp(5,1) = disp_matrix(u_m,1);
    local_disp(6,1) = disp_matrix(v_m,1);

    x_i = x_coor(i);
    x_j = x_coor(j);
    x_m = x_coor(m);
    y_i = y_coor(i);
    y_j = y_coor(j);
    y_m = y_coor(m);
    area = 0.5* (x_i *(y_j -y_m)+ x_j *(y_m -y_i)+ x_m *(y_i -y_j));

    beta_i  = y_j - y_m;
    beta_j  = y_m - y_i;
    beta_m  = y_i - y_j;
    gamma_i = x_m - x_j;
    gamma_j = x_i - x_m;
    gamma_m = x_j - x_i;

    b_const = 1/(2*area);

    b_matrix(1,1) = beta_i; 
    b_matrix(1,3) = beta_j;
    b_matrix(1,5) = beta_m;
    b_matrix(2,2) = gamma_i;
    b_matrix(2,4) = gamma_j;
    b_matrix(2,6) = gamma_m;
    b_matrix(3,1) = gamma_i ;
    b_matrix(3,2) = beta_i  ;
    b_matrix(3,3) = gamma_j ;
    b_matrix(3,4) = beta_j  ;
    b_matrix(3,5) = gamma_m ;
    b_matrix(3,6) = beta_m  ;

    b_matrix = b_const .* b_matrix;

    d_matrix = eye(3);
    d_matrix(1,2)= v;
    d_matrix(2,1) = v;
    d_matrix(3,3) = (1-v)/2;
    d_const = e/(1-(v^2));
    d_matrix = d_const .* d_matrix;

    %find element stress and record it
    stress = d_matrix * b_matrix * local_disp;
    sigma_x_all(a,1) = stress(1,1);
end
%find largest element stress of the matrix
[max,element] = max(abs(sigma_x_all));
disp('Absolute max normal stress magnitude in the X direction')
disp(max)
disp('It occurs in element:')
disp(element)

%to find element 70 uncomment the following in the for n=1:elements loop
%plot([x_i,x_j, x_m, x_i],[y_i,y_j,y_m, y_i],'bo-')
        %hold on 
        %x_c = 1/3*(x_i+x_j+x_m);
        %y_c = 1/3*(y_i+y_j+y_m);
        %text (x_c,y_c, num2str(n),'Color','r')
        %axis equal
    %axis off
    %axis equal
%% Problem 3 part 2 (done)

%find the largest element stress of the matrix but the force is varying in
%degree from 0 to 180 in 1 degree intervals

clear
clc
format longg
%constants
lx = 50;
ly = 10;
x_nodes = 21;
y_nodes = 11;
v = 0.3;
e = 30e6;
t = 1;

mssa = zeros(181,2); %mssa = max stress at specific angle, a matrix to record absolute max element stresses (in the x dir) for all angles
for theta = 0:1:180 %for 0 to 180 degrees with 1 degree intervals, run the loop from part a
    mssa(theta+1,1) = theta;
    element_y = y_nodes-1; %number of elements in between the y min and y max
    l_between = ly/element_y; %element length in the y axis

    %create triangular elements and element relations
    nodes = x_nodes * y_nodes;
    k_matrix = zeros(2*nodes,2*nodes);
    f_matrix = zeros(2*nodes,1);

    q = 1000;
    f_mag = q * l_between * t/2; %magnitude of distributed load per element
    f_mag_x = f_mag *cosd(theta);
    f_mag_y = f_mag *sind(theta);

    for N_y = 1:y_nodes;
        for N_x = 1:x_nodes;
            Node_number = N_x + (N_y -1)*x_nodes;
            x_coor (Node_number) = (N_x-1)*lx/(x_nodes- 1);
            y_coor (Node_number) = (N_y-1)*ly/(y_nodes- 1);
        end
    end
    elations = delaunay(x_coor,y_coor);
    %triplot(elation,x_coor,y_coor)
    [elements,verification_step] = size (elations);

    %using element information with corresponding coordinates, do the local
    %stiffness matrix stuff and add them to the global
    for n = 1:elements
        %extract node assignments from element relations and their
        %corresponding coordinates
        i= elations(n,1);
        j= elations(n,2);
        m= elations(n,3);
        x_i = x_coor(i);
        x_j = x_coor(j);
        x_m = x_coor(m);
        y_i = y_coor(i);
        y_j = y_coor(j);
        y_m = y_coor(m);

        %correspond each node direction to a corresponding spot on the global
        %stiffness matrix
        u_i = 2*i-1;
        u_j = 2*j-1;
        u_m = 2*m-1;
        v_i = 2*i;
        v_j = 2*j;
        v_m = 2*m;

        %plot([x_i,x_j, x_m, x_i],[y_i,y_j,y_m, y_i],'bo-')
            %hold on 
            %x_c = 1/3*(x_i+x_j+x_m);
            %y_c = 1/3*(y_i+y_j+y_m);
            %text (x_c,y_c, num2str(n),'Color','r')
            %axis equal
        %axis off
        %axis equal 

        %determine the area of the element
        area = 0.5* (x_i *(y_j -y_m)+ x_j *(y_m -y_i)+ x_m *(y_i -y_j));

        %figure out the greek letters
        beta_i  = y_j - y_m;
        beta_j  = y_m - y_i;
        beta_m  = y_i - y_j;
        gamma_i = x_m - x_j;
        gamma_j = x_i - x_m;
        gamma_m = x_j - x_i;
        %generate b matrix
        b_matrix(1,1) = beta_i; 
        b_matrix(1,3) = beta_j;
        b_matrix(1,5) = beta_m;
        b_matrix(2,2) = gamma_i;
        b_matrix(2,4) = gamma_j;
        b_matrix(2,6) = gamma_m;
        b_matrix(3,1) = gamma_i ;
        b_matrix(3,2) = beta_i  ;
        b_matrix(3,3) = gamma_j ;
        b_matrix(3,4) = beta_j  ;
        b_matrix(3,5) = gamma_m ;
        b_matrix(3,6) = beta_m  ;
        b_const = 1/(2*area);
        b_matrix = b_const .* b_matrix;
        %generate d matrix
        d_matrix = [1, v, 0;
                    v, 1, 0;
                    0, 0, (1-v)/2];
        d_const = e/(1-(v^2));
        d_matrix = d_const .* d_matrix;
        %determine local stiffness matrix using previous ones
        k_local = t * area * transpose(b_matrix) * d_matrix * b_matrix;
        % add them to their corresponding spot on the global stiffness matrix
        k_matrix(u_i,u_i) = k_matrix(u_i,u_i)+ k_local (1,1);
        k_matrix(u_i,v_i) = k_matrix(u_i,v_i)+ k_local (1,2);
        k_matrix(u_i,u_j) = k_matrix(u_i,u_j)+ k_local (1,3);
        k_matrix(u_i,v_j) = k_matrix(u_i,v_j)+ k_local (1,4);
        k_matrix(u_i,u_m) = k_matrix(u_i,u_m)+ k_local (1,5);
        k_matrix(u_i,v_m) = k_matrix(u_i,v_m)+ k_local (1,6);

        k_matrix(v_i,u_i) = k_matrix(v_i,u_i)+ k_local (2,1);
        k_matrix(v_i,v_i) = k_matrix(v_i,v_i)+ k_local (2,2);
        k_matrix(v_i,u_j) = k_matrix(v_i,u_j)+ k_local (2,3);
        k_matrix(v_i,v_j) = k_matrix(v_i,v_j)+ k_local (2,4);
        k_matrix(v_i,u_m) = k_matrix(v_i,u_m)+ k_local (2,5);
        k_matrix(v_i,v_m) = k_matrix(v_i,v_m)+ k_local (2,6);

        k_matrix(u_j,u_i) = k_matrix(u_j,u_i)+ k_local (3,1);
        k_matrix(u_j,v_i) = k_matrix(u_j,v_i)+ k_local (3,2);
        k_matrix(u_j,u_j) = k_matrix(u_j,u_j)+ k_local (3,3);
        k_matrix(u_j,v_j) = k_matrix(u_j,v_j)+ k_local (3,4);
        k_matrix(u_j,u_m) = k_matrix(u_j,u_m)+ k_local (3,5);
        k_matrix(u_j,v_m) = k_matrix(u_j,v_m)+ k_local (3,6);

        k_matrix(v_j,u_i) = k_matrix(v_j,u_i)+ k_local (4,1);
        k_matrix(v_j,v_i) = k_matrix(v_j,v_i)+ k_local (4,2);
        k_matrix(v_j,u_j) = k_matrix(v_j,u_j)+ k_local (4,3);
        k_matrix(v_j,v_j) = k_matrix(v_j,v_j)+ k_local (4,4);
        k_matrix(v_j,u_m) = k_matrix(v_j,u_m)+ k_local (4,5);
        k_matrix(v_j,v_m) = k_matrix(v_j,v_m)+ k_local (4,6);

        k_matrix(u_m,u_i) = k_matrix(u_m,u_i)+ k_local (5,1);
        k_matrix(u_m,v_i) = k_matrix(u_m,v_i)+ k_local (5,2);
        k_matrix(u_m,u_j) = k_matrix(u_m,u_j)+ k_local (5,3);
        k_matrix(u_m,v_j) = k_matrix(u_m,v_j)+ k_local (5,4);
        k_matrix(u_m,u_m) = k_matrix(u_m,u_m)+ k_local (5,5);
        k_matrix(u_m,v_m) = k_matrix(u_m,v_m)+ k_local (5,6);

        k_matrix(v_m,u_i) = k_matrix(v_m,u_i)+ k_local (6,1);
        k_matrix(v_m,v_i) = k_matrix(v_m,v_i)+ k_local (6,2);
        k_matrix(v_m,u_j) = k_matrix(v_m,u_j)+ k_local (6,3);
        k_matrix(v_m,v_j) = k_matrix(v_m,v_j)+ k_local (6,4);
        k_matrix(v_m,u_m) = k_matrix(v_m,u_m)+ k_local (6,5);
        k_matrix(v_m,v_m) = k_matrix(v_m,v_m)+ k_local (6,6);
    end

    %apply boundary conditions and solve for local displacements

    k_modified = k_matrix;
    for disp_checker =1:elements
        for node_check =1:3 %for every element, check every node
            node = elations(disp_checker,node_check);
            x_L = x_coor(node);
            y_L = y_coor(node);

            if abs(x_L-0)<1e-6 %if one node is at the left side, apply the fixed boundary condition
                u_LNN = node*2-1;
                v_LNN = node*2;
                k_modified(u_LNN,:)=0;
                k_modified(u_LNN,u_LNN)=1;
                k_modified(v_LNN,:)=0;
                k_modified(v_LNN,v_LNN)=1;       
            end
        end
    end

    for LNN = 1:nodes %check every node for boundary condition
        x_L = x_coor(LNN);
        y_L = y_coor(LNN);
            if abs(x_L-lx)<1e-6 % if the node lies on the end of the bar apply the distributed force boundary condition
                u_LNN = LNN*2-1;
                v_LNN = LNN*2;
                f_matrix(u_LNN) = f_mag_x;
                f_matrix(v_LNN) = f_mag_y;
            end
    end
    disp_matrix = k_modified\f_matrix;

    sigma_x_all = zeros(elements,1);
    for a = 1:elements
        i= elations(a,1);
        j= elations(a,2);
        m= elations(a,3);

        u_i = 2*i-1;
        u_j = 2*j-1;
        u_m = 2*m-1;
        v_i = 2*i;
        v_j = 2*j;
        v_m = 2*m;

        local_disp = zeros(6,1);
        local_disp(1,1) = disp_matrix(u_i,1);
        local_disp(2,1) = disp_matrix(v_i,1);
        local_disp(3,1) = disp_matrix(u_j,1);
        local_disp(4,1) = disp_matrix(v_j,1);
        local_disp(5,1) = disp_matrix(u_m,1);
        local_disp(6,1) = disp_matrix(v_m,1);

        x_i = x_coor(i);
        x_j = x_coor(j);
        x_m = x_coor(m);
        y_i = y_coor(i);
        y_j = y_coor(j);
        y_m = y_coor(m);
        area = 0.5* (x_i *(y_j -y_m)+ x_j *(y_m -y_i)+ x_m *(y_i -y_j));

        beta_i  = y_j - y_m;
        beta_j  = y_m - y_i;
        beta_m  = y_i - y_j;
        gamma_i = x_m - x_j;
        gamma_j = x_i - x_m;
        gamma_m = x_j - x_i;

        b_const = 1/(2*area);

        b_matrix(1,1) = beta_i; 
        b_matrix(1,3) = beta_j;
        b_matrix(1,5) = beta_m;
        b_matrix(2,2) = gamma_i;
        b_matrix(2,4) = gamma_j;
        b_matrix(2,6) = gamma_m;
        b_matrix(3,1) = gamma_i ;
        b_matrix(3,2) = beta_i  ;
        b_matrix(3,3) = gamma_j ;
        b_matrix(3,4) = beta_j  ;
        b_matrix(3,5) = gamma_m ;
        b_matrix(3,6) = beta_m  ;

        b_matrix = b_const .* b_matrix;

        d_matrix = eye(3);
        d_matrix(1,2)= v;
        d_matrix(2,1) = v;
        d_matrix(3,3) = (1-v)/2;
        d_const = e/(1-(v^2));
        d_matrix = d_const .* d_matrix;

        %find element stress and record it
        stress = d_matrix * b_matrix * local_disp;
        sigma_x_all(a,1) = stress(1,1);
    end
    %find largest element stress of the matrix record in angular
    %distribution
    mssa(theta+1,2) = max(abs(sigma_x_all));
end

plot(mssa(:,1),mssa(:,2))
xlabel('Angle (degrees)')
ylabel('Absolute Max Normal Stress in the x Direction (psi)')

[max,angle] = max(mssa(:,2));
disp('Absolute Max Normal Stress for All Angles (PSI)')
disp(max)
disp('Corresponding Angle (Degrees)')
disp(angle)