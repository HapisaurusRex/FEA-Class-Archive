%% Problem 1 (Plane Stress)

x = 0;
while(x==0)
    clear
    clc
    format longg

    nodes = 3;
    elements = 1;
    elation = [1,2,3]; %element nodal relations

    cc = 2*nodes; %creation constant
    b_matrix = zeros(3,cc);

    e = 105e6; %in MPa (N/mm^2)
    v = 0.25; %poisson ratio
    t = 1; %thickness (mm)
    x_coor = [20,80,50] ; %coordinates in mm
    y_coor = [30,30,120];

    x_disp = [0.5, 0.125, 0.75];
    y_disp = [0.25, 0, 0.25];
    disp_matrix = [x_disp(1);y_disp(1);x_disp(2);y_disp(2);x_disp(3);y_disp(3)];

    %figure out greek letters
    i= elation(1,1);
    j= elation(1,2);
    m= elation(1,3);
    x_i = x_coor(i);
    x_j = x_coor(j);
    x_m = x_coor(m);
    y_i = y_coor(i);
    y_j = y_coor(j);
    y_m = y_coor(m);

    beta_i  = y_j - y_m;
    gamma_i = x_m - x_j;

    beta_j  = y_m - y_i;
    gamma_j = x_i - x_m;

    beta_m  = y_i - y_j;
    gamma_m = x_j - x_i;

    %figure out area
    area = 0.5* (x_i *(y_j -y_m)+ x_j *(y_m -y_i)+ x_m *(y_i -y_j));

    %setup B matrix for plane matrix
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

    %setup D matrix for plane stress
    d_matrix = eye(3);
    d_matrix(1,2)= v;
    d_matrix(2,1) = v;
    d_matrix(3,3) = (1-v)/2;
    d_const = e/(1-(v^2));
    d_matrix = d_const .* d_matrix;

    %1. find element stiffness matrix
    k = t * area * transpose(b_matrix) * d_matrix * b_matrix %stiffness matrix in N/mm

    %2. find element stress components
    stress = d_matrix * b_matrix * disp_matrix

    %3. find element planar principal stresses
    stress_avg = (stress(1,1) + stress(2,1))/2;
    useful = (stress(1,1) - stress(2,1))/2;
    r = sqrt((useful^2)+(stress(3,1)^2));
    principal = [(stress_avg+r) ; (stress_avg-r)]
    %4. find associapted principal angles
    placeholder = 0.5*atan2d(stress(3,1),useful);
    principal_angles = [placeholder; placeholder+90]

    %5. find von mises stress
    s1 = stress(1,1)^2;
    s2 = stress(1,1) * stress(2,1);
    s3 = stress(2,1)^2;
    s4 = 3*stress(3,1)^2;
    vm = sqrt((s1-s2)+s3+s4)
    x=1;
end
%% Problem 2 (thin = plane stress)

x = 0;
while(x==0);
    clear
    clc
    format longg

    nodes = 4;
    elements = 2;
    elation = [1,3,2; 1,4,3];

    k_matrix = zeros(2*nodes,2*nodes);

    e = 30e6; %in psi)
    v = 0.3; %poisson ratio
    t = 1; %thickness (in)
    x_coor = [0,0,20,20] ; %coordinates in inches
    y_coor = [0,10,10,0];
    s = -5000

    for n = 1:elements
        i= elation(n,1);
        j= elation(n,2);
        m= elation(n,3);

        u_i = 2*i-1;
        u_j = 2*j-1;
        u_m = 2*m-1;
        v_i = 2*i;
        v_j = 2*j;
        v_m = 2*m;


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

        k_local = t * area * transpose(b_matrix) * d_matrix * b_matrix;

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
    
    k_modified = k_matrix;
    
    fixed_nodes = [1,2,3,4];
    for fixed = fixed_nodes
        k_modified (fixed,:) = 0;
        k_modified(fixed,fixed) = 1;
    end
    
    force_matrix = zeros(2*nodes,1);
    force_constant = s * t * max(y_coor)/2;
    force_matrix(8,1) = force_constant;
    force_matrix(6,1) = force_constant;
    
    disp_matrix = k_modified\force_matrix;
    
    for n = 1:elements
        i= elation(n,1);
        j= elation(n,2);
        m= elation(n,3);

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
        %1. find element stress for both elements
        stress = d_matrix * b_matrix * local_disp
        %2. find element planar stress for each element
        stress_avg = (stress(1,1) + stress(2,1))/2;
        useful = (stress(1,1) - stress(2,1))/2;
        r = sqrt((useful^2)+(stress(3,1)^2));
        principal = [(stress_avg+r) ; (stress_avg-r)]
        %3. find associated principal angles for each element
        placeholder = 0.5*atan2d(stress(3,1),useful);
        principal_angles = [placeholder; placeholder+90]
    end  
    x=1;
end
%% Problem 3 (thin = plane stress)

x=0;
while(x==0);

    %1. plane stress condition applies as the thickness is very small compared
    %to other dimensions (so the stress in the thickness direction is close to
    %0)

    clear
    clc
    format longg

    nodes = 10;
    elements = 9;
    elation = [1,3,2; 2,5,4; 2,3,5; 3,6,5; 4,8,7; 4,5,8; 5,9,8; 5,6,9; 6,10,9];

    k_matrix = zeros(2*nodes,2*nodes);
    x = linspace(0,60,4);
    y = linspace(0,30,4);

    e = 30000e3; %in psi)
    v = 0.3; %poisson ratio
    t = 1; %thickness (in)
    x_coor = [x(1),x(1),x(2),x(1),x(2),x(3),x(1),x(2),x(3),x(4)] ; %coordinates in inches
    y_coor = [y(1),y(2),y(2),y(3),y(3),y(3),y(4),y(4),y(4),y(4)];

    for n = 1:elements
        i= elation(n,1);
        j= elation(n,2);
        m= elation(n,3);

        u_i = 2*i-1;
        u_j = 2*j-1;
        u_m = 2*m-1;
        v_i = 2*i;
        v_j = 2*j;
        v_m = 2*m;


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

        k_local = t * area * transpose(b_matrix) * d_matrix * b_matrix;

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

    force_matrix = zeros(2*nodes,1);
    force_matrix(20,1) = - 50000;
    k_modified = k_matrix;

    fixed_nodes = [1,2,3,4,7,8,13,14];
    for fixed = fixed_nodes
        k_modified (fixed,:) = 0;
        k_modified(fixed,fixed) = 1;
    end

    disp_matrix = k_modified\force_matrix;
    %2. displacements at node 10
    disp('displacements at node 10 (horziontal first)')
    disp(disp_matrix(19,1))
    disp(disp_matrix(20,1))

    %3. element 5 principal and von msies stresses

    n = 5;
    i= elation(n,1);
    j= elation(n,2);
    m= elation(n,3);

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
    principal = [(stress_avg+r) ; (stress_avg-r)]

    %5. find von mises stress
    s1 = stress(1,1)^2;
    s2 = stress(1,1) * stress(2,1);
    s3 = stress(2,1)^2;
    s4 = 3*stress(3,1)^2;
    disp('von mises stress')
    vm = sqrt((s1-s2)+s3+s4)
    x=1;
end
%% Problem 4

clear
clc

p=0;
while(p==0);
    l_x = 20;
    l_y = 10;
    x_nodes = 2;
    y_nodes = 2;
    
    xx = linspace(0,l_x,x_nodes);
    yy = linspace(0,l_y,y_nodes);
    [xx,yy] = meshgrid(xx,yy);
    [row_s,col_s] = size(xx);
    
    nodes = row_s * col_s;
 
    x_coor = zeros(1,nodes);
    y_coor = zeros(1,nodes);
    
    for row_counter = 1:row_s;
        for col_counter = 1:col_s;
            chosennode = (row_counter-1)*col_s +col_counter;
            x_coor(chosennode) = xx(row_counter,col_counter);
            y_coor(chosennode) = yy(row_counter,col_counter);
        end
    end
    
    elation = delaunay(x_coor,y_coor)
    [elements,dummy] = size(elation);
    
    e = 30e6; %in psi
    v = 0.3; %poisson ratio
    t = 1;
    s = -5000;
    k_matrix = zeros(2*nodes,2*nodes);
    
    
    for n = 1:elements
        i= elation(n,1);
        j= elation(n,2);
        m= elation(n,3);
        
        x_i = x_coor(i);
        x_j = x_coor(j);
        x_m = x_coor(m);
        y_i = y_coor(i);
        y_j = y_coor(j);
        y_m = y_coor(m);
        area = 0.5* (x_i *(y_j -y_m)+ x_j *(y_m -y_i)+ x_m *(y_i -y_j));

        u_i = 2*i-1;
        u_j = 2*j-1;
        u_m = 2*m-1;
        v_i = 2*i;
        v_j = 2*j;
        v_m = 2*m;
        
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

        k_local = t * area * transpose(b_matrix) * d_matrix * b_matrix;

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
    
   
    
    
    contour_output_x = zeros(3,nodes);
    contour_output_y = zeros(3,nodes);
    
    sigmaEL_xx = rand(nodes,1);
    %{
    for el = 1:elements
        for local_n = 1:3
            lnn = elations(el,local_n)
            x_l = x_nodes(lnn)
            y_l = y_nodes(lnn)
            contour_output_x(local_n,el) = x_l
            contour_output_y(local_n,el) = y_l
        end
    end
    %}
    
    p=1;
end
%% Problem 4B

clear
clc
elation = readmatrix('NOC.csv');
coordinates = readmatrix('NODES.csv');
x_coor = coordinates (:,1);
y_coor = coordinates (:,2);

[nodes, nothelpful] = size(x_coor);
[elements, nothelpful] = size(elation);
k_matrix = zeros(2*nodes,2*nodes);

e = 30e6; %in psi
v = 0.3; %poisson ratio
t = 1;  
s = -5000;

    for n = 1:elements
        i= elation(n,1);
        j= elation(n,2);
        m= elation(n,3);
        
        x_i = x_coor(i);
        x_j = x_coor(j);
        x_m = x_coor(m);
        y_i = y_coor(i);
        y_j = y_coor(j);
        y_m = y_coor(m);
        area = 0.5* (x_i *(y_j -y_m)+ x_j *(y_m -y_i)+ x_m *(y_i -y_j));

        u_i = 2*i-1;
        u_j = 2*j-1;
        u_m = 2*m-1;
        v_i = 2*i;
        v_j = 2*j;
        v_m = 2*m;
        
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

        k_local = t * area * transpose(b_matrix) * d_matrix * b_matrix;

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