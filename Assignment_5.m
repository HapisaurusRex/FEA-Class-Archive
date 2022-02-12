%% Problem 1

clc
clear
nodes = 4;
elements = 3;
%constants
l = 2;
e = 120e9;
a = 500e-6;
k = (e*a);

%element relations (nodes connected and locations)
element_relations = [1,2; 1,3 ; 1,4];
x_coor = [3/5 * l, 0,0, 3/5*l];
y_coor = [4/5 *l, 4/5*l, 0, 0] ;

K_matrix = zeros(2*nodes,2*nodes);
disp_matrix = zeros(2*nodes,1);

for element_counter = 1:elements %truss assembly matrix
    %setup
    n1 = element_relations(element_counter,1);
    x1 = x_coor(n1);
    y1 = y_coor(n1);
    n2 = element_relations(element_counter,2);
    x2 = x_coor(n2);
    y2 = y_coor(n2);
    m11 = 2*n1-1;
    m12 = 2*n1;
    m21 = 2*n2-1;
    m22 = 2*n2;
    
    %determine angle and actual length
    delta_x = x2-x1;
    delta_y = y2-y1;
    element_length = sqrt((delta_y)^2+(delta_x)^2);
    element_angle = atan2(delta_y,delta_x);
    c = cos(element_angle);
    s = sin(element_angle);
    e_const = (e*a)/element_length;
    
    K_matrix(m11,m11) = K_matrix(m11,m11) + e_const * c^2;
    K_matrix(m11,m12) = K_matrix(m11,m12) + e_const * c*s;
    K_matrix(m11,m21) = K_matrix(m11,m21) - e_const * c^2;
    K_matrix(m11,m22) = K_matrix(m11,m22) - e_const * c*s;
   
    K_matrix(m12,m11) = K_matrix(m12,m11) + e_const * c*s;
    K_matrix(m12,m12) = K_matrix(m12,m12) + e_const * s^2;
    K_matrix(m12,m21) = K_matrix(m12,m21) - e_const * c*s;
    K_matrix(m12,m22) = K_matrix(m12,m22) - e_const * s^2;
    
    K_matrix(m21,m11) = K_matrix(m21,m11) - e_const * c^2;
    K_matrix(m21,m12) = K_matrix(m21,m12) - e_const * c*s;
    K_matrix(m21,m21) = K_matrix(m21,m21) + e_const * c^2;
    K_matrix(m21,m22) = K_matrix(m21,m22) + e_const * c*s;
    
    K_matrix(m22,m11) = K_matrix(m22,m11) - e_const * c*s;
    K_matrix(m22,m12) = K_matrix(m22,m12) - e_const * s^2;
    K_matrix(m22,m21) = K_matrix(m22,m21) + e_const * c*s;
    K_matrix(m22,m22) = K_matrix(m22,m22) + e_const * s^2;
    
end

%determine displacements by applying force BC and fixture BC
force_matrix = zeros(2*nodes,1);
force_matrix(1,1) = - 15000;
force_matrix(2,1) = - 20000;
K_modified = K_matrix;

fixed_nodes = [3,4,5,6,7,8];
for fixed = fixed_nodes
    K_modified (fixed,:) = 0;
    K_modified(fixed,fixed) = 1;
    
end

disp_matrix = K_modified\force_matrix;
disp('nodal displacements')
disp(disp_matrix)


for element_counter = 1:elements %truss element forces
    %setup
    n1 = element_relations(element_counter,1);
    x1 = x_coor(n1);
    y1 = y_coor(n1);
    n2 = element_relations(element_counter,2);
    x2 = x_coor(n2);
    y2 = y_coor(n2);
    m11 = 2*n1-1;
    m12 = 2*n1;
    m21 = 2*n2-1;
    m22 = 2*n2;
    
    thing = [1,-1;-1,1];
    local_disp_matrix = [disp_matrix(m11,1);disp_matrix(m12,1);disp_matrix(m21,1);disp_matrix(m22,1)];
    
    %determine angle and actual length
    delta_x = x2-x1;
    delta_y = y2-y1;
    element_length = sqrt((delta_y)^2+(delta_x)^2);
    element_angle = atan2(delta_y,delta_x);
    c = cos(element_angle);
    s = sin(element_angle);
    e_const = (e*a)/element_length;
    trig_matrix = [c,s,0,0;0,0,c,s];
    
    
    element_force = e_const * thing * trig_matrix * local_disp_matrix;
    explanation = ['Elemental forces for nodes ', num2str(n1), ', ', num2str(n2)];
    
    stress_matrix = [-1/element_length, 1/element_length];
    element_stress = e * stress_matrix * trig_matrix *local_disp_matrix;
    otherexplanation = ['Stress for element ', num2str(element_counter)];
    
    disp(explanation)
    disp(element_force)
    disp(otherexplanation)
    disp(element_stress)
    
    
end
%% problem 2

format longg
clc
clear
nodes = 4;
elements = 4;
%constants
e = 30e6;
a = [27/625,0.05,24389/625000];
k = (e*a);

%element relations (nodes connected and locations)
element_relations = [1,2; 1,3 ; 2,3; 3,4];
x_coor = [0, (3+4.2),3] *12;
y_coor = [0, 0, 4] *12 ;

K_matrix = zeros(2*nodes,2*nodes);
disp_matrix = zeros(2*nodes,1);

for element_counter = 1:elements %truss assembly matrix
    if element_counter <4
        %setup
        n1 = element_relations(element_counter,1);
        x1 = x_coor(n1);
        y1 = y_coor(n1);
        n2 = element_relations(element_counter,2);
        x2 = x_coor(n2);
        y2 = y_coor(n2);
        m11 = 2*n1-1;
        m12 = 2*n1;
        m21 = 2*n2-1;
        m22 = 2*n2;

        %determine angle and actual length
        delta_x = x2-x1;
        delta_y = y2-y1;
        element_length = sqrt((delta_y)^2+(delta_x)^2);
        element_angle = atan2(delta_y,delta_x);
        c = cos(element_angle);
        s = sin(element_angle);
        e_const = k(1,element_counter)/element_length;

        K_matrix(m11,m11) = K_matrix(m11,m11) + e_const * c^2;
        K_matrix(m11,m12) = K_matrix(m11,m12) + e_const * c*s;
        K_matrix(m11,m21) = K_matrix(m11,m21) - e_const * c^2;
        K_matrix(m11,m22) = K_matrix(m11,m22) - e_const * c*s;

        K_matrix(m12,m11) = K_matrix(m12,m11) + e_const * c*s;
        K_matrix(m12,m12) = K_matrix(m12,m12) + e_const * s^2;
        K_matrix(m12,m21) = K_matrix(m12,m21) - e_const * c*s;
        K_matrix(m12,m22) = K_matrix(m12,m22) - e_const * s^2;

        K_matrix(m21,m11) = K_matrix(m21,m11) - e_const * c^2;
        K_matrix(m21,m12) = K_matrix(m21,m12) - e_const * c*s;
        K_matrix(m21,m21) = K_matrix(m21,m21) + e_const * c^2;
        K_matrix(m21,m22) = K_matrix(m21,m22) + e_const * c*s;

        K_matrix(m22,m11) = K_matrix(m22,m11) - e_const * c*s;
        K_matrix(m22,m12) = K_matrix(m22,m12) - e_const * s^2;
        K_matrix(m22,m21) = K_matrix(m22,m21) + e_const * c*s;
        K_matrix(m22,m22) = K_matrix(m22,m22) + e_const * s^2;
    elseif element_counter == 4
        spring_matrix = 2500 *[1,-1;-1,1];
        K_matrix(5,5)= K_matrix(5,5) + spring_matrix (1,1);
        K_matrix(5,7)= K_matrix(5,7) + spring_matrix (1,2);
        K_matrix(7,5)= K_matrix(7,5) + spring_matrix (2,1);
        K_matrix(7,7)= K_matrix(7,7) + spring_matrix (2,2);
    end
    
end

%determine displacements by applying force BC and fixture BC
force_matrix = zeros(2*nodes,1);
force_matrix((2*3-1),1) = 1200; %force applied horizontally
K_modified = K_matrix;

fixed_nodes = [1,2,4,7,8];
for fixed = fixed_nodes
    K_modified (fixed,:) = 0;
    K_modified(fixed,fixed) = 1;
end

disp_matrix = K_modified\force_matrix;
disp('nodal displacements at node 3')
disp(disp_matrix(2*3-1,1))
disp(disp_matrix(2*3,1))

for element_counter = 1:elements
    %setup
    n1 = element_relations(element_counter,1);
    n2 = element_relations(element_counter,2);
    m11 = 2*n1-1;
    m12 = 2*n1;
    m21 = 2*n2-1;
    m22 = 2*n2;
   
    
    thing = [1,-1;-1,1];
    local_disp_matrix = [disp_matrix(m11,1);disp_matrix(m12,1);disp_matrix(m21,1);disp_matrix(m22,1)];
    if element_counter <4
        x1 = x_coor(n1);
        y1 = y_coor(n1);
        x2 = x_coor(n2);
        y2 = y_coor(n2);
        delta_x = x2-x1;
        delta_y = y2-y1;
        element_length = sqrt((delta_y)^2+(delta_x)^2);
        element_angle = atan2(delta_y,delta_x);
        e_const = k(1,element_counter)/element_length;
    elseif element_counter == 4
        e_const = 2500;
        element_angle = deg2rad(180);
    end
    %determine angle and actual length
    
    
    c = cos(element_angle);
    s = sin(element_angle);
    trig_matrix = [c,s,0,0;0,0,c,s];
    
    element_force = e_const * thing * trig_matrix * local_disp_matrix;
    
    explanation = ['Elemental forces for nodes ', num2str(n1), ', ', num2str(n2)];
    disp(explanation)
    disp(element_force)
    if element_counter == 4
        disp('the forces above is the magnitude of the spring force and it is in tension ^^')
    end
    
    if element_counter <4
        stress_matrix = [-1/element_length, 1/element_length];
        element_stress = e * stress_matrix * trig_matrix *local_disp_matrix;
        otherexplanation = ['Stress for element ', num2str(element_counter)];
        disp(otherexplanation)
        disp(element_stress)
    end
end
%% Problem 3

clc
clear
format longg
nodes = 5;
elements = 7;
%constants
e = 30e6;
a = 1;
k = (e*a);

%element relations (nodes connected and locations)

element_relations = [2,1; 1,3 ; 2,3; 2,4; 3,4; 3,5; 4,5];
x_coor = [2*(3*sqrt(3)),3*sqrt(3),3*sqrt(3),0,0] *12;
y_coor = [3,0,3,3,6] * 12;

K_matrix = zeros(2*nodes,2*nodes);
disp_matrix = zeros(2*nodes,1);

for element_counter = 1:elements %truss assembly matrix
    %setup
    n1 = element_relations(element_counter,1);
    x1 = x_coor(n1);
    y1 = y_coor(n1);
    n2 = element_relations(element_counter,2);
    x2 = x_coor(n2);
    y2 = y_coor(n2);
    m11 = 2*n1-1;
    m12 = 2*n1;
    m21 = 2*n2-1;
    m22 = 2*n2;
    
    %determine angle and actual length
    delta_x = x2-x1;
    delta_y = y2-y1;
    element_length = sqrt((delta_y)^2+(delta_x)^2);
    element_angle = atan2(delta_y,delta_x);
    c = cos(element_angle);
    s = sin(element_angle);
    e_const = (e*a)/element_length;
    
    K_matrix(m11,m11) = K_matrix(m11,m11) + e_const * c^2;
    K_matrix(m11,m12) = K_matrix(m11,m12) + e_const * c*s;
    K_matrix(m11,m21) = K_matrix(m11,m21) - e_const * c^2;
    K_matrix(m11,m22) = K_matrix(m11,m22) - e_const * c*s;
   
    K_matrix(m12,m11) = K_matrix(m12,m11) + e_const * c*s;
    K_matrix(m12,m12) = K_matrix(m12,m12) + e_const * s^2;
    K_matrix(m12,m21) = K_matrix(m12,m21) - e_const * c*s;
    K_matrix(m12,m22) = K_matrix(m12,m22) - e_const * s^2;
    
    K_matrix(m21,m11) = K_matrix(m21,m11) - e_const * c^2;
    K_matrix(m21,m12) = K_matrix(m21,m12) - e_const * c*s;
    K_matrix(m21,m21) = K_matrix(m21,m21) + e_const * c^2;
    K_matrix(m21,m22) = K_matrix(m21,m22) + e_const * c*s;
    
    K_matrix(m22,m11) = K_matrix(m22,m11) - e_const * c*s;
    K_matrix(m22,m12) = K_matrix(m22,m12) - e_const * s^2;
    K_matrix(m22,m21) = K_matrix(m22,m21) + e_const * c*s;
    K_matrix(m22,m22) = K_matrix(m22,m22) + e_const * s^2;
    
end

%determine displacements by applying force BC and fixture BC
force_matrix = zeros(2*nodes,1);
force_matrix(2,1) = - 2500;
K_modified = K_matrix;

fixed_nodes = [7,8,9,10];
for fixed = fixed_nodes
    K_modified (fixed,:) = 0;
    K_modified(fixed,fixed) = 1;
    
end

disp_matrix = K_modified\force_matrix;
disp('nodal displacements')
disp(disp_matrix)

for element_counter = 1:elements %truss element forces
    %setup
    n1 = element_relations(element_counter,1);
    x1 = x_coor(n1);
    y1 = y_coor(n1);
    n2 = element_relations(element_counter,2);
    x2 = x_coor(n2);
    y2 = y_coor(n2);
    m11 = 2*n1-1;
    m12 = 2*n1;
    m21 = 2*n2-1;
    m22 = 2*n2;
    
    thing = [1,-1;-1,1];
    local_disp_matrix = [disp_matrix(m11,1);disp_matrix(m12,1);disp_matrix(m21,1);disp_matrix(m22,1)];
    
    %determine angle and actual length
    delta_x = x2-x1;
    delta_y = y2-y1;
    element_length = sqrt((delta_y)^2+(delta_x)^2);
    element_angle = atan2(delta_y,delta_x);
    c = cos(element_angle);
    s = sin(element_angle);
    e_const = (e*a)/element_length;
    trig_matrix = [c,s,0,0;0,0,c,s];
    
    
    element_force = e_const * thing * trig_matrix * local_disp_matrix;
    explanation = ['Elemental forces for nodes ', num2str(n1), ', ', num2str(n2)];
    
    stress_matrix = [-1/element_length, 1/element_length];
    element_stress = e * stress_matrix * trig_matrix *local_disp_matrix;
    otherexplanation = ['Stress for element ', num2str(element_counter)];
    
    disp(explanation)
    disp(element_force)
    disp(otherexplanation)
    disp(element_stress)
    
    
end
%% Problem 4

clc
clear
format longg
nodes = 9;
elements = 15;
%constants
e = 120e9;
a = 5e-4;
k = (e*a);

%element relations (nodes connected and locations)

element_relations = [1,2; 1,3 ; 1,4; 2,4; 3,4; 3,5; 3,6;4,6; 5,6; 5,7; 5,8; 6,8; 7,8; 7,9; 8,9];
x_coor = [0,0,3,3,6,6,9,9,12];
y_coor = [5,0,5,(9*5/12),5,(6*5/12),5,(3*5/12),5];

K_matrix = zeros(2*nodes,2*nodes);
disp_matrix = zeros(2*nodes,1);

for element_counter = 1:elements %truss assembly matrix
    %setup
    n1 = element_relations(element_counter,1);
    x1 = x_coor(n1);
    y1 = y_coor(n1);
    n2 = element_relations(element_counter,2);
    x2 = x_coor(n2);
    y2 = y_coor(n2);
    m11 = 2*n1-1;
    m12 = 2*n1;
    m21 = 2*n2-1;
    m22 = 2*n2;
    %determine angle and actual length
    delta_x = x2-x1;
    delta_y = y2-y1;
    element_length = sqrt((delta_y)^2+(delta_x)^2);
    element_angle = atan2(delta_y,delta_x);
    c = cos(element_angle);
    s = sin(element_angle);
    e_const = (e*a)/element_length;
    
    K_matrix(m11,m11) = K_matrix(m11,m11) + e_const * c^2;
    K_matrix(m11,m12) = K_matrix(m11,m12) + e_const * c*s;
    K_matrix(m11,m21) = K_matrix(m11,m21) - e_const * c^2;
    K_matrix(m11,m22) = K_matrix(m11,m22) - e_const * c*s;
   
    K_matrix(m12,m11) = K_matrix(m12,m11) + e_const * c*s;
    K_matrix(m12,m12) = K_matrix(m12,m12) + e_const * s^2;
    K_matrix(m12,m21) = K_matrix(m12,m21) - e_const * c*s;
    K_matrix(m12,m22) = K_matrix(m12,m22) - e_const * s^2;
    
    K_matrix(m21,m11) = K_matrix(m21,m11) - e_const * c^2;
    K_matrix(m21,m12) = K_matrix(m21,m12) - e_const * c*s;
    K_matrix(m21,m21) = K_matrix(m21,m21) + e_const * c^2;
    K_matrix(m21,m22) = K_matrix(m21,m22) + e_const * c*s;
    
    K_matrix(m22,m11) = K_matrix(m22,m11) - e_const * c*s;
    K_matrix(m22,m12) = K_matrix(m22,m12) - e_const * s^2;
    K_matrix(m22,m21) = K_matrix(m22,m21) + e_const * c*s;
    K_matrix(m22,m22) = K_matrix(m22,m22) + e_const * s^2;
end

%determine displacements by applying force BC and fixture BC
force_matrix = zeros(2*nodes,1);
force_matrix(6,1) = -10000;
force_matrix(10,1) = -10000;
force_matrix(14,1) = -10000;
force_matrix(18,1) = -10000;
K_modified = K_matrix;

fixed_nodes = [1,2,3];
for fixed = fixed_nodes
    K_modified (fixed,:) = 0;
    K_modified(fixed,fixed) = 1;  
end

disp_matrix = K_modified\force_matrix;
disp('nodal displacements')
disp(disp_matrix)

for element_counter = 1:elements %truss element forces
    %setup
    n1 = element_relations(element_counter,1);
    x1 = x_coor(n1);
    y1 = y_coor(n1);
    n2 = element_relations(element_counter,2);
    x2 = x_coor(n2);
    y2 = y_coor(n2);
    m11 = 2*n1-1;
    m12 = 2*n1;
    m21 = 2*n2-1;
    m22 = 2*n2;
    
    thing = [1,-1;-1,1];
    local_disp_matrix = [disp_matrix(m11,1);disp_matrix(m12,1);disp_matrix(m21,1);disp_matrix(m22,1)];
    
    %determine angle and actual length
    delta_x = x2-x1;
    delta_y = y2-y1;
    element_length = sqrt((delta_y)^2+(delta_x)^2);
    element_angle = atan2(delta_y,delta_x);
    c = cos(element_angle);
    s = sin(element_angle);
    e_const = (e*a)/element_length;
    trig_matrix = [c,s,0,0;0,0,c,s];
    
    
    element_force = e_const * thing * trig_matrix * local_disp_matrix;
    explanation = ['Elemental forces for nodes ', num2str(n1), ', ', num2str(n2)];
    
    stress_matrix = [-1/element_length, 1/element_length];
    element_stress = e * stress_matrix * trig_matrix *local_disp_matrix;
    otherexplanation = ['Stress for element ', num2str(element_counter)];
    
    disp(explanation)
    disp(element_force)
    disp(otherexplanation)
    disp(element_stress)
    
    
end
%% Problem 5

clc
clear
format longg
nodes = 7;
elements = 11;
%constants
e = 30e6;
a = 1;
k = (e*a);

%element relations (nodes connected and locations)

element_relations = [1,2 ;1,3 ;2,3 ;2,4 ;3,4 ;3,5 ;5,4 ;4,6 ;5,6 ;5,7 ;7,6];
x_coor = [0,1,1,2,2,3,3] * 12;
y_coor = [0,2.5,0,5,0,7.5,0] * 12 ;

K_matrix = zeros(2*nodes,2*nodes);
disp_matrix = zeros(2*nodes,1);

for element_counter = 1:elements %truss assembly matrix
    %setup
    n1 = element_relations(element_counter,1);
    x1 = x_coor(n1);
    y1 = y_coor(n1);
    n2 = element_relations(element_counter,2);
    x2 = x_coor(n2);
    y2 = y_coor(n2);
    m11 = 2*n1-1;
    m12 = 2*n1;
    m21 = 2*n2-1;
    m22 = 2*n2;
    
    %determine angle and actual length
    delta_x = x2-x1;
    delta_y = y2-y1;
    element_length = sqrt((delta_y)^2+(delta_x)^2);
    element_angle = atan2(delta_y,delta_x);
    c = cos(element_angle);
    s = sin(element_angle);
    e_const = (e*a)/element_length;
    
    K_matrix(m11,m11) = K_matrix(m11,m11) + e_const * c^2;
    K_matrix(m11,m12) = K_matrix(m11,m12) + e_const * c*s;
    K_matrix(m11,m21) = K_matrix(m11,m21) - e_const * c^2;
    K_matrix(m11,m22) = K_matrix(m11,m22) - e_const * c*s;
   
    K_matrix(m12,m11) = K_matrix(m12,m11) + e_const * c*s;
    K_matrix(m12,m12) = K_matrix(m12,m12) + e_const * s^2;
    K_matrix(m12,m21) = K_matrix(m12,m21) - e_const * c*s;
    K_matrix(m12,m22) = K_matrix(m12,m22) - e_const * s^2;
    
    K_matrix(m21,m11) = K_matrix(m21,m11) - e_const * c^2;
    K_matrix(m21,m12) = K_matrix(m21,m12) - e_const * c*s;
    K_matrix(m21,m21) = K_matrix(m21,m21) + e_const * c^2;
    K_matrix(m21,m22) = K_matrix(m21,m22) + e_const * c*s;
    
    K_matrix(m22,m11) = K_matrix(m22,m11) - e_const * c*s;
    K_matrix(m22,m12) = K_matrix(m22,m12) - e_const * s^2;
    K_matrix(m22,m21) = K_matrix(m22,m21) + e_const * c*s;
    K_matrix(m22,m22) = K_matrix(m22,m22) + e_const * s^2;
    
end

%determine displacements by applying force BC and fixture BC
force_matrix = zeros(2*nodes,1);
force_matrix(2,1) = -20000;
force_matrix(10,1) = -20000;
force_matrix(14,1) = -10000;
K_modified = K_matrix;

fixed_nodes = [5,6,11,13];
for fixed = fixed_nodes
    K_modified (fixed,:) = 0;
    K_modified(fixed,fixed) = 1;
    
end

disp_matrix = K_modified\force_matrix;
disp('nodal displacements')
disp(disp_matrix)

for element_counter = 1:elements %truss element forces
    %setup
    n1 = element_relations(element_counter,1);
    x1 = x_coor(n1);
    y1 = y_coor(n1);
    n2 = element_relations(element_counter,2);
    x2 = x_coor(n2);
    y2 = y_coor(n2);
    m11 = 2*n1-1;
    m12 = 2*n1;
    m21 = 2*n2-1;
    m22 = 2*n2;
    
    thing = [1,-1;-1,1];
    local_disp_matrix = [disp_matrix(m11,1);disp_matrix(m12,1);disp_matrix(m21,1);disp_matrix(m22,1)];
    
    %determine angle and actual length
    delta_x = x2-x1;
    delta_y = y2-y1;
    element_length = sqrt((delta_y)^2+(delta_x)^2);
    element_angle = atan2(delta_y,delta_x);
    c = cos(element_angle);
    s = sin(element_angle);
    e_const = (e*a)/element_length;
    trig_matrix = [c,s,0,0;0,0,c,s];
    
    
    element_force = e_const * thing * trig_matrix * local_disp_matrix;
    explanation = ['Elemental forces for nodes ', num2str(n1), ', ', num2str(n2)];
    
    stress_matrix = [-1/element_length, 1/element_length];
    element_stress = e * stress_matrix * trig_matrix *local_disp_matrix;
    otherexplanation = ['Stress for element ', num2str(element_counter)];
    
    disp(explanation)
    disp(element_force)
    disp(otherexplanation)
    disp(element_stress)
    
    
end