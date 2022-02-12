%% Quiz 2 - Daniel Shah
clc
clear
nodes = 9;
elements = 15;
format longg
%constants
A = 1.5;
b = 2;
e = 200e9;
a = 400e-6;
k = (e*a);

%element relations (nodes connected and locations)
element_relations = [1,2; 2,3; 3,4; 4,5; 5,6; 6,7; 7,8; 8,9; 9,1; 9,2 ;9,3; 8,3; 8,4; 8,5; 7,5];
x_coor = [2*b,3*b,3*b,2*b,b,0,b,2*b,2*b];
y_coor = [0,A,2*A,3*A,3*A,2*A,2*A,2*A,A] ;

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
    e_const = k/element_length;
    
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
force_matrix(3,1) = - 40000; %x-force at node 2
force_matrix(5,1) = - 30000; %x-force at node 3
force_matrix(8,1) = - 50000; %y-force at node 4
force_matrix(10,1) = - 40000; %y-force at node 5
K_modified = K_matrix;

fixed_nodes = [1,2,12];
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

