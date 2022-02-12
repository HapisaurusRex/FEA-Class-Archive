%% Problem 1

clc
clear
n = input('Please give me a number: '); %gives number, assume it's an integer

for i = 2:1:n
    if rem(n,i) == 0 % if the remainder is 0, it's a divisor
        prime = 0; %this will remain 0 if no divisors exist between 2 and (divisor - 1), aka a prime number
        for j = 2:1:(i-1)
            if rem(i,j) == 0
                prime = prime + 1;
            end
        end
        if prime == 0 %if there are no divisors that aren't 1 and i, it's prime = 0 and we can show the number
            disp(i)
        end
    end
end
%% Problem 2

clc
clear
format shortg
nodes = 4;
elements = 3;
%constants
e = [15e6; 30e6 ;30e6]; %young's in psi
a = (pi()/4) * [1^2; (1^2 - 0.75^2); 0.5^2]; %area in in^2
k = e.*a; %constant made to increase readabilitiy

%element relations (nodes connected and locations)
element_relations = [1,2; 2,3 ; 2,4];
x_coor = [0, 20, 40, 70];


%setup matricies to be filled
K_matrix = zeros(nodes,nodes);
disp_matrix = zeros(nodes,1);
force_matrix = zeros(nodes,1); 

for element_counter = 1:elements %compute bar element matricies and add them to the global matrix
    n1 = element_relations(element_counter,1); %figure out 1st node, based of element relations, and coordinates
    x1 = x_coor(n1);
    n2 = element_relations(element_counter,2);%figure out 2nd node, based of element relations, and coordinates
    x2 = x_coor(n2);
    
    %determine element length and correct stiffness constant per element
    element_length = abs(x2-x1); 
    e_const = k(element_counter,1)/element_length;
    
    %add element matrix elements to their corresponding global position
    K_matrix(n1,n1) = K_matrix(n1,n1) + e_const ;
    K_matrix(n1,n2) = K_matrix(n1,n2) - e_const ;
    K_matrix(n2,n1) = K_matrix(n2,n1) - e_const ;
    K_matrix(n2,n2) = K_matrix(n2,n2) + e_const ;
    
end

%determine displacements by applying force BC and support BC
%force BC
force_matrix(4,1) = 2000;
%support BCs
fixed_nodes = [1,3]; %node 1 and 3 are fixed in the x-dir and we are only looking at x-dir
K_modified = K_matrix;

%modify K_matrix based on both BCs
for fixed = fixed_nodes
    K_modified (fixed,:) = 0;
    K_modified(fixed,fixed) = 1;
end

disp_matrix = K_modified\force_matrix; %solve modified matrix for nodal displacements
disp('Node 2 displacement')
disp(disp_matrix(2,1))
disp('Node 4 displacement')
disp(disp_matrix(4,1))

% finding nodal forces and stress (part 2 and 3)
for element_counter = 1:elements %element forces
    %setup
    % we need some things from the previous for loop (node relations,
    % coordinates, constants) as well as the displacement matrix
    n1 = element_relations(element_counter,1);
    x1 = x_coor(n1);
    n2 = element_relations(element_counter,2);
    x2 = x_coor(n2);
    element_length = abs(x2-x1);
    e_const = k(element_counter,1)/element_length;

    local_disp_matrix = [disp_matrix(n1,1);disp_matrix(n2,1)]; %take values from the displacement matrix that match the nodes' x and y directions
    
    stress_matrix = [-1/element_length, 1/element_length]; %matrix to be multiplied to obtain stress
    element_stress = e(element_counter,1) * stress_matrix *local_disp_matrix; %element stress eq as taught
    otherexplanation = ['Stress for element (in psi) ', num2str(element_counter)]; %text to accompany
    
    %show the answers for each element
    disp(otherexplanation)
    disp(element_stress)
    if element_stress < 0 %if element stress is negative, then the element is under compression
        disp('The element is under compression')
    else
        disp('The element is under tension')
    end
    
    %un comment to figure out support forces
    %thing = [1,-1;-1,1];
    %element_force = e_const * thing * local_disp_matrix; % element force calculation as learned
    %explanation = ['Elemental forces for nodes ', num2str(n1), ', ', num2str(n2)];
    %disp(explanation)
    %disp(element_force)
end
%% Problem 3

clc
clear
format shortg
nodes = 4;
elements = 3;
%constants
e = 70e9;
a = 400e-6;
k = (e*a);
c = 4 * tand(60); %get unknown height so we can assign coordinates to each node

%element relations (nodes connected and locations)
element_relations = [1,2; 3,1 ; 4,1];
x_coor = [4,0,0,0];
y_coor = [c,(c+3),c,0] ;

%setup matricies to be filled
K_matrix = zeros(2*nodes,2*nodes);
disp_matrix = zeros(2*nodes,1);
force_matrix = zeros(2*nodes,1); 

for element_counter = 1:elements %compute truss element eqs and add to global matrix
    n1 = element_relations(element_counter,1); %figure out 1st node and coordinates
    x1 = x_coor(n1);
    y1 = y_coor(n1);
    n2 = element_relations(element_counter,2);%figure out 2nd node and coordinates
    x2 = x_coor(n2);
    y2 = y_coor(n2);
    
    m11 = 2*n1-1; %translate 1st node x into global matrix format
    m12 = 2*n1;   %translate 1st node y into global matrix format
    m21 = 2*n2-1; %translate 2nd node x into global matrix format
    m22 = 2*n2;   %translate 2nd node y into global matrix format
    
    %determine angle and element length
    delta_x = x2-x1;
    delta_y = y2-y1;
    element_length = sqrt((delta_y)^2+(delta_x)^2); %pythagoras 
    element_angle = atan2(delta_y,delta_x); %inverse tan in radians
    
    %some constants to keep future code sane.
    c = cos(element_angle); 
    s = sin(element_angle);
    e_const = k/element_length;
    
    %add element matrix elements to their corresponding global position
    
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

%determine displacements by applying force BC and support BC
%force BC
force_matrix(1,1) = 60000;
%support BCs
fixed_nodes = [2,3,4,5,6,7,8]; %node 1 is only fixed in the y direction and everything else is pinned (fixed in x and y)
K_modified = K_matrix;

%modify K_matrix based on both BCs
for fixed = fixed_nodes
    K_modified (fixed,:) = 0;
    K_modified(fixed,fixed) = 1;
end

disp_matrix = K_modified\force_matrix;
disp('nodal displacements for node 1 (in meters) (x first, y last)')
disp(disp_matrix(1,1))
disp(disp_matrix(2,1))

% finding nodal forces and stress (part 2 and 3)
for element_counter = 1:elements %element forces
    %setup
    % we need some things from the previous for loop (node relations,
    % coordinates, element angle and length) as well as the displacement matrix
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
    delta_x = x2-x1;
    delta_y = y2-y1;
    element_length = sqrt((delta_y)^2+(delta_x)^2);
    element_angle = atan2(delta_y,delta_x);
    c = cos(element_angle);
    s = sin(element_angle);
    e_const = (e*a)/element_length;
    
    trig_matrix = [c,s,0,0;0,0,c,s]; %trig matrix to simplify calculations in the future
    thing = [1,-1;-1,1]; %thing matrix for force calculation
    local_disp_matrix = [disp_matrix(m11,1);disp_matrix(m12,1);disp_matrix(m21,1);disp_matrix(m22,1)]; %take values from the displacement matrix that match the nodes' x and y directions
    
    element_force = e_const * thing * trig_matrix * local_disp_matrix; %element forces at nodes in local coordinates, as taught
    explanation = ['Elemental forces (in Newtons) for nodes ', num2str(n1), ', ', num2str(n2)];
    
    stress_matrix = [-1/element_length, 1/element_length]; %element axial stresses at nodes in local coordinates, as taught
    element_stress = e * stress_matrix * trig_matrix *local_disp_matrix;
    otherexplanation = ['Stress (in Pascals) for element ', num2str(element_counter)];
    
    %show the answers for each element
    disp(explanation)
    disp(element_force)
    disp(otherexplanation)
    disp(element_stress)
    
end
%% Problem 4

clc
clear
format shortg

nodes = 8;
elements = 16;
%constants
e = 70e9;
a = 3e-6; %convert 3 cm^3 to m^3
k = (e*a);

%element relations (nodes connected and locations)
element_relations = [1,2; 2,4; 4,6; 6,8; 7,8; 5,7; 3,5; 1,3; 1,4; 3,4; 3,6; 5,6; 5,8; 7,6; 5,4; 3,2];
x_coor = [0,0,3,3,6,6,9,9];
y_coor = [0,3,0,3,0,3,0,3] ;

K_matrix = zeros(2*nodes,2*nodes);
disp_matrix = zeros(2*nodes,1);
force_matrix = zeros(2*nodes,1);

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
    
    %add element matrix elements to correct spot in global stiffness matrix
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
%force BCs
force_matrix(4,1) = - 50000;
force_matrix(8,1) = - 100000;
force_matrix(12,1) = - 50000;
force_matrix(16,1) = - 50000;
K_modified = K_matrix;

%support BCs (node 1 is pinned, and node 7 is a roller)
fixed_nodes = [1,2,14];
for fixed = fixed_nodes    
    K_modified (fixed,:) = 0;
    K_modified(fixed,fixed) = 1;
end

disp_matrix = K_modified\force_matrix;
%disp('nodal displacements')
%disp(disp_matrix)
disp('x and y displacements for asked nodes (in meters), x first followed by y displacement')
disp('node 3')
disp(disp_matrix(5,1))
disp(disp_matrix(6,1))
disp('node 5')
disp(disp_matrix(9,1))
disp(disp_matrix(10,1))
disp('node 7')
disp(disp_matrix(13,1))
disp(disp_matrix(14,1))

element_stresses = zeros(elements,1);

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

    delta_x = x2-x1;
    delta_y = y2-y1;
    element_length = sqrt((delta_y)^2+(delta_x)^2);
    element_angle = atan2(delta_y,delta_x);
    c = cos(element_angle);
    s = sin(element_angle);
    e_const = (e*a)/element_length;
    trig_matrix = [c,s,0,0;0,0,c,s];
    
    local_disp_matrix = [disp_matrix(m11,1);disp_matrix(m12,1);disp_matrix(m21,1);disp_matrix(m22,1)];
    stress_matrix = [-1/element_length, 1/element_length];
    
    element_stress = e * stress_matrix * trig_matrix *local_disp_matrix;
    %element_stresses(element_counter,1) = element_counter; I don't need to
    %annouce the element with the highest stress right ?
    element_stresses(element_counter,1) = element_stress;
    
    otherexplanation = ['Stress for element ', num2str(element_counter)];
    %disp(otherexplanation)
    %disp(element_stress) 
   
end

disp('The max magnitude axial stress (in Pa) occuring in an element is:')
disp(max(abs(element_stresses)))