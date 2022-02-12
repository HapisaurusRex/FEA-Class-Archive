%% Quiz 4 (plane strain) (wrong Von Mises)

clear
clc
format longg

nodes = 3;
elements = 1;
elation = [1,2,3]; %element nodal relations

cc = 2*nodes; %creation constant
b_matrix = zeros(3,cc);

t = 1; %thickness (mm) %maybe not needed


e = 6e6; %in psi
v = 0.35; %poisson ratio
x_coor = [1,6,3] ; %coordinates in in
y_coor = [1,3,5];

x_disp = [-0.001, 0.004, 0];
y_disp = [-0.002, 0, 0.003];
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

%setup D matrix for plane strain (new matrix elements and d_const)
d_matrix = eye(3);
d_matrix(1,1) = 1-v;
d_matrix(2,2) = 1-v;
d_matrix(1,2)= v;
d_matrix(2,1) = v;
d_matrix(3,3) = (1-(2*v))/2;
d_const = e/((1+v)*(1-2*v));
d_matrix = d_const .* d_matrix;

%1. find element stress components
stress = d_matrix * b_matrix * disp_matrix;

disp('Element Stress Components')
disp(stress) %first row = sigma_x, second row = sigma_y, third = tau_xy

%2. find element planar principal stresses
stress_avg = (stress(1,1) + stress(2,1))/2;
useful = (stress(1,1) - stress(2,1))/2;
r = sqrt((useful^2)+(stress(3,1)^2));
principal = [(stress_avg+r) ; (stress_avg-r)]; %sigma_1 and sigma_2

disp('Planar Principal Stresses')
disp(principal)

%3. find associapted principal angles
placeholder = 0.5*atan2d(stress(3,1),useful);
principal_angles = [placeholder; placeholder+90];

disp('Associated Principal Angles')
disp(principal_angles)

%4. find von mises stress
s1 = stress(1,1)^2;
s2 = stress(1,1) * stress(2,1);
s3 = stress(2,1)^2;
s4 = 3*stress(3,1)^2;
vm = sqrt((s1-s2)+s3+s4);

disp('Von Mises Effective Stress')
disp(vm)

%Tresca
sigma_1 = principal(1,1);
sigma_2 = principal(2,1);
sigma_3 = v* (stress(1,1) + stress(2,1));

t1 = abs(sigma_1 - sigma_2);
t2 = abs(sigma_2 - sigma_3);
t3 = abs(sigma_1 - sigma_3);
tresca_criteria = [t1;t2;t3];

disp('Tresca Effective Stress')
disp(max(tresca_criteria))
