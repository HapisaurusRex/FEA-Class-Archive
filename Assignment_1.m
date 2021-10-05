format shortG %stop with the scientific notation

%problem 1 (done)
%define variables
length = 15.6;
height = 6.3;
width = 5;
%calculates what's stated
volume = width*length*height;
surface_area = 2 * length * width + 2*length*height + 2*height*width;

%problem 2 (done)
%define variables
x = 0.88;
y = 2.35;
z = -0.7;
%plug variables into equations below
a = x^2+y^2;
b = (x-1)/(2*x*y - z);
c = y^(x-1)/ sqrt(z^2+1);
d = 1/(x-1) + y/(z-1);
e = (2*y^2 + 9/7*x^2*z^2 - 7*x*z^3) /(8/3 * z^2 * y^3 + 5*x^4);

%problem 3 (done)
a = 3;
b = 4;
c = 5;
P = (a+b+c)/2;
Area = sqrt(P*(P-a)*(P-b)*(P-c));

a = 3.5; %note this variable changed
b = 4;
c = 5;
P = (a+b+c)/2;
Area = sqrt(P*(P-a)*(P-b)*(P-c));

%problem 4 (done)
lower_radius = 12.8;
upper_radius = 11.35;
height = 29.4;
volume = (pi()*height)/3 *(upper_radius^2 + upper_radius * lower_radius + lower_radius^2);

%problem 5 (done)
g=9.81;
initial_height=3200;
time = transpose(linspace(0,10,11));
elevation =  -0.5 * g *time.^2 + initial_height;

%problem 6 (done)
v = [1,2,3,4,5,6,7,8,9];
v = linspace(1,9,9);
v = 1:1:9; %first no:spacing:last no

A = [v; 2*v; 3*v; 4*v; 5*v; 6*v; 7*v; 8*v; 9*v];
v_transpose = transpose(v);
A = v_transpose * v;

%problem 7 (done)
A =[-1, 6, -8, 7, -11;
    0, 4, 9, -2, 9;
    1, 2, -1, -1, 8;
    7, 0, 1, 5, 5;
    3, -3, 0, 0, 1;]; %create matrix of coefficents for all x values
B = [36; -29; -15; 1; 1]; %create matrix of constant terms
solution = inv(A) * B; %values for x listed from 1 to 5 (note x2 is 0 eventhough it says -1.9 * 10^-15)

%problem 8 (done)
syms d
x = (2*d+1)/2;
y = d-1;
z = 3-4*d;
P = 7*x +4*y -7*z == 8;

d = vpa(solve(P,d));
x = (2*d+1)/2;
y = d-1;
z = 3-4*d;

%problem 9 (done)
v1 = [3, 0, 0];
v2 = [0, 2, 0];
angle = acosd(dot(v1,v2)/(norm(v1) *norm(v2))); %answer in degrees (remove d for radians)

u1 = [3, 1, -1];
u2 = [-1, 2, 1];
angle = acosd(dot(u1,u2)/(norm(u1) *norm(u2))); %answer in degrees (remove d for radians)

w1 = [0, 6, 8];
w2 = [1.3, 0, 1.3];
angle = acosd(dot(w1,w2)/(norm(w1) *norm(w2))); %answer in degrees (remove d for radians)

%problem 9b (done)

p1 = [1,1,1];
p2 = [-1,-1,1];
p3 = [-1,1,-1];
p4 = [1,-1,-1];
%tetrahedron has 3 faces connected to a point, which we use p4 as that
%point, meaning there are 3 angles to be found
%turn the points into vectors using p4 as the base (see below for angle
%calculation)

%for the angle between p4p1 and p4p2
v1 = p4 - p1;
v2 = p4 - p2;
angle = acosd(dot(v1,v2)/(norm(v1) *norm(v2))); %answer in degrees (remove d for radians)

%for the angle between p4p1 and p4p3
v1 = p4 - p1;
v2 = p4 - p3;
angle = acosd(dot(v1,v2)/(norm(v1) *norm(v2)));

%for the angle between p4p2 and p4p3
v1 = p4 - p2;
v2 = p4 - p3;
angle = acosd(dot(v1,v2)/(norm(v1) *norm(v2)));

%problem 10 (done)
n= 0:1:10;
test = (-1).^n ./ (2.*n +1);
pi_estimate = 4* sum (test);
error = 100 * abs ((pi_estimate - pi())/pi());

n= 0:1:100;
test = (-1).^n ./ (2.*n +1);
pi_estimate = 4* sum (test);
error = 100 * abs ((pi_estimate - pi())/pi());

n= 0:1:1000;
test = (-1).^n ./ (2.*n +1);
pi_estimate = 4* sum (test);
error = 100 * abs ((pi_estimate - pi())/pi());

n= 0:1:10000;
test = (-1).^n ./ (2.*n +1);
pi_estimate = 4* sum (test);
error = 100 * abs ((pi_estimate - pi())/pi());

n= 0:1:100000;
test = (-1).^n ./ (2.*n +1);
pi_estimate = 4* sum (test);
error = 100 * abs ((pi_estimate - pi())/pi());

