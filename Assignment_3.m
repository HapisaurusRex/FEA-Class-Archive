%% Problem 1a (done)
clear
x = linspace(0,pi()/2,200);
y = x.^2 .* sin((50*x))
plot (x,y)
ylabel('y')
xlabel('x')
title('A plot of Y versus X')

%% Problem 1b (done)
clear
M = load('TheoryExperiment.txt')
x = M(:,1)
y_theory = M(:,2)
y_experiment = M(:,3)

graph = figure;
plot (x,y_theory,'-r')
hold on
plot (x,y_experiment,'*b')
xlabel('Time [s]')
ylabel('Concentration [mM]')
legend('Theory','Experiment')


%% Problem 2a (done)
clear
clc

f = 1000;
t = linspace(0,(4/f),1000);
c = 2 * pi() * f * t
pie = 4/pi
check  = 1
x_square = zeros(size(t))
hold on

for k = 1:1:26
    i = 2*k -1;
    x_square = x_square + pie * (sin(i *c))/i;
    if k == check
        figure(k) %remove this comment to group all graphs together
        plot(t,x_square)
        check = check + 5;
    end
end

%% Problem 2B(don)
clear
clc

f = 1000;
t = linspace(0,(4/f),1000);
c = -2 * pi() * f * t
pie = 2/pi
check  = 1
x_square = zeros(size(t))
hold on

for k = 1:1:26
    x_square = x_square + pie * (sin(k *c))/k;
    if k == check
        figure(k) %remove this comment to group all graphs together
        plot(t,x_square)
        check = check + 5;
    end
end

%% Problem 3 (done)
clear
x = linspace(-2,2,20);
y = linspace(-2,2,20);

[X,Y] = meshgrid(x,y);
Z = 2*X.^2 - Y.^2;
surf(X,Y,Z);
xlabel('X')
ylabel('Y')
zlabel('Z')

%% Problem 4 (done)
clear
x_h = -0.038;
y_h = 0.022;
z_h = 0.29;
P_vector = [-0.7, 0, -0.72]; %figure out that vector
P_x = 0.0275 * P_vector(1,1)
P_y = 0.0275 * P_vector(1,2)
P_z = 0.0275 * P_vector(1,3)
H_torso = 0.45;
R_torso = 0.15;
sigma = 0.09;
%copy and pasted
r = linspace(0,0,100)+R_torso;
[X,Y,z] = cylinder(r,100);
Z = z*H_torso;
phi = (P_x*(X-x_h) +P_y *(Y-y_h) + P_z*(Z-z_h))./(4*pi*sigma*(sqrt((X-x_h).^2+(Y-y_h).^2 + (Z-z_h).^2)).^3);
surf(X,Y,Z)
axis equal;
surface(X,Y,Z,phi)
shading interp

%% Problem 5 (done)
clear
n = [1,2,5,10,15,20];
delta_x = 1./n;
actual = log(2);
for a = 1:1:6
    interval = linspace(1,2,n(1,a));
    test = size(interval);
    values = 1./interval;
    i = 0;
    for k = 1:1:(test(1,2)-1)
        i = i + delta_x(1,a) * (values(1,k)+values(1,(k+1)))/2;
    end
    er(1,a) = abs(i-actual)/actual;
end

plot(n,er,'--*')
title('Error on trapozoidal rule calculating ln(2) vs number of intervals')
xlabel('Number of intervals')
ylabel('Estimation Error')

%% Problem 6A (done)
clear
clc
x_0 = 1;
x  = 0;
counter = 1;
tolerance = 0.01;
sigma = 3;
while sigma > tolerance
    x = cos(x_0).^2 +(pi/3)
    sigma = abs(x-x_0)
    x_0 = x;
    counter = counter + 1

end

%% Problem 6B (done?)

%this method does not work for this equation because the values do not
%converge and oscillate between 1.0032 and 1.5882 meaning a tolerance below
%0.05851 (the difference between the 2 values) is impossible and that value
%is the only way to get the loop to end

% or there are two roots that solve this problem

%
clear
clc
x_0 = 1;
x  = 0;
counter = 1;
tolerance = 0.5851;
sigma = 3;
while sigma > tolerance
    x = x_0^2 *cos(x_0) +(pi/3)
    sigma = abs(x-x_0)
    x_0 = x;
    counter = counter + 1

end
%


%% Problem 7 (done)
clear
clc

x1=3;
x2=7;
tol = 0.0001;
er = 3;


if  x1- x1^2 *cos(x1) -(pi/3)> 0
    xp = x1; %xp = x+
    xm = x2; %xm = x-
else
    xp = x2;
    xm = x1;
end

while abs(er)>tol
    fxp =xp- xp^2 *cos(xp) -(pi/3); %f(x+) done to keep my eyeballs working
    fxm =xm- xm^2 *cos(xm) -(pi/3); %f(x-)
    x3 = xp - ((fxp *(xm - xp))/(fxp-fxm))
    er = x3 - x3^2 *cos(x3) -(pi/3) 
    if er >0
        xp = x3;
    else
        xm =x3;
    end
end
