%% Problem 1 (done)

clear
clc

syms f(x)
f(x) = log(x) * (2*x - 2*x^4 - 2*x^2) + x^3 + 1; % simplified version of eq given
J = diff(f(x),x);

syms j(x)
j(x) = J;

x1 = 2; %initial guess
tolerance = 10^-4;
error = 10;
while error>tolerance
    x2 = x1 - double((f(x1)/j(x1)));
    error = abs((x2-x1)/x2);
    disp('initial value')
    disp(x1)
    disp('function value')
    disp(double(f(x1)))
    disp('function derivative')
    disp(double(j(x1)))
    x1 = x2;
    disp('new initial value')
    disp(x1)
    disp(' ')
    disp('below is the next iteration')
    disp(' ')
end
disp('the final answer')
disp(x1)
%% Problem 2

clear
clc
format shortg

syms f1(x,y)
f1(x,y) = sin(2*x+y) - exp(x-y)
syms f2(x,y)
f2(x,y) = cos(x+5) - (x*y)^2;

J = [diff(f1(x,y),x), diff(f1(x,y),y);
     diff(f2(x,y),x), diff(f2(x,y),y)];
F = [f1(x,y);f2(x,y)];

syms r(x,y)
r(x,y) = J;

syms q(x,y)
q(x,y) = F;

x1 = [0;0]; %x1 = i, x2 = i+1
%j = r(x1(1),x1(2))
%s = q(x1(1),x1(2))

tolerance = 1*10^-4;
error = [10;10];
%counter = 1;
iterative_table = {};

while error>tolerance
   j = double(r(x1(1),x1(2)));
   s = double(q(x1(1),x1(2)));
   x2 = x1 - (inv(j)*s);
   error = abs((x2-x1)./x2);
   %{
   iterative_table{counter,1} = {x1}; technically a table, and you can
   technically read it, but its hard to read so i'm using the code below
   iterative_table{counter,2} = {s};
   iterative_table{counter,3} = {j};
   iterative_table{counter,4} = {x2}; 
   counter = counter + 1;
   %} 
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

disp('the final answer')
disp(x1)

x1 = [1;1];
tolerance = 1*10^-4;
error = [10;10];
while error>tolerance
   j = double(r(x1(1),x1(2)));
   s = double(q(x1(1),x1(2)));
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

disp('the final answer')
disp(x1)

x1 = [0.5;0.25];
tolerance = 1*10^-4;
error = [10;10];
while error>tolerance
   j = double(r(x1(1),x1(2)));
   s = double(q(x1(1),x1(2)));
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

disp('the final answer')
disp(x1)

x1 = [0;3];
tolerance = 1*10^-4;
error = [10;10];
while error>tolerance
   j = double(r(x1(1),x1(2)));
   s = double(q(x1(1),x1(2)));
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

disp('the final answer')
disp(x1)
%% Problem 3 (Done)

clear
clc
format shortg

syms f(x)
f(x) = x^4 - 2*x + 1;
a = 0;
b = 2;

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
%% Problem 4 (done)

clear
clc
format shortg

syms f(x)
f(x) = exp(-x^3);
a = 0;
b = 2;

%evaluating directly
direct = int(f(x),a,b);
disp('evaluating directly:')
disp(double(direct))

%Gausian Quadrature (GQ) - 3 point (REFER TO TABLE FOR E and W)
syms r(g)
h = (b-a)/2 ;
r(g) = vpa(h * g + (b+a)/2);
dx = diff(r(g),g);

gp = [sqrt(3/5); 0; -sqrt(3/5)]; %gauss point refer to table and reference to number of points needed
w = [5/9; 8/9; 5/9];  %weights refer to table and refernece to number of points requested

e1 = r(gp(1,1));
e2 = r(gp(2,1));
e3 = r(gp(3,1));

NC = vpa(h * ((w(1,1)*f(e1)) + (w(2,1)*f(e2))+ (w(3,1)*f(e3))));

disp('evaluating using 2 poitn Gaussian Quadrature')
disp(double(NC))
%% Problem 5

clear
clc
format shortg


pressure = [336, 294.4, 266.4, 260.8, 260.5, 249.6, 193.6, 165.6];
volume = [0.5, 2, 3, 4, 6, 8, 10, 11];

work  = zeros(7,3);

for counter = 1:7
    a = volume(counter);
    b = volume(counter +1);
    %Trapezoidal Rule
    f0 = pressure(counter);
    f1 = pressure(counter +1);
    trap = (b-a) * (f0 + f1)/2;
    work(counter,1) = trap;
    %Simpson's 1/3 rule
    h = (b-a)/2;
    f = linspace (f0,f1,3);
    onethird = (b-a) * (f(1) + (4*f(2))+ f(3))/6;
    work(counter,2) = onethird;
    %Newton-Cotes (NC) - Simpson's 3/8 rule
    h = (b-a)/3;
    f = linspace (f0,f1,4);
    three8th = (b-a) * (f(1) + (3*f(2))+ (3*f(3))+ f(4))/8;
    work(counter,3) = three8th;
end

%uncoment to show work calcualated at each step disp(work)

total_work = sum(work)

disp('total work calculated using trapezoidal rule (kJ):')
disp(total_work (1,1))

disp('total work calculated Simpson''s 1/3 rule (kJ):')
disp(total_work (1,2))
disp('total work calculated Simpson''s 3/8 rule (kJ):')
disp(total_work (1,3))
%{
%Newton-Cotes (NC) - Simpson's 1/3 rule
h = (b-a)/2;
f0 = double(f(a + 0*h)) ;
f1 = double(f(a + 1*h)) ;
f2 = double(f(a + 2*h)) ;
NC = (b-a) *(f0 + (4*f1) + f2)/6;

%Newton-Cotes (NC) - Simpson's 3/8 rule
h = (b-a)/3;
f0 = double(f(a + 0*h)) ;
f1 = double(f(a + 1*h)) ;
f2 = double(f(a + 2*h)) ;
f3 = double(f(a + 3*h)) ;
NC = (b-a) *(f0 + (3*f1) + (3*f2) + f3)/8;
%}
%% Problem 6 (done)

clear
clc
format shortg

syms f(x,y)
f(x,y) = (x^2-x)* sin(y);

a1 = 0;
b1 = 3;
a2 = 0;
b2 = pi;

first = int(f(x,y),a1,b1);
direct = int(first, a2,b2);

disp('evaluating directly:')
disp(double(direct))

%Gausian Quadrature (GQ) - 2 point (REFER TO TABLE FOR E and W)
q = (b1-a1)/2;
r = (b2-a2)/2;


gp = [1/sqrt(3); -1/sqrt(3)]; %gauss point refer to table and reference to number of points needed
w = [1; 1];  %weights refer to table and refernece to number of points requested

syms s(g)
s(g) = vpa(q * g + (b1+a1)/2);
dx = diff(s(g),g);

syms t(n)
t(n) = vpa(r * n + (b2+a2)/2);
dy = diff(t(n),n);

syms h(g,n)
h(g,n) = f(s(g),t(n));

c1 = w(1,1) * w(1,1);
c2 = w(1,1) * w(2,1);
c3 = w(2,1) * w(1,1);
c4 = w(2,1) * w(2,1);

NC = double((q*r) * ((c1 * h(gp(1,1),gp(1,1))) + (c2 * h(gp(1,1),gp(2,1))) + (c3 * h(gp(2,1),gp(1,1))) +  (c4 * h(gp(2,1),gp(2,1))) ));

disp('evaluating using 2 point Gaussian Quadrature')
disp(double(NC))