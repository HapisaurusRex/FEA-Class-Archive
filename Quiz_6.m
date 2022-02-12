clear
clc
format shortg


clear
clc
format shortg

syms f1(x,y)
f1(x,y) = x^2 - sqrt(3) * x*y + 2*y^2 -10;
syms f2(x,y)
f2(x,y) = 4*x^2 + 3*sqrt(3)*x*y + y^2 - 22;
J = [diff(f1(x,y),x), diff(f1(x,y),y);
     diff(f2(x,y),x), diff(f2(x,y),y)];
F = [f1(x,y);f2(x,y)];

syms r(x,y)
r(x,y) = J;

syms q(x,y)
q(x,y) = F;

x1 = [1;1]; %x1 = i, x2 = i+1
%j = r(x1(1),x1(2))
%s = q(x1(1),x1(2))

error = 10;
tolerance = 1*10^-4;

while error>tolerance
   j = double(r(x1(1),x1(2)));
   s = double(q(x1(1),x1(2)));
   x2 = x1 - (inv(j)*s);
   error =  norm(x1 -x2);
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