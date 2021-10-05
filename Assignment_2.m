format shortG %stop with the scientific notation

%problem 1
disp('This program gives roots for the equation: ax^2 + bx + c = 0')
disp('Please input the following with respect to the equation above')
a = input('Please enter a value for a: ')
b = input('Please enter a value for b: ')
c = input('Please enter a value for c: ')

determinant = b^2 - 4 * a * c;
if determinant < 0
    disp('Error: No Real Roots Mathematically Possible')
elseif  determinant > 0
    disp('Two Real Roots Exist')
    root_1 = (-b + sqrt(determinant))/(2*a)
    root_2 = (-b - sqrt(determinant))/(2*a)
else
    disp('One Real Root Exists')
    root = -b/(2*a)
end


%problem 2 (done)
%{
x = input('Input a value for x:')
if x ==0
    x = 1
elseif x>0 && x<2
    x = -2*x + 2
elseif x>=2 && x<3
    x = -x
elseif x>-1 && x<0
    x = -2*x
elseif x<= -1 || x>= 3
    x = (x^2) - 2*x
end
%}

%problem 3 (done)

%{
sanity_check = 0;
while (sanity_check ~= 1)
    n = input('Please input a positive integer factorial you want calculated: ');
    disp(' ')
    if (n>0 && fix(n) == n) %integer check via fix, positive check via n>0
        sanity_check = 1;
    else
        disp('The number you entered is not a positive integer')
    end
end

factorial_n = 1;
x = 1;
while (x < n+1)
    factorial_n = factorial_n * x;
    x = x+1;
end

disp('The factorial/product of all positive integers from 1 to your number is:')
disp(factorial_n)

%}

%problem 4 (done)

%{
dont_tell = 123
disp('I am thinking of a number between 0 and 1000')
guess = 0
guess_counter = 1
while(guess ~= dont_tell)
    guess = input('Can you guess it ?: ')
    if guess == dont_tell && guess_counter<10
        disp('Awesome ! you guessed the secret number with fewer than 10 attempts')
    elseif guess ~= dont_tell %this can also just be an else statement and line 77
        disp('Try again !')
        guess_counter = guess_counter + 1
    end
end
%}

%problem 5 (done)

%{
sanity_check = 0;
while (sanity_check ~= 1)
    n = input('What do you reckon is a prime integer ?: ');
    disp(' ')
    if (n>0 && fix(n) == n) %integer check via fix, positive check via n>0
        sanity_check = 1;
    else
        disp('Get outta here partner, your number is a lie !')
    end
end
x = 2;
r = 0;
while(x<= n-1)
    rem(n,x);
    if rem(n,x) == 0
        disp('That aint no prime number partner')
        x = n-1;
        r=1;
    elseif rem(n,x)~= 0
        x = x+1;
    end
end
if (r==0)
    disp('Well shoot I guess it is a prime number partner')
    disp('The prime number: ')
    disp(n)
end
%}

%problem 6 (done)

%{
disp('Here are all the perfect numbers between 1 and 10000')
n = 1;
while (n <= 10000)
    dummy_variable = divisors(n);
    answer = sum(dummy_variable,'all') - n;
    if answer == n
        disp(n)
    end
    n = n+1;
end

%}

%problem 7A (done)

%{
n = 1000
counter = 0;  
for i = 1:n
    x0 = rand(1, 1);  
    y0 = rand(1, 1);  
    r0 = sqrt(x0^2 + y0^2);
    if r0 <= 1
        counter = counter + 1; 
    end
end
pi_value = counter/n*4;
display(pi_value)
%}

%problem 7b (done)

%{
error = 1
n = 0;
counter = 0; 
while error > 10^-6 %change to 9, but warning it takes forever
    n = n + 1;
    x0 = rand(1, 1);  
    y0 = rand(1, 1);
    r0 = sqrt(x0^2 + y0^2);
    if r0 <= 1
        counter = counter + 1; 
    end
    pi_value = counter/n*4;
    error = abs(pi()-pi_value)
end

disp('Times the program has been repeated: ')
disp(n)
%}

%problem 8a (done)

%{
x = pi()/6
n = 0:1:9
factor = (2*n) + 1;
denominator = factorial(factor);
test = (-1).^n .* (x.^factor);
series = test./denominator
answer = sum(series,'all')
%}

%problem 8b (done)

%{
x = pi()/6;
series = 0;
for n = 0:1:9
    factor = (2*n) + 1;
    denominator = factorial(factor);
    test = (-1).^n .* (x.^factor);
    disp(test./denominator)
    %series = series + (test./denominator) %uncoment for the sum of these
    %terms
end
%}

%problem 8c (done)

%{
difference = 1;
x = pi()/6;
n = 0;
first_term = 0;
second_term = 0;

while difference > 10^-11
    factor = (2*n) + 1;
    denominator = factorial(factor);
    test = (-1).^n .* (x.^factor);
    first_term =  first_term + (test./denominator);
    difference = abs(first_term - second_term)
    second_term = second_term + (test./denominator);
    n = n + 1
end
%}