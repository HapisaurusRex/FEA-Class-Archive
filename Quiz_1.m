% Finding abundant numbers - sum of divisors including the number itself,
% n,must be greated than 2*n
% This can done the easy way or the hard way. See the second half of the
% code for the hard way. 

%% Part A
n=1
while (n <= 1000)
    test = sum(divisors(n));
    if test > 2*n
        disp(n)
    end
    n = n + 1;
end
%% Part B

for n = 1:1:1000
    test = sum(divisors(n));
    if test > 2*n
        disp(n)
    end
end

%% Part A - Hard Way
clear
clc

n=1;
while (n <= 1000)
    i = 1;
    sum = 0;
    while (i<=n)
        if rem(n,i) == 0
            sum = sum + i;
        end
        i = i + 1;
    end
    if sum > 2*n
            disp(n)
    end
    n = n + 1;
end
%% Part B - Hard way
clear
clc

for n = 1:1:1000
    sum = 0;
    for i = 1:1:n
        if rem(n,i) == 0
            sum = sum + i;
        end
    end
    if sum > 2*n
        disp(n)
    end
end