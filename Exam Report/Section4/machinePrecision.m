
epsilon = 1;

while 1 + epsilon > 1
    epsilon = epsilon/2;
end

fprintf("\nMachine precision: %d (base 10) \n", round(log10(epsilon),0));
fprintf("Machine precision: %d (base 2) \n", round(log10(epsilon)/log10(2),0));
