function [e] = eqrm(s_or, s_ric, t_poz)
N = length(s_or);
m = length(t_poz);

i=0;
num = 0;
for i = 1:m
    num = num + (s_ric(t_poz(i)) - s_or(t_poz(i)))^2;
end
num = num/m;

i=0;
den = 0;
for i=1:N
    den = den + (s_or(i))^2;
end

den = den/N;

e = num/den;

end

