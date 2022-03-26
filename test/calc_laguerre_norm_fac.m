target_fun = @(x, p, l) abs(x).^abs(l) .* exp(-x.^2) .* laguerreL(p, abs(l), 2.0*x.^2);
p_max = 2;
l_max = 2;
norm_fac = zeros(1+p_max, 1+l_max);

for p = 0:p_max
    for l = 0:l_max
        [~, norm_fac(p+1, l+1)] = fminsearch(@(x)-target_fun(x,p,l), 0);
    end
end