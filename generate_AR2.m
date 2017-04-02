function x = generate_AR2(a1, a2, sigma, len, figures)

v = sqrt(sigma)*randn(len,1);

x = filter(1, [1, a1, a2], v);

if nargin > 4
    plot(x)
end
    