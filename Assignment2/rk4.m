function [t, Y] = rk4(f, T, y0, h)
%[t, Y] = rk4(f, [t0, tk], y0, h) poisce priblizek resitve diferencialne
%enacbe y' = f(t, y) z zacetnim pogojem y(t0) = y0. Uporabi Runge-Kutta 
%metodo 4. reda s korakom dolzine h.

t = T(1):h:T(2);
Y(:, 1) = y0;
for k = 1:(length(t) - 1)
	k1 = h*f(t(k), Y(:, k));
	k2 = h*f(t(k) + h/2, Y(:, k) + k1/2);
	k3 = h*f(t(k) + h/2, Y(:, k) + k2/2);
	k4 = h*f(t(k) + h, Y(:, k) + k3);
	Y(:, k + 1) = Y(:, k) + (k1 + 2*k2 + 2*k3 + k4)/6;
end