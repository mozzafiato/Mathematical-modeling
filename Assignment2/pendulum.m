function pos = pendulum(Phi0, Theta, T)
%za zacetne pogoje Phi0 = [phi0, omega0]T in funkcijo Theta(t) = [theta(t),  theta*(t), theta**(t)]T 
%resuje diferencialno enacbo iz 3. tocke, do koncega casa T 
%in vrne polozaj [x, y]T , v katerem je tockasta masa ob casu T.  
%Phi0 je vektor: Phi0 = [phi0; omega0], Theta je funkcija spremenljivke t,
% theta(t), T je stevilo. Pri resevanju bomo privzeli, da za nase nihalo velja 
%g = m = l = L = 1,   -1 < x0 < 1.  

%[t, Y] = rk4(f, [t0, tk], y0, h) poisce priblizek resitve diferencialne
%enacbe y' = f(t, y) z zacetnim pogojem y(t0) = y0. Uporabi Runge-Kutta 
%metodo 4. reda s korakom dolzine h.

%inicijalizacija
h=0.1;
t = 0:h:T;
L=1;
l=1;
g=1;

omega0 = Phi0(2);
phi0 = Phi0(1);

%pripravimo zacetni vektor Y
Y = zeros(2, length(t));
Y(1, 1) = phi0;
Y(2,1) = omega0;

%definiramo funkcijo f(t,Y)
f=@(t, Y) [Y(2); (-g*sin(Y(1))-l*Theta(t)(3)*cos(Theta(t)(1)-Y(1))+l*(Theta(t)(2))^2*sin(Theta(t)(1)-Y(1)))/L];

%metoda rk4
for k = 1:(length(t) - 1)
	k1 = h*f(t(k), Y(:, k));
	k2 = h*f(t(k) + h/2, Y(:, k) + k1/2);
	k3 = h*f(t(k) + h/2, Y(:, k) + k2/2);
	k4 = h*f(t(k) + h, Y(:, k) + k3);
	Y(:, k + 1) = Y(:, k) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

%resitev za phi
phi = Y(1,end);
%theta v casu T
theta = Theta(t(k+1))(1);

%izracunamo polozaj
x=l*sin(theta)+L*sin(phi);
y=-l*cos(theta)-L*cos(phi);
%rezultat
pos = [x; y];
end
%!test
%! Theta = @(t) [ 1/10+cos(3/5*t), -3/5*sin(3/5*t), -3/5*cos(3/5*t) ];
%! result = pendulum( [pi/6;0], Theta, 5);
%! assert( result, [  -1.699924695833261; -1.014559024821094 ] , sqrt( sqrt(eps) )  );

