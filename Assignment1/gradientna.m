function [d, t, s] = gradientna(F, G, t0, s0)
%Funkcija gradientna.m ima kot vhodne argumente dve funkciji F in G
%ki sta bodisi krivulji bodisi ploskvi. 
%Tocka t0 (funkcija F) in tocka s0 (funkcija G) sta zacetni tocki
%na krivuljah ali pa sta dva dvodimenzionalna vektorja 
%(z dimenzijo 2x1) na ploskvami. 
   
  maxit=10000;
  tol=1e-10;
  
  for k=1:maxit
   
      %H(t,s) = ||F(t0) - G(s0)||^2 
 
      %vrednosti v tocki t0
      [f,df,ddf]=feval(F,t0);
      %vrednosti v tocki s0
      [g,dg,ddg]=feval(G,s0);
      
      h1=2*(f-g)'*df; %odvod funkcije po t0
      h2=-2*(f-g)'*dg; %odvod funkcije po s0
      
      if(rows(t0)==2 && rows(s0)==2)
      %ce gre za ploskve transponiramo (da se dimenizije ujemajo)
        h1=h1';
        h2=h2';
      endif
      %dolocimo matriko H ki vsebuje gradiente 
      H=[h1;h2];
      
    if(norm(H)<tol)
      break;
    endif
    %naredimo gradientni spust
    t0=t0-0.01*h1;
    s0=s0-0.01*h2;
  end 
  t=t0;
  s=s0;
  r=(feval(F,t)-feval(G,s));
  d=sqrt(sum(r.^2));
  
  %Izpisemo opozorilo v primeru, da zadnji priblizek ni znotraj tolerancnega obmocja.
  if(k == maxit)
	  disp("Warning: The method did not converge after maxit iterations.")
  endif
  
endfunction

%!test
%! tZac= 1;
%! sZac= -1;
%! [d,t,s]=gradientna( @kroznica1, @kroznica2 ,tZac , sZac);
%! assert( d, 5, eps  );
%! assert( t , pi/2, sqrt(eps) );
%! assert( s, -pi/2, sqrt(eps) );
