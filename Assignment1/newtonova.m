function [d, t, s] = newtonova(F, G, t0, s0)
%Newtonova metoda ki poisce najmanjso razdaljo med krivuljami in ploskvami
%F in G sta funkciji podani v obliki F=[ff,df,ddf] (G=[gg,dg,gg]) kjer je :
%ff parametricno podana funkcija
%df prvi odvod funkcije
%ddf drugi odvod funkcije
%t0 (funkcija F) in s0 (funkcija G) sta zacetni tocki na krivulji, ali
%t0 in s0 sta dva dvodimenzionalna vektorja na ploskvami
%default
maxit=1000;
tol=1e-10;

for k = 1:maxit
  %vrednosti v tocki t0
  [f,df,ddf]=feval(F,t0);
  %vrednosti v tocki s0
  [g,dg,ddg]=feval(G,s0);

  if(rows(t0)==1 && rows(s0)==1)
  %krivulje
    h1=(f-g)'*df;
    h2=(f-g)'*dg;
    H=[h1 ; h2];
    
    %Jacobieva matrika JH
    h1f=df'*df+((f-g)'*ddf);
    h2f=(-df'*dg);
    h1g=dg'*df;
    h2g=(-dg'*dg)+((f-g)'*ddg);
    JH=[h1f h2f; h1g h2g];  
  elseif(rows(t0)==2 && rows(s0)==2)
   %ploskve
    h1=[ ((f-g)'*df(:,1)) ; ((f-g)'*df(:,2)) ];
    h2=[ ((f-g)'*dg(:,1)) ; ((f-g)'*dg(:,2)) ];
    H=[ h1; h2 ];
  %Jakobieva matrika JH
    JH= [df(:,1)'*df(:,1)+(f-g)'*ddf(:,1,1), df(:,2)'*df(:,1)+(f-g)'*ddf(:,1,2), -dg(:,1)'*df(:,1),                   -dg(:,2)'*df(:,1);
         df(:,1)'*df(:,2)+(f-g)'*ddf(:,2,1), df(:,2)'*df(:,2)+(f-g)'*ddf(:,2,2), -dg(:,1)'*df(:,2),                   -dg(:,2)'*df(:,2);     
         df(:,1)'*dg(:,1),                   df(:,2)'*dg(:,1),                   -dg(:,1)'*dg(:,1)+(f-g)'*ddg(:,1,1), -dg(:,2)'*dg(:,1)+(f-g)'*ddg(:,1,2);
         df(:,1)'*dg(:,2),                   df(:,2)'*dg(:,2),                   -dg(:,1)'*dg(:,2)+(f-g)'*ddg(:,2,1), -dg(:,2)'*dg(:,2)+(f-g)'*ddg(:,2,2)];
  endif
    
  tocki0=[t0;s0];
	%Izvedemo en korak Newtonove metode...
	tocki=tocki0-JH\H;
	%... in testiramo, ce je metoda ze 'konvergirala'.
	if(norm(tocki0-tocki) < tol)
		break;
	endif
    if(rows(t0)==1 && rows(s0)==1)
    %krivulje
    t0=tocki(1);
    s0=tocki(2);
  elseif(rows(t0)==2 && rows(s0)==2)
  %ploskve
    t0=tocki(1:2);
    s0=tocki(3:4);
   endif
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
%! [d,t,s]=newtonova( @kroznica1, @kroznica2 ,tZac , sZac);
%! assert( d, 5, eps  );
%! assert( t , pi/2, sqrt(eps) );
%! assert( s, -pi/2, sqrt(eps) );
