function vizualizacija = visualise( thetaVektor, omegaVektor, l, L)
  % uporabi, ce zelis ( samo zmanjsaj deltaT, drugace bo rabil program prevec casa, ker bo tock prevec )
  clf;
  stTock = size( thetaVektor,2 );
  %ce je stevilo tock preveliko, bom moral deltaCas zmanjsati. 
  plot3( 0,0,0, "kx" );   % oznacim izhodisce.
  hold on
  
  koeficient = 1;
  if stTock > 200
    
    koeficient = ceil(stTock / 200);
  endif
  
  for k = (1:koeficient:stTock)
    x1 =  sin( thetaVektor(1,k) ) * l; % x1 in y1 = pozicija prve gredi z dolzino l ( mali L ).
    y1 = -cos( thetaVektor(1,k) ) * l;
    
    x = x1 + sin( omegaVektor(1,k) ) * L; % x in y = pozicija dejanske masne tocke na koncu nihala ( veliki L ) 
    y = y1 - cos( omegaVektor(1,k) ) * L;
    
    plot3( 0,0,k, "kx" );  % oznacim izhodisce.
    plot3( x1, y1, k, "ro" ); % drugi konec mehanske gredi.
    plot3( x , y,  k, "bo" ); % dejanska masna tocka na koncu nihala.
    
    plot3( [0;x1],[0;y1],[k;k] ); % oznacim se daljice, da dobim lepe 3d - oblike, ki plesejo ??
    plot3( [x;x1],[y;y1],[k;k] );
    hold on
  endfor
 
endfunction