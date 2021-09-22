function [ g, dg, ddg ]=kroznica2(s)
  g=[ sin(s)+ 5; cos(s); 0 ];
  dg=[ cos(s) ; -sin(s) ; 0 ];
  ddg= [ -sin(s) ; -cos(s); 0 ];
endfunction
