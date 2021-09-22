function [ f, df, ddf ]=kroznica1(t)
  f=[2.*sin(t) - 3; 2.*cos(t); 0];
  df=[ 2.*cos(t); -2.*sin(t); 0];
  ddf= [ -2.*sin(t); -2.*cos(t); 0];
endfunction
