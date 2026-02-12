function y=bet(T,t,s)
if s>0
  y=(1+s.*t).^(T./t);
else if s<0
 y=(1-s.*t).^(-T./t);
 else
  y=0;
end
end
