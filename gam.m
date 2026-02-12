function y=gam(T,t,s)
if s>0
 y=(bet(T,t,s).^2-1)./(s.*(2+s.*t));
else if s<0
 y=(bet(T,t,s).^2-1)./(s.*(2-s.*t));
 else 
 y=0;
 end
end
