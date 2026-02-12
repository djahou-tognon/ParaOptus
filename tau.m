function y=tau(T,t1,t2,s)
  if s>0
    y=bet(T,t1,s)-gam(T,t1,s).*(bet(T,t1,s)-bet(T,t2,s))./(gam(T,t1,s)-gam(T,t2,s));
  else if s<0
    y=bet(T,t1,s)+gam(T,t1,s).*(bet(T,t1,s)-bet(T,t2,s))./abs(gam(T,t1,s)-gam(T,t2,s));
 else
     y=0;
end
end
