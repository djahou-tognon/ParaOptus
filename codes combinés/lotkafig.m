function lotkafig()
T	        =1;
a	        =07;
b	        =1;
c	        =0.2;
alpha	        =5e-2;
yin	        =[100;10];
ytg	        =[3;1];%[0;2];
L	        =8;
limit	        =14;
figlinear       =24;  
fignonlinear    =52;
ratio           =2.^(-(2:2:10));
%linearlotka(T,a,b,c,alpha,L,limit,yin,ytg,figlinear);

nonlinearlotka(T,a,b,c,alpha,L,limit,ratio,yin,ytg,fignonlinear);
end
