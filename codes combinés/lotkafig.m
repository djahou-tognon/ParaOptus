function lotkafig()
T	        =1;
a	        =10;
b	        =0.2;
c	        =0.2;
d           =10;
alpha	        =5e-1;
yin	        =[10;8];
ytg	        =[1,1];%[0;2];
L	        =8;
limit	        =14;
figlinear       =24;  
fignonlinear    =52;
ratio           =2.^(-(2:2:10));
linearlotka(T,a,b,c,d,alpha,L,limit,yin,ytg,figlinear);

nonlinearlotka(T,a,b,c,d,alpha,L,limit,ratio,yin,ytg,fignonlinear);
end
