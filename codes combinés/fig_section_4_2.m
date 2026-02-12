function fig_section_4_2()
T	        =1;
a1	        =10;
a2	        =0.2;
b1	        =0.2;
b2           =10;
alpha	        =5e-1;
yin	        =[20;15];
ytg	        =[1,1];%[0;2];
L	        =8;
tol         =1e-9;
limit	        =14;
figlinear       =79;  
fignonlinear    =52;
ratio           =2.^(-(2:2:10));
linearlotka(T,a1,a2,b1,b2,alpha,L,limit,yin,ytg,tol,figlinear);
nonlinearlotka(T,a1,a2,b1,b2,alpha,L,limit,ratio,yin,ytg,tol,fignonlinear);
end
