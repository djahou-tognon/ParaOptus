function PararoptControl20sp()
size_marker=8;warning off
L		= 10;
x0		= 20;
xtarget 	= 2;
y0		= 10;
ytarget 	=1 ;
	
N0		= 2^10;factorial(8)*3;.5*1e2;
N		= N0/L;


fprintf(2,'\t -- Computation of the solution  --\n')
compt		= 0;
	
alpha		= .05;
T_court		= 1/3;
T		= T_court;1;
t		= linspace(0,T,N0+1);
dt		= t(2)-t(1);
%Fig		= 1;
Newton		= 1;


Sol		= [linspace(x0,xtarget,L+1)' ; ...%x(1:N0/L:end)';...%
		   linspace(y0,ytarget,L+1)' ; ...%y(1:N0/L:end)';...%
		   ones(L,1);...
		   ones(L,1)];%zeros(4*L+2,1);

Tol_Newt	= 1e-13;
Err_Newt	= 1 ;

%initialisation
x_temp		= linspace(x0,xtarget,L*N+1);
y_temp		= linspace(y0,ytarget,L*N+1);
tabXl		= reshape(x_temp(2:end),N,L);
tabYl		= reshape(y_temp(2:end),N,L);
tabLaXl		= ones(N,L);zeros(N,L);
tabLaYl		= ones(N,L);zeros(N,L); 

ratio		= 10.^(-0);%[log10(N)]);%Liste_ratio		= 10.^(-[0:log10(N)]);
Napprox		= N*ratio;

fprintf(2,'\t\t ============== \n Ratio=%e|L=1e%i|N0=1e%i|N=1e%i|Napprox=1e%i\n\t\t ==============\n\n',ratio,log10(L),log10(N0),log10(N),log10(Napprox))
iter		= 0;
tic;
		while iter<20%Err_Newt>Tol_Newt
			iter			= iter + 1 ;
			X			= Sol(1:L+1);
			Y			= Sol(L+2:2*L+2);
			LaX			= Sol(2*L+3:3*L+2);
			LaY			= Sol(3*L+3:4*L+2);
			if L==1
				Compute_grad=1;
				[F,G,~,~,~,~,dFapprox,dGapprox]=pb_NL2(N,L,alpha,x0,xtarget,X,LaX,tabXl,tabLaXl,...
                                                         		 	        	 y0,ytarget,Y,LaY,tabYl,tabLaYl,dt,Newton,Cgrad);
			else
				fprintf(2,'\t\t *********** Pbs Fins *********** \n ')
				Compute_grad=0;
				[F,G,tabXl,tabLaXl,tabYl,tabLaYl]=pb_NL2(N,L,alpha,x0,xtarget,X,LaX,tabXl,tabLaXl,...
                                                        		 					 y0,ytarget,Y,LaY,tabYl,tabLaYl,dt,Newton,Compute_grad);
				fprintf(2,'\t\t ***** Pbs Grossiers ***** \n ')
				Compute_grad=1;
				[~,~,~,~,~,~,dFapprox,dGapprox]=pb_NL2(Napprox,L,alpha,x0,xtarget,X,LaX,tabXl(1:N/Napprox:end,:),tabLaXl(1:N/Napprox:end,:),...
                                                     	     			           y0,ytarget,Y,LaY,tabYl(1:N/Napprox:end,:),tabLaYl(1:N/Napprox:end,:),dt/ratio,Newton,Compute_grad);
			end%if
			dSol			= -[dFapprox;dGapprox]\[F;G];
			Sol			= Sol + dSol;
			Err_Newt		= norm(dSol);
			cost			= (tabXl(end,end)-xtarget)^2+(tabYl(end,end)-ytarget)^2+1/alpha*sum(tabLaXl(:).^2+tabLaYl(:).^2)*dt;
					h = figure(81);
				      	for l=1:L
		  				plot(...%[x((l-1)/L*N0+1:l/L*N0)], [y((l-1)/L*N0+1:l/L*N0)], 'b--',...
						     tabXl(:,l),tabYl(:,l),'k',x0,y0,'g.',xtarget,ytarget,'r.','Markersize',size_marker); ...

		  				hold on; 
	      			      	end;
	     			      	hold off;
	     				xlabel('Prey','interpreter','latex');
	     			   	ylabel('Predator','interpreter','latex');
						drawnow
			fprintf(2,'\t Ratio=%e|Iter=%i|Err_Newt=%e|J=%f\n',ratio,iter,Err_Newt,cost)
		end%while
		toc
		X=Sol;
		save('solutionX.mat','Sol');
end%function

%-------------------------------------------------------------------------------------------
function      [F,G,tabXl,tabLaXl,tabYl,tabLaYl,dF,dG]=pb_NL2(N,L,alpha,x0,xtarget,X,LaX,tabXl,tabLaXl,...
                              		 					  y0,ytarget,Y,LaY,tabYl,tabLaYl,dt,Newton,Compute_grad)

a = 10; b = 0.2; c = b; d = a;%Felix
f	=@(x,y)( a*x                      - b*x.*y ); ...
dfx	=@(x,y)( a*ones(size(x))          - b*y    );%
dfy	=@(x,y)( 0*x          		  - b*x    );%
dfxx	=@(x,y)( 0*x            		   );%
dfyy	=@(x,y)( 0*x                               );%
dfxy	=@(x,y)( 	       - Newton*b*ones(size(x)));%

g	=@(x,y)( c*x.*y		   	  - d*y    );
dgx	=@(x,y)( c*y		    	  + 0*y    );
dgy	=@(x,y)( c*x 	       - d*ones(size(y))   );
dgxx	=@(x,y)( 0*y		    	           );
dgyy	=@(x,y)( 0*x 	       		  + 0*y    );
dgxy	=@(x,y)( Newton*c*ones(size(y))    	           );

dF	=zeros(2*L+1,4*L+2);
dG	=zeros(2*L+1,4*L+2 );

F	=zeros(2*L+1,1);
G	=zeros(2*L+1,1);

F(1,1)  = X(1,1)-x0;
G(1,1)  = Y(1,1)-y0;
compt_par_FG1	=zeros(L,1);
compt_par_FG2	=zeros(L,1);
compt_par_dFG2l	=zeros(L,1);
compt_par_dFG2c	=zeros(L,4);

%parfor l=1:L %l=number of the sub-interval, starting form l=1. Y, La = vectors of initial conditions for Yl and Lal. La(:,1) not used.
parfor l=1:L

	%Actual unknowns on the subinterval:
	Xl  =tabXl(:,l);% 0.5*ones(N,1) ; %represents (y2     ,...,yN+1   )^T on subinterval l
	Yl  =tabYl(:,l);% 0.5*ones(N,1) ; %represents (y2     ,...,yN+1   )^T on subinterval l
	LaXl =tabLaXl(:,l);% 0.5*ones(N,1) ; %represents (lambda1,...,lambdaN)^T
	LaYl =tabLaYl(:,l);% 0.5*ones(N,1) ; %represents (lambda1,...,lambdaN)^T
	Err_Newt = 1; inf ;
	scomp    = 0;
	while Err_Newt>1e-10 
		%Newton system
		PXl	= Xl    - [X(l);Xl(1:end-1,1)]  - dt * (f(Xl,Yl)-1/alpha*[LaXl(2:end,1);LaX(l)]);
		PYl	= Yl    - [Y(l);Yl(1:end-1,1)]  - dt * (g(Xl,Yl)-1/alpha*[LaYl(2:end,1);LaY(l)]);

		dPXlX	=   speye(N) -   spdiags([ones(N-1,1);1],-1,N,N) - dt * spdiags(dfx(Xl,Yl),0,N,N) ;
		dPXlY	= 0*speye(N) - 0*spdiags([ones(N-1,1);1],-1,N,N) - dt * spdiags(dfy(Xl,Yl),0,N,N) ;
		dPYlX	= 0*speye(N) - 0*spdiags([ones(N-1,1);1],-1,N,N) - dt * spdiags(dgx(Xl,Yl),0,N,N) ;
		dPYlY	=   speye(N) -   spdiags([ones(N-1,1);1],-1,N,N) - dt * spdiags(dgy(Xl,Yl),0,N,N) ;

		dPXlLaX	=   dt * 1/alpha * spdiags([1;ones(N-1,1)],1,N,N) ; 
		dPXlLaY	= 0*dt * 1/alpha * spdiags([1;ones(N-1,1)],1,N,N) ; 
		dPYlLaX	= 0*dt * 1/alpha * spdiags([1;ones(N-1,1)],1,N,N) ; 
		dPYlLaY	=   dt * 1/alpha * spdiags([1;ones(N-1,1)],1,N,N) ; 

		QXl	= [LaXl(2:end,1);LaX(l)] - LaXl    + dt * dfx([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]).*LaXl...
							   + dt * dgx([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]).*LaYl;
		QYl	= [LaYl(2:end,1);LaY(l)] - LaYl    + dt * dgy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]).*LaYl...
							   + dt * dfy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]).*LaXl;

		dQXlX	=  dt*spdiags(dfxx([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),-1,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N)...
			  +dt*spdiags(dgxx([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),-1,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N);
		dQXlY	=  dt*spdiags(dfxy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),-1,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N)...
			  +dt*spdiags(dgxy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),-1,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N);
		dQYlX	=  dt*spdiags(dgxy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),-1,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N)...
			  +dt*spdiags(dfxy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),-1,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N);
		dQYlY	=  dt*spdiags(dgyy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),-1,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N)...
			  +dt*spdiags(dfyy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),-1,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N);

		dQXlLaX	=    spdiags([1;ones(N-1,1)],1,N,N) - speye(N) + spdiags(dt*dfx([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),0,N,N);
		dQXlLaY	=    						 spdiags(dt*dgx([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),0,N,N);
		dQYlLaX	=   						 spdiags(dt*dfy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),0,N,N);
		dQYlLaY	=    spdiags([1;ones(N-1,1)],1,N,N) - speye(N) + spdiags(dt*dgy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)]),0,N,N);
 
		%Newton update

		dS  	= -[dPXlX , dPXlY , dPXlLaX, dPXlLaY ;...
			    dPYlX , dPYlY , dPYlLaX, dPYlLaY ;...
			    dQXlX , dQXlY , dQXlLaX, dQXlLaY ;...
			    dQYlX , dQYlY , dQYlLaX, dQYlLaY ]\[PXl;PYl;QXl;QYl];

		Xl  	=  dS(  1:N       , 1) + Xl ;
		Yl  	=  dS(  N+1:2*N   , 1) + Yl ;
		LaXl 	=  dS(2*N+1:3*N   , 1) + LaXl;
		LaYl 	=  dS(3*N+1:4*N   , 1) + LaYl;%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Err_Newt= norm(dS);
			fprintf(2,'\t\t\t Sub-Newton proc. l=%i|sub-iter=%i|Err_Newt=%e|dt=%f|N.L.dt=%f\n',l,scomp,Err_Newt,dt,N*L*dt);
			scomp    = 1+scomp;
	end
	tabXl  (:,l) 	= Xl;
	tabYl  (:,l) 	= Yl;
	tabLaXl(:,l)	= LaXl;
	tabLaYl(:,l)	= LaYl;

%fprintf(2,'\t\t Assemblage F, l=%i\n',l)
	compt_par_FG1(l)=[l+1];
	F_par1(l) = X(l+1) - Xl(end,1);		%This is the equation on the last  time step of the sub-int l , for Xl ;
	G_par1(l) = Y(l+1) - Yl(end,1);		%This is the equation on the last  time step of the sub-int l , for Yl ;
	if l>1
		compt_par_FG2(l)=[L+l];
		F_par2(l)=LaX(l-1) - LaXl(1,1);	%This the equation on the first time step of the sub-int l (with l>1), for LaXl; La(:,1) not used.
		G_par2(l)=LaY(l-1) - LaYl(1,1);	%This the equation on the first time step of the sub-int l (with l>1), for LaYl; La(:,1) not used.
	end


%fprintf(2,'\t\t Fin Assemblage F, l=%i, Calcul inverse fonc. impl.\n',l)
if Compute_grad
	Z	 = -([dPXlX , dPXlY , dPXlLaX, dPXlLaY ;...
   		      dPYlX , dPYlY , dPYlLaX, dPYlLaY ;...
		      dQXlX , dQXlY , dQXlLaX, dQXlLaY ;...
		      dQYlX , dQYlY , dQYlLaX, dQYlLaY ])...
				\...
			[ [-1;zeros(N-1,1)] ,    zeros(N  ,1)   ,  dt/alpha*[zeros(N-1,1);1],  zeros(N,1)               ;...
			      zeros(N  ,1)  ,[-1;zeros(N-1,1)]  ,            zeros(N  ,1)   ,  dt/alpha*[zeros(N-1,1);1];...
dt*[dfxx(X(l),Y(l))*LaXl(1);zeros(N-1,1)]+dt*[dgxx(X(l),Y(l))*LaYl(1);zeros(N-1,1)],...
dt*[dfxy(X(l),Y(l))*LaXl(1);zeros(N-1,1)]+dt*[dgxy(X(l),Y(l))*LaYl(1);zeros(N-1,1)],       [zeros(N-1,1);1] ,   zeros(N  ,1);...
dt*[dgxy(X(l),Y(l))*LaYl(1);zeros(N-1,1)]+dt*[dfxy(X(l),Y(l))*LaXl(1);zeros(N-1,1)],...
dt*[dgyy(X(l),Y(l))*LaYl(1);zeros(N-1,1)]+dt*[dfyy(X(l),Y(l))*LaXl(1);zeros(N-1,1)],        zeros(N  ,1)    ,  [zeros(N-1,1);1];...
];
		% d(P,Q)/d(Xl,Yl,LaX,LaY) * d(Xl,Yl,LaXl,LaYl)/d(X,Y,LaX,LaY) + d(P,Q)/d(X,Y,LaX,LaY) = 0 And Z= d(Xl,Yl,LaXl,LaYl)/d(X,Y,LaX,LaY)

%fprintf(2,'\t\t Extract derivatives.\n')
        dxdX   = Z(N    ,1); dxdY   = Z(N    ,2); dxdLaX   = Z(N    ,3); dxdLaY   = Z(N    ,4);
        dydX   = Z(2*N  ,1); dydY   = Z(2*N  ,2); dydLaX   = Z(2*N  ,3); dydLaY   = Z(2*N  ,4);
        dlaxdX = Z(2*N+1,1); dlaxdY = Z(2*N+1,2); dlaxdLaX = Z(2*N+1,3); dlaxdLaY = Z(2*N+1,4);
        dlaydX = Z(3*N+1,1); dlaydY = Z(3*N+1,2); dlaydLaX = Z(3*N+1,3); dlaydLaY = Z(3*N+1,4);
%ZZ(16*(l-1)+1:16*l)=[dxdX  ;dydX  ;dlaxdX  ;dlaydX;...
%		     dxdY  ;dydY  ;dlaxdY  ;dlaydY;...
%		     dxdLaX;dydLaX;dlaxdLaX;dlaydLaX;...
%		     dxdLaY;dydLaY;dlaxdLaY;dlaydLaY;];

%fprintf(2,'\t\t Assemblage dF, l=%i\n',l)
%(dX(l+1)- (A*  dxdX*dX(l) + B*  dxdY*dY(l) + C*  dxdLaX*dLaX(l) + D*  dxdLaY*dLaY(l))
	%dF(l+1 ,       l)	= -dxdX  ;
	%dF(l+1 ,   L+1+l)	= -dxdY  ;%verifier la taille de Y et La !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	%dF(l+1 , 2*L+2+l)	= -dxdLaX;%verifier la taille de Y et La !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	%dF(l+1 , 3*L+2+l)	= -dxdLaY;%verifier la taille de Y et La !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	dF_par1(l,:)=[-dxdX -dxdY -dxdLaX -dxdLaY];
	compt_par_dFG1l(l,:)=l+1;
	compt_par_dFG1c(l,:)=[l L+1+l 2*L+2+l 3*L+2+l];

	if l>1
		%dF(L+l ,       l)	= -dlaxdX ;
		%dF(L+l ,   L+1+l)	= -dlaxdY ;
		%dF(L+l , 2*L+2+l)	= -dlaxdLaX;
		%dF(L+l , 3*L+2+l)	= -dlaxdLaY;
		dF_par2(l,:)=[-dlaxdX -dlaxdY -dlaxdLaX -dlaxdLaY];
		compt_par_dFG2l(l,:)=L+l ;
		compt_par_dFG2c(l,:)=[l L+1+l 2*L+2+l 3*L+2+l];
	end
	%dG(l+1 ,       l)	= -dydX  ;
	%dG(l+1 ,   L+1+l)	= -dydY  ;%verifier la taille de Y et La !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	%dG(l+1 , 2*L+2+l)	= -dydLaX;%verifier la taille de Y et La !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	%dG(l+1 , 3*L+2+l)	= -dydLaY;%verifier la taille de Y et La !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	dG_par1(l,:)=[-dydX -dydY -dydLaX -dydLaY];

	if l>1
		%dG(L+l , l      )	= -dlaydX ;
		%dG(L+l , L+l+1  )	= -dlaydY ;
		%dG(L+l , 2*L+l+2)	= -dlaydLaX;
		%dG(L+l , 3*L+l+2)	= -dlaydLaY;
		dG_par2(l,:)=[-dlaydX -dlaydY -dlaydLaX -dlaydLaY];
	end
end%if nargout
end%parfor

	%F(l+1,1) = X(l+1) - Xl(end,1);		%This is the equation on the last  time step of the sub-int l , for Yl ;
	%F(L+l,1)=LaX(l-1) - LaXl(1,1);		%This the equation on the first time step of the sub-int l (with l>1), for Lal; La(:,1) not used.


for l=1:L
	F(compt_par_FG1(l))=F_par1(l);		%This is the equation on the last  time step of the sub-int l , for Xl ;
	G(compt_par_FG1(l))=G_par1(l);		%This is the equation on the last  time step of the sub-int l , for Xl ;
	if l>1
		F(compt_par_FG2(l))=F_par2(l);		%This is the equation on the last  time step of the sub-int l , for Xl ;
		G(compt_par_FG2(l))=G_par2(l);		%This is the equation on the last  time step of the sub-int l , for Xl ;
	end
end

F(2*L+1,1)	= LaX(L) - X(L+1)+xtarget;
G(2*L+1,1)	= LaY(L) - Y(L+1)+ytarget;

if Compute_grad
	for l=1:L
		dF(compt_par_dFG1l(l),compt_par_dFG1c(l,:))=dF_par1(l,:)';
		dG(compt_par_dFG1l(l),compt_par_dFG1c(l,:))=dG_par1(l,:)';
		if l>1
			dF(compt_par_dFG2l(l),compt_par_dFG2c(l,:))=dF_par2(l,:)';
			dG(compt_par_dFG2l(l),compt_par_dFG2c(l,:))=dG_par2(l,:)';
		end
	end

	dF(2*L+1,L+1)	=        -1;
	dF		= dF + [  eye(L+1)  , zeros(L+1)   , zeros(L+1,L) , zeros(L+1,L) ;...
	 			zeros(L,L+1), zeros(L,L+1) ,   eye(L)     , zeros(L,L) ];%%%Problematique...
	dG(2*L+1,2*L+2)	=        -1;
	dG		= dG + [ zeros(L+1) , eye(L+1)     , zeros(L+1,L) , zeros(L+1,L) ;...
	 			zeros(L,L+1), zeros(L,L+1) , zeros(L,L)   , eye(L)];%%%Problematique...
end

%fprintf(2,'\t\t Fin assemblage F et dF\n')

end
