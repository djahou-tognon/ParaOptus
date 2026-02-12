function [Sol,err]=pendulum_test(ratio,ll,ExactSol)
size_marker=8;warning off;
L		=8;
x0		= pi/4;
xtarget 	= 0.0;
y0		= pi/6;
ytarget 	=0.0;
    %%%%%%%%%%%%%%%%%%%%%
N0		=2e5;factorial(8)*3;.5*1e2;
N		= N0/L;
fprintf(2,'\t -- Computation of the solution  --\n')
compt		= 0;
alpha		= 0.01;
T_court		=1;
T		= T_court;1/3;
t		= linspace(0,T,N0+1);
dt		= t(2)-t(1);
%Fig		= 1;
%Newton		= 1;


Sol		= [linspace(x0,xtarget,L+1)' ; ...%x(1:N0/L:end)';...%
		   linspace(y0,ytarget,L+1)' ; ...%y(1:N0/L:end)';...%
		   ones(L,1);...
		   ones(L,1)];%zeros(4*L+2,1);

Tol_Newt	= 1e-12;
Err_Newt	= 1 ;

%initialisation
x_temp		=[x0,ones(1,L*N-1),xtarget];
y_temp		=[y0,ones(1,L*N-1),ytarget];
tabXl		= reshape(x_temp(2:end),N,L);
tabYl		= reshape(y_temp(2:end),N,L);
tabLaXl		= zeros(N,L);ones(N,L);
tabLaYl		= zeros(N,L);ones(N,L);
%ratio		= 10.^(-r);
Napprox		= N*ratio;
fprintf(2,'\t\t ============== \n Ratio=%e|L=1e%i|N0=1e%i|N=1e%i|Napprox=1e%i\n\t\t ==============\n\n',ratio,log10(L),log10(N0),log10(N),log10(Napprox))
iter		= 0;
err=[];
		while Err_Newt>Tol_Newt
			iter			= iter + 1 ;
			X			= Sol(1:L+1);
			Y			= Sol(L+2:2*L+2);
			LaX			= Sol(2*L+3:3*L+2);
			LaY			= Sol(3*L+3:4*L+2);
			if L==1
				Compute_grad=1;
				[F,G,~,~,~,~,dFapprox,dGapprox]=pb_NL2(N,L,alpha,x0,xtarget,X,LaX,tabXl,tabLaXl,...
                                                         		 	        	 y0,ytarget,Y,LaY,tabYl,tabLaYl,dt,Compute_grad,ll);
			else
				fprintf(2,'\t\t *********** Pbs Fins *********** \n ')
				Compute_grad=0;
				
				[F,G,tabXl,tabLaXl,tabYl,tabLaYl]=pb_NL2(N,L,alpha,x0,xtarget,X,LaX,tabXl,tabLaXl,...
                                                        		 					 y0,ytarget,Y,LaY,tabYl,tabLaYl,dt,Compute_grad,ll);
				fprintf(2,'\t\t ***** Pbs Grossiers ***** \n ')
				Compute_grad=1;
				[~,~,~,~,~,~,dFapprox,dGapprox]=pb_NL2(Napprox,L,alpha,x0,xtarget,X,LaX,tabXl(1:N/Napprox:end,:),tabLaXl(1:N/Napprox:end,:),...
                                                     	     			           y0,ytarget,Y,LaY,tabYl(1:N/Napprox:end,:),tabLaYl(1:N/Napprox:end,:),dt/ratio,Compute_grad,ll);
			end%if
			dSol			= -[dFapprox;dGapprox]\[F;G];
			Sol			= Sol + dSol;
			Err_Newt		= norm(dSol);
           tabx=tabXl;
            taby=tabYl;
			cost= (tabXl(end,end)-xtarget)^2+(tabYl(end,end)-ytarget)^2+1/alpha*sum(tabLaXl(:).^2+tabLaYl(:).^2)*dt;
					h = figure(81);
			for l=1:L
                %plot(tabXl(:,l),tabYl(:,l),'b',x0,y0,'g.',xtarget,ytarget,'r.','Markersize',size_marker);
                hold on;
            end
            
	     			      	hold off;
	     				xlabel('Position','interpreter','latex');
	     		   	ylabel('Speed','interpreter','latex');
				%		drawnow
			fprintf(2,'\t Ratio=%e|Iter=%i|Err_Newt=%e|J=%f\n',ratio,iter,Err_Newt,cost)
           %err(iter)=norm(SoL-Sol);
           if length(ExactSol)>1
                err(iter)=norm(Sol-ExactSol,2);
                it=iter;
           end
          %Cost(iter)=cost;
      %z=Sol;
		end%while

end%function

%----------------------------- --------------------------------------------------------------
function      [F,G,tabXl,tabLaXl,tabYl,tabLaYl,dF,dG]=pb_NL2(N,L,alpha,x0,xtarget,X,LaX,tabXl,tabLaXl,...
                              		 					  y0,ytarget,Y,LaY,tabYl,tabLaYl,dt,Compute_grad,ll)
gg=9.81;
omega=gg/ll;

f	=@(x,y,Lax,Lay) ( y ); ...
dfx	=@(x,y,Lax,Lay)( 0*x);%
dfy	=@(x,y,Lax,Lay)(ones(size(y)) );%
%dfxx	=@(x,y)( 0*x           		   );%
%dfyy	=@(x,y)( 0*y                              );%
%dfxy	=@(x,y)( 	     0*x);%
dfxLa	=@(x,y,Lax,Lay) (0*y);
dfyLa  =@(x,y,Lax,Lay) (0*y);
%%%%%%%%

g	=@(x,y,Lax,Lay)( omega*sin(x)-Lay.*(cos(x).^2)/alpha);
dgx	=@(x,y,Lax,Lay)(  omega*cos(x)+Lay.*sin(2*x)/alpha);
dgy	=@(x,y,Lax,Lay)( 0*x  );
dgxLa	=@(x,y,Lax,Lay)(0*Lax		    	           );
dgyLa	=@(x,y,Lax,Lay)(-(cos(x).^2)/alpha);
%dgxy	=@(x,y)( 0*x   );
%%%%
%%%%%%%%%%%%
f1	=@(x,y,Lax,Lay) (0*x);
%f3	=@(x,y,Lax,Lay)(0*y);
df1x	=@(x,y,Lax,Lay) (0*x);
df1y	=@(x,y,Lax,Lay) (0*y);
df1Lax	=@(x,y,Lax,Lay) (0*Lax);
df1Lay	=@(x,y,Lax,Lay) (0*Lay);

f3	=@(x,y,Lax,Lay) (omega*cos(x)+Lay.*sin(2*x)/(2*alpha));
df3x	=@(x,y,Lax,Lay) (-omega*sin(x)+Lay.*cos(2*x)/alpha);
df3y	=@(x,y,Lax,Lay) (0*y);
df3Lax	=@(x,y,Lax,Lay) (0*Lax);
df3Lay	=@(x,y,Lax,Lay) (sin(2*x)/(2*alpha));

f2	=@(x,y,Lax,Lay) (ones(size(x)));
df2x	=@(x,y,Lax,Lay) (0*x);
df2y	=@(x,y,Lax,Lay) (0*y);
df2Lax	=@(x,y,Lax,Lay) (0*Lax);
df2Lay	=@(x,y,Lax,Lay) (0*Lay);

f4	=@(x,y,Lax,Lay) (0*x);
df4x	=@(x,y,Lax,Lay) (0*x);
df4y	=@(x,y,Lax,Lay) (0*y);
df4Lax	=@(x,y,Lax,Lay) (0*Lax);
df4Lay	=@(x,y,Lax,Lay) (0*Lay);


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

for l=1:L%
%for l=1:L
	%Actual unknowns on the subinterval:
	Xl  =tabXl(:,l);% 0.5*ones(N,1) ; %represents (y2     ,...,yN+1   )^T on subinterval l
	Yl  =tabYl(:,l);% 0.5*ones(N,1) ; %represents (y2     ,...,yN+1   )^T on subinterval l
	LaXl =tabLaXl(:,l);% 0.5*ones(N,1) ; %represents (lambda1,...,lambdaN)^T
	LaYl =tabLaYl(:,l);% 0.5*ones(N,1) ; %represents (lambda1,...,lambdaN)^T
	Err_Newt = 1; inf ;
	scomp    = 0;
	while Err_Newt>1e-10
		%Newton system
		%PXl	= Xl    - [X(l);Xl(1:end-1,1)]  - dt * (f(Xl,Yl)-1/alpha*[LaXl(2:end,1);LaX(l)]);
		%PYl	= Yl    - [Y(l );Yl(1:end-1,1)]  - dt * (g(Xl,Yl)-1/alpha*[LaYl(2:end,1);LaY(l)]);

		PXl	= Xl- [X(l);Xl(1:end-1,1)]  - dt*f([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],[LaXl(2:end,1);LaX(l)],[LaYl(2:end,1);LaY(l)]);
		PYl	= Yl- [Y(l);Yl(1:end-1,1)]  - dt*g([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],[LaXl(2:end,1);LaX(l)],[LaYl(2:end,1);LaY(l)]);



		dPXlX	=   speye(N) -   spdiags([ones(N-1,1);1],-1,N,N) - dt * spdiags(dfx([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],[LaXl(2:end,1);LaX(l)],[LaYl(2:end,1);LaY(l)]),-1,N,N) ;
		dPXlY	= 0*speye(N) - 0*spdiags([ones(N-1,1);1],-1,N,N) - dt * spdiags(dfy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],[LaXl(2:end,1);LaX(l)],[LaYl(2:end,1);LaY(l)]),-1,N,N) ;
		dPYlX	= 0*speye(N) - 0*spdiags([ones(N-1,1);1],-1,N,N) - dt * spdiags(dgx([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],[LaXl(2:end,1);LaX(l)],[LaYl(2:end,1);LaY(l)]),-1,N,N) ;
		dPYlY	=   speye(N) -   spdiags([ones(N-1,1);1],-1,N,N) - dt * spdiags(dgy([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],[LaXl(2:end,1);LaX(l)],[LaYl(2:end,1);LaY(l)]),-1,N,N) ;

	
		DfxLa=dfxLa([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],[LaXl(2:end,1);LaX(l)],[LaYl(2:end,1);LaY(l)]);
		DfyLa=dfyLa([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],[LaXl(2:end,1);LaX(l)],[LaYl(2:end,1);LaY(l)]);
		DgxLa=dgxLa([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],[LaXl(2:end,1);LaX(l)],[LaYl(2:end,1);LaY(l)]);
		DgyLa=dgyLa([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],[LaXl(2:end,1);LaX(l)],[LaYl(2:end,1);LaY(l)]);
		dPXlLaX	=   0*speye(N) -   0*spdiags([ones(N-1,1);1],1,N,N) - dt * spdiags([0;DfxLa(2:end)],0,N,N) ;
		dPXlLaY	= 0*speye(N) - 0*spdiags([ones(N-1,1);1],1,N,N) - dt * spdiags([0;DfyLa(2:end)],0,N,N) ;
		dPYlLaX	= 0*speye(N) - 0*spdiags([ones(N-1,1);1],1,N,N) - dt * spdiags([0;DgxLa(2:end)],0,N,N) ;
		dPYlLaY	=   0*speye(N) -   0*spdiags([ones(N-1,1);1],1,N,N) - dt * spdiags([0;DgyLa(2:end)],0,N,N) ;


		QXl	= [LaXl(2:end,1);LaX(l)] - LaXl    + dt * f1([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl).*LaXl...
							   + dt * f3([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl).*LaYl;
		QYl	= [LaYl(2:end,1);LaY(l)] - LaYl    + dt * f4([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl).*LaYl...
							   + dt * f2([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl).*LaXl;



		dQXlX	=  dt*spdiags(df1x([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),-1,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N)...
			  +dt*spdiags(df3x([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),-1,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N);
		dQXlY	=  dt*spdiags(df1y([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),-1,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N)...
			  +dt*spdiags(df3y([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),-1,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N);
			  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		dQYlX	=  dt*spdiags(df4x([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),-1,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N)...
			  +dt*spdiags(df2x([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),-1,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N);
		dQYlY	=  dt*spdiags(df4y([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),-1,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N)...
			  +dt*spdiags(df2y([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),-1,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N);



		dQXlLaX	= spdiags([1;ones(N-1,1)],1,N,N) - speye(N)...
		           		 + dt*spdiags(df1Lax([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N)...
		          		 +dt*spdiags(f1([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)...
			  		 +dt*spdiags(df3Lax([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N);

		dQXlLaY	=  dt*spdiags(df1Lay([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N)...
					+dt*spdiags(f3([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)...
			  		+dt*spdiags(df3Lay([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N);
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		dQYlLaX	=  dt*spdiags(df4Lax([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N)...
					+dt*spdiags(f2([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)...
			 		 +dt*spdiags(df2Lax([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N);

		dQYlLaY	=  spdiags([1;ones(N-1,1)],1,N,N) - speye(N)+ dt*spdiags(df4Lay([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)*spdiags([LaYl(2:end,1); LaY(l) ],0,N,N)...
				+dt*spdiags(f4([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)...
			  	+dt*spdiags(df2Lay([X(l);Xl(1:end-1,1)],[Y(l);Yl(1:end-1,1)],LaXl,LaYl),0,N,N)*spdiags([LaXl(2:end,1); LaX(l) ],0,N,N);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z	 = -([dPXlX , dPXlY , dPXlLaX, dPXlLaY ;...
   	      dPYlX , dPYlY , dPYlLaX, dPYlLaY ;...
	      dQXlX , dQXlY , dQXlLaX, dQXlLaY ;...
	      dQYlX , dQYlY , dQYlLaX, dQYlLaY ])...
				\...
			[ [-1;zeros(N-1,1)]+dt*[dfx(X(l),Y(l),LaXl(l),LaYl(l));zeros(N-1,1)] ,    zeros(N  ,1)+dt*[dfy(X(l),Y(l),LaXl(l),LaYl(l));zeros(N-1,1)]  ,  dt*[zeros(N-1,1);dfxLa(X(l),Y(l),LaXl(l),LaYl(l))],  dt*[zeros(N-1,1);dfyLa(X(l),Y(l),LaXl(l),LaYl(l))]              ;...
			      zeros(N  ,1)+dt*[dgx(X(l),Y(l),LaXl(l),LaYl(l));zeros(N-1,1)]  ,[-1;zeros(N-1,1)]+dt*[dgy(X(l),Y(l),LaXl(l),LaYl(l));zeros(N-1,1)] ,           dt*[zeros(N-1,1);dgxLa(X(l),Y(l),LaXl(l),LaYl(l))]   , dt*[zeros(N-1,1);dgyLa(X(l),Y(l),LaXl(l),LaYl(l))] ;...
dt*[df1x(X(l),Y(l),LaXl(l),LaYl(l))*LaXl(l);zeros(N-1,1)]+dt*[df3x(X(l),Y(l),LaXl(l),LaYl(l))*LaYl(l);zeros(N-1,1)],...
dt*[df1y(X(l),Y(l),LaXl(l),LaYl(l))*LaXl(l);zeros(N-1,1)]+dt*[df3y(X(l),Y(l),LaXl(l),LaYl(l))*LaYl(l);zeros(N-1,1)],       [zeros(N-1,1);1] ,   zeros(N  ,1);...
dt*[df4x(X(l),Y(l),LaXl(l),LaYl(l))*LaYl(l);zeros(N-1,1)]+dt*[df2x(X(l),Y(l),LaXl(l),LaYl(l))*LaXl(l);zeros(N-1,1)],...
dt*[df4y(X(l),Y(l),LaXl(l),LaYl(l))*LaYl(l);zeros(N-1,1)]+dt*[df2y(X(l),Y(l),LaXl(l),LaYl(l))*LaXl(l);zeros(N-1,1)],        zeros(N  ,1)    ,  [zeros(N-1,1);1];...
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% d(P,Q)/d(Xl,Yl,LaX,LaY) * d(Xl,Yl,LaX,LaY)/d(X,Y,LaX,LaY) + d(P,Q)/d(X,Y,LaX,LaY) = 0 And Z= d(Xl,Yl,LaX,LaY)/d(X,Y,LaX,LaY)

%fprintf(2,'\t\t Extract derivatives.\n')
        dxdX   = Z(N    ,1); dxdY   = Z(N    ,2); dxdLaX   = Z(N    ,3); dxdLaY   = Z(N    ,4);
        dydX   = Z(2*N  ,1); dydY   = Z(2*N  ,2); dydLaX   = Z(2*N  ,3); dydLaY   = Z(2*N  ,4);
        dlaxdX = Z(2*N+1,1); dlaxdY = Z(2*N+1,2); dlaxdLaX = Z(2*N+1,3); dlaxdLaY = Z(2*N+1,4);
        dlaydX = Z(3*N+1,1); dlaydY = Z(3*N+1,2); dlaydLaX = Z(3*N+1,3); dlaydLaY = Z(3*N+1,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

