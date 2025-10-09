% 01/06/2023
% This function works for differents scheme on Cgrid and Fgird
 
function linearlotka(T,a,b,c,alpha,L,limit,yin,ytg,figlinear)%[errDt,sol]=paraoptNumLotkaEq()
    %%%%% Espace data
    %T=1;
    line=1.6; Fontsize=14;
    marker=4;
    colors=[1,0,0;0,0,1;0,0.6667,0];
    labelsize=12;
    sigma=a;%8;
    %b2=0.2;
%alpha=5e-2;
    %yin=[20;5];
    %ytg=[0;1];%[0;2];
    B=eye(2);%[0;1];%eye(2);
    r=2;
    %Delta=-(r+1)^2*spdiags([-ones(r,1), 2*ones(r,1),-ones(r,1)],-1:1,r,r);
    %%%%Time Data
    %T=1e-1; 
    %L=10;
    DT=T/L;
    M=2^limit;
    calL=[sigma,0;0,-b];%[8,0;0,-1];
    %if scheme==3
    %M=2^10;
    %%% Euler explicit 
    A=0;
    b=1;
    c=0;
    d=1;
    s=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%d
    
    dt=DT/M;
    %%% Linear system %%%%%%%%%%%%%%%%%%%%
     Mdt=matrix(calL,B,r,A,b,s,L,M,dt,alpha,d);
     bf=zeros(r*(2*L+1),1);
     bf(1:2)=yin(:);
     bf(end-1:end)=ytg(:);
     tol = 1e-13;
     restart = 100;
     %[sol, flag, relres, iter] = gmres(Mdt, bf, restart, tol);
     sol=Mdt\bf;
     %disp(relres)
     %iter=0;
     %limit=10;
     k=1:limit-1;
     N=2.^k;
     itermax=15;
     errDt=zeros(limit-1,itermax);
     %raidus=zeros(10,1);
     Dt=DT./N;
     for i=1:limit-1
	     MDt=matrix(calL,B,r,A,b,s,L,N(i),Dt(i),alpha,d);
	     Mat=itermatrix(MDt,Mdt);
	     [x,err]=solveriterative(Mat,MDt\bf,sol,itermax);
	     errDt(i,:)=err;
	end
    %%% Linear system %%%%%%%%%%%%%%%%%%%%
     Mdt=matrix(calL,B,r,A,b,s,L,M,dt,alpha,d);
     raidus=0*k;
     EstimateBound=0*k;
     convFactor=0*k;
     for i=1:limit-1
	     MDt=matrix(calL,B,r,A,b,s,L,N(i),Dt(i),alpha,d);
	     raidus(i)=spectralradius(itermatrix(MDt,Mdt));
	     EstimateBound(i)=sigma*(Dt(i)-dt)*(0.5+...
    (sigma*(Dt(i)-dt)/2+1)*exp(2*sigma*DT)); 
    convFactor(i)=errDt(i,2)/errDt(i,1);     
	end
	figure(figlinear)
   loglog(Dt,EstimateBound, ...
    'Color', [colors(1,:)], ...      
    'LineStyle', '-', ...           
    'LineWidth', line, ...             
    'Marker', 'd', ...              
    'MarkerSize', marker,...
    'DisplayName','Upper Bound'...
    );
	hold on
    loglog(Dt,raidus, ...
    'Color', [colors(2,:)], ...      
    'LineStyle', '-', ...           
    'LineWidth', line, ...             
    'Marker', 'o', ...              
    'MarkerSize', marker,...
    'DisplayName','$\rho$'...
    );

	loglog(Dt,convFactor, ...
    'Color', [colors(3,:)], ...      
    'LineStyle', '--', ...           
    'LineWidth', line, ...             
    'Marker', 's', ...              
    'MarkerSize', marker,...
    'DisplayName','$\hat\rho$'...
    );

	 xlabel('$\Delta t$','interpreter','latex','FontSize',Fontsize);
     ylabel('$\rho$','interpreter','latex','FontSize',Fontsize);
     legend('interpreter','latex','FontWeight','bold','Location','northwest','FontSize',12);
     set(gca,'FontSize',labelsize);
     grid on
     hold off

    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [x,err]=solveriterative(M,b,sol,itermax)
    x=b;
    err=zeros(itermax,1);
	    for i=1:itermax
	    	x=M*x+b;
	    	err(i)=norm(x-sol,2);
	    end
    %err=norm(x-sol,2);
    end
    
    function rho_test(L,B,r,A1,b1,s1,A2,b2,s2,Len,DT,alpha,d1,d2,N_test,M_test,cgrid,fgrid,figname)
    %%%%% Dicretization Parameters Test
    Fontsize=19;
    Size_line=1.5;
    labelsize=16;
    %N_test=1:1:10; %SDIRK
   
    %N_test=2;
    %M_test=3:1:16;
    %cgrid=DT./2.^(N_test);
    %fgrid=DT./2.^M_test;
    %%%%%
    fMatrix=matrix(L,B,r,A1,b1,s1,Len,DT/fgrid,fgrid,alpha,d1);
    %cMatrix=matrix(L,B,r,A,b,s,Len,fgrid,cgrid,alpha,d);
    convergence_factor=[];zeros(1,length(fgrid));
    %A1=0;b1=1;c1=0;d1=1;
        for p=1:length(cgrid)
            cMatrix=matrix(L,B,r,A2,b2,s2,Len,DT/cgrid(p),cgrid(p),alpha,d2);
            %fMatrix=matrix(L,B,r,A,b,s,Len,DT/fgrid(p),fgrid(p),alpha,d);
            convergence_factor(p)=spectral(cMatrix,fMatrix);
        end
    
     regression=polyfit(log10(cgrid),log10(convergence_factor),1);
	reg=polyval(regression,log10(cgrid));
     figure(figname)
     loglog(cgrid,convergence_factor,'-ob',...
                cgrid,10.^reg,'d--r','LineWidth',1.4,'MarkerSize',4);
     grid on
     xlabel('$\Delta t$','interpreter','latex','FontSize',Fontsize);
     ylabel('$\rho$','interpreter','latex','FontSize',Fontsize);
     set(gca,'FontSize',labelsize);

     %figure(12)
     %loglog(cgrid,convergence_factor,'-ob',...
     %           cgrid,(polyval(regression1,cgrid)),'s-r');
     %xlabel("$\log10(\Delta t)$","interpreter","latex")
     %ylabel("$\log10(\rho)$","interpreter","latex")
     %grid on
     %save("SDIRK");
     
     %plot(Lambda,J_ini,Lambda,J_opt(5,:),Lambda,J_NL,'LineWidth',Size_line)
     %xlabel('$$\Delta t$$','interpreter','latex','FontSize',Fontsize);
     %ylabel('$\rho$','interpreter','latex','FontSize',Fontsize);
     %legend({'$$J_\lambda(\gamma_\lambda^\ast,c_\lambda^\ast)$$','$$J_\lambda(\gamma_\lambda^c,c_\lambda^c)$$',...
     %    '$$J_\lambda(\gamma_\lambda^{true},c_\lambda^{true})$$'},...
     % 'interpreter','latex','FontWeight','bold','Location','Southeast','FontSize',Fontsize);
     %a = get(gca,'XTickLabel');
     %set(gca,'TickLabelInterpreter','latex')
     %set(gca,'XTickLabel',a,'FontName','Times','fontsize',Fontsize)
     %xticks(Lambda0(1 : 12:  end));xticklabels(num2str(Lambda0(1 : 12:  end)));
     %grid on
     %figname        = sprintf('J_l');
     %print(h,'-dpng',strcat(figname,'.png'));
     %print(h,'-deps',strcat(figname,'.eps'),'-r3000'); 
     
    %plot(log10(cgrid),log10(convergence_factor),'-ob');
        grid on	
     fprintf(2,'\t\t Regression %.2f, %.2f\n\n',regression(1),regression(2))
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     function [Mat]=matrix(calL,B,r,A,b,s,L,M,t,alpha,d)
     %Computation of the RK operator
     Z=(eye(r*s)-t*kron(A,calL))\kron(b',eye(r));
     W=zeros(r,r);
     W1=zeros(r,r,s);
     Qd=zeros(r,r); % Quadrature operator inside R
        for j=1:s 
                W1(:,:,j)=Z(:,r*(j-1)+1:r*j);
                W=W+W1(:,:,j);
                Qd=Qd+W1(:,:,j)*B*B'*W1(:,:,j)'/d(j);
        end
     RK=eye(r)+t*W*calL; %RK operator
     %% Operator P
     P=RK^M;
     R=zeros(r,r);
         for k=1:M
          R=R+t*RK^(k-1)*Qd*(RK')^(k-1);
         end
      Mat=eye(r*(2*L+1));
          %%Fill the Firt bloc row
          for k=1:L
             Mat(k*r+1:(k+1)*r,(k-1)*r+1:k*r)=-P;
             Mat(k*r+1:(k+1)*r,(L+k)*r+1:(L+k+1)*r)=R/alpha;
          end
          %%Fill the second bloc row
          Mat(2*L*r+1:(2*L+1)*r,L*r+1:(L+1)*r)=-eye(r,r);
          for k=1:L-1
          Mat((L+k)*r+1:(L+1)*r+k*r,(L+1)*r+k*r+1:(L+1)*r+(k+1)*r)=-P';
          end
     end
      
    function Mat=itermatrix(cMat,fMat)
      Ik=eye(length(cMat));
       Mat=Ik-cMat\fMat;
       %rho=max(abs(eig(Mat)));
    end
        function rho=spectralradius(M)
       rho=max(abs(eig(M)));
    end
    
