% 01/06/2023
% This function works for differents scheme on Cgrid and Fgird
 
function linearheateq()
    Fontsize=20;
    Fontsizelegende=14.5;
    Fontsizescale=13;
    %%%%% Espace data
    T=1;
    sigma=16;
    alpha=1;
    yin=1;
    ytg=2;%[0;2];
    B=1;%[0;1];%eye(2);
    r=1;
    %Delta=-(r+1)^2*spdiags([-ones(r,1), 2*ones(r,1),-ones(r,1)],-1:1,r,r);
    %%%%Time Data
    %T=1e-1; 
    L=16;
    DT=T/L;
    calL=-sigma;%[8,0;0,-1];
    %if scheme==3
    M=2^14;
    %%% Euler explicit 
    A=0;
    b=1;
    c=0;
    d=1;
    s=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%d
 
     limit=14;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  M=2^limit;
  dt=DT/M;
  k=2:limit-1;
  N=2.^k;
    %%% Linear system %%%%%%%%%%%%%%%%%%%%
     Mdt=matrix(calL,B,r,A,b,s,L,M,dt,alpha,d);
     radius=0*k;
     EstimateBound=0*k;
     for i=1:limit-2
     %N=5;
	     Dt=DT/N(i);
	     MDt=matrix(calL,B,r,A,b,s,L,N(i),Dt,alpha,d);
	     radius(i)=spectralradius(itermatrix(MDt,Mdt));
	     EstimateBound(i)=2.32*sigma*(Dt-dt);
	 end
	%disp(size(radius))
	figure(028)
	% loglog(DT./N,radius,'o-',DT./N,EstimateBound,'d-');
    loglog(DT./N,EstimateBound,'-or','LineWidth', 1.2, 'MarkerSize', 4);
	hold on
	%disp(length(convFactor));
	%loglog(DT./N,convFactor,'-sg');
    loglog(DT./N,radius,'-sb','LineWidth', 1.2, 'MarkerSize', 4);
	 xlabel('$\Delta t$','interpreter','latex','FontSize',Fontsize);
     ylabel('$\hat \rho$','interpreter','latex','FontSize',Fontsize);
     % legend({'$\rho$','Estimate','$\hat\rho$'},'interpreter','latex','FontWeight','bold','Location','Southeast','FontSize',Fontsize);
     legend({'Estimate(44)','$\hat \rho$'},'interpreter','latex','FontWeight','bold','Location','Southeast','FontSize',Fontsizelegende);
     set(gca,'FontSize',Fontsizescale); % <-- met les graduations (ticks) en fontsize 12
     grid on
     hold off
     print(gcf, 'fig_Dt_negative_sigma.eps', '-depsc', '-r300');  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  M=2^3;
  dt=DT/M;
  N=2^2;
  Dt=DT/N;
    %%% Linear system %%%%%%%%%%%%%%%%%%%%
    Alpha=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3];
     
     radius=0*Alpha;
     EstimateBound=0*Alpha;
     k=length(Alpha);
     for i=1:k
     %N=5
         Mdt=matrix(calL,B,r,A,b,s,L,M,dt,Alpha(i),d);
	     MDt=matrix(calL,B,r,A,b,s,L,N,Dt,Alpha(i),d);
	     radius(i)=spectralradius(itermatrix(MDt,Mdt));
	     EstimateBound(i)=2.32*sigma*(Dt-dt);
	 end
	figure(085)
    loglog(Alpha,EstimateBound,'-or','LineWidth', 1.2, 'MarkerSize', 4);
	hold on
    loglog(Alpha,radius,'-sb','LineWidth', 1.2, 'MarkerSize', 4);
    ylim([1e-2,1]);
	 xlabel('$\alpha$','interpreter','latex','FontSize',Fontsize);
     ylabel('$\hat \rho$','interpreter','latex','FontSize',Fontsize);
     legend({'Estimate (44)','$\hat \rho$'},'interpreter','latex','FontWeight','bold','Location','Southwest','FontSize',Fontsizelegende);
     set(gca,'FontSize',Fontsizescale);
     grid on
     hold off
     print(gcf, 'fig_alpha_negative_sigma.eps', '-depsc', '-r300');  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      M=2^3;
  L=100;
  DT=T/L;
  dt=DT/M;
  N=2^2;
  Dt=DT/N;
    %%% Linear system %%%%%%%%%%%%%%%%%%%%
    calL=-[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2];
     
     radius=0*calL;
     EstimateBound=0*calL;
     k=length(calL);
     for i=1:k
     %N=5
         Mdt=matrix(calL(i),B,r,A,b,s,L,M,dt,alpha,d);
	     MDt=matrix(calL(i),B,r,A,b,s,L,N,Dt,alpha,d);
	     radius(i)=spectralradius(itermatrix(MDt,Mdt));
	     EstimateBound(i)=2.32*abs(calL(i))*(Dt-dt);%/(1+sigma*alpha));
	 end
	figure(052)
	% loglog(DT./N,radius,'o-',DT./N,EstimateBound,'d-');
    loglog(abs(calL),EstimateBound,'-or','LineWidth', 1.2, 'MarkerSize', 4);
	hold on
	%disp(length(convFactor));
	%loglog(DT./N,convFactor,'-sg');
    loglog(abs(calL),radius,'-sb','LineWidth', 1.2, 'MarkerSize', 4);
	xlabel('$|\sigma|$','interpreter','latex','FontSize',Fontsize);
    ylabel('$\hat \rho$','interpreter','latex','FontSize',Fontsize);
     % legend({'$\rho$','Estimate','$\hat\rho$'},'interpreter','latex','FontWeight','bold','Location','Southeast','FontSize',Fontsize);
    legend({'Estimate (44)','$\hat \rho$'},'interpreter','latex','FontWeight','bold','Location','Southeast','FontSize',Fontsizelegende);
    set(gca,'FontSize',Fontsizescale);
     grid on
     hold off
     print(gcf, 'fig_sigma_negative_sigma.eps', '-depsc', '-r300');  


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
    
