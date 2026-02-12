% 01/06/2023
% This function works for differents scheme on Cgrid and Fgird
 
function linearlotka(T,a,b,c,d,alpha,L,limit,yin,ytg,tol,figlinear)
    % --- Graphical Configuration ---
    line = 1.6; 
    Fontsize = 20;
    marker = 6;
    labelsize = 17;
    % colors: Red (Estimate), Blue (Spectral Radius), Green (Empirical Factor)
    colors = [1,0,0; 0,0,1; 0,0.6667,0]; 
    
    % --- Physical and Numerical Parameters ---
    sigma = a;              % Growth/Decay rate (often the spectral bound)
    B = eye(2);             % Control/Input matrix
    r = 2;                  % System dimension
    DT = T/L;               % Time interval per subdomain
    M = 2^limit;            % Fine discretization factor
    calL = [a,0; 0,-d];     % Linear operator (Jacobian matrix)
    
    % --- Time Step Initialization ---
    dt = DT/M;              % Fine time step size
    
    % --- Step 1: Reference Fine Solution ---
    % Construct the global system matrix for the fine discretization
    % Numercial scheme 
    A=0;
    b=1;
    s=1;
    M_fine = matrix(calL,B,r,A,b,s,L,M,dt,alpha,d);
    
    % Construct the right-hand side vector (Boundary/Initial conditions)
    rhs = zeros(r*(2*L+1),1);
    rhs(1:2) = yin(:);
    rhs(end-1:end) = ytg(:);
    
    % Direct solve for the reference solution (High-fidelity)
    sol_ref = M_fine \ rhs;
    
    % --- Step 2: Multi-level Coarse Discretization Analysis ---
    k_indices = 1:2:limit-1;
    N_coarse = 2.^k_indices;           % Number of steps for coarse levels
    Dt_vals = DT ./ N_coarse;          % Coarse time step sizes
    
    num_tests = length(k_indices);
    errDt = cell(num_tests);           % To store error trajectories
    convFactor = zeros(1, num_tests);  % Empirical convergence factor (rho_hat)
    rho_spectral = zeros(1, num_tests);% Theoretical spectral radius (rho)
    EstimateBound = zeros(1, num_tests);% Analytical bound (Estimate 38)
    
    Tau_vals = zeros(1, num_tests);    % Tau stability parameter
    L0_vals = zeros(1, num_tests);     % L0 operator constant
    
    itermax = 30; % Maximum iterations for the iterative solver

    % --- Step 3: Loop over Coarse Time Steps ---
    for i = 1:num_tests
        % Construct Coarse system matrix
        M_coarse = matrix(calL,B,r,A,b,s,L,N_coarse(i),Dt_vals(i),alpha,d);
        
        % Construct the iteration matrix (G = I - M_coarse^-1 * M_fine)
        G_iter = itermatrix(M_coarse, M_fine);
        
        % Perform iterative solving to measure empirical convergence
        [~, err_history] = solveriterative(G_iter, M_coarse\rhs, sol_ref, itermax, tol);
        errDt{i} = err_history;
        
        % Calculate  Convergence Factor (Ratio of last two residuals)
        if length(err_history) > 1
            convFactor(i) = err_history(end-1) / err_history(end-2);
        end
        
        % --- Step 4: Spectral and Analytical Estimates ---
        % Compute the spectral radius of the iteration matrix
        rho_spectral(i) = spectralradius(G_iter);
        
        % Compute the analytical upper bound (Estimate 38)
        diff_dt = Dt_vals(i) - dt;
        EstimateBound(i) = max(sigma * diff_dt * (0.5 + ...
            (sigma * diff_dt / 2 + 1) * exp(2 * sigma * DT)), 2.32 * diff_dt * sigma); 
        
        % Compute auxiliary stability constants
        Tau_vals(i) = tau(DT, Dt_vals(i), dt, sigma);
        L0_vals(i) = computeL0(T, L, Dt_vals(i), dt, sigma);
    end
    
    % --- Step 5: Visualization (Convergence Factors) ---
    figure(figlinear)
    
    % Plot Theoretical Estimate (38)
    loglog(Dt_vals, EstimateBound, 'Color', colors(1,:), 'LineStyle', '-', ...
        'LineWidth', line, 'Marker', 'd', 'MarkerSize', marker, 'DisplayName', 'Estimate $(38)$');
    hold on
    
    % Plot Spectral Radius (rho)
    loglog(Dt_vals, rho_spectral, 'Color', colors(2,:), 'LineStyle', '-', ...
        'LineWidth', line, 'Marker', 'o', 'MarkerSize', marker, 'DisplayName', '$\tilde \rho$');
    
    % Plot Convergence Factor (rho_hat)
    loglog(Dt_vals, convFactor, 'Color', colors(3,:), 'LineStyle', '--', ...
        'LineWidth', line, 'Marker', 's', 'MarkerSize', marker, 'DisplayName', '$\bar \rho$');
    
    % Axis and Legend Formatting
    xlabel('$\Delta t$', 'interpreter', 'latex', 'FontSize', Fontsize);
    legend('interpreter', 'latex', 'FontWeight', 'bold', 'Location', 'northwest', 'FontSize', 20);
    set(gca, 'FontSize', labelsize);
    grid on
    hold off
    
    % --- Step 6: Console Output and Error Histories ---
    fprintf('\nTau range: [%.2f, %.2f] | alpha*L0 range: [%.2f, %.2f]\n\n', ...
        min(Tau_vals), max(Tau_vals), min(alpha*abs(L0_vals)), max(alpha*abs(L0_vals)));
    
    figure(80)
    for i = 1:num_tests 
        semilogy(1:length(errDt{i}), errDt{i}, "o-", 'LineWidth', 1);
        hold on
    end
    xlabel('Iteration count');
    ylabel('Residual Error');
    title('Convergence History for Linearized Lotka-Volterra');
    grid on
    hold off
    close(80);
end    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [x,err]=solveriterative(M,b,sol,itermax,tol)
    x=b;
    err=[];
	    for i=1:itermax
	    	x=M*x+b;
	    	err(i)=norm(x-sol,2);
            if err(i)<tol
                break
            end 
	    end
    %err=norm(x-sol,2);
    end
    function L0 = computeL0(T_final, numSteps, Dt_coarse, dt_fine, sigma_val)
tau_val = tau(T_final/numSteps, Dt_coarse, dt_fine, sigma_val);
beta_val = bet(T_final/numSteps, Dt_coarse, sigma_val);
gamma_val = gam(T_final/numSteps, Dt_coarse, sigma_val);
R = zeros(2*numSteps, 1);
R(1) = 1; R(2) = -1; R(end-1) = -1; R(end) = 1 / beta_val;
root_vals = roots(R);
tau_0 = max(real(root_vals(abs(imag(root_vals)) < 1e-12 & real(root_vals) > 1 & real(root_vals) < 1.5)));
L0 = (beta_val - tau_val) / (gamma_val * (tau_val - tau_0));
end
    function rho_test(L,B,r,A1,b1,s1,A2,b2,s2,Len,DT,alpha,d1,d2,N_test,M_test,cgrid,fgrid,figname)
    %%%%% Dicretization Parameters Test
    Fontsize=19;
    Size_line=1.5;
    labelsize=16;
    fMatrix=matrix(L,B,r,A1,b1,s1,Len,DT/fgrid,fgrid,alpha,d1);
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
    
