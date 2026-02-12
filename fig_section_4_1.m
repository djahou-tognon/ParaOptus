function fig_section_4_1()
% Simulates an optimal control problem where alpha is the regularization parameter.
% Performs several convergence tests (Test1–Test5) based on time discretization.
line=1.5;Fontsize=18;
alpha_reg = 1.0;      % Regularization parameter α
sigma_val = 9.0; 1.0;      % System parameter σ
T_final = 1.0;        % Final simulation time
numSteps = 10;        % (L) Number of coarse time intervals
DT = T_final / numSteps;  % Coarse time step

%%%%%%%%%%%%% Test 1 %%%%%%%%%%%
% Convergence with varying fine grid refinement
Dt_coarse = DT / 4;                    % Coarse substep
dt_fine = Dt_coarse ./ (2.^(1:20));    % Sequence of fine time steps
run_test1(alpha_reg, sigma_val, T_final, numSteps, Dt_coarse, dt_fine,line,Fontsize);
clear dt_fine Dt_coarse

%%%%%%%%%%%%% Test 2 %%%%%%%%%%%
% Convergence with varying coarse grid refinement
Dt_coarse = DT ./ 2.^(1:19);
dt_fine = Dt_coarse(1) * 2^(-20);
run_test2(alpha_reg, sigma_val, T_final, numSteps, Dt_coarse, dt_fine,line,Fontsize);

%%%%%%%%%%%%% Test 4 %%%%%%%%%%%
% Convergence with varying regularization alpha refinement
alpha_range = logspace(-3, 4, 30);  % Define a range of alpha values for testing
Dt_coarse = DT / 2;                  % Set coarse time step for this test
dt_fine = DT / 2^4;                  % Set fine time step for this test
run_test4(alpha_range, sigma_val, T_final, numSteps, Dt_coarse, dt_fine, line, Fontsize);

%%%%%%%%%%%%% Test 5 %%%%%%%%%%%
% Convergence with varying sigma refinement
% This test is to show the gap between the estimate (34) and (35).
sigma_val=[linspace(1e-2,1e-1,5),linspace(1e-1,1e0,5),linspace(1e0,1e1,5),linspace(1e1,1e2,5)];
Dt_coarse = DT / 2;
dt_fine = DT / 2^4;
run_test5(alpha_reg, sigma_val, T_final, numSteps, Dt_coarse, dt_fine,line,Fontsize);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 1: Convergence with decreasing fine step dt (fixed Dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_test1(alpha_reg, sigma_val, T_final, numSteps, Dt_coarse, dt_fine,line,Fontsize)
num_dt = length(dt_fine);
Tau = zeros(num_dt,1);
Results = zeros(3,num_dt);
L0_vals = zeros(num_dt,1);
Beta_vals = zeros(num_dt,1);

for k = 1:num_dt
    [Results(1,k), Results(2,k), Results(3,k)] = computeSpectralRadius(alpha_reg, sigma_val, T_final, numSteps, Dt_coarse, dt_fine(k));
    Tau(k) = computeTau(T_final/numSteps, Dt_coarse, dt_fine(k), sigma_val);
    L0_vals(k) = computeL0(T_final, numSteps, Dt_coarse, dt_fine(k), sigma_val);
    Beta_vals(k) = computeBeta(T_final/numSteps, Dt_coarse, sigma_val);
end

fprintf('\nTau is in =[%.2f,%.2f]| Beta in [%.2f,%.2f]| alpha*L0 in [%.2f,%.2f]\n\n', min(Tau),max(Tau), min(Beta_vals),max(Beta_vals), min(alpha_reg*abs(L0_vals)),max(alpha_reg*abs(L0_vals)));

% Log-log plot of spectral radius and estimates
figure(111)
loglog(dt_fine, Results(1,:), '--ob', 'LineWidth', line);
hold on
loglog(dt_fine, Results(2,:), '-sg', 'LineWidth', line);
loglog(dt_fine, Results(3,:), '-dr', 'LineWidth', line);
grid on
xlabel('$\delta t$', 'Interpreter', 'latex', 'FontSize', Fontsize);
ylabel('$\rho$', 'Interpreter', 'latex','FontSize', Fontsize);
legend('$\rho$','Estimate $(35)$', 'Estimate $(34)$', 'Interpreter', 'latex','FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
%ylim([5e-3,1e-1]);
hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 2: Convergence with decreasing coarse step Dt (fixed dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_test2(alpha_reg, sigma_val, T_final, numSteps, Dt_coarse, dt_fine,line,Fontsize)
num_Dt = length(Dt_coarse);
Tau = zeros(num_Dt,1);
Results = zeros(3,num_Dt);
L0_vals = zeros(num_Dt,1);
Beta_vals = zeros(num_Dt,1);

for k = 1:num_Dt
    [Results(1,k), Results(2,k), Results(3,k)] = computeSpectralRadius(alpha_reg, sigma_val, T_final, numSteps, Dt_coarse(k), dt_fine);
    Tau(k) = computeTau(T_final/numSteps, Dt_coarse(k), dt_fine, sigma_val);
    L0_vals(k) = computeL0(T_final, numSteps, Dt_coarse(k), dt_fine, sigma_val);
    Beta_vals(k) = computeBeta(T_final/numSteps, Dt_coarse(k), sigma_val);
end

fprintf('\nTau is in =[%.2f,%.2f]| Beta in [%.2f,%.2f]| alpha*L0 in [%.2f,%.2f]\n\n', min(Tau),max(Tau), min(Beta_vals),max(Beta_vals), min(alpha_reg*abs(L0_vals)),max(alpha_reg*abs(L0_vals)));

% Plot spectral radius vs Dt
figure(222)
loglog(Dt_coarse, Results(1,:), '--ob', 'LineWidth', line);
hold on

loglog(Dt_coarse, Results(2,:), '-sg', 'LineWidth', line);
loglog(Dt_coarse, Results(3,:), '-dr', 'LineWidth', line);
grid on
xlabel('$\Delta t$', 'Interpreter', 'latex', 'FontSize', Fontsize);
ylabel('$\rho$', 'Interpreter', 'latex','FontSize', Fontsize);
legend('$\rho$','Estimate $(35)$', 'Estimate $(34)$', 'Interpreter', 'latex','FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 3: Convergence with vectorized fine steps dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_test3(alpha_reg, sigma_val, T_final, numSteps, Dt_coarse, dt_fine)
num_dt = length(dt_fine);
Tau = zeros(num_dt,1);
Results = zeros(3, length(Dt_coarse));
L0_vals = zeros(length(Dt_coarse),1);

for k = 1:length(Dt_coarse)
    [Results(1,k), Results(2,k), Results(3,k)] = computeSpectralRadius(alpha_reg, sigma_val, T_final, numSteps, Dt_coarse(k), dt_fine(k));
    Tau(k) = computeTau(T_final/numSteps, Dt_coarse(k), dt_fine(k), sigma_val);
    L0_vals(k) = computeL0(T_final, numSteps, Dt_coarse(k), dt_fine(k), sigma_val);
end

% Identify tau threshold crossings
tau_idx = 0;
for k = 1:length(Tau)-1
    if Tau(k)-1.5>0 && Tau(k)*Tau(k+1)-9/4<0
        tau_idx = k;
    elseif Tau(k)-1.5<0 && Tau(k)*Tau(k+1)-9/4<0
        tau_idx = k+1;
    end
end
disp(tau_idx)

% Plot results
figure(333)
loglog(Dt_coarse, Results(1,:), '-bo', 'MarkerSize', 4, 'LineWidth', 2.2);
hold on
%loglog(Dt_coarse, Results(2,:), '-sg', 'LineWidth', 1);
loglog(Dt_coarse, Results(3,:), '-dr', 'LineWidth', 2);
grid on
xlabel('$\Delta t$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\rho$', 'Interpreter', 'latex');
legend('$\rho$', 'Estimate $(34)$', 'Interpreter', 'latex');
hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 4: Convergence vs regularization parameter alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_test4(alpha_range, sigma_val, T_final, numSteps, Dt_coarse, dt_fine,line,Fontsize)
Tau_val = computeTau(T_final/numSteps, Dt_coarse, dt_fine, sigma_val);
L0_val = computeL0(T_final, numSteps, Dt_coarse, dt_fine, sigma_val);
L0_alpha = alpha_range .* L0_val;
Beta_vals = bet(T_final/numSteps,Dt_coarse,sigma_val);

Results = zeros(3,length(alpha_range));
for k = 1:length(alpha_range)
    [Results(1,k), Results(2,k), Results(3,k)] = computeSpectralRadius(alpha_range(k), sigma_val, T_final, numSteps, Dt_coarse, dt_fine);
end
fprintf('\n Tau =%.2f| Beta =%.2f| alpha*L0 in [%.2f,%.2f]\n\n', Tau_val, Beta_vals, min(abs(L0_alpha)),max(abs(L0_alpha)));

% Plot results
figure(444)
loglog(alpha_range, Results(1,:), '--ob', 'LineWidth', line);
hold on
plot(alpha_range, Results(3,:), '-dr',  'LineWidth', line);
grid on; grid minor
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', Fontsize);
ylabel('$\rho$', 'Interpreter', 'latex', 'FontSize', Fontsize);
legend('$\rho$', 'Estimate $(34)$', 'Interpreter', 'latex','FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 5: Convergence vs number of coarse steps L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_test5(alpha_reg, sigma_val, T_final, L_vec, Dt_vec, dt_vec,line,Fontsize)
Results = zeros(3,length(sigma_val));

for k = 1:length(sigma_val)
    [Results(1,k), Results(2,k), Results(3,k)] = computeSpectralRadius(alpha_reg, sigma_val(k), T_final, L_vec, Dt_vec, dt_vec);
end

% Plot results
figure(555)
loglog(sigma_val, Results(1,:), '--bo', 'LineWidth', line);
hold on
loglog(sigma_val, Results(2,:), '-sg', 'LineWidth', line);
loglog(sigma_val, Results(3,:), '-dr', 'LineWidth', line);
grid on
xlabel('$\sigma$', 'Interpreter', 'latex','FontSize', Fontsize);
ylabel('$\rho$', 'Interpreter', 'latex','FontSize', Fontsize);
legend('$\rho$','Estimate $(35)$', 'Estimate $(34)$', 'Interpreter', 'latex','FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute spectral radius and estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho_val, est_spectrum, est_bound] = computeSpectralRadius(alpha_reg, sigma_val, T_final, numSteps, Dt_coarse, dt_fine)
DT = T_final / numSteps;
beta_Dt = beta_fun(DT, Dt_coarse, sigma_val);
dbeta = beta_Dt - beta_fun(DT, dt_fine, sigma_val);
gamma_Dt = gamma_fun(DT, Dt_coarse, sigma_val);
dgamma = gamma_Dt - gamma_fun(DT, dt_fine, sigma_val);
eta_val = dbeta / dgamma;
tau_val = beta_fun(DT, Dt_coarse, sigma_val) - gamma_fun(DT, Dt_coarse, sigma_val) * eta_val;

% Characteristic polynomial for eigenvalues
poly_coef = zeros(2*numSteps, 1);
poly_coef(1:2:end) = 1;
poly_coef(2:2:end-1) = -tau_val;
poly_coef(end) = alpha_reg * eta_val - tau_val;

roots_val = roots(poly_coef);
eigen_vals = dbeta ./ (beta_Dt - roots_val);
rho_val = max(abs(eigen_vals)); % Spectral radius

% Estimates
est_spectrum = abs(dgamma) / gamma_Dt;
est_bound = sigma_val * (Dt_coarse - dt_fine) * (0.5 + (sigma_val * Dt_coarse + 1) * exp(2*sigma_val*DT));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = beta_fun(T, t, sigma_val)
val = (1 + sigma_val * t)^(T / t);
end

function val = gamma_fun(T, t, sigma_val)
val = ((beta_fun(T, t, sigma_val)^2) - 1) / (sigma_val * (2 + sigma_val * t));
end

function val = computeBeta(T, Dt, sigma_val)
val = beta_fun(T, Dt, sigma_val);
end

function val = computeTau(T, Dt, dt, sigma_val)
val = beta_fun(T, Dt, sigma_val) - gamma_fun(T, Dt, sigma_val) * ((beta_fun(T, Dt, sigma_val) - beta_fun(T, dt, sigma_val)) / ...
    (gamma_fun(T, Dt, sigma_val) - gamma_fun(T, dt, sigma_val)));
end

function L0 = computeL0(T_final, numSteps, Dt_coarse, dt_fine, sigma_val)
tau_val = computeTau(T_final/numSteps, Dt_coarse, dt_fine, sigma_val);
beta_val = beta_fun(T_final/numSteps, Dt_coarse, sigma_val);
gamma_val = gamma_fun(T_final/numSteps, Dt_coarse, sigma_val);
R = zeros(2*numSteps, 1);
R(1) = 1; R(2) = -1; R(end-1) = -1; R(end) = 1 / beta_val;
root_vals = roots(R);
tau_0 = max(real(root_vals(abs(imag(root_vals)) < 1e-12 & real(root_vals) > 1 & real(root_vals) < 1.5)));
L0 = (beta_val - tau_val) / (gamma_val * (tau_val - tau_0));
end

