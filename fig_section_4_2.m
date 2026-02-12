function fig_section_4_2()
% =========================================================================
% Numerical Simulation Script for ParaOpt on Lotka-Volterra Systems
% This script initializes the parameters and executes both linear and 
% non-linear convergence analyses.
% =========================================================================

% --- Time and Model Parameters ---
T            = 1;          % Final simulation time
a1           = 10;         % Growth/Interaction rate coefficient (sigma)
a2           = 0.2;        % Interaction/Decay coefficient
b1           = 0.2;        % Interaction coefficient
b2           = 10;         % Interaction/Decay coefficient
alpha        = 5e-1;       % Regularization/Coupling parameter (alpha)

% --- Initial and Target States ---
yin          = [20; 15];   % Initial condition vector z(0)
ytg          = [1, 1];     % Target state vector z(T) for the boundary value problem

% --- Discretization and Solver Settings ---
L            = 8;          % Number of temporal subdomains (time slices)
tol          = 1e-9;       % Tolerance for the iterative solver convergence
limit        = 14;         % Refinement level for the fine-grid reference solution

% --- Visualization and Output Configuration ---
figlinear    = 79;         % Figure handle for linear system analysis
fignonlinear = 52;         % Figure handle for non-linear system analysis

% --- Sensitivity Analysis Parameters ---
% Defines the set of discretization ratios r = delta_t / Delta_t
% This generates values: 2^-2, 2^-4, ..., 2^-10
ratio        = 2.^(-(2:2:10));

% =========================================================================
% 1. Execution of the Linearized Convergence Study
% Investigates spectral radius (rho) and empirical convergence for the 
% linearized Lotka-Volterra operator.
% =========================================================================
linearlotka(T, a1, a2, b1, b2, alpha, L, limit, yin, ytg, tol, figlinear);

% =========================================================================
% 2. Execution of the Non-linear Sensitivity Study
% Evaluates the robustness of the convergence estimates for the full 
% non-linear system across varying time-step ratios.
% =========================================================================
nonlinearlotka(T, a1, a2, b1, b2, alpha, L, limit, ratio, yin, ytg, tol, fignonlinear);
end
