function fig_section_4_3()
    % --- Configuration des paramètres ---
    ratio = 10.^(-[0,1,2,3]);
    convergence_factors = zeros(1, length(ratio));
    all_errors = cell(1, length(ratio)); 
    EstimateBound=zeros(1,length(ratio));
    gg=9.81; ll=0.25;
    sigma=sqrt(gg/ll);
    
    % --- Simulation et Collecte des données ---
    fprintf('\t Starting simulations for different ratios...\n');
    %load("errorexplicitparaopt1.mat");
    dt=5e-6;
    iter=0;
    ll=0.5;
    DT=1/8;
    fprintf(2,'\t Exact solution computation using Newtion iteration\n')
    [ExactSol,~]=pendulum_test(1,ll,0);
    fprintf(2,'\t Error computation using ParaOpt\n')
    for idx = 1:length(ratio)
        r = ratio(idx);
        % On suppose que explicitparaopt(r) retourne le vecteur d'erreur par itération
        [~,errors_vector] = pendulum_test(r,ll,ExactSol); 
        all_errors{idx} = errors_vector;
        %errors_vector=all_errors{idx};
        % Calcul du facteur de convergence empirique (pente de l'erreur)
        if length(errors_vector) > 1
            convergence_factors(idx) = errors_vector(end-1) / errors_vector(end-2);
            diff_dt =dt/r;
            EstimateBound(idx) = max(sigma * diff_dt * (0.5 + ...
                (sigma * diff_dt / 2 + 1) * exp(2 * sigma * DT)), 2.32 * diff_dt * sigma);
        end
    end

    % Sauvegarde des résultats
    save("errorexplicitparaopt2.mat", "all_errors", "convergence_factors", "ratio","EstimateBound");

    % --- Graphique 1 : Erreur en fonction des itérations ---
    figure(70);
    for idx = 1:length(ratio)
        current_err = all_errors{idx};
        semilogy(1:length(current_err)-1, current_err(1:end-1), '-o', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Ratio $r = 10^{-%d}$', idx-1));
        hold on;
    end
    grid on;
    xlabel('\# Iters', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\|\cdot\|_2$-{error}', 'Interpreter', 'latex', 'FontSize', 14);
    %title('Convergence History for Pendulum Problem', 'FontSize', 16);
    legend('Interpreter', 'latex', 'Location', 'northeast');
    hold off;

    % --- Graphique 2 : Facteur de convergence en fonction du Ratio ---
    figure(267);
    loglog(ratio, convergence_factors, 'd--g', 'LineWidth', 2, ...
        'MarkerSize', 6);
    hold on
    loglog(ratio, EstimateBound, 's-r', 'LineWidth', 2, ...
        'MarkerSize', 6);
    grid on;
    xlabel('Ratio $r = \frac{\delta t}{\Delta t}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Convergence Factor $\hat{\rho}$', 'Interpreter', 'latex', 'FontSize', 14);
    legend("$\bar \rho$","Estimate (38)","interpreter","latex");
    %title('Sensitivity of Convergence Factor to Step Size Ratio', 'FontSize', 16);
    
    % Inverser l'axe X pour montrer la réduction du ratio vers la droite
   % set(gca, 'XDir', 'reverse'); 
end