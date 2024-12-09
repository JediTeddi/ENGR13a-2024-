function sirvd_cohorts
    % Parameters for younger cohort (cohort 1)
    beta1 = 0.4;  % Infection rate for younger
    gamma1 = 0.1; % Recovery rate for younger
    mu1 = 0.005;  % Death rate for younger
    nu1 = 0.05;   % Vaccination rate for younger
    rho1 = 0.02;  % Reinfection rate for younger

    % Parameters for older cohort (cohort 2)
    beta2 = 0.3;  % Infection rate for older
    gamma2 = 0.08; % Recovery rate for older
    mu2 = 0.02;   % Death rate for older
    nu2 = 0.04;   % Vaccination rate for older
    rho2 = 0.01;  % Reinfection rate for older

    % Interaction between cohorts (contact rates)
    contact_matrix = [0.8, 0.5;  % Younger interacting with younger and older
                      0.5, 0.7]; % Older interacting with younger and older

    % Total population
    N1 = 600; % Population size of younger cohort
    N2 = 400; % Population size of older cohort

    % Initial conditions
    S1_0 = 590; I1_0 = 10; R1_0 = 0; D1_0 = 0; % Younger cohort
    S2_0 = 390; I2_0 = 10; R2_0 = 0; D2_0 = 0; % Older cohort
    initial_conditions = [S1_0, I1_0, R1_0, D1_0, S2_0, I2_0, R2_0, D2_0];

    % Time vector
    tspan = [0 100];

    % Solve the system of ODEs
    params = struct('beta1', beta1, 'gamma1', gamma1, 'mu1', mu1, 'nu1', nu1, 'rho1', rho1, ...
                    'beta2', beta2, 'gamma2', gamma2, 'mu2', mu2, 'nu2', nu2, 'rho2', rho2, ...
                    'contact_matrix', contact_matrix, 'N1', N1, 'N2', N2);

    [t, Pop] = ode45(@(t, Pop) sirvd_ode(t, Pop, params), tspan, initial_conditions);

    % Extract results
    S1 = Pop(:, 1); I1 = Pop(:, 2); R1 = Pop(:, 3); D1 = Pop(:, 4);
    S2 = Pop(:, 5); I2 = Pop(:, 6); R2 = Pop(:, 7); D2 = Pop(:, 8);

    % Plot results
    figure;
    subplot(2, 1, 1);
    plot(t, S1, 'b', 'LineWidth', 2); hold on;
    plot(t, I1, 'r', 'LineWidth', 2);
    plot(t, R1, 'g', 'LineWidth', 2);
    plot(t, D1, 'k', 'LineWidth', 2);
    xlabel('Time');
    ylabel('Population');
    legend('S1 (Younger)', 'I1 (Younger)', 'R1 (Younger)', 'D1 (Younger)');
    title('Younger Cohort');
    grid on;

    subplot(2, 1, 2);
    plot(t, S2, 'b--', 'LineWidth', 2); hold on;
    plot(t, I2, 'r--', 'LineWidth', 2);
    plot(t, R2, 'g--', 'LineWidth', 2);
    plot(t, D2, 'k--', 'LineWidth', 2);
    xlabel('Time');
    ylabel('Population');
    legend('S2 (Older)', 'I2 (Older)', 'R2 (Older)', 'D2 (Older)');
    title('Older Cohort');
    grid on;
end

% ODE Function
function dPop = sirvd_ode(~, Pop, params)
    % Extract parameters
    beta1 = params.beta1; gamma1 = params.gamma1; mu1 = params.mu1; nu1 = params.nu1; rho1 = params.rho1;
    beta2 = params.beta2; gamma2 = params.gamma2; mu2 = params.mu2; nu2 = params.nu2; rho2 = params.rho2;
    contact_matrix = params.contact_matrix; N1 = params.N1; N2 = params.N2;

    % Population compartments
    S1 = Pop(1); I1 = Pop(2); R1 = Pop(3); D1 = Pop(4);
    S2 = Pop(5); I2 = Pop(6); R2 = Pop(7); D2 = Pop(8);

    % Force of infection (interaction between cohorts)
    lambda1 = beta1 * (contact_matrix(1, 1) * I1 / N1 + contact_matrix(1, 2) * I2 / N2);
    lambda2 = beta2 * (contact_matrix(2, 1) * I1 / N1 + contact_matrix(2, 2) * I2 / N2);

    % ODEs for younger cohort
    dS1 = -lambda1 * S1 - nu1 * S1 + rho1 * R1;
    dI1 = lambda1 * S1 - gamma1 * I1 - mu1 * I1;
    dR1 = gamma1 * I1 + nu1 * S1 - rho1 * R1;
    dD1 = mu1 * I1;

    % ODEs for older cohort
    dS2 = -lambda2 * S2 - nu2 * S2 + rho2 * R2;
    dI2 = lambda2 * S2 - gamma2 * I2 - mu2 * I2;
    dR2 = gamma2 * I2 + nu2 * S2 - rho2 * R2;
    dD2 = mu2 * I2;

    % Combine results
    dPop = [dS1; dI1; dR1; dD1; dS2; dI2; dR2; dD2];
end
