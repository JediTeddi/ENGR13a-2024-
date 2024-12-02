% SIRD Model with Age Cohorts
clc; clear; close all;

% Parameters for each cohort
beta = [0.4, 0.3, 0.2];   % Infection rates for 0-30, 30-60, 60+ years
gamma = [0.15, 0.1, 0.05]; % Recovery rates for 0-30, 30-60, 60+ years
mu = [0.005, 0.02, 0.1];  % Death rates for 0-30, 30-60, 60+ years
N = [500, 300, 200];      % Total population for each age cohort

% Interaction matrix (contact rates between cohorts)
contact_matrix = [1, 0.5, 0.2; 
                  0.5, 1, 0.4; 
                  0.2, 0.4, 1];

% Initial conditions
I0 = [1, 1, 1];          % Initial number of infectious individuals per cohort
R0 = [0, 0, 0];          % Initial number of recovered individuals per cohort
D0 = [0, 0, 0];          % Initial number of deceased individuals per cohort
S0 = N - I0;             % Initial number of susceptible individuals per cohort

% Time span
tspan = [0 100]; % Days
initial_conditions = [S0, I0, R0, D0];

% Solve the ODEs
[t, SIRD] = ode45(@(t, SIRD) sird_ode(t, SIRD, beta, gamma, mu, N, contact_matrix), tspan, initial_conditions);

% Extract results for each cohort
S = SIRD(:, 1:3);
I = SIRD(:, 4:6);
R = SIRD(:, 7:9);
D = SIRD(:, 10:12);

% Plot results
figure;
for i = 1:3
    subplot(3, 1, i);
    plot(t, S(:, i), '-b', 'LineWidth', 2); hold on;
    plot(t, I(:, i), '-r', 'LineWidth', 2);
    plot(t, R(:, i), '-g', 'LineWidth', 2);
    plot(t, D(:, i), '-k', 'LineWidth', 2);
    xlabel('Time (days)');
    ylabel('Population');
    legend('Susceptible', 'Infectious', 'Recovered', 'Deceased');
    title(['Age Cohort ', num2str((i-1)*30), '-', num2str(i*30 - 1), ' years']);
    grid on;
end

% Function definition at the end
function dSIRD = sird_ode(~, SIRD, beta, gamma, mu, N, contact_matrix)
    % Unpack state variables
    S = SIRD(1:3);
    I = SIRD(4:6);
    R = SIRD(7:9);
    D = SIRD(10:12);

    % Initialize derivatives
    dS = zeros(3, 1);
    dI = zeros(3, 1);
    dR = zeros(3, 1);
    dD = zeros(3, 1);
    
    % Calculate interaction terms
    for i = 1:3
        % Force of infection from all cohorts
        lambda = 0;
        for j = 1:3
            lambda = lambda + beta(i) * contact_matrix(i, j) * I(j) / N(j);
        end
        
        % SIRD equations for cohort i
        dS(i) = -lambda * S(i);
        dI(i) = lambda * S(i) - gamma(i) * I(i) - mu(i) * I(i);
        dR(i) = gamma(i) * I(i);
        dD(i) = mu(i) * I(i);
    end
    
    % Combine derivatives into a single vector
    dSIRD = [dS; dI; dR; dD];
end
