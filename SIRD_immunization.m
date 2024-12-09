% Parameters
% beta = 0.8;   % Infection rate
% gamma = 0.001;  % Recovery rate
% mu = 0.1;    % Death rate
% nu = 0.9;    % Vaccination rate
% rho = 0.8;   % Reinfection rate
% N = 10000;     % Total population

% Initial conditions
S0 = 5000;     % Initial susceptible population
I0 = 100;      % Initial infectious population
R0 = 10;       % Initial recovered population
D0 = 100;       % Initial deaths
initial_conditions = [S0, I0, R0, D0];

% Time vector
tspan = [0 1000];

% Solve the system of ODEs
[t, Pop] = ode45(@sirvd_ode, tspan, initial_conditions);

% Extract results
S = Pop(:, 1);
I = Pop(:, 2);
R = Pop(:, 3);
D = Pop(:, 4);

% Plot results
figure;
plot(t, S, 'b', 'LineWidth', 2); hold on;
plot(t, I, 'r', 'LineWidth', 2);
plot(t, R, 'g', 'LineWidth', 2);
plot(t, D, 'k', 'LineWidth', 2);
xlabel('Time');
ylabel('Population');
legend('Susceptible', 'Infectious', 'Recovered', 'Deaths');
title('SIRD Model with Vaccination and Reinfection');
grid on;

% ODE function
function dPop = sirvd_ode(~, Pop, tspan)
    % Parameters
    beta = 0.8 * (1 + 0.2 * sin(0.01 * tspan));   % Infection rate
    gamma = 0.1;  % Recovery rate
    mu = 0.1;    % Death rate
    nu = 0.6 + 0.2 * exp(-0.005 * tspan);    % Vaccination rate
    rho = 0.8;   % Reinfection rate
    N = 10000;     % Total population

    % Population compartments
    S = Pop(1);
    I = Pop(2);
    R = Pop(3);
    D = Pop(4);

    % ODE system
    dS = -beta * S * I / N - nu * S + rho * R;
    dI = beta * S * I / N - gamma * I - mu * I;
    dR = gamma * I + nu * S - rho * R;
    dD = mu * I;

    dPop = [dS; dI; dR; dD];
end
