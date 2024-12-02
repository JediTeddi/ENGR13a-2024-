% SIRD Model (Susceptible-Infectious-Recovered-Deceased)
clc; clear; close all;

% Parameters
beta = 0.9;   % Infection rate
gamma = 0.1;  % Recovery rate
mu = 0.5;    % Death rate
N = 10000;     % Total population

% Initial conditions
I0 = 1;       % Initial number of infectious individuals
R0 = 0;       % Initial number of recovered individuals
D0 = 0;       % Initial number of deceased individuals
S0 = N - I0;  % Initial number of susceptible individuals

% Time span
tspan = [0 1000]; % Days
initial_conditions = [S0 I0 R0 D0];

% Solve the ODEs (pass parameters as additional arguments)
[t, SIRD] = ode45(@(t, SIRD) sird_ode(t, SIRD, beta, gamma, mu, N), tspan, initial_conditions);

% Extract results
S = SIRD(:, 1);
I = SIRD(:, 2);
R = SIRD(:, 3);
D = SIRD(:, 4);

% Plot results
figure;
plot(t, S, '-b', 'LineWidth', 2); hold on;
plot(t, I, '-r', 'LineWidth', 2);
plot(t, R, '-g', 'LineWidth', 2);
plot(t, D, '-k', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Population');
legend('Susceptible', 'Infectious', 'Recovered', 'Deceased');
title('SIRD Model with Deaths');
grid on;

% Function definition at the end
function dSIRD = sird_ode(~, SIRD, beta, gamma, mu, N)
    S = SIRD(1);
    I = SIRD(2);
    R = SIRD(3);
    D = SIRD(4);
    
    dS = -beta * S * I / N;
    dI = beta * S * I / N - gamma * I - mu * I;
    dR = gamma * I;
    dD = mu * I;
    
    dSIRD = [dS; dI; dR; dD];
end
