% SIR Model Simulation
%Params (CHANGE THESE) 
beta = 0.3;    % Transmission rate (probability of contact rates and transmission 
gamma = 0.1;   % Recovery rate
N = 1000000;   % Total population
I0 = 1000;     % Initial number of infected individuals
R0 = 0;        % Initial number of recovered individuals
S0 = N - I0 - R0; % Initial number of susceptible individuals
days = 100;    % Duration of simulation (days)
tspan = [0 days]; % Time Span

% Initial Conditions
y0 = [S0; I0; R0];

% Differential equations of SIR Model 
function dydt = sir_model(t, y, beta, gamma, N)
    S = y(1); % Susceptible
    I = y(2); % Infected
    R = y(3); % Recovered

    dS = -beta * S * I / N;            % Rate of change of Susceptible
    dI = beta * S * I / N - gamma * I; % Rate of change of Infected
    dR = gamma * I;                    % Rate of change of Recovered

    dydt = [dS; dI; dR];
end

% Solve the ODEs
[t, y] = ode45(@(t, y) sir_model(t, y, beta, gamma, N), tspan, y0);

% Extract Results
S = y(:, 1); % Susceptible
I = y(:, 2); % Infected
R = y(:, 3); % Recovered

% Plot the Results
figure;
hold on;
plot(t, S, 'b', 'LineWidth', 2); % Susceptible
plot(t, I, 'r', 'LineWidth', 2); % Infected
plot(t, R, 'g', 'LineWidth', 2); % Recovered
hold off;

% Customize Plot
title('SIR Model Simulation');
xlabel('Time (days)');
ylabel('Population');
legend('Susceptible', 'Infected', 'Recovered');
grid on;
