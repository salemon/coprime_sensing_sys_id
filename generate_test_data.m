% Generate simulated DC motor data with PRBS input
% This creates a test dataset for the multirate system identification
% Bob Xiaohai Hu, UW ME, 2024

clear; clc; close all;

%% Simulation parameters
Fs = 1000; % Sampling frequency (Hz)
Ts = 1/Fs; % Sampling time
duration = 100; % Total duration (seconds)
N = duration * Fs; % Total number of samples
t = (0:N-1)' * Ts; % Time vector

%% Create a realistic DC motor model (4th order)
% DC motor transfer function: G(s) = K / ((Js+b)(Ls+R) + K^2)
% Parameters for a typical DC motor
J = 0.01;   % Moment of inertia (kg.m^2)
b = 0.1;    % Damping coefficient (N.m.s)
K = 0.01;   % Electromotive force constant
R = 1;      % Electrical resistance (Ohm)
L = 0.5;    % Electrical inductance (H)

% State-space representation
A_true = [-R/L, -K/L, 0, 0;
          K/J, -b/J, 0, 0;
          0, 1, 0, 0;
          0, 0, 1, 0];
B_true = [1/L; 0; 0; 0];
C_true = [0, 0, 1, 0]; % Position output (encoder)
D_true = 0;

% Create continuous system
sys_cont = ss(A_true, B_true, C_true, D_true);

% Discretize the system
sys_discrete = c2d(sys_cont, Ts, 'zoh');

%% Generate PRBS input signal
% Using the create_prbs function
ValUinit = 0;      % Initial value
ValAmpli = 2;      % Amplitude (Volts)
ValDecal = 0;      % DC offset
ValLgReg = 10;     % Register length (affects PRBS period)
ValDivi = 5;       % Frequency divider (affects switching rate)
Nsamp = N;         % Number of samples
Tappli = 1000;     % Application instant (samples)

% Generate PRBS
input_voltage = create_prbs(ValUinit, ValAmpli, ValDecal, ValLgReg, ValDivi, Nsamp, Tappli);
input_voltage = input_voltage(:); % Ensure column vector

%% Simulate the system response
% Initialize state
x = zeros(4, 1);
output_encoder = zeros(N, 1);

% Add some measurement noise
noise_level = 0.01; % Small noise

for i = 1:N
    % State update
    x = sys_discrete.A * x + sys_discrete.B * input_voltage(i);
    
    % Output with noise
    output_encoder(i) = sys_discrete.C * x + sys_discrete.D * input_voltage(i) + ...
                        noise_level * randn();
end

%% Create the data structure matching the expected format
% Structure: m_data.in_voltage.signals.values
%            m_data.in_voltage.time
%            m_data.out_encoder.signals.values

m_data = struct();

% Input voltage structure
m_data.in_voltage = struct();
m_data.in_voltage.time = t;
m_data.in_voltage.signals = struct();
m_data.in_voltage.signals.values = input_voltage;

% Output encoder structure  
m_data.out_encoder = struct();
m_data.out_encoder.signals = struct();
m_data.out_encoder.signals.values = output_encoder;

%% Save the data
% Create directory if it doesn't exist
if ~exist('Model Data', 'dir')
    mkdir('Model Data');
end

save('Model Data/prbs_run_1.mat', 'm_data');

fprintf('✓ Test data generated successfully!\n');
fprintf('  - Sampling frequency: %.0f Hz\n', Fs);
fprintf('  - Duration: %.0f seconds\n', duration);
fprintf('  - Total samples: %d\n', N);
fprintf('  - File saved: Model Data/prbs_run_1.mat\n');

%% Visualize the generated data (optional)
figure('Position', [100, 100, 1000, 600]);

subplot(2,1,1);
plot(t, input_voltage, 'b', 'LineWidth', 1);
grid on;
xlabel('Time (s)');
ylabel('Input Voltage (V)');
title('PRBS Input Signal');
xlim([0, min(20, duration)]);

subplot(2,1,2);
plot(t, output_encoder, 'r', 'LineWidth', 1);
grid on;
xlabel('Time (s)');
ylabel('Encoder Output');
title('Simulated Motor Response');
xlim([0, min(20, duration)]);

%% Display system information
fprintf('\n--- Simulated System Info ---\n');
fprintf('True system poles (continuous):\n');
disp(pole(sys_cont));
fprintf('True system poles (discrete, Ts=%.4f):\n', Ts);
disp(pole(sys_discrete));
