%% Simplified Multirate System Identification Test
% This script tests the multirate N4SID algorithm without requiring
% external data files or Simulink models
% 
% Bob Xiaohai Hu, UW ME, 2024
% Based on: "State-Space System Identification beyond the Nyquist Frequency 
% with Collaborative Non-Uniform Sensing Data"

clear; clc; close all;

fprintf('\n');
fprintf('========================================================\n');
fprintf('  Multirate System Identification - Simplified Test\n');
fprintf('========================================================\n\n');

%% 1. SYSTEM SETUP
fprintf('Step 1: Setting up test system...\n');

% Sampling parameters
Fs = 1024;              % Base sampling frequency (Hz)
Ts = 1/Fs;              % Sampling time
m = 2;                  % Zero-order hold period multiplier
n1 = 3;                 % Downsampling factor for sensor 1
n2 = 5;                 % Downsampling factor for sensor 2

fprintf('  - Base sampling freq: %d Hz\n', Fs);
fprintf('  - Zero-order hold: %dT\n', m);
fprintf('  - Sensor 1 period: %dT (Nyquist: %.1f Hz)\n', n1, Fs/(2*n1));
fprintf('  - Sensor 2 period: %dT (Nyquist: %.1f Hz)\n', n2, Fs/(2*n2));
fprintf('  - Collaborative Nyquist: %.1f Hz\n', Fs/(2*m));

% Create a 4th-order test system (DC motor-like)
A_true = [-10, -5, 0, 0;
          5, -8, 0, 0;
          0, 1, -2, 0;
          0, 0, 1, -1];
B_true = [1; 0; 0; 0];
C_true = [0, 0, 1, 0.5];
D_true = 0;

sys_cont = ss(A_true, B_true, C_true, D_true);
sys_discrete = c2d(sys_cont, Ts, 'zoh');

fprintf('\n  True system poles (continuous):\n');
disp(pole(sys_cont)');
fprintf('  True system order: %d\n', size(A_true, 1));

%% 2. GENERATE PRBS INPUT
fprintf('\nStep 2: Generating PRBS input signal...\n');

duration = 60;          % seconds
N = duration * Fs;

% PRBS parameters
ValUinit = 0;
ValAmpli = 1.5;
ValDecal = 0;
ValLgReg = 10;
ValDivi = 5;
Nsamp = N;
Tappli = 1000;

input_signal = create_prbs(ValUinit, ValAmpli, ValDecal, ValLgReg, ...
                           ValDivi, Nsamp, Tappli);
input_signal = input_signal(:);

fprintf('  - Signal duration: %d seconds\n', duration);
fprintf('  - Total samples: %d\n', N);
fprintf('  - PRBS register length: %d\n', ValLgReg);

%% 3. SIMULATE SYSTEM RESPONSE
fprintf('\nStep 3: Simulating system response...\n');

% Simulate full-rate output
x = zeros(4, 1);
output_full = zeros(N, 1);
noise_level = 0.01;

for i = 1:N
    x = sys_discrete.A * x + sys_discrete.B * input_signal(i);
    output_full(i) = sys_discrete.C * x + noise_level * randn();
end

% Apply zero-order hold to input
input_zoh = input_signal(1:m:end);

% Downsample outputs for sensors
output_sensor1 = output_full(1:n1:end);
output_sensor2 = output_full(1:n2:end);

fprintf('  - Full-rate samples: %d\n', length(output_full));
fprintf('  - ZOH input samples: %d\n', length(input_zoh));
fprintf('  - Sensor 1 samples: %d\n', length(output_sensor1));
fprintf('  - Sensor 2 samples: %d\n', length(output_sensor2));

%% 4. LIFTING PROCESS
fprintf('\nStep 4: Applying lifting transformation...\n');

n1n2 = n1 * n2;
n1m = n1 * m;
n2m = n2 * m;

% Align dimensions
Q1 = fix(length(input_zoh) / n1n2);
R1 = mod(length(input_zoh), n1n2);
Q2 = fix(length(output_sensor1) / n2m);
R2 = mod(length(output_sensor1), n2m);
Q3 = fix(length(output_sensor2) / n1m);
R3 = mod(length(output_sensor2), n1m);

% Reshape for lifting
U_lifted = reshape(input_zoh(1:end-R1), 1, n1n2, []);
Y1_lifted = reshape(output_sensor1(1:end-R2), 1, n2m, []);
Y2_lifted = reshape(output_sensor2(1:end-R3), 1, n1m, []);

% Stack outputs
Y_stack = cat(2, Y1_lifted, Y2_lifted);
U_lifted2(:,:) = U_lifted(1,:,:);
Y_stack2(:,:) = Y_stack(1,:,:);

uu = U_lifted2';
yy = Y_stack2';

% Remove mean
uu = uu - mean(uu);
yy = yy - mean(yy);

fprintf('  - Lifted input size: %dx%d\n', size(uu));
fprintf('  - Lifted output size: %dx%d\n', size(yy));
fprintf('  - Lifted period: %dT = %.4f seconds\n', m*n1n2, m*n1n2*Ts);

%% 5. N4SID IDENTIFICATION (LIFTED SYSTEM)
fprintf('\nStep 5: Running N4SID on lifted system...\n');

k = 6;  % Hankel block order
try
    [A_lift, B_lift, C_lift, D_lift, xf, n_est, Uf, Yf, U, Y] = ...
        n4sidkatamodar(uu, yy, k);
    
    fprintf('  ✓ N4SID completed successfully\n');
    fprintf('  - Estimated order: %d\n', n_est);
    fprintf('  - Lifted A matrix size: %dx%d\n', size(A_lift));
    
catch ME
    fprintf('  ✗ N4SID failed: %s\n', ME.message);
    fprintf('  Trying alternative N4SID version...\n');
    
    [A_lift, B_lift, C_lift, D_lift, xf, n_est, Uf, Yf, U, Y] = ...
        n4sidkatamodar_TC(uu, yy, k);
    
    fprintf('  ✓ Alternative N4SID completed\n');
    fprintf('  - Estimated order: %d\n', n_est);
end

%% 6. RECOVER FAST SINGLE-RATE MODEL
fprintf('\nStep 6: Recovering fast single-rate model...\n');

% Extract last two B columns
B_m = B_lift(:, end);
B_m2 = B_lift(:, end-1);

% Eigenvalue decomposition approach
[V, Lambda] = eig(A_lift);
V_inv = inv(V);

% Compute lambda ratios
V_inv_B_n1 = V_inv * B_m;
V_inv_B_n2 = V_inv * B_m2;
lambda_ratios = V_inv_B_n2 ./ V_inv_B_n1;

% Reconstruct fast-rate A matrix
A_fast = real(V * diag(lambda_ratios) * V_inv);

% Reconstruct B matrix
B_fast = real(inv(eye(n_est) + A_fast) * B_m);

% C and D matrices
C_fast = C_lift(1, :);
D_fast = D_lift(1, 1);

% Create identified system
sys_identified = ss(A_fast, B_fast, C_fast, D_fast, Ts);

fprintf('  ✓ Fast model recovered\n');
fprintf('  - Fast A matrix size: %dx%d\n', size(A_fast));
fprintf('  Identified poles:\n');
disp(pole(sys_identified)');

%% 7. FREQUENCY RESPONSE COMPARISON
fprintf('\nStep 7: Comparing frequency responses...\n');

% Frequency vector
freq_Hz = logspace(-1, log10(Fs/2), 1000);
freq_rad = freq_Hz * 2 * pi;

% Compute Bode responses
[mag_true, ~, ~] = bode(sys_discrete, freq_rad);
[mag_id, phase_id, ~] = bode(sys_identified, freq_rad);

mag_true = 20*log10(squeeze(mag_true));
mag_id = 20*log10(squeeze(mag_id));
phase_id = squeeze(phase_id);

% Compute fit percentage
fit_percentage = 100 * (1 - norm(mag_true - mag_id) / norm(mag_true));

fprintf('  ✓ Bode plot generated\n');
fprintf('  - Magnitude fit: %.2f%%\n', fit_percentage);

%% 8. PLOTTING
fprintf('\nStep 8: Creating visualizations...\n');

figure('Position', [100, 100, 1200, 800], 'Name', 'Multirate System ID Results');

% Plot 1: Input signal
subplot(3,2,1);
t = (0:length(input_signal)-1) * Ts;
plot(t(1:5000), input_signal(1:5000), 'b', 'LineWidth', 1);
grid on;
xlabel('Time (s)');
ylabel('Input');
title('PRBS Input Signal (first 5s)');

% Plot 2: Sensor outputs
subplot(3,2,2);
t_s1 = (0:length(output_sensor1)-1) * n1 * Ts;
t_s2 = (0:length(output_sensor2)-1) * n2 * Ts;
plot(t_s1(1:500), output_sensor1(1:500), 'ro-', 'MarkerSize', 4); hold on;
plot(t_s2(1:500), output_sensor2(1:500), 'bs-', 'MarkerSize', 4);
grid on;
xlabel('Time (s)');
ylabel('Output');
title('Downsampled Sensor Outputs');
legend('Sensor 1 (n_1=3)', 'Sensor 2 (n_2=5)', 'Location', 'best');

% Plot 3 & 4: Bode magnitude
subplot(3,2,3);
semilogx(freq_Hz, mag_true, 'k-', 'LineWidth', 2); hold on;
semilogx(freq_Hz, mag_id, 'r--', 'LineWidth', 1.5);
xline(Fs/(2*n1), 'b:', 'LineWidth', 1.5, 'Label', sprintf('Nyq_{s1} %.0fHz', Fs/(2*n1)));
xline(Fs/(2*n2), 'g:', 'LineWidth', 1.5, 'Label', sprintf('Nyq_{s2} %.0fHz', Fs/(2*n2)));
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Response - Magnitude');
legend('True System', 'Identified', 'Location', 'best');
xlim([freq_Hz(1), freq_Hz(end)]);

subplot(3,2,4);
semilogx(freq_Hz, phase_id, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
title('Frequency Response - Phase');
xlim([freq_Hz(1), freq_Hz(end)]);

% Plot 5: Pole-zero map
subplot(3,2,5);
pzmap(sys_discrete, 'b'); hold on;
pzmap(sys_identified, 'r');
title('Pole-Zero Map');
legend('True', 'Identified', 'Location', 'best');

% Plot 6: Error analysis
subplot(3,2,6);
mag_error = abs(mag_true - mag_id);
semilogx(freq_Hz, mag_error, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude Error (dB)');
title(sprintf('Identification Error (Fit: %.1f%%)', fit_percentage));
xlim([freq_Hz(1), freq_Hz(end)]);

fprintf('  ✓ Plots created\n');

%% 9. SUMMARY
fprintf('\n');
fprintf('========================================================\n');
fprintf('  SUMMARY\n');
fprintf('========================================================\n');
fprintf('True System:\n');
fprintf('  - Order: %d\n', size(A_true, 1));
fprintf('  - Continuous poles: '); fprintf('%.2f ', real(pole(sys_cont)')); fprintf('\n');
fprintf('\nIdentified System:\n');
fprintf('  - Order: %d\n', n_est);
fprintf('  - Discrete poles: '); fprintf('%.4f ', abs(pole(sys_identified)')); fprintf('\n');
fprintf('\nPerformance:\n');
fprintf('  - Magnitude fit: %.2f%%\n', fit_percentage);
fprintf('  - Sensor 1 Nyquist: %.1f Hz\n', Fs/(2*n1));
fprintf('  - Sensor 2 Nyquist: %.1f Hz\n', Fs/(2*n2));
fprintf('  - Collaborative range: up to %.1f Hz\n', Fs/(2*m));
fprintf('\n✓ Test completed successfully!\n');
fprintf('========================================================\n\n');
