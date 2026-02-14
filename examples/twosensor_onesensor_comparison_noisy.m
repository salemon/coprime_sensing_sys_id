%Collabrative System ID by X.H Hu T.C Chu, X. Chen @UW 
%This demo use two sensors
clf, clear all, close all, clc
Fs = 1024;% Basic frequency for this system
% Define parameters for zero-order hold and sampling periods
m = 2;  % Zero-order hold for input is 2
n1 = 5; % Sampling period for output channel 1 is 3
n2 = 7; % Sampling period for output channel 2 is 5
n3 = 3;
n1n2 = n1*n2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what are these for?
n1m = n1*m;
n2m = n2*m;
n3m=n3*m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mn1n2=n1*n2*m;
%% Setup the system
% Define the continuous-time transfer function and state-space representation
GG = tf(1, [1 3 1 1]); % Real system without input delay

% n_GG =8; % system order %% change this for different order test
% %%%%%
% 
% n_k = n_GG+2; % order for n4sid
% n_GG_im = round(n_GG/6); % determine number of complex conjugate pairs
% n_GG_re = n_GG-n_GG_im*2; % number of real poles
% 
% % create random vector ranging from 0 - 1 to determine value of poles
% n_rand_re = -rand(n_GG_re,1); % random values for real poles
% n_rand_im1 = -rand(n_GG_im,1); % random value for re(conjugate poles)
% n_rand_im2 = rand(n_GG_im,1); % random value for im(conjugate poles)
% n_rand_zero = -rand(round(n_GG/2),1); % random value for zeros
% pole_re_dom = 30; % range of pole will be -20 to -50 (range = 30)
% pole_im_dom = 30; % range of zeros from 0 to -30 (range = 30)
% pole_start = -20; % starting maximum value of the pole, -20
% pole_im_start = 5; % starting point for imaginary part
% % assigning pole and zero values
% zero_re = n_rand_zero*pole_re_dom; % generating zeros
% pole_re = n_rand_re*pole_re_dom+pole_start; % generating poles
% pole_im = n_rand_im1*pole_im_dom+(pole_im_start+n_rand_im2*pole_im_dom)*j; % generaing complex conjugates
% pole_im = [pole_im; conj(pole_im)]; % find the conjugates
% den = poly([pole_re', pole_im']); % coefficient of denoninator
% num = poly([zero_re']); % coefficient of numerator
% 
% GG = tf(num, den); % Real system without input delay      
% a_dt = c2d(GG, 1/Fs); % Discretize the system
% fprintf('Is discretized sys stable? %.2f',isstable(a_dt))

%% === Generate PRBS signal as input ===
% Set display format for more precision in output
format long g

% Initialization of PRBS (Pseudo Random Binary Sequence) parameters
ValUinit = 0; % Initial value for the PRBS generation
ValAmpli = 1; % Amplitude of the PRBS signal
ValDecal = 0; % DC offset for the PRBS signal
ValLgReg = 10; % Length of the shift register used in PRBS generation
ValDivi = 1; % Frequency divider for controlling the PRBS generation rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% added m for nsamp
% Nsamp = 20*n_mult*2^ValLgReg; % Total number of samples to generate, based on the register length
Nsamp = 200*m*2^ValLgReg; % Total number of samples to generate, based on the register length
Tappli = 0; % when PRBS starts

% Generate the PRBS signal with specified parameters
prbs = create_prbs(ValUinit, ValAmpli, ValDecal, ValLgReg, ValDivi, Nsamp, Tappli)'; % transposed as well

%% === Simulation Preparation ===
x = prbs;             % simulation input
x_noisy=awgn(x,20,'measured');
t_sim = 0:1/Fs:(length(x_noisy)-1)/Fs;    % simulation time vector
t_sim = t_sim(:);
Ts_zoh = m * 1 / Fs;  % ZOH time period
StopTime = t_sim(end);
simdata = sim('models_multirate_data_example_2022b_noisy');% make sure that simulink is well named and is in the same path
figure(1)
plot(x_noisy)
hold on
plot(x)
legend('x_noisy','x')

% === Zero-Order Hold (ZOH) Processing === for input x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% for t_zoh, should the time be m/Fs instead? %%%%%%%
t_zoh = 0:1/Fs:t_sim(end);  % ZOH time vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_zoh = simdata.x_zoh_out; %Obtain the ZOH-processed input signal from simulation results.

% === Signal Trimming and Decimation ===

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signal trimming mean we pick a range?
div_mn1n2 = round(Nsamp/(4*mn1n2)); % pick 4 as the divisor
% t_strt = 400*n_mult;
% t_end = 600*n_mult;
t_strt = round(2*div_mn1n2); % guarentee this will always fall in the middle
t_end = round(3*div_mn1n2);
% Adjust the ZOH time vector to focus on a specific interval of interest, defined by the variables 'n1m', '400', and '600'.
% t_zoh = t_zoh(:,(mn1n2*t_strt+1):(mn1n2*t_end)); % Trim the ZOH time vector to a specific interval.
% x_zoh = x_zoh((mn1n2*t_strt+1):m:(mn1n2*t_end)); % Decimate the ZOH-processed signal within the specified interval.
t_zoh_twosensor = t_zoh(:,(mn1n2*t_strt+1):(mn1n2*t_end)); 
x_zoh_twosensor = x_zoh((mn1n2*t_strt+1):m:(mn1n2*t_end)); 

x_zoh_twosensor = x_zoh_twosensor - mean(x_zoh_twosensor); % Subtract the mean from the decimated ZOH signal to normalize it.

% Extract simulated output data
y_prbs_8_1_scal1 = simdata.y;
y_prbs_8_1_scal1 = y_prbs_8_1_scal1 -mean(y_prbs_8_1_scal1);%subtract sample means
y1 = simdata.y1;% Extract secondary output data (y1).
y1=y1((n2m*t_strt+1):n2m*t_end);
y1 = y1-mean(y_prbs_8_1_scal1);
y2 = simdata.y2;% Extract secondary output data (y2).
y2 =y2((n1m*t_strt+1):n1m*t_end);
y2 = y2-mean(y_prbs_8_1_scal1);

% Refine the time vector 't' to focus on a specific range of interest. This range is determined by the variable 'n1m' and the specific indices calculated.
t = t_sim((mn1n2*t_strt+1):mn1n2*t_end);

% Generate refined time vectors for 'y1' and 'y2' outputs, taking into account their respective sampling rates.
t_y1 = t(1:n1:end);  % Create a time vector for 'y1', downsampling by factor 'n1' to match its sampling rate.
t_y2 = t(1:n2:end);  % Create a time vector for 'y2', downsampling by factor 'n2' to match its sampling rate.


% === Align Dimensions for Lifting Processing ===
% Determine the quotient and remainder for resizing operations, ensuring that the data aligns correctly for matrix operations.
Q1 = fix(length(x_zoh_twosensor) / n1n2);
R1 = mod(length(x_zoh_twosensor), n1n2);
Q2 = fix(length(y1) / n2m);
R2 = mod(length(y1), n2m);
Q3 = fix(length(y2) / n1m);
R3 = mod(length(y2), n1m);

% Reshape the input and output vectors to match the dimensions required for further analysis or processing, using MATLAB's reshape function.
U_lifted = reshape(x_zoh_twosensor(1:end-R1,:), 1, n1n2, []); % Reshape 'x_zoh', excluding the remainder to fit the 'n1' dimension.
Y1_lifted = reshape(y1(1:end-R2,:), 1, n2m, []); % Reshape 'y1', excluding the remainder to fit the 'm' dimension.
Y2_lifted = reshape(y2(1:end-R3,:), 1, n1m,[]);
clear U_lifted2 Y_stack2
% Flatten the first dimension of 'U_lifted' and 'Y_stack' to simplify their structure for analysis.
Y_stack = cat(2,Y1_lifted,Y2_lifted);
U_lifted2(:,:) = U_lifted(1,:,:); % Flatten 'U_lifted' by removing the singleton first dimension.
Y_stack2(:,:) = Y_stack(1,:,:); % Flatten 'Y_stack' similarly.

% Convert the matrices 'Y_stack2' and 'U_lifted2' into cell arrays, allowing for more flexible data manipulation and access.
yy = mat2cell(Y_stack, size(Y_stack,1), size(Y_stack,2), ones(1, size(Y_stack,3))); 
uu = mat2cell(U_lifted, size(U_lifted,1), size(U_lifted,2), ones(1, size(U_lifted,3)));

% Transpose 'Y_stack2' and 'U_lifted2' to match expected dimensions for subsequent processing.
yy = Y_stack2'; 
uu = U_lifted2';


% Execute the N4SID (Numerical Algorithm for Subspace State Space System Identification) method with the provided input and output data.
% The function returns state-space matrices (A, B, C, D) among other outputs, for a system order defined by k=5.
[aaa, bbb, ccc, ddd, xf, nn, Uf, Yf] = n4sid_test(uu, yy, 5);

% Set the estimated system order.
nx = nn;

% Trim the input (uu) and output (yy) data by removing the first and last 5 samples to align with the system dynamics explored.
uu = uu(5:end-5,:);
yy = yy(5:end-5,:);

% Prepare the initial state vector for simulation by using the final states estimated from N4SID.
% This is transposed to match the expected dimensionality for state-space functions.
x_lsim = [ xf(:,1:end-1) ]';

% Determine the total number of integration steps based on the trimmed state vector length.
number_intotal = length(x_lsim);

% Prepare matrices for constructing the estimated state-space model.
A_B_left = zeros((nx), (nx + n1n2));
C_0 = zeros(1, nx); % Initial guess for C matrix, assuming no prior information.

% Extract the last and second to last columns of B matrix from N4SID results to use in eigenvector approach.
B_m = bbb(1:nx, end);
B_m2 = bbb(1:nx, end-1);

% Compute the A matrix via an eigenvector approach for system identification.
[u, lamda] = eig(aaa); % Eigendecomposition of A matrix from N4SID.
u_inv = inv(u); % Inverse of eigenvector matrix.
u_inv_B_n_minus_1 = u_inv * B_m;
u_inv_B_n_minus_2 = u_inv * B_m2;
lamda_matrix = u_inv_B_n_minus_2 ./ u_inv_B_n_minus_1; % Calculate the lambda matrix.
A_m = u * diag(lamda_matrix) * u_inv; % Reconstruct A matrix using lambda matrix.
A = real(sqrtm(A_m)); % Ensure A matrix is real and apply square root for stabilization.
B = real(inv(eye(nx) + A) * B_m); % Calculate B matrix and ensure it's real.


C = C_0;

% Construct the estimated continuous-time state-space system model with identified A, B, and C matrices.
sys_est = ss(A, B, ccc(1, 1:end), [], 1/Fs);

% Convert the estimated state-space model to a transfer function for easier comparison and analysis.
[Num, Den] = ss2tf(sys_est.A, sys_est.B, sys_est.C, sys_est.D);
G_est = tf(Num, Den, 1/Fs); % The estimated transfer function model.

% Discretize the actual continuous-time system model 'GG' using zero-order hold (ZOH) for a fair comparison.
GGd = c2d(GG, 1/Fs, 'zoh');






% %% === Simulation Preparation ===
% % 
% % t_sim = 0:1/Fs:(length(x_noisy)-1)/Fs;    % simulation time vector
% % t_sim = t_sim(:);
% % Ts_zoh = m * 1 / Fs;  % ZOH time period
% % StopTime = t_sim(end);
% % simdata = sim('models_multirate_data_example_2022b_noisy');% make sure that simulink is well named and is in the same path
% 
% 
% % === Zero-Order Hold (ZOH) Processing === for input x
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%% for t_zoh, should the time be m/Fs instead? %%%%%%%
% t_zoh = 0:1/Fs:t_sim(end);  % ZOH time vector
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% x_zoh = simdata.x_zoh_out; %Obtain the ZOH-processed input signal from simulation results.
% 
% % === Signal Trimming and Decimation ===

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signal trimming mean we pick a range?
div_mn3 = round(Nsamp/(4*n3m)); % pick 4 as the divisor
% t_strt = 400*n_mult;
% t_end = 600*n_mult;
t_strt = round(2*div_mn3); % guarentee this will always fall in the middle
t_end = round(3*div_mn3);
% Adjust the ZOH time vector to focus on a specific interval of interest, defined by the variables 'n1m', '400', and '600'.
% t_zoh = t_zoh(:,(mn1n2*t_strt+1):(mn1n2*t_end)); % Trim the ZOH time vector to a specific interval.
% x_zoh = x_zoh((mn1n2*t_strt+1):m:(mn1n2*t_end)); % Decimate the ZOH-processed signal within the specified interval.
t_zoh_onesensor = t_zoh(:,(n3m*t_strt+1):(n3m*t_end)); 
x_zoh_onesensor = x_zoh((n3m*t_strt+1):m:(n3m*t_end)); 

x_zoh_onesensor = x_zoh_onesensor - mean(x_zoh_onesensor); % Subtract the mean from the decimated ZOH signal to normalize it.

% % Extract simulated output data
% y_prbs_8_1_scal1 = simdata.y;
% y_prbs_8_1_scal1 = y_prbs_8_1_scal1 -mean(y_prbs_8_1_scal1);%subtract sample means
y3 = simdata.y3;% Extract secondary output data (y1).
y3=y3((m*t_strt+1):m*t_end);
y3 = y3-mean(y_prbs_8_1_scal1);


% Refine the time vector 't' to focus on a specific range of interest. This range is determined by the variable 'n1m' and the specific indices calculated.
t = t_sim((n3m*t_strt+1):n3m*t_end);

% Generate refined time vectors for 'y1' and 'y2' outputs, taking into account their respective sampling rates.
t_y1 = t(1:n1:end);  % Create a time vector for 'y1', downsampling by factor 'n1' to match its sampling rate.
t_y2 = t(1:n2:end);  % Create a time vector for 'y2', downsampling by factor 'n2' to match its sampling rate.
t_y3 = t(1:n3:end);

% === Align Dimensions for Lifting Processing ===
% Determine the quotient and remainder for resizing operations, ensuring that the data aligns correctly for matrix operations.
Q3 = fix(length(x_zoh_onesensor) / n3);
R3 = mod(length(x_zoh_onesensor), n3);
Q2 = fix(length(y3) / m);
R2 = mod(length(y3), m);
% Q3 = fix(length(y2) / n1m);
% R3 = mod(length(y2), n1m);

% Reshape the input and output vectors to match the dimensions required for further analysis or processing, using MATLAB's reshape function.
U_lifted = reshape(x_zoh_onesensor(1:end-R3,:), 1, n3, []); % Reshape 'x_zoh', excluding the remainder to fit the 'n1' dimension.
Y1_lifted = reshape(y3(1:end-R2,:), 1, m, []); % Reshape 'y1', excluding the remainder to fit the 'm' dimension.
% Y2_lifted = reshape(y2(1:end-R3,:), 1, n1m,[]);
clear U_lifted2 Y_stack2

Y_stack = Y1_lifted;
% Flatten the first dimension of 'U_lifted' and 'Y_stack' to simplify their structure for analysis.
U_lifted2(:,:) = U_lifted(1,:,:); % Flatten 'U_lifted' by removing the singleton first dimension.
Y_stack2(:,:) = Y_stack(1,:,:); % Flatten 'Y_stack' similarly.

% Convert the matrices 'Y_stack2' and 'U_lifted2' into cell arrays, allowing for more flexible data manipulation and access.
yy = mat2cell(Y_stack, size(Y_stack,1), size(Y_stack,2), ones(1, size(Y_stack,3))); 
uu = mat2cell(U_lifted, size(U_lifted,1), size(U_lifted,2), ones(1, size(U_lifted,3)));

% Transpose 'Y_stack2' and 'U_lifted2' to match expected dimensions for subsequent processing.
yy = Y_stack2'; 
uu = U_lifted2';


% Execute the N4SID (Numerical Algorithm for Subspace State Space System Identification) method with the provided input and output data.
% The function returns state-space matrices (A, B, C, D) among other outputs, for a system order defined by k=5.
[aaa, bbb, ccc, ddd, xf, nn, Uf, Yf] = n4sid_test(uu, yy, 5);

% Set the estimated system order.
nx = nn;

% Trim the input (uu) and output (yy) data by removing the first and last 5 samples to align with the system dynamics explored.
uu = uu(5:end-5,:);
yy = yy(5:end-5,:);

% Prepare the initial state vector for simulation by using the final states estimated from N4SID.
% This is transposed to match the expected dimensionality for state-space functions.
x_lsim = [ xf(:,1:end-1) ]';

% Determine the total number of integration steps based on the trimmed state vector length.
number_intotal = length(x_lsim);

% Prepare matrices for constructing the estimated state-space model.
A_B_left = zeros((nx), (nx + n1n2));
C_0 = zeros(1, nx); % Initial guess for C matrix, assuming no prior information.

% Extract the last and second to last columns of B matrix from N4SID results to use in eigenvector approach.
B_m = bbb(1:nx, end);
B_m2 = bbb(1:nx, end-1);

% Compute the A matrix via an eigenvector approach for system identification.
[u, lamda] = eig(aaa); % Eigendecomposition of A matrix from N4SID.
u_inv = inv(u); % Inverse of eigenvector matrix.
u_inv_B_n_minus_1 = u_inv * B_m;
u_inv_B_n_minus_2 = u_inv * B_m2;
lamda_matrix = u_inv_B_n_minus_2 ./ u_inv_B_n_minus_1; % Calculate the lambda matrix.
A_m = u * diag(lamda_matrix) * u_inv; % Reconstruct A matrix using lambda matrix.
A = real(sqrtm(A_m)); % Ensure A matrix is real and apply square root for stabilization.
B = real(inv(eye(nx) + A) * B_m); % Calculate B matrix and ensure it's real.


C = C_0;

% Construct the estimated continuous-time state-space system model with identified A, B, and C matrices.
sys_est_onesensor = ss(A, B, ccc(1, 1:end), [], 1/Fs);

% Convert the estimated state-space model to a transfer function for easier comparison and analysis.
[Num, Den] = ss2tf(sys_est_onesensor.A, sys_est_onesensor.B, sys_est_onesensor.C, sys_est_onesensor.D);
G_est_onesensor = tf(Num, Den, 1/Fs); % The estimated transfer function model.



% Plot the Bode plot for both the actual discretized system and the estimated system for comparison.
figure;
bode(GGd, 'r', sys_est_onesensor,'g',sys_est, 'b--');
legend('G_actual','G_onesensor', 'G_estimated');
