

% Collaborative Sensing Example
% This code is associated with the file "models_multirate_data_example_2022b.slx", which contains Simulink blocks.
% To use this code, make sure to add the "macs-matlab-toolbox-master" to the Matlab path.
% Authors: Bob Xiaohai Hu (huxh@uw.edu), Thomas Chu(tchu@uw.edu),X. Chen (chx@uw.edu), UW Mechanical Eng department

clf, clear all, close all
Fs = 1024;
% Define parameters for zero-order hold and sampling periods
m = 2;  % Zero-order hold for input is 2
n1 = 3; % Sampling period for output channel 1 is 3
n2 = 5; % Sampling period for output channel 2 is 5
n1n2 = n1*n2;
n1m = n1*m;
n2m = n2*m;
mn1n2=n1*n2*m;
%% Setup the system
% Define the continuous-time transfer function and state-space representation
GG = tf(1, [1 3 1 1]); % Real system without input delay
[aa, bb, cc, dd] = tf2ss(1, [1 3 1 1]); 
a_dt = c2d(GG, 1/Fs); % Discretize the system
pole(a_dt);
%sys1 = ss(aa, bb, cc, dd,1/Fs);

%% Generate PRBS signal as input
format long g
flag_spec = 0;
flag_fprint = 1;

ValUinit = 0;
ValAmpli = 1;
ValDecal = 0;
%ValDecal = 0;
ValLgReg = 10;
ValDivi = 1;
Nsamp = 20*2^ValLgReg;
Tappli = 0;


%  "Entry parameters" are :
% 	  ValUinit  : Initial steady state
%     ValAmpli  : Magnitude
%     ValDecal  : Add-on DC component
%     ValLgReg  : Register length
%     ValDivi   : Frequency divider
%     Nsamp     : Number of samples
%     Tappli    : Application instant

prbs = create_prbs(ValUinit, ValAmpli, ValDecal, ValLgReg, ValDivi, Nsamp, Tappli);
prbs = prbs';

figure, plot(prbs)

if 0
    %%
    specPRBs = specCale(prbs,Fs);           % calculate the signal spectrum 
    figure, plot(specPRBs.f,specPRBs.amp)
end

% different scaling of the random signal
prbs_8_1 = prbs;
prbs_8_1_scal1 = prbs;
% return
% start with a small amplitude for safety, gradually incease the amplitude
% in practice to find the best signal-to-noise ratio
%% Simulation
x = prbs_8_1_scal1;             % simulation input (real)
t = 0:1/Fs:(length(x)-1)/Fs;    % simulation time
t = t(:);
StopTime = t(end);
simdata = sim('models_multirate_data_example_2022b');
% Simulation
% Apply zero-order hold (ZOH) to the input signal x
Ts_zoh = m * 1 / Fs;  % ZOH time period
t_zoh = 0:1/Fs:t(end);  % ZOH time vector
% x_zoh = interp1(t, x, t_zoh, 'previous');  % Apply ZOH to the input signal
x_zoh = simdata.x_zoh_out;
% Simulate the continuous-time system using ZOH input
%ytest = lsim(GG, x_zoh, t_zoh);
%return
t_zoh = t_zoh(:,(n1m*400+1):(n1m*600));
x_zoh = x_zoh((n1m*400+1):2:(n1m*600));
x_zoh = x_zoh -mean(x_zoh);%subtract sample means
%return
% ytest=ytest(12801:end);
% ytest=ytest-mean(ytest);
% Simulate the system using Simulink and collect the output data


% Extract simulated output data
y_prbs_8_1_scal1 = simdata.y;
%y_prbs_8_1_scal1 = y_prbs_8_1_scal1(25601:end);

y_prbs_8_1_scal1 = y_prbs_8_1_scal1 -mean(y_prbs_8_1_scal1);%subtract sample means
y1 = simdata.y1;
y1=y1((m*400+1):m*600);
y1 = y1-mean(y_prbs_8_1_scal1);
y2 = simdata.y2;
y2 =y2((m*400+1):m*600);
y2 = y2-mean(y_prbs_8_1_scal1);
% return
% Plot the output to check how much should be segmented

% Plot ytest and y_prbs_8_1_scal1
% subplot(2,1,1);
% %plot(t_zoh, ytest, 'b', 'DisplayName', 'ytest');
% hold on;
t=t((n1m*400+1):n1m*600);
% plot(t, y_prbs_8_1_scal1, 'r', 'DisplayName', 'y\_prbs\_8\_1\_scal1');
% hold off;
% xlabel('Time');
% ylabel('Output');
% title('Continuous-Time System Output Comparison');
% legend;
% grid on;

% Plot y1 and y2
% subplot(2,1,2);
figure,
t_y1 = t(1:n1:end);  % Generate time vector for y1 based on its sampling rate
t_y2 = t(1:n2:end);  % Generate time vector for y2 based on its sampling rate

plot(t_y1, y1, 'b', 'DisplayName', 'y1');
% hold on;
% plot(t_y2, y2, 'r', 'DisplayName', 'y2');
% hold off;
xlabel('Time');
ylabel('Output');
title('Sampled Output Channels Comparison');
legend;
grid on;

% sgtitle('Output Comparison');
% return


% align the dimensions
Q1 = fix(length(x_zoh) / n1);
R1 = mod(length(x_zoh), n1);
Q2 = fix(length(y1) / m);
R2 = mod(length(y1), m);
% Q3 = fix(length(y2) / n1m);
% R3 = mod(length(y2), n1m);
%use matlab reshape function for lift operator representation
% x_zoh=x_zoh';
U_lifted = reshape(x_zoh(1:end-R1,:),1,n1,[]);
Y1_lifted = reshape(y1(1:end-R2,:),1,m,[]);
% Y2_lifted = reshape(y2(1:end-R3,:), 1, n1m,[]);

% stack Y matrix,  Concatenate along the second dimension
%  Y_stack = cat(2,Y1_lifted,Y2_lifted);
Y_stack = Y1_lifted;
 U_lifted2(:,:)= U_lifted(1,:,:);
 Y_stack2(:,:) = Y_stack(1,:,:);

yy =  mat2cell(Y_stack, size(Y_stack,1), size(Y_stack,2), ones(1,size(Y_stack,3))); 
uu = mat2cell(U_lifted, size(U_lifted,1), size(U_lifted,2), ones(1,size(U_lifted,3)));
yy = Y_stack2'; uu = U_lifted2';

% return
[aaa bbb ccc ddd xf nn Uf Yf] = n4sidkatamodar(uu,yy,5);
nx = nn;
%[y_lsim, tt_lsim, x_lsim] = lsim(sys,uu,t_index_M,x0);
x_lsim = [ xf(:,1:end-1) ]';
%x_lsim = [x0  xf(:,1:end-1)]';
number_intotal  = length(x_lsim);
uu = uu(5:end-5,:);
yy = yy(5:end-5,:);


A_B_left = zeros((nx),(nx+n1n2));
C_0 = zeros(1,nx);

B_m = bbb(1:nx,end);
B_m2= bbb(1:nx,end-1);
 %compute A_stack via eigenvector approach
[u, lamda] = eig(aaa);
u_inv = inv(u);
u_inv_B_n_minus_1 = u_inv*B_m;
u_inv_B_n_minus_2 = u_inv*B_m2;
lamda_matrix = u_inv_B_n_minus_2 ./ u_inv_B_n_minus_1;
A_m = u*diag(lamda_matrix)*u_inv;
% A = round(sqrtm(A_m));
A = real(sqrtm(A_m));
B = inv(eye(nx)+A)*B_m;
B = real(B);
C = C_0;
sys_est = ss(A,B,ccc(1,1:end),[],1/Fs);
% sys_est = ss(aaa,bbb(1,1),ccc(1,1),[],1/Fs)
[Num Den] = ss2tf(sys_est.A,sys_est.B,sys_est.C,sys_est.D);
G_est = tf(Num,Den,1/Fs)
GGd = c2d(GG,1/Fs,'zoh');
figure
bode(GGd,'r',sys_est,'b--')
%bode(GGd,'r',G_est,'b--')

legend('G_actual','G_estimated')

if 0
    %% Validate the equation: x(k+1) = Ax(k) + Bu(k)
    x_k_1_kk = zeros(size(x_lsim));
    time_vector = 1:length(x_lsim); % Generate a time vector (assuming unit time steps)
    for i = 1:length(x_lsim)
        x_k_1_kk(i,:) = A_stack * x_lsim(i,:) + B_stack * uu(i,:)';
        predict_error_x = x_kplus1 - x_k_1_kk;
    end
    
 % Plot the prediction error and x_k_1_kk against time
    figure;
%     subplot(2, 1, 1);
    plot(time_vector, predict_error_x, 'b', 'DisplayName', 'Prediction Error');
    hold on;
    plot(time_vector, x_k_1_kk, 'r', 'DisplayName', 'x(k+1)');
    xlabel('Time');
    ylabel('Values');
    title('Prediction Error and x(k+1) vs. Time');
    legend('Location', 'best');
    grid on;
    hold off;

    %Validate y = CX + D

%     y_hat = sys_est.C * x_lsim + sys_est.D * uu; % Compute y_hat using the estimated system
%     y_residual = yy - y_hat; % Compute the residual (measurement - prediction)
% 
%     subplot(2, 1, 2);
%     plot(time_vector, y_residual, 'g', 'DisplayName', 'Residual (y - y_{hat})');
%     xlabel('Time');
%     ylabel('Residual');
%     title('Residual (Measurement - Prediction) vs. Time');
%     legend('Location', 'best');
%     grid on;
end