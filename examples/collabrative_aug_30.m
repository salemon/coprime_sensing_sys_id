
% Collaborative Sensing Example
% This code is associated with the file "models_multirate_data_example_2022b.slx", which contains Simulink blocks.
% To use this code, make sure to add the "macs-matlab-toolbox-master" to the Matlab path.
% Authors: Bob Xiaohai Hu (huxh@uw.edu), X. Chen (chx@uw.edu), UW Mechanical Eng department

clf, clear all, close all
Fs = 1024;
% Define parameters for zero-order hold and sampling periods
m = 2;  % Zero-order hold for input is 2
n1 = 3; % Sampling period for output channel 1 is 3
n2 = 5; % Sampling period for output channel 2 is 5

%% Setup the system
% Define the continuous-time transfer function and state-space representation
GG = tf(1, [20 1]); % Real system without input delay
[aa, bb, cc, dd] = tf2ss(1, [20 1]); 
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
ValLgReg = 8;
ValDivi = 2;
Nsamp = 200*2^ValLgReg;
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
return
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

t_zoh = t_zoh(:,25601:end);
x_zoh = x_zoh(25601:2:end);
x_zoh = x_zoh -mean(x_zoh);%subtract sample means
%return
% ytest=ytest(12801:end);
% ytest=ytest-mean(ytest);
% Simulate the system using Simulink and collect the output data


% Extract simulated output data
y_prbs_8_1_scal1 = simdata.y;
y_prbs_8_1_scal1 = y_prbs_8_1_scal1(25601:end);
% return
y_prbs_8_1_scal1 = y_prbs_8_1_scal1 -mean(y_prbs_8_1_scal1);%subtract sample means
y1 = simdata.y1;

y1=y1(8534:end);
y1 = y1-mean(y_prbs_8_1_scal1);
y2 = simdata.y2;
y2 =y2(5121:end);
y2 = y2-mean(y_prbs_8_1_scal1);
%return
% Plot the output to check how much should be segmented
figure;

% Plot ytest and y_prbs_8_1_scal1
subplot(2,1,1);
%plot(t_zoh, ytest, 'b', 'DisplayName', 'ytest');
hold on;
t=t(25601:end);
plot(t, y_prbs_8_1_scal1, 'r', 'DisplayName', 'y\_prbs\_8\_1\_scal1');
hold off;
xlabel('Time');
ylabel('Output');
title('Continuous-Time System Output Comparison');
legend;
grid on;

% Plot y1 and y2
subplot(2,1,2);
t_y1 = t(1:n1:end);  % Generate time vector for y1 based on its sampling rate
t_y2 = t(1:n2:end);  % Generate time vector for y2 based on its sampling rate

plot(t_y1, y1, 'b', 'DisplayName', 'y1');
hold on;
plot(t_y2, y2, 'r', 'DisplayName', 'y2');
hold off;
xlabel('Time');
ylabel('Output');
title('Sampled Output Channels Comparison');
legend;
grid on;

sgtitle('Output Comparison');
% return

n1n2 = n1*n2;
n1m = n1*m;
n2m = n2*m;
% align the dimensions
Q1 = fix(length(x_zoh) / n1n2);
R1 = mod(length(x_zoh), n1n2);
Q2 = fix(length(y1) / n2m);
R2 = mod(length(y1), n2m);
Q3 = fix(length(y2) / n1m);
R3 = mod(length(y2), n1m);

%use matlab reshape function for lift operator representation
% x_zoh=x_zoh';
U_lifted = reshape(x_zoh(1:end-R1,:),1,n1n2,[]);
Y1_lifted = reshape(y1(1:end-R2,:),1,n2m,[]);
Y2_lifted = reshape(y2(1:end-R3,:), 1, n1m,[]);

% stack Y matrix,  Concatenate along the second dimension
 Y_stack = cat(2,Y1_lifted,Y2_lifted);
 U_lifted2(:,:)= U_lifted(1,:,:);
 Y_stack2(:,:) = Y_stack(1,:,:);

yy =  mat2cell(Y_stack, size(Y_stack,1), size(Y_stack,2), ones(1,size(Y_stack,3))); 
uu = mat2cell(U_lifted, size(U_lifted,1), size(U_lifted,2), ones(1,size(U_lifted,3)));
yy = Y_stack2'; uu = U_lifted2';

%return
%%
num_samples_uu = n1n2;
num_samples_yy = m*(n1+n2);
%return
% Create the timetable 'tt2' using the 'timetable' function with 'uu1' column
% tt2 = timetable(uu(:, 1), 'SampleRate', Fs, 'VariableNames', {'uu1'});
tt2 = timetable(uu(:, 1), 'SampleRate', Fs/(n1n2*m), 'VariableNames', {'uu1'},'StartTime',seconds(25));
%tt2 = timetable(uu(:, 1), 'SampleRate', Fs/(n1n2*m), 'VariableNames', {'uu1'});
% Loop to add 'uu2' to 'uu140' columns to the timetable
for i = 2:num_samples_uu
    column_name = sprintf('uu%d', i);
    tt2 = addvars(tt2, uu(:, i), 'NewVariableNames', {column_name});
end
column_names = cell(1, num_samples_uu); % Initialize a cell array to store the column names
for i = 1:num_samples_uu
    column_names_uu{i} = sprintf('uu%d', i); % Assign column names to each element of the cell array
end

for i = 1:num_samples_yy
    column_name2 = sprintf('yy%d', i);
    tt2 = addvars(tt2, yy(:, i), 'NewVariableNames', {column_name2});
end
for i = 1:num_samples_yy
    column_names_yy{i} = sprintf('yy%d', i); % Assign column names to each element of the cell array
end

%% Use n4sid to estimate the system model

% Estimate the system model using 'n4sid' with all variable names
inputNames = {};
outputNames = {};

% Loop from 1 to 140
for i = 1:m*n1n2
    % Construct the string "uuX" and "yyX" for each iteration
    uu_str = strcat('uu', num2str(i));
    yy_str = strcat('yy', num2str(i));
    
    % Add the co structed strings to the respective arrays
    inputNames = [inputNames, uu_str];
    outputNames = [outputNames, yy_str];
end
% Define the numeric ranges for uu and yy
uu_range = 1:n1n2;
yy_range = 1:m*(n1+n2);

% Create string arrays for input and output names
input_names = "uu" + string(uu_range);
output_names = "yy" + string(yy_range);
nx =3
% opt = n4sidOptions('N4Weight','MOESP','Display','on');
% [sys,x0] = n4sid(tt2,nx,'InputName',[input_names ],'OutputName',[output_names],opt)

 % simulate to extract  the time evolution of the state values in response to the input signal
number_intotal = floor(Fs/(n1n2*m) * ((length(x)-1)/Fs)/2);
t_index_M = 0:1/Fs*(n1n2*m):(length(x)-1)/Fs;
t_index_M = StopTime/2:1/Fs*(n1n2*m):StopTime;
t_index_M = t_index_M(1:number_intotal)';
[aaa bbb ccc ddd xf nn Uf Yf] = n4sidkatamodar(uu,yy,7);
nx = nn;
%[y_lsim, tt_lsim, x_lsim] = lsim(sys,uu,t_index_M,x0);
x_lsim = [ xf(:,1:end-1) ]';
%x_lsim = [x0  xf(:,1:end-1)]';
number_intotal  = length(x_lsim);
%return
% prepare the least square method parameters
for i  = 1:(number_intotal)
    x_u_stack(i,:) =  [x_lsim(i,:) uu(i,:)];
end

%phi
phi = zeros((nx+n1n2),(nx+n1n2));
phi_0 = zeros(nx,nx);
for i =1:(number_intotal)
%     phi = phi+x_u_stack(i,:)*x_u_stack(i,:)';
     phi = phi + (x_u_stack(i,:)')*x_u_stack(i,:);
     phi_0 = phi_0 + (x_lsim(i,:)')*x_lsim(i,:);
end
phi = (1/(number_intotal)).*phi;
phi_0 = (1/(number_intotal)).*phi_0;

A_B_left = zeros((nx),(nx+n1n2));
C_0 = zeros(1,nx);
% x_kplus1= x_lsim(2:end,:);
%x_kplus1 = [ x_lsim(1:8,:)' xf(:,1:end) ]';
x_kplus1 = xf(:,2:end)';
uu = uu(5:end-5,:);
yy = yy(5:end-5,:);
%return
for i = 1:(number_intotal)
    A_B_left = A_B_left + [x_kplus1(i,:)'*x_lsim(i,:) x_kplus1(i,:)'*uu(i,:) ];
    C_0 = C_0+ yy(i,1)*x_lsim(i,:);
end
display(C_0)
A_B_left = 1/(number_intotal).*A_B_left*inv(phi);
C_0 = 1/(number_intotal) .*C_0*(inv(phi_0)) ;
A_stack= A_B_left(1:nx,1:nx);
B_m = A_B_left(1:nx,end);

 %compute A_stack via eigenvector approach
[u, lamda] = eig(A_stack);
u_inv = inv(u);
u_inv_B_n_minus_1 = u_inv*B_m;
u_inv_B_n_minus_2 = u_inv*(A_B_left(1:nx,end-1));
lamda_matrix = u_inv_B_n_minus_2 ./ u_inv_B_n_minus_1;
A_m = u*diag(lamda_matrix)*u_inv;
A = A_m;
A = round(sqrtm(A));
A = real(A_m);
B = inv(eye(nx)+A)*B_m;
B = real(B_m);
C = C_0;
sys_est = ss(A,B,C,[],-1);
[Num Den] = ss2tf(sys_est.A,sys_est.B,sys_est.C,sys_est.D);
G_est = tf(Num,Den)
GGd = c2d(GG,1/Fs,'zoh');
figure
bode(GGd,'r',sys_est,'b--')
%bode(GGd,'r',G_est,'b--')
legend('G_actual','G_estimated')
return
% another method to recover fast single rate model
% tau depends on the dimension of A_m
B_m = A_B_left(1:nx,end);
tau_c = [B_m A_B_left(1:nx,end-1) A_B_left(1:nx,end-2)];
tau = [A_B_left(1:nx,end-1) A_B_left(1:nx,end-2) A_B_left(1:nx,end-3)];
A_m_2 = tau*tau_c'*inv(tau_c*tau_c');
A2 = real(A_m_2);
A2 = real(A2);
B = inv(eye(nx)+A2)*B_m;
B = real(B_m);
C = C_0;
sys_est2 = ss(A2,B,C,[],1/Fs)
GGd = c2d(GG,1/Fs,'zoh');
w = linspace(0,10*pi,128);
figure
bode(GG,'r',sys_est2,'b--')
legend('G_actual','G_estimated')
return