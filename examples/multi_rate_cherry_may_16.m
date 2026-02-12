% Clear the command window, close all figures, and clear the workspace
clc; close all; clear;

% DC motor test using minimum segment length
% Load the sweep signal data
m_data = load('./Model Data/prbs_run_1.mat');

% ==============user change variables ============================
% n4SID orders
nx = 4; % order guess for singlerate n4SID
n_mult = 4; % order guess for the multirate SID
i_single = 1; % increase to make data sparser for single-rate n4sid

% PRBS time evaluation for n4SID
time_prbs = 50; % time span for using SID on prbs

% multirate parameters to change
m = 1;%zero-order hold
n1 = 2;%downsampling factor
n2 = 3;%downsampling factor channel 2
% ==============user change variables ============================

% collect data
in_V = m_data.in_voltage.signals.values;
out_E = m_data.out_encoder.signals.values;
in_V = squeeze(in_V);
out_E = squeeze(out_E);
in_V = double(in_V);
out_E = double(out_E);
% apply filter
out_E_f = filter([1,-1],1,out_E);
in_V_f = in_V;
% in_V_f = filter([1,-1],1,in_V);
% in_V_f = filter([1 -1],[1 0],in_V_f);
% Extract the time vector from the loaded data
tt = m_data.in_voltage.time;

% Calculate the sampling time (assuming uniform sampling)
Ts = tt(2) - tt(1);
Fs = 1/Ts; % Sampling frequency
tt_n = length(tt); % number of iterations
time_n = tt_n - time_prbs*Fs; % number of iterations within this time span
prbs_strt = ceil(time_n/2); % find the starting point
prbs_end = tt_n - floor(time_n/2); % find the ending point
% middle time span of the input and output prbs signal are extracted
input_prbs = in_V_f(prbs_strt:prbs_end);
output_prbs = out_E_f(prbs_strt:prbs_end);
% process the input and output signals (PRBS), system identification will use prbs signal
% remove the mean value
input_prbs = input_prbs - mean(input_prbs);
output_prbs = output_prbs - mean(output_prbs);

%% ============= SYSTEM IDENTIFICATION USING buildin N4SID ===============
DAT_prbs = iddata(output_prbs(1:i_single:end), input_prbs(1:i_single:end), i_single*Ts);
sys_identified_buildin_prbs = n4sid(DAT_prbs, nx);

%% ============= SYSTEM IDENTIFICATION USING multirate N4SID ===============
uu = input_prbs;
yy = output_prbs;
n1n2=n1*n2;n1m=n1*m;n2m=n2*m; mn1n2=m*n1n2;

%from channel 1
y1 = yy(1:n1:end);
y2 = yy(1:n2:end);

uu_zoh = uu(1:m:end);

% aligning dimensions for lifting process
Q1 = fix(length(uu_zoh) / n1n2);
R1 = mod(length(uu_zoh), n1n2);
Q2 = fix(length(y1) / n2m);
R2 = mod(length(y1), n2m);
Q3 = fix(length(y2) / n1m);
R3 = mod(length(y2), n1m);


U_lifted = reshape(uu_zoh(1:end-R1,:), 1, n1n2, []); % Reshape 'x_zoh', excluding the remainder to fit the 'n1' dimension.
Y1_lifted = reshape(y1(1:end-R2,:), 1, n2m, []); % Reshape 'y1', excluding the remainder to fit the 'm' dimension.
Y2_lifted = reshape(y2(1:end-R3,:), 1, n1m,[]);
% Flatten the first dimension of 'U_lifted' and 'Y_stack' to simplify their structure for analysis.
Y_stack = cat(2,Y1_lifted,Y2_lifted);
U_lifted2(:,:) = U_lifted(1,:,:); % Flatten 'U_lifted' by removing the singleton first dimension.
Y_stack2(:,:) = Y_stack(1,:,:); % Flatten 'Y_stack' similarly.

yyy = Y_stack2';
uuu = U_lifted2';
% [aaa, bbb, ccc, ddd, xf, n] = n4sid_test(uuu, yyy, 6);
% [aaa,bbb,ccc,ddd,nx]=ort_test(uuu,yyy,7);

% make uuu yyy the same size
[sizea,sizeb] = size(uuu);
[sizec,sized] = size(yyy);
if sizea>sizec
    uuu = uuu(1:sizec,:);
elseif sizea<sizec
    yyy = yyy(1:sizea,:);
end

idata_multi = iddata(yyy,uuu,Ts);
sys_identified_buildin_multi = n4sid(idata_multi, n_mult);

% Extract the A, B, and C matrices from the N4SID output.
aaa = sys_identified_buildin_multi.A;
bbb = sys_identified_buildin_multi.B;
ccc = sys_identified_buildin_multi.C;
ddd = sys_identified_buildin_multi.D;
B_m = bbb(1:n_mult,end);
B_m2 = bbb(1:n_mult,end-1);


% Compute the A matrix via an eigenvector approach for system identification.
[u, lamda] = eig(aaa); % Eigendecomposition of A matrix from N4SID.
u_inv = inv(u); % Inverse of eigenvector matrix.
u_inv_B_n_minus_1 = u_inv * B_m;
u_inv_B_n_minus_2 = u_inv * B_m2;
lamda_matrix = u_inv_B_n_minus_2 ./ u_inv_B_n_minus_1; % Calculate the lambda matrix.
A_m = u * diag(lamda_matrix) * u_inv; % Reconstruct A matrix using lambda matrix.
% A = real(sqrtm(A_m)); % Ensure A matrix is real and apply square root for stabilization.
A =real(A_m);

%B = real(inv(eye(nx) + A) * B_m); % Calculate B matrix and ensure it's real.
B = real(inv(A+A^0) * B_m);
B= B_m;
% Construct the estimated continuous-time state-space system model with identified A, B, and C matrices.
sys_est = ss(A, B, ccc(1, 1:end),ddd(1,1), Ts);


%% FREQUENCY RESPONSE PLOT
color_all = {[1 0 0], [0 1 0], [0 0 1]};
line_style = {'-','--',':'};
line_width = {2,1.5,1.5};
% actual system
[mag,freq,pha] = m_freq_resp_cal(output_prbs,input_prbs,Fs);
if size(mag,1)<size(mag,2)
    mag = mag';
end
if size(freq,1)<size(freq,2)
    freq = freq';
end
% multirate system
[mag_prbs, phase_prbs, w_prbs] = bode(sys_est,freq*2*pi);
mag_prbs = squeeze(mag_prbs);
mag_prbs = (20*log10(mag_prbs));
phase_prbs = squeeze(phase_prbs);

% single rate system
[mag_single, phase_single, w_single] = bode(sys_identified_buildin_prbs,freq*2*pi);
mag_single = squeeze(mag_single);
mag_single = (20*log10(mag_single));
phase_single = squeeze(phase_single);

% phase matching
if (phase_prbs(2)<0) ~= (pha(2)<0) && phase_prbs(2)<0
    phase_prbs = phase_prbs+360;
elseif (phase_prbs(2)<0) ~= (pha(2)<0) && phase_prbs(2)>0
    phase_prbs = phase_prbs-360;
else
    phase_prbs = phase_prbs;
end
if (phase_single(2)<0) ~= (pha(2)<0) && phase_single(2)<0
    phase_single = phase_single+360;
elseif (phase_single(2)<0) ~= (pha(2)<0) && phase_single(2)>0
    phase_single = phase_single-360;
else
    phase_single = phase_single;
end
figure
subplot(2,1,1)
h = semilogx(freq,mag,freq,mag_prbs,freq,mag_single); % actual
for i = 1:3
    h(i).Color = color_all{i};
    h(i).LineStyle = line_style{i};
    h(i).LineWidth = line_width{i};
end
xlabel('Hz')
ylabel('Magnitude (dB)')
xlim([0.5 Fs/2])
title(['Time span: ',num2str(time_prbs),' sec; Model Order: ',num2str(nx),'; Data Points: ',num2str(prbs_end-prbs_strt)])
subplot(2,1,2)
h = semilogx(freq,pha,freq,phase_prbs,freq,phase_single); % actual
for i = 1:3
    h(i).Color = color_all{i};
    h(i).LineStyle = line_style{i};
    h(i).LineWidth = line_width{i};
end
xlim([0.5 Fs/2])
legend('actual system ',' multirate n4sid','singlerate n4sid')
xlabel('Hz')
ylabel('Magnitude (dB)')

%% plotting spectrum response
spec_in_f = specCale(input_prbs,Fs);
spec_out_f = specCale(output_prbs,Fs);
spec_in = specCale(in_V(prbs_strt:prbs_end),Fs);
spec_out = specCale(out_E(prbs_strt:prbs_end),Fs);
% spec_in_f = freq_resp_cal(input_prbs,Fs);
% spec_out_f = freq_resp_cal(output_prbs,Fs);
% spec_in = freq_resp_cal(in_V(prbs_strt:prbs_end),Fs);
% spec_out = freq_resp_cal(out_E(prbs_strt:prbs_end),Fs);
figure()
subplot(2,1,1)
h = semilogy(spec_in.f,spec_in.amp,spec_in_f.f,spec_in_f.amp);
title('Input Spectrum')
xlabel('Hz')
ylabel('Magnitude')
subplot(2,1,2)
h = semilogy(spec_out.f,spec_out.amp,spec_out_f.f,spec_out_f.amp);
title('Output Spectrum')
xlabel('Hz')
ylabel('Magnitude')
legend('Unfiltered','Filtered')

% plot filtered and unfiltered inputs/outputs
figure
subplot(2,1,1)
plot(tt(prbs_strt:prbs_end),out_E(prbs_strt:prbs_end))
title('Unfiltered Output')
xlim([tt(prbs_strt) tt(prbs_end)])
xlabel('Time (sec)')
ylabel('Encoder Count')
subplot(2,1,2)
plot(tt(prbs_strt:prbs_end),out_E_f(prbs_strt:prbs_end))
title('Filtered Output')
xlim([tt(prbs_strt) tt(prbs_end)])
xlabel('Time (sec)')