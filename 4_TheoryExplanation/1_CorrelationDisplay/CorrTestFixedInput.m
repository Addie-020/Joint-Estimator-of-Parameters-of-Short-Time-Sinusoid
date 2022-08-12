% Description:  Test Program for Correlation Method
%               Parameters from keyboard input
% Projet:       Short Sequence Parameter Estimation
% Date:         May 3, 2022
% Author:       Zhiyu Shen

clear
clc

fprintf('-------- Short Sequence Signal Frequency and Initial Phase Co-Estimation --------\n\n');

row = 3;
line = 2;

%%% Fetch test signal parameters from keyboard input

fprintf('Please enter parameter values of test signal...\n');
a_t = input('Amplitude (1V): ');
f_t = input('Frequency (0 ~ 1 Hz): ');
p_t = input('Initial phase (-pi ~ pi rad): ');
fprintf('Please input sampling parameters...\n');
Fs = input('Sampling rate (Hz): ');
t_s = input('Sampling time (s): ');
n_e = input('Adding noise? (1: Yes, 0: No): ');
fprintf('The test signal is %.0f*sin(2pi*%.2f+%.2f)\n', a_t, f_t, p_t);
fprintf('Signal sampling rate is %.2fHz and sample for %.2fs\n', Fs, t_s);

tic

%%% Generate test signal and preparation

fprintf('\n------------------------------- Estimation Start -------------------------------\n');
% % Variable definition
% A_t = 1;            % Original signal's amplitude (V)
% f_t = 0.1;          % Original signal's frequency (Hz)
% p_t = -3/7*pi;      % Original signal's initial phase (rad)
% Fs = 50;            % Sampling rate (Hz)
% t_s = 5;            % Sampling time (s)

% Generate a sine wave for test
[fn_r, L, index_t] = test_wave(a_t, f_t, p_t, Fs, t_s);

% Add gaussian white noise
if n_e == 0
    fn = fn_r;
else
    SNR = 15;
    PWR = 0;
    fn = awgn(fn_r, SNR, PWR);
end

% Plot time domain signal
x_f = index_t/Fs;
subplot(row, line, 1)
hold on
plot(x_f, fn, 'b-', 'LineWidth', 0.8);
axis([0 t_s -1.2*a_t 1.2*a_t]);
title('Time domain signal');
xlabel('\it t /\rm s');
ylabel('\it f \rm(\itt\rm)');
set(gca, 'FontSize', 10 , 'FontName', 'Times New Roman');
grid on;
hold off;

%% Estimation process 1
% Frequency scale: 0.01 ~ 1 Hz
% Frequency precision: 0.01 Hz
% Phase scale: -pi ~ pi rad
% Phase precision: pi/100 rad
est_time = 1;
% Sweep freqeuncy and initial phase setting
A_e = 1;                                    % Amplitude (V)
f_head = 0.01;                              % Starting frequency (Hz)
f_end = 1;                                  % Ending frequency (Hz)
f_inc = 0.01;                               % Frequency increment (Hz)
p_head = -pi;                               % Starting phase (rad)
p_end = pi;                                 % Ending phase (rad)
p_inc = pi/100;                             % Phase increment (rad)

[R, i_max, j_max] = correlation_sweep_time(A_e, f_head, f_end, f_inc, p_head, p_end, p_inc, ...
    index_t, fn, Fs);

% Plot correlation coefficient diagram
x_n = 0 : i_max;                            % Index matrix of frequency
y_n = 0 : j_max;                            % Index matrix of phase
x_pos = f_head + x_n * f_inc;               % X coordinate: frequency axis
y_pos = p_head + y_n * p_inc;               % Y coordinate: phase axis
[X, Y] = meshgrid(x_pos, y_pos);            % Generate coordinate grid
Z = (R(x_n + 1, y_n + 1)).';                % Z coordinate: cross correlation coefficient axis
subplot(row, line, est_time + 1);
s = surf(X, Y, Z);
s.FaceAlpha = 1;
s.EdgeColor = 'flat';
s.Marker = 'none';
title('Correlation Coefficient');
xlabel('\it f /\rm Hz');
ylabel('\it \phi /\rm rad');
zlabel('\it R \rm(\itf,\phi\rm)');
set(gca, 'FontSize', 10 , 'FontName', 'Times New Roman');
grid on;

% Find maximal value of correlation coefficient
r_t = R;
[m, i_m] = max(r_t);
[m2, i_m2] = max(m);
i_hat = i_m(i_m2) - 1;
j_hat = i_m2 - 1;
R_max = m2;
f_hat = f_head + i_hat * f_inc;
p_hat = p_head + j_hat * p_inc;
fprintf('Estimation No.%d:\n', est_time);
fprintf('Precision: %.5f Hz, %.5f pi rad\n', f_inc, p_inc/pi);
fprintf('Maximal correlation coefficient: %.5f\n', R_max);
fprintf('Frequency: %.5f Hz\n', f_hat);
fprintf('Phase: %.5f rad (%.5fpi)\n', p_hat, p_hat/pi);

%% Estimation process 2
% Frequency scale: 0.(x-1) ~ 0.(x+1) Hz
% Precision: 0.002 Hz
% Phase scale: 0.(x-1)*pi ~ 0.(x+1)*pi rad
% Phase precision: pi/500 rad
est_time = est_time + 1;
est_f = round(f_hat * 10);
est_p = round(p_hat / pi * 10);

% Sweep freqeuncy and initial phase setting
A_e = 1;                                    % Amplitude (V)
f_head = max(0, (est_f - 1)) / 10;          % Starting frequency (Hz)
f_end = (est_f + 1) / 10;                   % Ending frequency (Hz)
f_inc = 0.002;                              % Frequency increment (Hz)
p_head = (est_p - 1) / 10 * pi;             % Starting phase (rad)
p_end = (est_p + 1) / 10 * pi;              % Ending phase (rad)
p_inc = pi/500;                             % Phase increment (rad)

[R, i_max, j_max] = correlation_sweep_time(A_e, f_head, f_end, f_inc, p_head, p_end, p_inc, ...
    index_t, fn, Fs);

% Plot correlation coefficient diagram
x_n = 0 : i_max;                            % Index matrix of frequency
y_n = 0 : j_max;                            % Index matrix of phase
x_pos = f_head + x_n * f_inc;               % X coordinate: frequency axis
y_pos = p_head + y_n * p_inc;               % Y coordinate: phase axis
[X, Y] = meshgrid(x_pos, y_pos);            % Generate coordinate grid
Z = (R(x_n + 1, y_n + 1)).';                % Z coordinate: cross correlation coefficient axis
subplot(row, line, est_time + 1);
s = surf(X, Y, Z);
s.FaceAlpha = 1;
s.EdgeColor = 'flat';
s.Marker = 'none';
title('Correlation Coefficient');
xlabel('\it f /\rm Hz');
ylabel('\it \phi /\rm rad');
zlabel('\it R \rm(\itf,\phi\rm)');
set(gca, 'FontSize', 10 , 'FontName', 'Times New Roman');
grid on;

% Find maximal value of correlation coefficient
r_t = R;
[m, i_m] = max(r_t);
[m2, i_m2] = max(m);
i_hat = i_m(i_m2) - 1;
j_hat = i_m2 - 1;
R_max = m2;
f_hat = f_head + i_hat * f_inc;
p_hat = p_head + j_hat * p_inc;
fprintf('Estimation No.%d:\n', est_time);
fprintf('Precision: %.5f Hz, %.5f pi rad\n', f_inc, p_inc/pi);
fprintf('Maximal correlation coefficient: %.5f\n', R_max);
fprintf('Frequency: %f Hz\n', f_hat);
fprintf('Phase: %.5f rad (%.5fpi)\n', p_hat, p_hat/pi);

%% Estimation process 3
% Frequency scale: 0.x(y-1) ~ 0.x(y+1) Hz
% Precision: 0.0002 Hz
% Phase scale: 0.x(y-1)*pi ~ 0.x(y+1)*pi rad
% Phase precision: pi/5000 rad
est_time = est_time + 1;
est_f = round(f_hat * 100);
est_p = round(p_hat / pi * 100);

% Sweep freqeuncy and initial phase setting
A_e = 1;                                    % Amplitude (V)
f_head = max(0, (est_f - 1)) / 100;         % Starting frequency (Hz)
f_end = (est_f + 1) / 100;                  % Ending frequency (Hz)
f_inc = 0.0002;                             % Frequency increment (Hz)
p_head = (est_p - 1) / 100 * pi;            % Starting phase (rad)
p_end = (est_p + 1) / 100 * pi;             % Ending phase (rad)
p_inc = pi/5000;                            % Phase increment (rad)

[R, i_max, j_max] = correlation_sweep_time(A_e, f_head, f_end, f_inc, p_head, p_end, p_inc, ...
    index_t, fn, Fs);

% Plot correlation coefficient diagram
x_n = 0 : i_max;                            % Index matrix of frequency
y_n = 0 : j_max;                            % Index matrix of phase
x_pos = f_head + x_n * f_inc;               % X coordinate: frequency axis
y_pos = p_head + y_n * p_inc;               % Y coordinate: phase axis
[X, Y] = meshgrid(x_pos, y_pos);            % Generate coordinate grid
Z = (R(x_n + 1, y_n + 1)).';                % Z coordinate: cross correlation coefficient axis
subplot(row, line, est_time + 1);
s = surf(X, Y, Z);
s.FaceAlpha = 1;
s.EdgeColor = 'flat';
s.Marker = 'none';
title('Correlation Coefficient');
xlabel('\it f /\rm Hz');
ylabel('\it \phi /\rm rad');
zlabel('\it R \rm(\itf,\phi\rm)');
set(gca, 'FontSize', 10 , 'FontName', 'Times New Roman');
grid on;

% Find maximal value of correlation coefficient
r_t = R;
[m, i_m] = max(r_t);
[m2, i_m2] = max(m);
i_hat = i_m(i_m2) - 1;
j_hat = i_m2 - 1;
R_max = m2;
f_hat = f_head + i_hat * f_inc;
p_hat = p_head + j_hat * p_inc;
fprintf('Estimation No.%d:\n', est_time);
fprintf('Precision: %.5f Hz, %.5f pi rad\n', f_inc, p_inc/pi);
fprintf('Maximal correlation coefficient: %.5f\n', R_max);
fprintf('Frequency: %.5f Hz\n', f_hat);
fprintf('Phase: %.5f rad (%.5fpi)\n', p_hat, p_hat/pi);

%% Estimation process 4
% Frequency scale: 0.xy(z-1) ~ 0.xy(z+1) Hz
% Precision: 0.00002 Hz
% Phase scale: 0.xy(z-1)*pi ~ 0.xy(z+1)*pi rad
% Phase precision: pi/500000 rad
est_time = est_time + 1;
est_f = round(f_hat * 1000);
est_p = round(p_hat / pi * 1000);

% Sweep freqeuncy and initial phase setting
A_e = 1;                                    % Amplitude (V)
f_head = max(0, (est_f - 1)) / 1000;        % Starting frequency (Hz)
f_end = (est_f + 1) / 1000;                 % Ending frequency (Hz)
f_inc = 0.00002;                            % Frequency increment (Hz)
p_head = (est_p - 1) / 1000 * pi;           % Starting phase (rad)
p_end = (est_p + 1) / 1000 * pi;            % Ending phase (rad)
p_inc = pi/50000;                           % Phase increment (rad)

[R, i_max, j_max] = correlation_sweep_time(A_e, f_head, f_end, f_inc, p_head, p_end, p_inc, ...
    index_t, fn, Fs);

% Plot correlation coefficient diagram
x_n = 0 : i_max;                            % Index matrix of frequency
y_n = 0 : j_max;                            % Index matrix of phase
x_pos = f_head + x_n * f_inc;               % X coordinate: frequency axis
y_pos = p_head + y_n * p_inc;               % Y coordinate: phase axis
[X, Y] = meshgrid(x_pos, y_pos);            % Generate coordinate grid
Z = (R(x_n + 1, y_n + 1)).';                % Z coordinate: cross correlation coefficient axis
subplot(row, line, est_time + 1);
s = surf(X, Y, Z);
s.FaceAlpha = 1;
s.EdgeColor = 'flat';
s.Marker = 'none';
title('Correlation Coefficient');
xlabel('\it f /\rm Hz');
ylabel('\it \phi /\rm rad');
zlabel('\it R \rm(\itf,\phi\rm)');
set(gca, 'FontSize', 10 , 'FontName', 'Times New Roman');
grid on;

% Find maximal value of correlation coefficient
r_t = R;
[m, i_m] = max(r_t);
[m2, i_m2] = max(m);
i_hat = i_m(i_m2) - 1;
j_hat = i_m2 - 1;
R_max = m2;
f_hat = f_head + i_hat * f_inc;
p_hat = p_head + j_hat * p_inc;
fprintf('Estimation No.%d:\n', est_time);
fprintf('Precision: %.5f Hz, %.5f pi rad\n', f_inc, p_inc/pi);
fprintf('Maximal correlation coefficient: %.5f\n', R_max);
fprintf('Frequency: %.5f Hz\n', f_hat);
fprintf('Phase: %.5f rad (%.5fpi)\n', p_hat, p_hat/pi);

%% Calculate efficiency indicators
fprintf('\n-------------------------- Display estimation results --------------------------');
tot_time = toc;
f_prec = abs(f_hat - f_t) / f_t;
p_prec = abs(p_hat - p_t) / abs(p_t);

fprintf("\n");
fprintf('True value: %.5f Hz, %.5f rad (%.5fpi) s\n', f_t, p_t, p_t/pi);
fprintf('Estimation result: %.5f Hz, %.5f rad (%.5fpi) s\n', f_hat, p_hat, p_hat/pi);
fprintf('Sampling time: %.5f s\n', t_s);
fprintf('Estimation time: %.5f s\n', tot_time);
fprintf('Total time: %.5f s\n', tot_time + t_s);
fprintf('Frequency estimation precision: %.5f%%\n', f_prec * 100);
fprintf('Initial phase estimation precision: %.5f%%\n', p_prec * 100);