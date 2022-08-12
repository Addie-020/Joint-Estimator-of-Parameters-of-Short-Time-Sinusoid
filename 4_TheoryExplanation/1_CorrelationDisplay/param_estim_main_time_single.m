% Test script of cross-correlation
clear
clc

tic
% Variable definition
A_t = 1;            % Original signal's amplitude (V)
f_t = 0.05;          % Original signal's frequency (Hz)
p_t = pi/3;         % Original signal's initial phase (rad)
Fs = 50;            % Sampling rate (Hz)
t_s = 5;            % Sampling time (s)

% Generate a sine wave for test
[fn_r, L, index_t] = test_wave(A_t, f_t, p_t, Fs, t_s);

% % Add gaussian white noise
% SNR = 40;
% PWR = 0;
% fn = awgn(fn_r, SNR,PWR);
fn =fn_r;

% Generate estimated signal (Example)
A_e = 1;                            % Estimated signal's amplitude (V)
f_e = 10;                           % Estimated signal's frequency (Hz)
p_e = 0;                            % Estimated signal's initial phase (rad)
[gn, ~, ~] = test_wave(A_e, f_e, p_e, Fs, t_s);

% Plot time domain signal
x_f = index_t/Fs;
subplot(2, 1, 1)
hold on
plot(x_f, fn, 'b-', 'LineWidth', 0.8);
plot(x_f, gn, 'r-', 'LineWidth', 0.8);
axis([0 t_s -1.2*A_t 1.2*A_t]);
title('Time domain signal');
xlabel('\it t /\rm s');
ylabel('\it f \rm(\itt\rm), \it g \rm(\itt\rm)');
set(gca, 'FontSize', 10 , 'FontName', 'Times New Roman');
grid on;
hold off;

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
subplot(2, 1, 2);
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
fprintf('Maximal correlation coefficient: %f\n', R_max);
fprintf('Frequency: %f Hz\n', f_hat);
fprintf('Phase: %f (%fpi)\n', p_hat, p_hat/pi);

toc