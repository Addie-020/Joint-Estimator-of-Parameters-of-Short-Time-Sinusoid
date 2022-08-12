% Sweep frequency and initial phase
% Estimatie frequency and initial phase among the given scale
% Input amplitude, frequency (head, end, increment), phase (head, end, increment)
% Input sampling rate, sequence index, length, original sequence spectrum
% Return correlation coefficient matrix, maximal i index, maximal j index
function [R, i_max, j_max] = correlation_sweep_time(A_e, f_head, f_end, f_inc, p_head, p_end, p_inc, ...
                                index, fn, Fs)
i_max = ceil((f_end - f_head) / f_inc);     % Maximum value of frequency index
j_max = ceil((p_end - p_head) / p_inc);     % Maximum value of phase index
f_num = uint8(i_max + 1);                   % Number of frequency sweeping points
p_num = uint8(j_max + 1);                   % Number of phase sweeping points
R = zeros(f_num, p_num);                    % Generate a matrix to store correlation coefficients

% Cross correlation with for loop
for i = 0 : i_max                               % Outer loop: frequency sweeping
    f_e = f_head + i.*f_inc;
    w_e = 2 * pi * f_e;
    for j = 0 : j_max-1                         % Inner loop: phase sweeping
        p_e = p_head + j*p_inc;
        ge = @(x) A_e * sin(w_e*x + p_e);
        gn = ge(index/Fs);
        r = corrcoef(fn, gn);
        R(i + 1, j + 1) = r(1, 2);
    end
end

end