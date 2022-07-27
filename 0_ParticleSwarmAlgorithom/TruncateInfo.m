function info = TruncateInfo(info,maxIter,iter)
%
% Removes the empty entries in the info struct
%

names = fieldnames(info);
for i = 1 : length(names)
    if (isnumeric(info.(names{i})))   % Check if it's a matrix
        if size(info.(names{i}), 2) == maxIter    % Check if it is iteration data
            info.(names{i}) = info.(names{i})(:, 1 : iter);
        end
    end
end

end