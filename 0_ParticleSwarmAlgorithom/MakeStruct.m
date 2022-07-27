function S = MakeStruct(varargin)
%
% A struct is created with the property that each field corresponds to one
% of the arguments passed to this function.
%
% Example:
%
%   If defines:
%       a = 1;
%       b = 2;
%       c = 0;
%       S = makeStruct(a,b,c);
%   Then
%       S.a = 1;
%       S.b = 2;
%       S.c = 0;
%
% Notes:
%
%   Input names should be unique.
%

N_Inputs = length(varargin);

for i = 1 : N_Inputs
    name = inputname(i);
    S.(name) = varargin{i};
end

end