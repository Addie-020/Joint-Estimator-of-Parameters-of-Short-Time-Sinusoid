function S = MakeStruct(varargin)
%
% A struct is created with the property that each field corresponds to one
% of the arguments passed to this function.
%
% Example:
%
%   If defines:
%       S = makeStruct(a,b,c);
%   Then
%       S.a == a;
%       S.b == b;
%       S.c == c;
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