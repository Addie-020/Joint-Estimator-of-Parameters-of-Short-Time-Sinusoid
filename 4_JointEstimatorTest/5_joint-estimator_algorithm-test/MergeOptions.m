function output = MergeOptions(default, user, name)
%
% Merge a default options struct with a user-defined options struct. Works
% recursively, and will issue warning messages if the user attempts to
% define a field that is not in the default options.
%
% DESCRIPTION:
%
% - All fields in DEFAULT will be present in OUTPUT
% - If a field is in both DEFAULT and USER, then the value from USER is
% present in OUTPUT
% - If a field is present in USER, but not DEFAULT, then issue a warning.
% - Applies recursively
%
% NOTES:
%
%   The argument "name" is optional, and contains a string specifying the
%   name of the options struct. This is primarily used for printing
%   warnings to the user.
%
%   This function works recursively. For example, if there is a struct
%   inside of a struct, then it will recursively apply this merge.
%

% Start by assuming that the OUTPUT is just the DEFAULT
output = default;

% Check if user define option name
if nargin == 2
    structName = '';
else
    structName = [name '.'];
end % end: if

% Merge user-define options with default ones
if ~isempty(user)
    % Check for any overriding fields in the USER-defined struct
    defaultFields = fieldnames(default);
    for i = 1 : length(defaultFields)
        if isfield(user, defaultFields{i})
            C0 = isstruct(default.(defaultFields{i}));
            C1 = isstruct(user.(defaultFields{i}));
            if C0 && C1         % Both are structs
                output.(defaultFields{i}) = MergeOptions(...
                    default.(defaultFields{i}), ...
                    user.(defaultFields{i}), ...
                    [structName defaultFields{i}]);
            elseif ~C0 && ~C1   % Both are fields
                output.(defaultFields{i}) = user.(defaultFields{i});
            elseif C0 && ~C1    %default is struct, user is a field
                disp(['WARNING: ' structName defaultFields{i} ' should be a struct!']);
            elseif ~C0 && C1    %default is struct, user is a field
                disp(['WARNING: ' structName defaultFields{i} ' should not be a struct!']);
            end % end: if
        end % end: if
    end % end: for
    % Check for any fields in USER that are not in DEFAULT
    userFields = fieldnames(user);
    for i = 1 : length(userFields)
        if ~isfield(default, userFields{i})
            disp(['WARNING: unrecognized option: ' structName userFields{i}]);
        end % end: if
    end % end: for
end % end: if

end % end: function MergeOptions

