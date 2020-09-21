function varargout = mmv2struct(varargin)
%MMV2STRUCT Pack/Unpack Variables to/from a Scalar Structure.
% MMV2STRUCT(X,Y,Z,. . .) returns a structure having fields X,Y,Z,. . .
% containing the corresponding data stored in X,Y,Z,. . .
% Inputs that are not variables are stored in fields named ansN
% where N is an integer identifying the Nth unnamed input.
%
% MMV2STRUCT(S)assigns the contents of the fields of the scalar structure
% S to variables in the calling workspace having names equal to the
% corresponding field names.
%
% [A,B,C,. . .] = MMV2STRUCT(S) assigns the contents of the fields of the
% scalar structure S to the variables A,B,C,. . . rather than overwriting
% variables in the caller. If there are fewer output variables than
% there are fields in S, the remaining fields are not extracted. Variables
% are assigned in the order given by fieldnames(S).

if nargin==0
    error('Input Arguments Required.')
elseif nargin==1 % Unpack Unpack Unpack Unpack Unpack Unpack Unpack Unpack
    argin = varargin{1};
    if ~isstruct(argin) || length(argin)~=1
        error('Single Input Must be a Scalar Structure.')
    end
    names = fieldnames(argin);
    if nargout==0 % assign in caller
        for i = 1:length(names)
            assignin('caller',names{i},argin.(names{i}))
        end
    else % deal fields into variables in caller
        varargout = cell(1,nargout); % preallocate output
        for i = 1:nargout
            varargout{i} = argin.(names{i});
        end
    end
else % Pack Pack Pack Pack Pack Pack Pack Pack Pack Pack Pack Pack Pack Pack
    args = cell(2,nargin);
    num = 1;
    for i = 1:nargin % build cells for call to struct
        args(:,i) = {inputname(i); varargin{i}};
        if isempty(args{1,i})
            args{1,i} = sprintf('ans%d',num);
            num = num + 1;
        end
    end
    varargout{1} = struct(args{:}); % create struct using comma-separated list
end


