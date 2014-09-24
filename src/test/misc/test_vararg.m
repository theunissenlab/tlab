function test_vararg(a, b, varargin)

    fprintf('length of extra vars=%d\n', length(varargin));
    
    proc_vararg(varargin{:});
end

function proc_vararg(d, e)

    if nargin == 0
        d = 99;
        e = 100;
    end

    fprintf('d=%d, e=%d\n', d, e);


end