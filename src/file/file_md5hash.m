function md5hash = file_md5hash(filename)

% Compute MD5 hash of file specified by filename
% n.b. it's best to pass full file paths to this function, since it
% executes in part on the command line

switch computer
    case {'PCWIN'}
    	% error('To enable support on windows machines, you must first install an MD5 binary')
        % If you want this to run on your PC, you need to install a binary
        % that returns the MD5 hash digest for a file. Do that, and put the
        % path to that binary here.
        md5command = 'md5';
    case {'MACI64'}
        md5command = 'md5 -r';
    case {'GLNX86' 'GLNXA64'}
        md5command = 'md5sum';
    otherwise
    	error('Operating system not supported')
end

syscommand = sprintf('%s "%s"',md5command,filename);
[status,result] = system(syscommand);

%Store MD5 hash only if program successfully generated one
if status == 0
    md5hash = upper(result(1:32));
else
    error('MyFcns:CommandLineError',...
        'Command line utility returned error message:\n%s',result);
end
