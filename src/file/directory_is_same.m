function result = directory_is_same(directory1,directory2)

%% Return true if the two directories are identical; false otherwise

% Get file list for both directories
d1_files = dir(directory1);
d2_files = dir(directory2);
nfiles = length(d1_files);
d2_names = {d2.name};

% Return 0 unless the file lists are the same length
if length(d2_files) != nfiles
	result = 0;
	return
end

% Initialize results vector
results = zeros(1,nfiles);

% Loop over all files in directory 1, checking each against directory 2
for jj = 1:length(d1_files)
	% Assign file object from d1 list
	file = d1_files(jj)
	
	% Create full file names
	fname1 = fullfile(directory1,file.name)
	fname2 = fullfile(directory2,file.name)
	
	% Check whether the files are identical
	if exist(fname2,'file')
		result(jj) = file_is_same(fname1,fname2);
	else
		result(jj) = 0;
	end
end

result = all(result>0)