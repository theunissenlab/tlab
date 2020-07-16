function result = file_is_same(fname1,fname2)

%% Return 1 if files are identical; return 0 if not identical;
%% return -1 if either file does not exist
% Currently, this routine uses the md5 hash to test whether the files are identical.

if exist(fname1,'file') && exist(fname2,'file')
	return filemd5hash(fname1) == filemd5hash(fname2)
else
	return -1
end