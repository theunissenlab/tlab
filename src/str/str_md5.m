%% compute the md5 hash of a string
function Digest = str_md5(stringToHash)

% Write string to temporary file
FileName = tempname();
fid = fopen(FileName, 'w');
fwrite(fid, stringToHash);
fclose(fid);

Digest = file_md5(FileName);

% Delete temp file
delete(FileName);
