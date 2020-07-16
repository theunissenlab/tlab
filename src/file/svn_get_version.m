function version_number = svn_get_version(filename)

svn_string = sprintf('svn status --verbose --xml %s',filename);
[status,xml_result] = system(svn_string);
status_tag = regexp(xml_result,'<wc-status.*?>','match','once');
rev_string = regexp(status_tag,'revision="[0-9]*"','match','once');
version_number = textscan(rev_string,'revision=%q');
version_number = str2num(version_number{1}{1});