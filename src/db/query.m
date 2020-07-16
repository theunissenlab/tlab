function result = query(class_name,object_id,property,server_addr)

%% Get query from database

if nargin < 4
	server_addr = 'http://app.fet.berkeley.edu/query';
end

url = sprintf('%s/%s/%d/%s', server_addr, class_name, object_id, property);
result = urlread(url);