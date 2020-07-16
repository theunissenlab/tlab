
function f = h5utils()
    
    f = struct;
    
    f.get_subgroups = @(x, y) h5_get_subgroups(x, y);
    
    f.get_ds = @(x, y) h5_get_ds(x, y);
    f.set_ds = @(w, x, y, z) h5_set_ds(w, x, y, z);
    
    f.get_attr = @(x, y, z) h5_get_attr(x, y, z);
    f.set_attr = @(a, b, c, d) h5_set_attr(a, b, c, d);
    
    f.create_group = @(x,y) h5_create_group_rec(x, y);
    f.create = @(x) h5_create_file(x);
    
    f.open = @(x,y) h5_open_file(x,y);
    f.close = @(x) h5_close_file(x);
    
    f.get_attrs = @(x, y) h5_get_attrs(x, y);
end

function fid = h5_open_file(fileName, perm)

    h5perm = 'H5F_ACC_RDWR';
    if nargin == 2
        if strncmp(perm,'a', 1)
            h5perm = 'H5F_ACC_RDWR';
        elseif strncmp(perm, 'r', 1)
            h5perm = 'H5F_ACC_RDONLY';
        else
            warning('Unknown permission option in h5_open_file. Will use RW');
        end
    end
    
    fid = H5F.open(fileName, h5perm,'H5P_DEFAULT');
end

function fid = h5_close_file(fid)
    H5F.close(fid);
end


function fid = h5_create_file(fileName)
    fcpl = H5P.create('H5P_FILE_CREATE');
    fapl = H5P.create('H5P_FILE_ACCESS');
    fid = H5F.create(fileName, 'H5F_ACC_TRUNC', fcpl, fapl);
end

function h5_create_group(fid, groupPath)
    try
        gid = H5G.open(fid, groupPath);    
    catch exception
        gid = H5G.create(fid, groupPath, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');        
    end    
    H5G.close(gid);
end

function h5_create_group_rec(fid, groupPath)

    if groupPath(1) ~= '/'
        error('HD5F group paths should start with leading slash!');
    end

    gpath = '';
    grps = regexp(groupPath, '/', 'split');
    for k = 2:length(grps)        
        gpath = [gpath '/' grps{k}];
        h5_create_group(fid, gpath);
    end
end

function h5_set_ds(fid, groupPath, dsName, dsVal)
    
    typeId = H5T.copy('H5T_NATIVE_DOUBLE');
    % dims = fliplr(size(dsVal));
    dims = size(dsVal);              % We are not flipping the order of the array...
    maxdims = dims;
    spaceId = H5S.create_simple(length(dims), dims, maxdims);    
    
    h5_create_group_rec(fid, groupPath);
    gid = H5G.open(fid, groupPath);

    try
        dsid = H5D.open(gid, dsName);
    catch exception
        dsid = H5D.create(gid, dsName, typeId, spaceId, ...
                          'H5P_DEFAULT');
    end
         
    % But we are writing the transpose - works for 2 d arrays only
    % The transpose is needed either in the read or the write
    H5D.write(dsid, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', dsVal');
                  
    H5D.close(dsid);
    H5G.close(gid);
end


function h5_set_attr(fid, groupPath, attrName, attrValue)

    if isnumeric(attrValue)
        typeId = H5T.copy('H5T_NATIVE_DOUBLE');
        spaceId = H5S.create('H5S_SCALAR');
        attrValue = double(attrValue);
    elseif ischar(attrValue)        
        typeId = H5T.copy('H5T_C_S1');
        if ~isempty(attrValue)
            % Don't do this when working with empty strings.
            H5T.set_size(typeId,numel(attrValue));
        else
            return;
        end
		H5T.set_strpad(typeId,'H5T_STR_NULLTERM');
        spaceId = H5S.create('H5S_SCALAR');
    end
    
    h5_create_group_rec(fid, groupPath);
    gid = H5G.open(fid, groupPath);
        
    aid = H5A.create(gid, attrName, typeId, spaceId, ...
                     H5P.create('H5P_ATTRIBUTE_CREATE'), 'H5P_DEFAULT');
    
    H5A.write(aid, typeId, attrValue);
                 
    H5A.close(aid);
    H5G.close(gid);
end

function clist = h5_get_subgroups(fid, groupPath)

    gid = H5G.open(fid, groupPath);
    
    clist = {};    
    [status, clist] = H5O.visit(gid, 'H5_INDEX_NAME',  'H5_ITER_INC', ...
                                @h5_append_subgroup, clist);
    H5G.close(gid);                            
end


function [status, clist] = h5_append_subgroup(gid, name, clist)
    if isempty(strfind(name, '/'))
        if ~strcmp('.', name)
            clist{end+1} = name;
        end
    end        
    status = 0;
end


function val = h5_get_ds(fid, dsPath)    
    did = H5D.open(fid, dsPath);
    val = H5D.read(did, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
    H5D.close(did);
end


function val = h5_get_attr(fid, groupPath, attrName)

    gid = H5G.open(fid, groupPath);

    aid = H5A.open(gid, attrName, 'H5P_DEFAULT');
    val = H5A.read(aid,'H5ML_DEFAULT');
    
    if ischar(val)
        if size(val, 1) ~= 1
            val = val';
        end
    end
    
    if iscell(val)
        val = val{1};
    end
    
    H5A.close(aid);
    H5G.close(gid);
end


function attrs = h5_get_attrs(fid, groupPath)
    
    attrs = struct;
    gid = H5G.open(fid, groupPath);
    
    [status, idx_stop, cdataOut] = H5A.iterate(gid, 'H5_INDEX_NAME', 'H5_ITER_NATIVE', 0, @h5_iter_attr, attrs);
    
    attrs = cdataOut;
    H5G.close(gid);   

end


function [status, cdataOut] = h5_iter_attr(locId, attrName, info, cdataIn)

    attrId = H5A.open(locId, attrName, 'H5P_DEFAULT');
    attrVal = H5A.read(attrId, 'H5ML_DEFAULT');
    if ischar(attrVal)
        if size(attrVal, 1) ~= 1
            attrVal = attrVal';
        end
    end
    H5A.close(attrId);
    %fprintf('h5_iter_attr: locId=%d, attrName=%s, val=%s\n', locId, attrName, attrVal);    
    cdataOut = cdataIn;
    cdataOut.(attrName) = attrVal;
    status = 0;
end



