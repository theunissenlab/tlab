% Testing h5 write functions

testdata = [1; 2; 3; 4; 5];

fname = '/Users/frederictheunissen/Desktop/test.h5';
h5 = h5utils();
    
% Create file
fid = h5.create(fname);


h5.set_attr(fid, '/', 'testattr', 'This is test data');
h5.set_attr(fid, '/', 'testval', 0.456);

h5.set_ds(fid, '/', 'testdata', testdata);


h5.close(fid);