function saveTableToHDF5(fname, lname, tdata)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   fname: File name
%   lname: Location name in HDF5
%   tdata: table data
sname = ['/' lname];

% Convert table to serialized byte vector
toSave = hlp_serialize(tdata);

% Create hdf5 container and allocate resources
h5create(fname, sname, size(toSave));

% Write data into hdf5 container
h5write(fname, sname, toSave);

disp('Done');

end

