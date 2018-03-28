function [hdr img] = read_nii(fn, read_img_flag)

% read_nii Read NIfTI header and image

% Author: Ali Tabesh
% Version: 2.5
% Last modified: 07/25/12

% set header only flag to zero (to read both header and image)

if nargin == 1
    read_img_flag = 1;
end

% -------------------------------------------------------------------------
% read header
% -------------------------------------------------------------------------

% initialize

hdr = [];

% check if gzipped

[pathstr, name, ext] = fileparts(fn);
if strcmpi(ext, '.gz')
    gzipped_flag = 1;
    rand('seed', sum(100 * clock));
    dir_new = fullfile(pathstr, ['tmp' num2str(round(rand * 1e9))]);
    gunzip(fn, dir_new);
    fn = fullfile(dir_new, name);
else
    gzipped_flag = 0;
end
hdr.fn = fn;

% open as little endian

hdr.endian = 'l';
fid = fopen(fn, 'r', hdr.endian);
if fid == -1
    error('File does not exist: %s!\n', fn);
end

hdr.sizeof_hdr = fread(fid, 1, 'int');

% if sizeof_hdr is incorrect, try opening as big endian

if(hdr.sizeof_hdr ~= 348)
    fclose(fid);
    hdr.endian = 'b';
    fid = fopen(niftifile, 'r', hdr.endian);
    hdr.sizeof_hdr = fread(fid, 1, 'int');
    if(hdr.sizeof_hdr ~= 348)
        fclose(fid);
        fidrintf('Invalid NIfTI header size = %d! Should be 348.\n', hdr.sizeof_hdr);
    end
end

hdr.data_type = fread(fid, 10, 'char');
hdr.db_name = fread(fid, 18, 'char');
hdr.extents = fread(fid, 1, 'int32');
hdr.session_error = fread(fid, 1, 'int16');
hdr.regular = fread(fid, 1, 'char');
hdr.dim_info = fread(fid, 1, 'char');
hdr.dim = fread(fid, 8, 'int16')';
hdr.intent_p1 = fread(fid, 1, 'float32');
hdr.intent_p2 = fread(fid, 1, 'float32');
hdr.intent_p3 = fread(fid, 1, 'float32');
hdr.intent_code = fread(fid, 1, 'int16');
hdr.datatype = fread(fid, 1, 'int16');
hdr.bitpix = fread(fid, 1, 'int16');
hdr.slice_start = fread(fid, 1, 'int16');
hdr.pixdim = fread(fid, 8, 'float32')';
hdr.vox_offset = fread(fid, 1, 'float32');
hdr.scl_slope = fread(fid, 1, 'float32');
hdr.scl_inter = fread(fid, 1, 'float32');
hdr.slice_end = fread(fid, 1, 'int16');
hdr.slice_code = fread(fid, 1, 'char');
hdr.xyzt_units = fread(fid, 1, 'char');
hdr.cal_max = fread(fid, 1, 'float32');
hdr.cal_min = fread(fid, 1, 'float32');
hdr.slice_duration = fread(fid, 1, 'float32');
hdr.toffset = fread(fid, 1, 'float32');
hdr.glmax = fread(fid, 1, 'int32');
hdr.glmin = fread(fid, 1, 'int32');
hdr.descrip = fread(fid, 80, 'char');
hdr.aux_file = fread(fid, 24, 'char');
hdr.qform_code = fread(fid, 1, 'int16');
hdr.sform_code = fread(fid, 1, 'int16');
hdr.quatern_b = fread(fid, 1, 'float32');
hdr.quatern_c = fread(fid, 1, 'float32');
hdr.quatern_d = fread(fid, 1, 'float32');
hdr.quatern_x = fread(fid, 1, 'float32');
hdr.quatern_y = fread(fid, 1, 'float32');
hdr.quatern_z = fread(fid, 1, 'float32');
hdr.srow_x = fread(fid, 4, 'float32');
hdr.srow_y = fread(fid, 4, 'float32');
hdr.srow_z = fread(fid, 4, 'float32');
hdr.intent_name = fread(fid, 16, 'char');
hdr.magic = fread(fid, 4, 'char');

% -------------------------------------------------------------------------
% read image data
% -------------------------------------------------------------------------

if read_img_flag == 1

    % skip the header

    fseek(fid, round(hdr.vox_offset), 'bof');

    % determine matlab data type
    
    switch hdr.datatype

        case 2
            precision = 'uint8';
        case 4
            precision = 'int16';
        case 8
            precision = 'int32';
        case 16
            precision = 'float32';
        case 64
            precision = 'float64';
        case 128
            precision = 'uint8';
        case 256
            precision = 'int8';
        case 512
            precision = 'uint16';
        case 768
            precision = 'uint32';
        case 1024
            precision = 'int64';
        case 1280
            precision = 'uint64';
        otherwise
            error('Unsupported NIfTI data type %s!', num2str(hdr.datatype));

    end

    [img nvoxels_read] = fread(fid, inf, precision);

    % check if number of read voxels was same as expected voxels

    dim = hdr.dim(2:end);
    dim(dim == 0) = 1;
    nvoxels_expected = prod(dim);

    if nvoxels_expected ~= nvoxels_read
        error('Number of image voxels read was different than number of expected voxels in %s!', fn);
    end

    % reshape array and scale according to header

    img = double(reshape(img, dim));
    if hdr.scl_slope ~= 0
        img = img * hdr.scl_slope + hdr.scl_inter;
    end

end

fclose(fid);

% if gzipped delete tmp file

if gzipped_flag == 1
    delete(fn);
	rmdir(dir_new);
end
