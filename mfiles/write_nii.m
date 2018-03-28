function write_nii(hdr, img)

% write_nii Write NIfTI header and image

% Author: Ali Tabesh
% Version: 2.5
% Last modified: 07/24/12

% -------------------------------------------------------------------------
% write header
% -------------------------------------------------------------------------

fid = fopen(hdr.fn, 'w', hdr.endian);
if fid == -1
    error('Cannot open output file: %s!\n', fn);
    return
end

fwrite(fid, hdr.sizeof_hdr, 'int32');
fwrite(fid, hdr.data_type, 'char');
fwrite(fid, hdr.db_name, 'char');
fwrite(fid, hdr.extents, 'int32');
fwrite(fid, hdr.session_error, 'int16');
fwrite(fid, hdr.regular, 'char');
fwrite(fid, hdr.dim_info, 'char');
fwrite(fid, hdr.dim, 'int16');
fwrite(fid, hdr.intent_p1, 'float32');
fwrite(fid, hdr.intent_p2, 'float32');
fwrite(fid, hdr.intent_p3, 'float32');
fwrite(fid, hdr.intent_code, 'int16');
fwrite(fid, hdr.datatype, 'int16');
fwrite(fid, hdr.bitpix, 'int16');
fwrite(fid, hdr.slice_start, 'int16');
fwrite(fid, hdr.pixdim, 'float32');
fwrite(fid, hdr.vox_offset, 'float32');
fwrite(fid, hdr.scl_slope, 'float32');
fwrite(fid, hdr.scl_inter, 'float32');
fwrite(fid, hdr.slice_end, 'int16');
fwrite(fid, hdr.slice_code, 'char');
fwrite(fid, hdr.xyzt_units, 'char');
fwrite(fid, hdr.cal_max, 'float32');
fwrite(fid, hdr.cal_min, 'float32');
fwrite(fid, hdr.slice_duration, 'float32');
fwrite(fid, hdr.toffset, 'float32');
fwrite(fid, hdr.glmax, 'int32');
fwrite(fid, hdr.glmin, 'int32');
fwrite(fid, hdr.descrip, 'char');
fwrite(fid, hdr.aux_file, 'char');
fwrite(fid, hdr.qform_code, 'int16');
fwrite(fid, hdr.sform_code, 'int16');
fwrite(fid, hdr.quatern_b, 'float32');
fwrite(fid, hdr.quatern_c, 'float32');
fwrite(fid, hdr.quatern_d, 'float32');
fwrite(fid, hdr.quatern_x, 'float32');
fwrite(fid, hdr.quatern_y, 'float32');
fwrite(fid, hdr.quatern_z, 'float32');
fwrite(fid, hdr.srow_x, 'float32');
fwrite(fid, hdr.srow_y, 'float32');
fwrite(fid, hdr.srow_z, 'float32');
fwrite(fid, hdr.intent_name, 'char');
fwrite(fid, hdr.magic, 'char');

% zero-pad to vox_offset
fwrite(fid, zeros(hdr.vox_offset - ftell(fid), 1), 'uchar');

% -------------------------------------------------------------------------
% write image data
% -------------------------------------------------------------------------

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

if hdr.scl_slope ~= 0
    img = (img - hdr.scl_inter) / hdr.scl_slope;
end

nvoxels_written = fwrite(fid, img, precision);

% check if number of written voxels was same as expected voxels

dim = hdr.dim(2:end);
dim(dim == 0) = 1;
nvoxels_expected = prod(dim);

if nvoxels_expected ~= nvoxels_written
    error('Number of image voxels written was different than number of expected voxels in %s!', hdr.fn);
    return
end

fclose(fid);
