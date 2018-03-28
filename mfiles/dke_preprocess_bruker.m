function dke_preprocess_bruker(basedir, bval, options)

%dke_preprocess Preprocess diffusional kurtosis imaging data
%       Data must be structured into folders DKI1, DKI2, ..., DKIn, and
%       DKI_B0 [optional], where n is the number of averages; all extra b = 0 images
%       must be copied to DKI_B0 folder
%
%       The function carries out the following processing steps
%
%       1. Conversion of all DICOM images into NIfTI format (using
%          spm_dicom_convert) [optional]
%
%       2. Coregistration of b = 0 images in DKI2, ..., DKIn, and first b =
%          0 images of all extra b = 0 series in DKI_B0 to b = 0 image in
%          DKI1 (using spm_coreg)
%
%       3. Averaging of coregistered images
%
%Syntax
%
%       dke_preprocess(basedir, options)
%
%Inputs
%
%       basedir Root folder for data
%
%       options Preprocessing options
%               format      Input image format ('dicom' or 'nifti') (default: 'dicom')
%               navg        Number of DKI averages (series) (i.e., DKI1, DKI2, ..., DKInavg)
%               extra_b0    Whether (1) or not (0) extra b = 0 images are available (default: 1)
%                           If extra_b0 = 1, number of extra b = 0 averages is automatically
%                           determined from the images in DKI_B0 folder
%               coreg_flag  Whether or not to coregister diffusion-weighted images across DKI averages (default: 1)
%
%Outputs
%
%       Folder intermediate_processing will be created which will contain the following folders
%
%       - DKI1_nii, ..., DKIn_nii, DKI_B0_nii contain original NIfTI images
%
%       - DKI1_coreg, ..., DKIn_coreg, DKI_B0_coreg contain NIfTI images coregistered to b = 0 in DKI1 series
%
%       - DKI_comb contains averaged coregistered NIfTI images
%
%NOTES
%
%       1. DICOM images must be in Siemens mosaic format.
%       2. SPM8 with update revision number 4667 must be used.

% Author: Edward S. Hui
% Modified: Ali Tabesh
% Version: 2.5
% Last modified: 07/17/12

%--------------------------------------------------------------------------
% convert DICOM images to NIfTI format
%--------------------------------------------------------------------------

fprintf('Converting input Bruker images to NIfTI... ')

% make intermediate_processing folder
mkdir(fullfile(basedir, 'intermediate_processing'))

bruker2nii(basedir);

fprintf('complete.\n')

%--------------------------------------------------------------------------
% co-register all b0 images to b0 image in DKI1 series
%--------------------------------------------------------------------------

if options.coreg_flag == 1
    
    fprintf('Co-registering images...\n')
    
    % coregister all b = 0 images

    targetdir = fullfile(basedir, 'intermediate_processing', 'dki_b0_nii');
    list = dir(fullfile(targetdir, 'dki_0*.nii'));
    fn_target = fullfile(targetdir, list(1).name);
    
    sourcedir = fullfile(basedir, 'intermediate_processing', 'dki_b0_nii');
    fn_source = spm_select('fplist', sourcedir, '^dki_0_*.*');

    coregister(fn_target, fn_source);

    coregdir = fullfile(basedir, 'intermediate_processing', 'dki_b0_coreg');
    mkdir(coregdir)
    movefile(fullfile(sourcedir, 'rdki*.*'), coregdir)
    
    % 1. average DWIs for each nonzero b-value
    % 2. coregister the mean DWI image to the mean DWI image for previous
    %    b-value (for the first nonzero b-value, the target is the first b = 0 image)
    % 3. coregister individual DWIs for that b-value to the coregistered
    %    mean DWI
    
    for iavg = 1:options.navg
    
        sourcedir = fullfile(basedir, 'intermediate_processing', ['dki_avg', num2str(iavg), '_nii']);
        coregdir  = fullfile(basedir, 'intermediate_processing', ['dki_avg', num2str(iavg), '_coreg']);
    
        for Bval = 2:length(bval)
            
            % calculate mean DWI for each non-zero b-value
            list = dir(fullfile(sourcedir, ['dki_', num2str(bval(Bval)), '_*.nii']));
            for i = 1:size(list, 1)
                fn_array{i} = fullfile(sourcedir, list(i).name);
            end
            image_avg(fn_array, fullfile(sourcedir, ['adki_', num2str(bval(Bval)), '.nii']));
            
            if Bval == 2
                % co-register mean DWI for the first nonzero b-value to first b = 0 image
                coregister(fn_target, fullfile(sourcedir, ['adki_', num2str(bval(Bval)), '.nii']));
            else
                % co-register mean DWI for b-values other than the first b-value at a smaller bvalue to a
                coregister(fullfile(sourcedir, ['radki_', num2str(bval(Bval - 1)), '.nii']), ...
                    fullfile(sourcedir, ['adki_', num2str(bval(Bval)), '.nii']));
            end
        end
        parfor Bval = 2:length(bval)

            % co-register all DWIs at a given b-value to mean coregistered
            % DWI at the same b-value
            fn_source = spm_select('fplist', sourcedir, ['^dki_', num2str(bval(Bval)), '_[0-9]+.nii']);
            coregister(fullfile(sourcedir, ['radki_', num2str(bval(Bval)), '.nii']), fn_source);

        end
        
        mkdir(coregdir)
        movefile(fullfile(sourcedir, 'rdki*.*'), coregdir)

    end
    
    fprintf('Co-registration complete.\n')
    
elseif options.coreg_flag == 0
    
    sourcedir = fullfile(basedir, 'intermediate_processing', 'dki_b0_nii');
    list = dir(fullfile(sourcedir, 'dki*.nii'));
    coregdir = fullfile(basedir, 'intermediate_processing', 'dki_b0_coreg');
    mkdir(coregdir)
    for i = 1:size(list, 1)
        fn_source = fullfile(sourcedir, list(i).name);
        fn_coreg = fullfile(coregdir, ['r', list(i).name]);
        copyfile(fn_source, fn_coreg)
    end
    
    for iavg = 1:options.navg
        sourcedir = fullfile(basedir, 'intermediate_processing', ['dki_avg', num2str(iavg), '_nii']);
        list = dir(fullfile(sourcedir, 'dki*.nii'));
        coregdir = fullfile(basedir, 'intermediate_processing', ['dki_avg', num2str(iavg), '_coreg']);
        mkdir(coregdir)
        for i = 1:size(list, 1)
            fn_source = fullfile(sourcedir, list(i).name);
            fn_coreg = fullfile(coregdir, ['r', list(i).name]);
            copyfile(fn_source, fn_coreg)
        end
    end
    
else
    
    error('Invalid ''coreg_flag''! ''coreg_flag'' must be 0 or 1.')
    
end

%--------------------------------------------------------------------------
% average images
%--------------------------------------------------------------------------

fprintf('Averaging images... ')

folder_avg1 = fullfile(basedir, 'intermediate_processing', 'dki_avg1_coreg');
folder_output = fullfile(basedir, 'intermediate_processing', 'combined');

mkdir(folder_output)

% average b0 images
clear fn_array
folder_b0 = fullfile(basedir, 'intermediate_processing', 'DKI_B0_coreg');
list = dir(fullfile(folder_b0, 'rdki_*.nii'));
for ib0 = 1:size(list, 1)
    fn_array{ib0} = fullfile(basedir, 'intermediate_processing', 'dki_b0_coreg', list(ib0).name);
end
image_avg(fn_array, fullfile(folder_output, 'rdki_0.nii'))


% average non-b0 images
clear fn_array
list = dir(fullfile(folder_avg1, 'rdki_*_*.nii'));
for i = 1:size(list, 1)
    fn_array{1} = fullfile(folder_avg1, list(i).name);
    if options.navg>1
        for iavg = 2:options.navg
            fn_array{iavg} = fullfile(basedir, 'intermediate_processing', ['dki_avg', num2str(iavg), '_coreg'], list(i).name);
        end
        image_avg(fn_array, fullfile(folder_output, list(i).name))
    else
        copyfile(fn_array{1}, fullfile(folder_output, list(i).name))
    end
end

fprintf('complete.\n')



%--------------------------------------------------------------------------
% convert bruker raw data to nifti format
%--------------------------------------------------------------------------

function bruker2nii(basedir)

[DIM, VOX, NEX, NBo, bval, reco_wordtype] = brukeracq(basedir);
NDiffExp = NBo + DIM(4) * DIM(5);

switch reco_wordtype
    case 'int16'
        datatype = 4;
    case 'int32'
        datatype = 8;
    case 'single'
        datatype = 16;
end

mkdir(fullfile(basedir, 'intermediate_processing', 'dki_b0_nii'));

f=fopen(fullfile(basedir, '2dseq'), 'r');
for WhichNEX=1:NEX
    mkdir(fullfile(basedir, 'intermediate_processing', ['dki_avg', num2str(WhichNEX), '_nii']));
    for WhichB0 = 1:NBo
        im = reshape(fread(f, DIM(1) * DIM(2) * DIM(3), reco_wordtype), DIM(1:3));
        img_nii = make_nii(im, 10 * VOX, [0, 0, 0], datatype); % x10 the original voxels for SPM registration
        outputdir = fullfile(basedir, 'intermediate_processing', 'dki_b0_nii', ...
            ['dki_0_', num2str(WhichB0 + NBo * (WhichNEX - 1)), '.nii']);
        save_nii(img_nii, outputdir);
    end
    
    for WhichB0 = 1:NDiffExp - NBo
        im = reshape(fread(f, DIM(1) * DIM(2) * DIM(3), reco_wordtype), DIM(1:3));
        img_nii = make_nii(im, 10 * VOX, [0, 0, 0], datatype); % x10 the original voxels for SPM registration
        outputdir = fullfile(basedir, 'intermediate_processing', ['dki_avg', num2str(WhichNEX), '_nii'], ...
            ['dki_', num2str(bval(mod(WhichB0 - 1, DIM(4)) + 1)), '_' ,num2str(floor((WhichB0 - 1) / DIM(4)) + 1), '.nii']);
        save_nii(img_nii, outputdir);
    end
end
fclose(f);


%--------------------------------------------------------------------------
% co-register two b0 images
%--------------------------------------------------------------------------

function coregister(fn_target, fn_source)

% fn_other = spm_select('fplist', folder_other, fn_other_filt);

% coregistration and reslicing parameters
estflg.cost_fun = 'nmi';
estflg.sep      = [4, 2];
estflg.tol      = [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001];
estflg.fwhm     = [7, 7];
wrtflg.interp   = 1;
wrtflg.wrap     = [0, 0, 0];
wrtflg.mask     = 0;

hdr_trg = spm_vol(fn_target);
hdr_src = spm_vol(fn_source);

x  = spm_coreg(hdr_trg, hdr_src, estflg);
for j = 1:size(fn_source,1)
    M  = inv(spm_matrix(x(j, :)));
    MM = spm_get_space(deblank(fn_source(j, :)));
    spm_get_space(deblank(fn_source(j,:)), M * MM);
end
spm_reslice(char(fn_target, fn_source),wrtflg);



%--------------------------------------------------------------------------
% average images
%--------------------------------------------------------------------------

function image_avg(fn_array, fn_out)

hdr = spm_vol(fn_array{1});
imgavg = spm_read_vols(hdr);
for i = 2:length(fn_array)
    hdr = spm_vol(fn_array{i});
    img = spm_read_vols(hdr);
    imgavg = imgavg + img;
end
imgavg = imgavg / length(fn_array);

hdr.fname = fn_out;
spm_write_vol(hdr, imgavg);


%--------------------------------------------------------------------------
% Read Bruker's header files
%--------------------------------------------------------------------------

function [DIM VOX NEX NBo bval reco_wordtype] = brukeracq(path)

file = textread(fullfile(path, 'reco'),'%s','delimiter','\n','whitespace','');
for i=1:length(file)
    tmp = cell2mat(file(i));
    
    if length(tmp)>=11 && ~isempty(strfind(tmp, '##$RECO_fov'))
        fov_xy=str2num(cell2mat(file(i+1)));
    end
    
    if length(tmp)>=12 && ~isempty(strfind(tmp, '##$RECO_size'))
        xy=str2num(cell2mat(file(i+1)));
    end
    
    if length(tmp)>=16 && ~isempty(strfind(tmp, '##$RECO_wordtype'))
        pos1=strfind(tmp, '=');
        reco_wordtype = tmp(pos1+1:end);
        switch reco_wordtype
            case '_16BIT_SGN_INT'
                reco_wordtype='int16';
            case '_32BIT_SGN_INT'
                reco_wordtype='int32';
            case '_8BIT_UNSGN_INT'
                reco_wordtype='uint8';
            case '_32BIT_FLOAT'
                reco_wordtype='single';
            otherwise
                error('Unknown RECO_wordtype!')
        end
    end
end
DIM = [xy 1 1];

file = textread(fullfile(path, 'method'),'%s','delimiter','\n','whitespace','');
for i=1:length(file)
    tmp = cell2mat(file(i));
    
    if length(tmp)>=23 && ~isempty(strfind(tmp, '##$PVM_DwUsedSliceThick'))
        pos1=strfind(tmp, '=');
        vox_z = str2num(tmp(pos1+1:end));
    end
    
    if length(tmp)>=17 && ~isempty(strfind(tmp, '##$PVM_DwNDiffDir'))
        pos1=strfind(tmp, '=');
        DIM(5) = str2num(tmp(pos1+1:end));
    end
    
    if length(tmp)>=17 && ~isempty(strfind(tmp, '##$PVM_DwAoImages'))
        pos1=strfind(tmp, '=');
        NBo = str2num(tmp(pos1+1:end));
    end
    
    if length(tmp)>=17 && ~isempty(strfind(tmp, '##$PVM_DwBvalEach'))
        pos1=strfind(tmp, '(');
        pos2=strfind(tmp, ')');
        DIM(4) = str2num(tmp(pos1+2:pos2-2));
        bval=str2num(cell2mat(file(i+1)));
        offset=1;
        while length(bval)~=DIM(4),
            bval=[bval str2num(cell2mat(file(i+1+offset)));];
            offset=offset+1;
        end
    end
    
    if length(tmp)>=22 && ~isempty(strfind(tmp, '##$PVM_SPackArrNSlices'))
        DIM(3)=str2num(cell2mat(file(i+1)));
    end
    
    if length(tmp)>=19 && ~isempty(strfind(tmp, '##$PVM_NRepetitions'))
        pos1=strfind(tmp, '=');
        NEX = str2num(tmp(pos1+1:end));
    end
end
VOX = [fov_xy./xy*10 vox_z];

%  Make nii structure specified by 3D matrix [x y z]. It also takes
%  4D matrix like [x y z t]. Optional features can also be included,
%  such as: voxel_size, origin, datatype, and description.
%
%  Usage: nii = make_nii(img, [voxel_size], [origin], [datatype], ...
%		[description])
%
%  Where:
%
%	img:			3D matrix [x y z], or 4D matrix that
%				includes all the images along the
%				time course [x y z t]
%
%	voxel_size (optional):	Voxel size in millimeter for each
%				dimension. Default is [1 1 1].
%
%	origin (optional):	The AC origin. Default is [0 0 0].
%
%	datatype (optional):	Storage data type. Default is float32 [16]:
%		2 - uint8,  4 - int16,  8 - int32,  16 - float32,
%		32 - complex64,  64 - float64,  128 - RGB24,
%		256 - int8,  512 - uint16,  768 - uint32, 1792 - complex128
%
%	description (optional):	Description of data. Default is ''.
%
%  e.g.:
%     origin = [33 44 13]; datatype = 64;
%     nii = make_nii(img, [], origin, datatype);    % default voxel_size
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%
function nii = make_nii(varargin)

nii.img = varargin{1};
dims = size(nii.img);
dims = [length(dims) dims ones(1,8)];
dims = dims(1:8);

voxel_size = [0 ones(1,7)];
origin = zeros(1,5);
datatype = 16;
descrip = '';

if nargin > 1 & ~isempty(varargin{2})
    voxel_size(2:4) = double(varargin{2});
end

if nargin > 2 & ~isempty(varargin{3})
    origin(1:3) = double(varargin{3});
end

if nargin > 3 & ~isempty(varargin{4})
    datatype = double(varargin{4});
end

if nargin > 4 & ~isempty(varargin{5})
    descrip = varargin{5};
end

%MPH: Fix dims in the case of RGB datatype
%
if datatype == 128
    dims(1) = dims(1)-1;
    dims(5:8) = [dims(6:8) 1];
end

maxval = round(double(max(nii.img(:))));
minval = round(double(min(nii.img(:))));

nii.hdr = make_header(dims, voxel_size, origin, datatype, ...
    descrip, maxval, minval);

switch nii.hdr.dime.datatype
    case 2
        nii.img = uint8(nii.img);
    case 4
        nii.img = int16(nii.img);
    case 8
        nii.img = int32(nii.img);
    case 16
        nii.img = single(nii.img);
    case 32
        nii.img = single(nii.img);
    case 64
        nii.img = double(nii.img);
    case 128
        nii.img = uint8(nii.img);
    case 256
        nii.img = int8(nii.img);
    case 512
        nii.img = uint16(nii.img);
    case 768
        nii.img = uint32(nii.img);
    case 1792
        nii.img = double(nii.img);
    otherwise
        error('Datatype is not supported by make_nii.');
end

return;					% make_nii


%---------------------------------------------------------------------
function hdr = make_header(dims, voxel_size, origin, datatype, ...
    descrip, maxval, minval)

hdr.hk   = make_header_header_key;
hdr.dime = make_header_image_dimension(dims, voxel_size, datatype, maxval, minval);
hdr.hist = make_header_data_history(origin, descrip);

return;					% make_header


%---------------------------------------------------------------------
function hk = make_header_header_key

hk.sizeof_hdr       = 348;			% must be 348!
hk.data_type        = '';
hk.db_name          = '';
hk.extents          = 0;
hk.session_error    = 0;
hk.regular          = 'r';
hk.dim_info         = 0;

return;					% header_key


%---------------------------------------------------------------------
function dime = make_header_image_dimension(dims, voxel_size, datatype, maxval, minval)

dime.dim = dims;
dime.intent_p1 = 0;
dime.intent_p2 = 0;
dime.intent_p3 = 0;
dime.intent_code = 0;
dime.datatype = datatype;

switch dime.datatype
    case   2,
        dime.bitpix = 8;  
        precision = 'uint8';
    case   4,
        dime.bitpix = 16; 
        precision = 'int16';
    case   8,
        dime.bitpix = 32; 
        precision = 'int32';
    case  16,
        dime.bitpix = 32; 
        precision = 'float32';
    case  32,
        dime.bitpix = 64; 
        precision = 'float32';
    case  64,
        dime.bitpix = 64; 
        precision = 'float64';
    case  128,
        dime.bitpix = 24; 
        precision = 'uint8';
    case 256
        dime.bitpix = 8;  
        precision = 'int8';
    case 512
        dime.bitpix = 16; 
        precision = 'uint16';
    case 768
        dime.bitpix = 32; 
        precision = 'uint32';
    case  1792,
        dime.bitpix = 128; 
        precision = 'float64';
    otherwise
        error('Datatype is not supported by make_nii.');
end

dime.slice_start = 0;
dime.pixdim = voxel_size;
dime.vox_offset = 0;
dime.scl_slope = 0;
dime.scl_inter = 0;
dime.slice_end = 0;
dime.slice_code = 0;
dime.xyzt_units = 0;
dime.cal_max = 0;
dime.cal_min = 0;
dime.slice_duration = 0;
dime.toffset = 0;
dime.glmax = maxval;
dime.glmin = minval;

return;					% image_dimension



%---------------------------------------------------------------------
function hist = make_header_data_history(origin, descrip)

hist.descrip = descrip;
hist.aux_file = 'none';
hist.qform_code = 0;
hist.sform_code = 0;
hist.quatern_b = 0;
hist.quatern_c = 0;
hist.quatern_d = 0;
hist.qoffset_x = 0;
hist.qoffset_y = 0;
hist.qoffset_z = 0;
hist.srow_x = zeros(1,4);
hist.srow_y = zeros(1,4);
hist.srow_z = zeros(1,4);
hist.intent_name = '';
hist.magic = '';
hist.originator = origin;

return;					% data_history

%  Save NIFTI dataset. Support both *.nii and *.hdr/*.img file extension.
%  If file extension is not provided, *.hdr/*.img will be used as default.
%
%  Usage: save_nii(nii, filename, [old_RGB])
%
%  nii.hdr - struct with NIFTI header fields.
%  nii.img - 3D (or 4D) matrix of NIFTI data.
%  filename - NIFTI file name.
%
%  old_RGB    - an optional boolean variable to handle special RGB data
%       sequence [R1 R2 ... G1 G2 ... B1 B2 ...] that is used only by
%       AnalyzeDirect (Analyze Software). Since both NIfTI and Analyze
%       file format use RGB triple [R1 G1 B1 R2 G2 B2 ...] sequentially
%       for each voxel, this variable is set to FALSE by default. If you
%       would like the saved image only to be opened by AnalyzeDirect
%       Software, set old_RGB to TRUE (or 1).
%
%  Tip: to change the data type, set nii.hdr.dime.datatype,
%	and nii.hdr.dime.bitpix to:
%
%     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN
%     1 Binary                         (ubit1, bitpix=1) % DT_BINARY
%     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8
%     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16
%     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32
%    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32
%    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
%    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64
%   128 Red-Green-Blue            (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24
%   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8
%   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16
%   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32
%  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
%  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64
%  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128
%  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128
%  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128
%
%  Part of this file is copied and modified under GNU license from
%  MRI_TOOLBOX developed by CNSP in Flinders University, Australia
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%  - "old_RGB" related codes in "save_nii.m" are added by Mike Harms (2006.06.28)
%
function save_nii(nii, fileprefix, old_RGB)

if ~exist('nii','var') | isempty(nii) | ~isfield(nii,'hdr') | ...
        ~isfield(nii,'img') | ~exist('fileprefix','var') | isempty(fileprefix)
    
    error('Usage: save_nii(nii, filename, [old_RGB])');
end

if ~exist('old_RGB','var'), old_RGB = 0; end

filetype = 1;

%  Note: fileprefix is actually the filename you want to save
%
if findstr('.nii', fileprefix)
    filetype = 2;
    fileprefix = strrep(fileprefix,'.nii','');
end

if findstr('.hdr', fileprefix)
    fileprefix = strrep(fileprefix,'.hdr','');
end

if findstr('.img', fileprefix)
    fileprefix = strrep(fileprefix,'.img','');
end

write_nii(nii, filetype, fileprefix, old_RGB);

if filetype == 1
    
    %  So earlier versions of SPM can also open it with correct originator
    %
    M=[[diag(nii.hdr.dime.pixdim(2:4)), -[nii.hdr.hist.originator(1:3) .* nii.hdr.dime.pixdim(2:4)]']; [0, 0, 0, 1]];
    save(fileprefix, 'M');
end

return					% save_nii


%-----------------------------------------------------------------------------------
function write_nii(nii, filetype, fileprefix, old_RGB)

hdr = nii.hdr;

if isfield(nii,'ext')
    ext = nii.ext;
    [ext, esize_total] = verify_nii_ext(ext);
else
    ext = [];
end

switch double(hdr.dime.datatype),
    case   1,
        hdr.dime.bitpix = int16(1 ); 
        precision = 'ubit1';
    case   2,
        hdr.dime.bitpix = int16(8 ); 
        precision = 'uint8';
    case   4,
        hdr.dime.bitpix = int16(16); 
        precision = 'int16';
    case   8,
        hdr.dime.bitpix = int16(32); 
        precision = 'int32';
    case  16,
        hdr.dime.bitpix = int16(32); 
        precision = 'float32';
    case  32,
        hdr.dime.bitpix = int16(64); 
        precision = 'float32';
    case  64,
        hdr.dime.bitpix = int16(64); 
        precision = 'float64';
    case 128,
        hdr.dime.bitpix = int16(24); 
        precision = 'uint8';
    case 256
        hdr.dime.bitpix = int16(8 ); 
        precision = 'int8';
    case 512
        hdr.dime.bitpix = int16(16); 
        precision = 'uint16';
    case 768
        hdr.dime.bitpix = int16(32); 
        precision = 'uint32';
    case 1024
        hdr.dime.bitpix = int16(64); 
        precision = 'int64';
    case 1280
        hdr.dime.bitpix = int16(64); 
        precision = 'uint64';
    case 1792,
        hdr.dime.bitpix = int16(128); 
        precision = 'float64';
    otherwise
        error('This datatype is not supported');
end

hdr.dime.glmax = round(double(max(nii.img(:))));
hdr.dime.glmin = round(double(min(nii.img(:))));

if filetype == 2
    fid = fopen(sprintf('%s.nii',fileprefix),'w');
    
    if fid < 0,
        msg = sprintf('Cannot open file %s.nii.',fileprefix);
        error(msg);
    end
    
    hdr.dime.vox_offset = 352;
    
    if ~isempty(ext)
        hdr.dime.vox_offset = hdr.dime.vox_offset + esize_total;
    end
    
    hdr.hist.magic = 'n+1';
    save_nii_hdr(hdr, fid);
    
    if ~isempty(ext)
        save_nii_ext(ext, fid);
    end
else
    fid = fopen(sprintf('%s.hdr',fileprefix),'w');
    
    if fid < 0,
        msg = sprintf('Cannot open file %s.hdr.',fileprefix);
        error(msg);
    end
    
    hdr.dime.vox_offset = 0;
    hdr.hist.magic = 'ni1';
    save_nii_hdr(hdr, fid);
    
    if ~isempty(ext)
        save_nii_ext(ext, fid);
    end
    
    fclose(fid);
    fid = fopen(sprintf('%s.img',fileprefix),'w');
end

ScanDim  = double(hdr.dime.dim(5));		% t
SliceDim = double(hdr.dime.dim(4));		% z
RowDim   = double(hdr.dime.dim(3));		% y
PixelDim = double(hdr.dime.dim(2));		% x
SliceSz  = double(hdr.dime.pixdim(4));
RowSz    = double(hdr.dime.pixdim(3));
PixelSz  = double(hdr.dime.pixdim(2));

x = 1:PixelDim;

if filetype == 2 & isempty(ext)
    skip_bytes = double(hdr.dime.vox_offset) - 348;
else
    skip_bytes = 0;
end

if double(hdr.dime.datatype) == 128
    
    %  RGB planes are expected to be in the 4th dimension of nii.img
    %
    if(size(nii.img,4)~=3)
        error(['The NII structure does not appear to have 3 RGB color planes in the 4th dimension']);
    end
    
    if old_RGB
        nii.img = permute(nii.img, [1, 2, 4, 3, 5]);
    else
        nii.img = permute(nii.img, [4, 1, 2, 3, 5]);
    end
end

%  For complex float32 or complex float64, voxel values
%  include [real, imag]
%
if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
    real_img = real(nii.img(:))';
    nii.img = imag(nii.img(:))';
    nii.img = [real_img; nii.img];
end

if skip_bytes
    fwrite(fid, ones(1,skip_bytes), 'uint8');
end

fwrite(fid, nii.img, precision);
%   fwrite(fid, nii.img, precision, skip_bytes);        % error using skip
fclose(fid);

return;					% write_nii

%  internal function

%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)

function save_nii_hdr(hdr, fid)

if ~exist('hdr','var') | ~exist('fid','var')
    error('Usage: save_nii_hdr(hdr, fid)');
end

if ~isequal(hdr.hk.sizeof_hdr,348),
    error('hdr.hk.sizeof_hdr must be 348.');
end

if hdr.hist.qform_code == 0 & hdr.hist.sform_code == 0
    hdr.hist.sform_code = 1;
    hdr.hist.srow_x(1) = hdr.dime.pixdim(2);
    hdr.hist.srow_x(2) = 0;
    hdr.hist.srow_x(3) = 0;
    hdr.hist.srow_y(1) = 0;
    hdr.hist.srow_y(2) = hdr.dime.pixdim(3);
    hdr.hist.srow_y(3) = 0;
    hdr.hist.srow_z(1) = 0;
    hdr.hist.srow_z(2) = 0;
    hdr.hist.srow_z(3) = hdr.dime.pixdim(4);
    hdr.hist.srow_x(4) = (1-hdr.hist.originator(1))*hdr.dime.pixdim(2);
    hdr.hist.srow_y(4) = (1-hdr.hist.originator(2))*hdr.dime.pixdim(3);
    hdr.hist.srow_z(4) = (1-hdr.hist.originator(3))*hdr.dime.pixdim(4);
end

write_header(hdr, fid);

return;					% save_nii_hdr


%---------------------------------------------------------------------
function write_header(hdr, fid)

%  Original header structures
%  struct dsr				/* dsr = hdr */
%       {
%       struct header_key hk;            /*   0 +  40       */
%       struct image_dimension dime;     /*  40 + 108       */
%       struct data_history hist;        /* 148 + 200       */
%       };                               /* total= 348 bytes*/

write_header_header_key(fid, hdr.hk);
write_header_image_dimension(fid, hdr.dime);
write_header_data_history(fid, hdr.hist);

%  check the file size is 348 bytes
%
fbytes = ftell(fid);

if ~isequal(fbytes,348),
    msg = sprintf('Header size is not 348 bytes.');
    warning(msg);
end

return;					% write_header


%---------------------------------------------------------------------
function write_header_header_key(fid, hk)

fseek(fid,0,'bof');

%  Original header structures
%  struct header_key                      /* header key      */
%       {                                /* off + size      */
%       int sizeof_hdr                   /*  0 +  4         */
%       char data_type[10];              /*  4 + 10         */
%       char db_name[18];                /* 14 + 18         */
%       int extents;                     /* 32 +  4         */
%       short int session_error;         /* 36 +  2         */
%       char regular;                    /* 38 +  1         */
%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
%       };                               /* total=40 bytes  */

fwrite(fid, hk.sizeof_hdr(1),    'int32');	% must be 348.

% data_type = sprintf('%-10s',hk.data_type);	% ensure it is 10 chars from left
% fwrite(fid, data_type(1:10), 'uchar');
pad = zeros(1, 10-length(hk.data_type));
hk.data_type = [hk.data_type  char(pad)];
fwrite(fid, hk.data_type(1:10), 'uchar');

% db_name   = sprintf('%-18s', hk.db_name);	% ensure it is 18 chars from left
% fwrite(fid, db_name(1:18), 'uchar');
pad = zeros(1, 18-length(hk.db_name));
hk.db_name = [hk.db_name  char(pad)];
fwrite(fid, hk.db_name(1:18), 'uchar');

fwrite(fid, hk.extents(1),       'int32');
fwrite(fid, hk.session_error(1), 'int16');
fwrite(fid, hk.regular(1),       'uchar');	% might be uint8

% fwrite(fid, hk.hkey_un0(1),    'uchar');
% fwrite(fid, hk.hkey_un0(1),    'uint8');
fwrite(fid, hk.dim_info(1),      'uchar');

return;					% header_key


%---------------------------------------------------------------------
function write_header_image_dimension(fid, dime)

%  Original header structures
%  struct image_dimension
%       {                                /* off + size      */
%       short int dim[8];                /* 0 + 16          */
%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
%       short int intent_code;   % short int unused1;   /* 28 + 2 */
%       short int datatype;              /* 30 + 2          */
%       short int bitpix;                /* 32 + 2          */
%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
%       float pixdim[8];                 /* 36 + 32         */
%			/*
%				pixdim[] specifies the voxel dimensions:
%				pixdim[1] - voxel width
%				pixdim[2] - voxel height
%				pixdim[3] - interslice distance
%				pixdim[4] - volume timing, in msec
%					..etc
%			*/
%       float vox_offset;                /* 68 + 4          */
%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
%       float scl_inter;   % float funused1;      /* 76 + 4 */
%       short slice_end;   % float funused2;      /* 80 + 2 */
%       char slice_code;   % float funused2;      /* 82 + 1 */
%       char xyzt_units;   % float funused2;      /* 83 + 1 */
%       float cal_max;                   /* 84 + 4          */
%       float cal_min;                   /* 88 + 4          */
%       float slice_duration;   % int compressed; /* 92 + 4 */
%       float toffset;   % int verified;          /* 96 + 4 */
%       int glmax;                       /* 100 + 4         */
%       int glmin;                       /* 104 + 4         */
%       };                               /* total=108 bytes */

fwrite(fid, dime.dim(1:8),        'int16');
fwrite(fid, dime.intent_p1(1),  'float32');
fwrite(fid, dime.intent_p2(1),  'float32');
fwrite(fid, dime.intent_p3(1),  'float32');
fwrite(fid, dime.intent_code(1),  'int16');
fwrite(fid, dime.datatype(1),     'int16');
fwrite(fid, dime.bitpix(1),       'int16');
fwrite(fid, dime.slice_start(1),  'int16');
fwrite(fid, dime.pixdim(1:8),   'float32');
fwrite(fid, dime.vox_offset(1), 'float32');
fwrite(fid, dime.scl_slope(1),  'float32');
fwrite(fid, dime.scl_inter(1),  'float32');
fwrite(fid, dime.slice_end(1),    'int16');
fwrite(fid, dime.slice_code(1),   'uchar');
fwrite(fid, dime.xyzt_units(1),   'uchar');
fwrite(fid, dime.cal_max(1),    'float32');
fwrite(fid, dime.cal_min(1),    'float32');
fwrite(fid, dime.slice_duration(1), 'float32');
fwrite(fid, dime.toffset(1),    'float32');
fwrite(fid, dime.glmax(1),        'int32');
fwrite(fid, dime.glmin(1),        'int32');

return;					% image_dimension


%---------------------------------------------------------------------
function write_header_data_history(fid, hist)

% Original header structures
%struct data_history
%       {                                /* off + size      */
%       char descrip[80];                /* 0 + 80          */
%       char aux_file[24];               /* 80 + 24         */
%       short int qform_code;            /* 104 + 2         */
%       short int sform_code;            /* 106 + 2         */
%       float quatern_b;                 /* 108 + 4         */
%       float quatern_c;                 /* 112 + 4         */
%       float quatern_d;                 /* 116 + 4         */
%       float qoffset_x;                 /* 120 + 4         */
%       float qoffset_y;                 /* 124 + 4         */
%       float qoffset_z;                 /* 128 + 4         */
%       float srow_x[4];                 /* 132 + 16        */
%       float srow_y[4];                 /* 148 + 16        */
%       float srow_z[4];                 /* 164 + 16        */
%       char intent_name[16];            /* 180 + 16        */
%       char magic[4];   % int smin;     /* 196 + 4         */
%       };                               /* total=200 bytes */

% descrip     = sprintf('%-80s', hist.descrip);     % 80 chars from left
% fwrite(fid, descrip(1:80),    'uchar');
pad = zeros(1, 80-length(hist.descrip));
hist.descrip = [hist.descrip  char(pad)];
fwrite(fid, hist.descrip(1:80), 'uchar');

% aux_file    = sprintf('%-24s', hist.aux_file);    % 24 chars from left
% fwrite(fid, aux_file(1:24),   'uchar');
pad = zeros(1, 24-length(hist.aux_file));
hist.aux_file = [hist.aux_file  char(pad)];
fwrite(fid, hist.aux_file(1:24), 'uchar');

fwrite(fid, hist.qform_code,    'int16');
fwrite(fid, hist.sform_code,    'int16');
fwrite(fid, hist.quatern_b,   'float32');
fwrite(fid, hist.quatern_c,   'float32');
fwrite(fid, hist.quatern_d,   'float32');
fwrite(fid, hist.qoffset_x,   'float32');
fwrite(fid, hist.qoffset_y,   'float32');
fwrite(fid, hist.qoffset_z,   'float32');
fwrite(fid, hist.srow_x(1:4), 'float32');
fwrite(fid, hist.srow_y(1:4), 'float32');
fwrite(fid, hist.srow_z(1:4), 'float32');

% intent_name = sprintf('%-16s', hist.intent_name);	% 16 chars from left
% fwrite(fid, intent_name(1:16),    'uchar');
pad = zeros(1, 16-length(hist.intent_name));
hist.intent_name = [hist.intent_name  char(pad)];
fwrite(fid, hist.intent_name(1:16), 'uchar');

% magic	= sprintf('%-4s', hist.magic);		% 4 chars from left
% fwrite(fid, magic(1:4),           'uchar');
pad = zeros(1, 4-length(hist.magic));
hist.magic = [hist.magic  char(pad)];
fwrite(fid, hist.magic(1:4),        'uchar');

return;					% data_history
