function dke_preprocess_prisma(basedir, fn_params)

% dke_preprocess_prisma Preprocess diffusional kurtosis imaging data
% (convert DICOM to NIfTI; optionally denoise, correct for Rician noise bias,
% and correct for Gibbs ringing artifact; coregister diffusion-weighted images;
% average images)

%--------------------------------------------------------------------------
% Check inputs
%--------------------------------------------------------------------------

if nargin ~= 2
    fprintf('\nUsage: dke_preprocess_prisma basedir paramsfile\n')
    fprintf('\nbasedir  input Prisma data folder')
    fprintf('\nparamsfile  DKE processing parameters file; see included example\n\n')
    return
end

if ~exist(fn_params, 'file')
    fprintf('\n')
    fprintf('Input parameter file %s does not exist!\n\n', fn_params)
    return
end

warning('off','all')

%--------------------------------------------------------------------------
% Read parameter file, check settings of flags
%--------------------------------------------------------------------------

fid=fopen(fn_params);
file = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
for i = 1:length(file{1})
    eval(file{1}{i})
end
fclose(fid);

options = preprocess_options;

if options.denoise_flag ~= 0 && options.denoise_flag ~= 1
    error('Invalid ''denoise_flag'' parameter! ''denoise_flag'' must be 0 or 1.')
end

if options.rician_corr_flag ~= 0 && options.rician_corr_flag ~= 1
    error('Invalid ''rician_corr_flag'' parameter! ''rician_corr_flag'' must be 0 or 1.')
end

if options.gibbs_corr_flag ~= 0 && options.gibbs_corr_flag ~= 1
    error('Invalid ''gibbs_corr_flag'' parameter! ''gibbs_corr_flag'' must be 0 or 1.')
end

if options.rician_corr_flag == 1 && options.denoise_flag == 0
    fprintf('Will not correct for Rician noise bias (denoise_flag = 0)\n')
    fprintf('- Setting Rician bias correction flag to 0\n')
    options.rician_corr_flag = 0;
end

expected_num_series = options.navg;      % Number of DKI series
if options.extra_b0
    expected_num_series = expected_num_series+1;  % Number of DKI + b=0 series
end

series_description_length = length(options.series_description);
if series_description_length ~= expected_num_series
    error('Length of series_description parameter does not match expected number of series.\nLength of series_description parameter = %d\nnavg = %d\nextra_b0 = %d\nExpected number of series = %d', series_description_length, options.navg, options.extra_b0, expected_num_series)
end

%--------------------------------------------------------------------------
% Save current working directory
%--------------------------------------------------------------------------

orig_dir = pwd;

%--------------------------------------------------------------------------
% Convert DICOM images to NIfTI format
%--------------------------------------------------------------------------

b12root = fullfile(basedir, 'intermediate_processing');
mkdir(b12root);

nifti_dir = fullfile(b12root, 'nifti');
mkdir(nifti_dir);

fprintf('Converting from DICOM to NIfTI with dcm2niix (MRIcroGL)...  ')
command=['/Applications/MRIcroGL/dcm2niix -f %d ' basedir];
[status,cmdout] = system(command);
if status == 0
    fprintf('done.\n')
    list=dir(fullfile(basedir,'*.nii'));
    fprintf('\t- NIfTI output:\n');
    for i=1:length(list)
        fprintf('\t  %s\n', list(i).name)
    end
else
    fprintf('\nAn error occurred. Output of dcm2niix command:\n');
    cmdout
    error('Error running dcm2niix.')
end

% move new NIfTI files to subdirectories of nifti_dir

% this part is needed when you have an additional separate b0 sequence
%     in=fullfile(basedir, '*B0*.nii');
%     copyfile(in, b0_dir);

%--------------------------------------------------------------------------
% Copy extra b=0 file into a separate subdirectory and split into 3D images
%--------------------------------------------------------------------------

if options.extra_b0
    extra_b0_dir = fullfile(nifti_dir, 'extra_b0');
    mkdir(extra_b0_dir);
    extra_b0_file = [options.series_description{end} '.nii'];
    extra_b0_file_path = fullfile(basedir, extra_b0_file);
    fprintf('Copying %s to %s\n', extra_b0_file_path, extra_b0_dir)
    copyfile(extra_b0_file_path, extra_b0_dir);
    fprintf('and splitting into 3D volumes\n')
    extra_b0_3D_files = spm_file_split(extra_b0_file_path, extra_b0_dir);
end

%--------------------------------------------------------------------------
% Copy file to preprocess to a new subdirectory
%--------------------------------------------------------------------------

iavg = 1;

dki_dir = fullfile(nifti_dir, ['dki_avg' num2str(iavg)]);
mkdir(dki_dir);

b0_dir = fullfile(nifti_dir, ['dki_avg' num2str(iavg) '_b0']);
mkdir(b0_dir);

% current_4D_file is the name of the file being preprocessed with denoising,
% Rician bias correction, and Gibbs ringing artifact correction (depending
% on settings of flags in the preprocess_options struct) and changes
% after each preprocessing step
current_4D_file = [options.series_description{iavg} '.nii'];
in=fullfile(basedir, current_4D_file);
if ~exist(in, 'file')
    error_msg='The series_description parameter must match the series descriptions in the DICOM headers';
    error('Input file %s does not exist.\n%s', in, error_msg)
end
out=fullfile(dki_dir, current_4D_file);
copyfile(in,out);

fprintf('Preprocessing image file %s\n', out);

%--------------------------------------------------------------------------
% Read b values from the .bval file
%--------------------------------------------------------------------------

bval_file = [options.series_description{iavg} '.bval'];
in=fullfile(basedir, bval_file);
if ~exist(in, 'file')
    error_msg='Could not read b values.';
    error('Input file %s does not exist.\n%s', in, error_msg)
end
bval_file_id = fopen(in);
bvals = textscan(bval_file_id, '%s');
fclose(bval_file_id);

%--------------------------------------------------------------------------
% Switch to dki_dir
%--------------------------------------------------------------------------

% for separate b0
% Vdir = dir(fullfile(b0_dir, '*.nii'));
% V=fullfile(b0_dir, Vdir(1).name);
% Vo = spm_file_split(V, b0_dir);

cd(dki_dir);

%   Vdir = dir('*.nii');
%   V=fullfile(dki_dir, Vdir(1).name);
%   Vo = spm_file_split(V,dki_dir);

%--------------------------------------------------------------------------
% Denoise, correct for Rician noise bias, correct for Gibbs ringing artifact
%--------------------------------------------------------------------------

if options.denoise_flag
    [current_4D_file, noise_file] = denoise(current_4D_file);
end

if options.rician_corr_flag
    current_4D_file = rician_bias_correct(current_4D_file, noise_file);
end

if options.gibbs_corr_flag
    current_4D_file = gibbs_ringing_correct(current_4D_file);
end

%--------------------------------------------------------------------------
% Split preprocessed 4D file into several 3D dki_vol#.nii files
%--------------------------------------------------------------------------

fprintf('- Splitting 4D file into 3D files\n')

dki_file = 'dki.nii';
copyfile(current_4D_file, dki_file);

% for separate b0
% Vdir = dir(fullfile(b0_dir, '*.nii'));
% V=fullfile(b0_dir, Vdir(1).name);
% Vo = spm_file_split(V, b0_dir);

V=fullfile(dki_dir, dki_file);
Vo = spm_file_split(V, dki_dir);

%--------------------------------------------------------------------------
% Rename NIfTI images -- append b values to file names (before .nii)
%--------------------------------------------------------------------------

fprintf('- Renaming image files to include b values\n')

num_images = length(Vo);
num_bvals = length(bvals{1});
if num_images ~= num_bvals
    error('The number of 3D NIfTI images in %s (%d) does not match the number of b values in the .bval file (%d)', dki_dir, num_images, num_bvals)
end

for k=1:length(Vo)
    app_string = ['_b' bvals{1}{k}];
    new_name = append_to_name(Vo(k).fname, app_string);
    movefile(Vo(k).fname, new_name);
end

%--------------------------------------------------------------------------
% Move b0s to a separate folder
%--------------------------------------------------------------------------

list=dir(fullfile(dki_dir,'*b0*.nii'));

for l=1:length(list)
    in=fullfile(dki_dir,list(l).name);
    out=fullfile(b0_dir, list(l).name);
    movefile(in,out);
end

%--------------------------------------------------------------------------
% Coregister b0s to b0
%--------------------------------------------------------------------------

%ONLY USE WITH INTERLEAVED B0S
% dirb0=dir(fullfile(b0_dir, '*.nii'));
% dirdki=dir(fullfile(dki_dir, '*.nii'));
% 
% fprintf('Co-registering images...\n')
% fn_source = fullfile(b0_dir, dirb0(1).name); % source file is the first b = 0 image in the series returned by the operating system
% 
%             fn_target = fullfile(dki_dir, dirdki(1).name);
%             M=coregister(fn_target, fn_source, b0_dir, '.nii');
% 
% delete(fullfile(dki_dir, 'r*.nii'))
% movefile(fullfile(b0_dir, 'r*.nii'), fullfile(b12root, 'nifti/B0_coreg'));
% fprintf('Co-registration complete.\n')

%--------------------------------------------------------------------------
% Average b0's
%--------------------------------------------------------------------------

fprintf('- Averaging b=0 files\n')

list = dir(fullfile(b0_dir, '*.nii'));
hdr = spm_vol(fullfile(b0_dir, list(1).name));

imgavg = spm_read_vols(hdr);
for j = 2:length(list)
    hdr = spm_vol(fullfile(b0_dir, list(j).name));
    img = spm_read_vols(hdr);
    imgavg = imgavg + img;
end

imgavg = imgavg / (length(list));

combined_dir = fullfile(nifti_dir, 'combined');
mkdir(combined_dir);

hdr.dt=[64 0];
hdr.fname = fullfile(combined_dir, 'b0_avg.nii');
imgavg(isnan(imgavg))=0;
spm_write_vol(hdr, imgavg);

%--------------------------------------------------------------------------
% Move DKI images and make 4D NIfTI images
%--------------------------------------------------------------------------

copyfile(fullfile(dki_dir, '*00*'), combined_dir);
cd(combined_dir);
files = dir('*.nii');

fprintf('Assembling 3D volumes into final 4D image\n')

final_dir = fullfile(b12root, 'dke');
mkdir(final_dir);
final_image = fullfile(final_dir, '4D.nii');

for j = 1:length(files)
    hdr = spm_vol(files(j).name);
    img = spm_read_vols(hdr);
    img(isnan(img)) = 0;
    hdr.fname = final_image;
    hdr.n = [j, 1];
    hdr.dt=[64 0];
    spm_write_vol(hdr, img);
end

fprintf('\t- Image is %s\n', final_image)

%--------------------------------------------------------------------------
% Make gradient file
%--------------------------------------------------------------------------

fprintf('Creating gradient vector file\n')

cd(basedir);
bvec_file = [options.series_description{1} '.bvec'];
A=importdata(bvec_file);
B=A(:,any(A));
Gradient=B';
Gradient1=Gradient(1:(round(end/2)),:);
gradient_file = fullfile(final_dir, 'gradient_dke.txt');
save(gradient_file, 'Gradient1', '-ASCII')
fprintf('\t- Gradient vector file is %s\n', gradient_file)

%--------------------------------------------------------------------------
% Return to original working directory
%--------------------------------------------------------------------------

cd(orig_dir);


%--------------------------------------------------------------------------
% Co-register two b = 0 images
%--------------------------------------------------------------------------

function M=coregister(fn_target, fn_source, folder_other, fn_other_filt)

fn_other = spm_select('fplist', folder_other, fn_other_filt);

% coregistration and reslicing parameters
estflg.cost_fun = 'nmi';
estflg.sep      = [4 2];
estflg.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
estflg.fwhm     = [7 7];
wrtflg        = spm_get_defaults('realign.write');
wrtflg.interp   = 1;
wrtflg.which    = [2 0];

hdr_trg = spm_vol(fn_target);
hdr_src = spm_vol(fn_source);

x  = spm_coreg(hdr_trg, hdr_src, estflg);
M  = inv(spm_matrix(x));
MM = zeros(4, 4, size(fn_other, 1));
for j=1:size(fn_other, 1)
    MM(:,:,j) = spm_get_space(deblank(fn_other(j,:)));
end
for j = 1:size(fn_other, 1)
    spm_get_space(deblank(fn_other(j,:)), M*MM(:,:,j));
end;

fn_other = char(fn_target, fn_other);
spm_reslice(fn_other, wrtflg);


%--------------------------------------------------------------------------
% Generate a file name from appending a string to the input file name
% before the extension
%--------------------------------------------------------------------------
function newname = append_to_name(oldname, appstring)

[old_path, old_fn, old_ext] = fileparts(oldname);
newname = fullfile(old_path, [old_fn appstring old_ext]);


%--------------------------------------------------------------------------
% Denoise an image, estimate the noise map, and calculate the residual image
%--------------------------------------------------------------------------
function [denoised_file, noise_file] = denoise(image_file)

denoised_file = append_to_name(image_file, '_dn');
noise_file = append_to_name(image_file, '_noise');
residual_file = append_to_name(image_file, '_res');

fprintf('- Denoising data with dwidenoise (MRtrix)...  ')
command=['/usr/local/mrtrix3/bin/dwidenoise ' image_file  ' ' denoised_file ' -noise ' noise_file];
[status,cmdout] = system(command);
if status == 0
    fprintf('done.\n')
    fprintf('\t- denoised output file is %s.\n', denoised_file)
    fprintf('\t- noise estimate output file is %s.\n', noise_file)
else
    fprintf('\nAn error occurred. Output of dwidenoise command:\n');
    cmdout
    error('Error running dwidenoise.')
end

fprintf('- Calculating residuals with mrcalc (MRtrix)...  ')
command=['/usr/local/mrtrix3/bin/mrcalc ' image_file  ' ' denoised_file ' -subtract ' residual_file];
[status,cmdout] = system(command);
if status == 0
    fprintf('done.\n')
    fprintf('\t- residual output file is %s.\n', residual_file)
else
    fprintf('\nAn error occurred. Output of mrcalc command:\n');
    cmdout
    error('Error running mrcalc.')
end


%--------------------------------------------------------------------------
% Correct for Rician noise bias
% Set voxels where noise > signal to 0
%--------------------------------------------------------------------------
function rician_corrected_file = rician_bias_correct(image_file, noise_file)

rician_corrected_file = append_to_name(image_file, '_rc');

fprintf('- Correcting for Rician noise bias...  ')
hdr = spm_vol(image_file);
signal = spm_read_vols(hdr);
noise = spm_read_vols(spm_vol(noise_file));
signal = sqrt(signal.^2 - noise.^2);
signal = real(signal);
signal(isnan(signal)) = 0;

hdr(1).fname = rician_corrected_file;
hdr(1).dt = [64 0];
for j = 1:length(hdr)
    hdr(1).n = [j, 1];
    spm_write_vol(hdr(1), signal(:, :, :, j));
end

fprintf('done.\n')
fprintf('\t- Rician noise bias corrected output file is %s.\n', rician_corrected_file)


%--------------------------------------------------------------------------
% Correct for Gibbs ringing artifact
%--------------------------------------------------------------------------
function gibbs_corrected_file = gibbs_ringing_correct(image_file)

gibbs_corrected_file = append_to_name(image_file, '_gr');

fprintf('- NOTE: You should not use mrdegibbs on partial Fourier data\n')
fprintf('- Correcting for Gibbs ringing artifact with mrdegibbs (MRtrix)...  ')
command=['/usr/local/mrtrix3/bin/mrdegibbs -datatype float64 ' image_file ' ' gibbs_corrected_file];
[status,cmdout] = system(command);
if status == 0
    fprintf('done.\n')
    fprintf('\t- Gibbs ringing artifact corrected output file is %s.\n', gibbs_corrected_file)
else
    fprintf('\nAn error occurred. Output of mrdegibbs command:\n');
    cmdout
    error('Error running mrdegibbs.')
end
