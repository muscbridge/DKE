function dke_preprocess_prisma(basedir, fn_params)

% dke_preprocess_prisma Preprocess diffusional kurtosis imaging data (convert Siemens Prisma DICOM to NIfTI, coregister diffusion-weighted images, average images)

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
% Read parameter file
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
% Set file names
%
% current_4D_file is the name of the file being preprocessed with denoising,
% Rician bias correction, and Gibbs ringing artifact correction (depending
% on settings of flags in the preprocess_options struct) and changes
% after each preprocessing step
%--------------------------------------------------------------------------

dki_file = 'dki.nii';
current_4D_file = '4D.nii';

%--------------------------------------------------------------------------
% Convert DICOM images to NIfTI format
%--------------------------------------------------------------------------

b12root = fullfile(basedir, 'intermediate_processing');
mkdir(b12root);

nifti_dir = fullfile(b12root, 'nifti');
mkdir(nifti_dir);

b0_dir = fullfile(nifti_dir, 'B0');
mkdir(b0_dir);

eval(['!/Applications/MRIcroGL/dcm2niix -f %d ' basedir])

% move new NIfTI files to subdirectories of nifti_dir

% this part is needed when you have an additional separate b0 sequence
%     in=fullfile(basedir, '*B0*.nii');
%     copyfile(in, b0_dir);

dki1_dir = fullfile(nifti_dir, 'DKI1');
mkdir(dki1_dir);
file1 = [options.series_description{1} '.nii'];
in=fullfile(basedir, file1);
if ~exist(in, 'file')
    error_msg='The series_description parameter must match the series descriptions in the DICOM headers';
    error('Input file %s does not exist.\n%s', in, error_msg)
end
out=fullfile(dki1_dir, dki_file);
copyfile(in,out);

%--------------------------------------------------------------------------
% Read b values from the .bval file
%--------------------------------------------------------------------------

bval_file1 = [options.series_description{1} '.bval'];
in=fullfile(basedir, bval_file1);
if ~exist(in, 'file')
    error_msg='Could not read b values.';
    error('Input file %s does not exist.\n%s', in, error_msg)
end
bval_file_id = fopen(in);
bvals = textscan(bval_file_id, '%s');
fclose(bval_file_id);

%--------------------------------------------------------------------------
% Switch to dki1_dir and rename dki_file to current_4D_file
%--------------------------------------------------------------------------

% for separate b0
% Vdir = dir(fullfile(b0_dir, '*.nii'));
% V=fullfile(b0_dir, Vdir(1).name);
% Vo = spm_file_split(V, b0_dir);

cd(dki1_dir);
%   Vdir = dir('*.nii');
%   V=fullfile(dki1_dir, Vdir(1).name);
%   Vo = spm_file_split(V,dki1_dir);

in=fullfile(dki1_dir, dki_file);
out=fullfile(dki1_dir, current_4D_file);
movefile(in,out);

%--------------------------------------------------------------------------
% Denoise if denoise_flag = 1
%--------------------------------------------------------------------------

if options.denoise_flag == 1
    fprintf('Denoising data with dwidenoise (MRtrix)...  ')

    denoised_file = append_to_name(current_4D_file, '_dn');
    noise_file = append_to_name(current_4D_file, '_noise');
    residual_file = append_to_name(current_4D_file, '_res');

    command=['/usr/local/mrtrix3/bin/dwidenoise ' current_4D_file  ' ' denoised_file ' -noise ' noise_file];
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
    fprintf('Calculating residuals with mrcalc (MRtrix)...  ')
    command=['/usr/local/mrtrix3/bin/mrcalc ' current_4D_file  ' ' denoised_file ' -subtract ' residual_file];
    [status,cmdout] = system(command);
    if status == 0
        fprintf('done.\n')
        fprintf('\t- residual output file is %s.\n', residual_file)
    else
        fprintf('\nAn error occurred. Output of mrcalc command:\n');
        cmdout
        error('Error running mrcalc.')
    end
    current_4D_file = denoised_file;  % The file being processed is now the denoised file
else
    fprintf('Not denoising data\n')
end

%--------------------------------------------------------------------------
% Correct for Rician noise bias if rician_corr_flag = 1 and denoise_flag = 1
%--------------------------------------------------------------------------

if options.denoise_flag == 1
    if options.rician_corr_flag == 1
        fprintf('Correcting for Rician noise bias...  ')

        rician_corrected_file = append_to_name(current_4D_file, '_rc');

        hdr_DN = spm_vol(current_4D_file);
        DN = spm_read_vols(hdr_DN);
        noise = spm_read_vols(spm_vol(noise_file));
        DN = sqrt(DN.^2 - noise.^2);
        DN = real(DN);
        DN(isnan(DN)) = 0;
        make_4D_nii(hdr_DN, DN, rician_corrected_file);
        fprintf('done.\n')
        fprintf('\t- Rician noise bias corrected output file is %s.\n', rician_corrected_file)
        current_4D_file = rician_corrected_file;  % The file being processed is now the Rician bias corrected file
    else
        fprintf('Not correcting for Rician noise bias\n')
    end
else
    if options.rician_corr_flag == 1
        fprintf('Not correcting for Rician noise bias (denoise_flag = 0)\n')
    end
end

%--------------------------------------------------------------------------
% Correct for Gibbs ringing artifact
%--------------------------------------------------------------------------

if options.gibbs_corr_flag == 1
    fprintf('NOTE: You should not use mrdegibbs on partial Fourier data\n')
    fprintf('Correcting for Gibbs ringing artifact with mrdegibbs (MRtrix)...  ')

    gibbs_corrected_file = append_to_name(current_4D_file, '_gr');

    command=['/usr/local/mrtrix3/bin/mrdegibbs ' current_4D_file ' ' gibbs_corrected_file];
    [status,cmdout] = system(command);
    if status == 0
        fprintf('done.\n')
        fprintf('\t- Gibbs ringing artifact corrected output file is %s.\n', gibbs_corrected_file)
    else
        fprintf('\nAn error occurred. Output of mrdegibbs command:\n');
        cmdout
        error('Error running mrdegibbs.')
    end
    current_4D_file = gibbs_corrected_file;  % The file being processed is now the Gibbs corrected file

else
    fprintf('Not correcting for Gibbs artifact\n')
end

%--------------------------------------------------------------------------
% Split preprocessed 4D file into several 3D dki_vol#.nii files
%--------------------------------------------------------------------------

copyfile(current_4D_file, dki_file);

% for separate b0
% Vdir = dir(fullfile(b0_dir, '*.nii'));
% V=fullfile(b0_dir, Vdir(1).name);
% Vo = spm_file_split(V, b0_dir);

V=fullfile(dki1_dir, dki_file);
Vo = spm_file_split(V, dki1_dir);

%--------------------------------------------------------------------------
% Rename NIfTI images -- append b values to file names (before .nii)
%--------------------------------------------------------------------------

list=dir(fullfile(dki1_dir,'*00*.nii'));

num_images = length(list);
num_bvals = length(bvals{1});
if num_images ~= num_bvals
    error('The number of 3D NIfTI images in %s (%d) does not match the number of b values in the .bval file (%d)', dki1_dir, num_images, num_bvals)
end

for k=1:length(list)
    app_string = ['_b' bvals{1}{k}];
    new_name = append_to_name(list(k).name, app_string);
    movefile(list(k).name, new_name);
end

% move b0s to folder 'intermediate_processing/nifti/b0'

clear list
list=dir(fullfile(dki1_dir,'*b0*.nii'));

for l=1:length(list)
    in=fullfile(dki1_dir,list(l).name);
    out=fullfile(b0_dir, list(l).name);
    movefile(in,out);
end

%--------------------------------------------------------------------------
% Make gradient file
%--------------------------------------------------------------------------

dke_dir = fullfile(b12root, 'dke');
mkdir(dke_dir);

cd(basedir);
name=dir('*.bvec');
A=importdata([name(1).name]);
B=A(:,any(A));
Gradient=B';
Gradient1=Gradient(1:(round(end/2)),:);
save(fullfile(dke_dir, 'gradient_dke.txt'),'Gradient1','-ASCII')

%--------------------------------------------------------------------------
% Coregister b0s to b0
%--------------------------------------------------------------------------

%ONLY USE WITH INTERLEAVED B0S
% dirb0=dir(fullfile(b0_dir, '*.nii'));
% dirdki=dir(fullfile(dki1_dir, '*.nii'));
% 
% fprintf('Co-registering images...\n')
% fn_source = fullfile(b0_dir, dirb0(1).name); % source file is the first b = 0 image in the series returned by the operating system
% 
%             fn_target = fullfile(dki1_dir, dirdki(1).name);
%             M=coregister(fn_target, fn_source, b0_dir, '.nii');
% 
% delete(fullfile(dki1_dir, 'r*.nii'))
% movefile(fullfile(b0_dir, 'r*.nii'), fullfile(b12root, 'nifti/B0_coreg'));
% fprintf('Co-registration complete.\n')

%--------------------------------------------------------------------------
% Average b0's
%--------------------------------------------------------------------------
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

hdr.dt=[16 0];
hdr.fname = fullfile(combined_dir, 'b0_avg.nii');
imgavg(isnan(imgavg))=0;
spm_write_vol(hdr, imgavg);

%--------------------------------------------------------------------------
% Move DKI images and make 4D NIfTI images
%--------------------------------------------------------------------------
copyfile(fullfile(dki1_dir, '*00*'), combined_dir);
files = dir(fullfile(combined_dir, '*.nii'));
make_4D_nii(combined_dir, {files.name}, '4D.nii');
movefile(fullfile(combined_dir, '4D.nii'), fullfile(dke_dir, '4D.nii'))

img=spm_read_vols(spm_vol(fullfile(dke_dir, '4D.nii')));
img(isnan(img))=0;
make_4D_nii(spm_vol(fullfile(dke_dir, '4D.nii')), img, '4D.nii');

%--------------------------------------------------------------------------
% Return to original working directory
%--------------------------------------------------------------------------

cd(orig_dir);


fprintf('complete.\n')


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

