function dke_preprocess_dicom(basedir, fn_params)

% dke_preprocess Preprocess diffusional kurtosis imaging data (convert DICOM to NIfTI, coregister diffusion-weighted images, average images)

% Author: Ali Tabesh
% Version: 2.5.1
% Last modified: 12/26/12

% Portions of this code are based on John's Gems:
%
% http://www-personal.umich.edu/~nichols/JohnsGems.html

if nargin ~= 2
    fprintf('\nUsage: dke_preprocess_dicom basedir paramsfile\n')
    fprintf('\nbasedir  input dicom folder')
    fprintf('\nparamsfile  DKE processing parameters file; see included example\n\n')
    return
end

warning('off','all')

%--------------------------------------------------------------------------
% read parameters file
%--------------------------------------------------------------------------


fid=fopen(fn_params); %EM
file = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', ''); %EM
for i = 1:length(file{1}) %EM
    eval(file{1}{i})%EM
end
fclose(fid);

options = preprocess_options;

if options.extra_b0 ~= 0 && options.extra_b0 ~= 1
    error('Invalid ''extra_b0'' parameter! ''extra_b0'' must be 0 or 1.')
end
    
%--------------------------------------------------------------------------
% convert DICOM images to NIfTI format
%--------------------------------------------------------------------------

fprintf('Converting input DICOM images to NIfTI... ')

% create output folders

mkdir(fullfile(basedir, 'intermediate_processing'))

for iavg = 1:options.navg
    outputdir = fullfile(basedir, 'intermediate_processing', ['dki_avg' num2str(iavg) '_nii']);
    mkdir(outputdir)
end

if options.extra_b0 == 1
    outputdir = fullfile(basedir, 'intermediate_processing', 'dki_b0_nii');
    mkdir(outputdir)
end

dcm2nii(basedir, options.series_description, options.extra_b0)

fprintf('complete.\n')

%--------------------------------------------------------------------------
% co-register all b = 0 images to b = 0 image in DKI1 series
%--------------------------------------------------------------------------

if options.coreg_flag == 1

    fprintf('Co-registering images...\n')

    fn_target = fullfile(basedir, 'intermediate_processing', 'dki_avg1_nii', 'dki_0.nii');

    for iavg = 1:options.navg
        fn_source = fullfile(basedir, 'intermediate_processing', ['dki_avg' num2str(iavg) '_nii'], 'dki_0.nii');
        folder_other = fullfile(basedir, 'intermediate_processing', ['dki_avg' num2str(iavg) '_nii']);
        fn_other_filt = '^dki*.*';
        coregister(fn_target, fn_source, folder_other, fn_other_filt);
        coregdir = fullfile(basedir, 'intermediate_processing', ['dki_avg' num2str(iavg) '_coreg']);
        mkdir(coregdir)
        movefile(fullfile(folder_other, 'rdki*.*'), coregdir, 'f')
    end
    
    if options.extra_b0 == 1

        folder_other = fullfile(basedir, 'intermediate_processing', 'dki_b0_nii');
        list = dir(fullfile(folder_other, 'dki_0_*_1.nii'));
        for iseries = 1:size(list, 1)
            idx = strfind(list(iseries).name, '_');
            series_str = list(iseries).name(idx(2)+1:idx(3)-1);
            fn_source = fullfile(folder_other, list(iseries).name);   % source file is the first b = 0 image in the series returned by the operating system
            fn_other_filt = ['^dki_0_' series_str '_*.*'];
            coregister(fn_target, fn_source, folder_other, fn_other_filt);
        end
        coregdir = fullfile(basedir, 'intermediate_processing', 'dki_b0_coreg');
        mkdir(coregdir)
        movefile(fullfile(folder_other, 'rdki*.*'), coregdir, 'f')

    end
    
    fprintf('Co-registration complete.\n')

elseif options.coreg_flag == 0

    for iavg = 1:options.navg
        sourcedir = fullfile(basedir, 'intermediate_processing', ['dki_avg' num2str(iavg) '_nii']);
        list = dir(fullfile(sourcedir, 'dki*.nii'));
        coregdir = fullfile(basedir, 'intermediate_processing', ['dki_avg' num2str(iavg) '_coreg']);
        mkdir(coregdir)
        for i = 1:size(list, 1)
            fn_source = fullfile(sourcedir, list(i).name);
            fn_coreg = fullfile(coregdir, ['r' list(i).name]);
            copyfile(fn_source, fn_coreg)
        end
    end

    if options.extra_b0 == 1
        
        sourcedir = fullfile(basedir, 'intermediate_processing', 'dki_b0_nii');
        list = dir(fullfile(sourcedir, 'dki*.nii'));
        coregdir = fullfile(basedir, 'intermediate_processing', 'dki_b0_coreg');
        mkdir(coregdir)
        for i = 1:size(list, 1)
            fn_source = fullfile(sourcedir, list(i).name);
            fn_coreg = fullfile(coregdir, ['r' list(i).name]);
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

% average non-b = 0 images
list = dir(fullfile(folder_avg1, 'rdki_*_*.nii'));
fn_array=cell(1,options.navg);
for i = 1:size(list, 1)
    fn_array{1} = fullfile(folder_avg1, list(i).name);
    for iavg = 2:options.navg
        fn_array{iavg} = fullfile(basedir, 'intermediate_processing', ['dki_avg' num2str(iavg) '_coreg'], list(i).name);
    end
    image_avg(fn_array, fullfile(folder_output, list(i).name))
end

% average b = 0 images
for iavg = 1:options.navg
    fn_array{iavg} = fullfile(basedir, 'intermediate_processing', ['dki_avg' num2str(iavg) '_coreg'], 'rdki_0.nii');
end

if options.extra_b0 == 1

    folder_b0 = fullfile(basedir, 'intermediate_processing', 'dki_b0_coreg');
    list = dir(fullfile(folder_b0, 'rdki_*.nii'));
    for ib0 = 1:size(list, 1)
        fn_array{options.navg+ib0} = fullfile(basedir, 'intermediate_processing', 'dki_b0_coreg', list(ib0).name);
    end

end

image_avg(fn_array, fullfile(folder_output, 'rdki_0.nii'))

fprintf('complete.\n')

function dcm2nii(folder_dcm, series_description, extra_b0)

P = spm_select('fplist', folder_dcm, '.*');  % fplist: list w/ full path
hdr = spm_dicom_headers(P);

current_path = pwd;
cd(folder_dcm)
spm_dicom_convert_series_description(hdr, 'all', 'flat', 'nii', series_description, extra_b0);
cd(current_path)

%--------------------------------------------------------------------------
% co-register two b = 0 images
%--------------------------------------------------------------------------

function coregister(fn_target, fn_source, folder_other, fn_other_filt)

fn_other = spm_select('fplist', folder_other, fn_other_filt);

% coregistration and reslicing parameters
estflg.cost_fun = 'nmi';
estflg.sep      = [4 2];
estflg.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
estflg.fwhm     = [7 7];
wrtflg.interp   = 1;
wrtflg.mean    = 0;

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
