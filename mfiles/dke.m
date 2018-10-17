function dke(fn_params)

% dke Diffusional Kurtosis Estimator
%   inputs are dicom or nifti diffusion-weighted images 
%   outputs are diffusion and diffusional kurtosis tensors and tensor-derived maps

% Author: Ali Tabesh
% Version: 2.6.0
% Last modified: 11/5/2014 by EM

% -------------------------------------------
% Constants
% -------------------------------------------

% if set to 1, dke will output tensors, e'vecs, and e'vals, 
% as well as d2/k2, d3/k3, color fa, and white matter model-derived maps
internal_flag = 1;

% -------------------------------------------
% Check inputs
% -------------------------------------------

[dkeVersion, dkeDate] = GetDKEVersion;
fprintf('%s, %s\n', dkeVersion, dkeDate)


if nargin ~= 1
    fprintf('\n')
    fprintf('Usage: dke paramsfile\n')
    fprintf('paramsfile  DKE processing parameters file\n\n')
    return
end

if ~exist(fn_params, 'file')
    fprintf('\n')
    fprintf('Input parameters file %s does not exist!\n\n', fn_params)
    return
end

% -------------------------------------------
% Read parameters file
% -------------------------------------------

fid=fopen(fn_params); %EM
file = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', ''); %EM
for i = 1:length(file{1}) %EM
    eval(file{1}{i})%EM
end
fclose(fid);

if exist('subject_list','var') %EM
else
subject_list={''};
end

% -------------------------------------------
% Output file names
% -------------------------------------------

fn_out_struc.dtype_out = 'float32';                 % output data format

% DKI outputs
fn_out_struc.kmean = 'kmean';                       % Mean kurtosis
fn_out_struc.k2    = 'k2';                          % Kurtosis along direction of medium diffusion
fn_out_struc.k3    = 'k3';                          % Kurtosis along direction of minimum diffusion
fn_out_struc.kax   = 'kax';                         % Axial kurtosis
fn_out_struc.krad  = 'krad';                        % Radial kurtosis

fn_out_struc.dmean        = 'dmean';                % Mean diffusivity
fn_out_struc.d2           = 'd2';                   % Diffusivity along direction of medium diffusion
fn_out_struc.d3           = 'd3';                   % diffusivity along direction of minimum diffusion
fn_out_struc.dax          = 'dax';                  % axial diffusivity
fn_out_struc.drad         = 'drad';                 % radial diffusivity
fn_out_struc.fa           = 'fa';                   % fractional anisotropy (FA)
fn_out_struc.fa_color_map = 'fa_color_map';         % FA color map

fn_out_struc.nexcl           = 'noutlier';          % number of outliers with the robust method (with method.robust_option ~= 0)
fn_out_struc.fit_err         = 'fit_err';           % fraction of unexplained signal variance by DKI model
fn_out_struc.fit_err_avg     = 'fit_err_avg';       % fraction of unexplained signal variance by DKI model
fn_out_struc.abs_fit_err_map = 'abs_fit_err_map';   % fraction of unexplained signal variance by DKI model
fn_out_struc.fit_err_map     = 'fit_err_map';       % fraction of unexplained signal variance by DKI model
fn_out_struc.abs_fit_err_avg = 'abs_fit_err_avg';   % fraction of unexplained signal variance by DKI model

fn_out_struc.d_viol    = 'd_viol';                  % fraction of constraint violations on directional diffusivity
fn_out_struc.kmin_viol = 'kmin_viol';               % fraction of constraint violations on minimum directional kurtosis
fn_out_struc.kmax_viol = 'kmax_viol';               % fraction of constraint violations on maximum directional kurtosis

fn_out_struc.kt  = 'KT.mat';                        % kurtosis tensor
fn_out_struc.dt  = 'DT.mat';                        % DKI diffusion tensor
fn_out_struc.kfa = 'kfa';                           % Kurtosis FA
fn_out_struc.meankt='mkt';                          % Mean Kurtosis Tensor 

fn_out_struc.dt_eigval = 'dt_eigval';               % DKI diffusion tensor eigenvalues
fn_out_struc.dt_eigvec = 'dt_eigvec';               % DKI diffusion tensor eigenvectors

% DTI outputs
fn_out_struc.dmean_dti = 'dmean_dti';               % mean diffusivity from DTI computation
fn_out_struc.d2_dti    = 'd2_dti';                  % diffusivity along direction of medium diffusion from DTI computation
fn_out_struc.d3_dti    = 'd3_dti';                  % diffusivity along direction of minimum diffusion from DTI computation
fn_out_struc.dax_dti   = 'dax_dti';                 % axial diffusivity from DTI computation
fn_out_struc.drad_dti  = 'drad_dti';                % radial diffusivity from DTI computation
fn_out_struc.fa_dti    = 'fa_dti';                  % FA from DTI computation

fn_out_struc.fa_color_map_dti    = 'fa_color_map_dti';      % FA color map
fn_out_struc.nexcl_dti           = 'noutlier_dti';          % number of outliers with the robust method (with method.robust_option ~= 0) from DTI computation
fn_out_struc.fit_err_dti         = 'fit_err_dti';           % fraction of unexplained signal variance by DTI model
fn_out_struc.fit_err_avg_dti     = 'fit_err_avg_dti';       % fraction of unexplained signal variance by DKI model
fn_out_struc.abs_fit_err_map_dti = 'abs_fit_err_map_dti';   % fraction of unexplained signal variance by DKI model
fn_out_struc.fit_err_map_dti     = 'fit_err_map_dti';       % fraction of unexplained signal variance by DKI model
fn_out_struc.abs_fit_err_avg_dti = 'abs_fit_err_avg_dti';   % fraction of unexplained signal variance by DKI model

fn_out_struc.dt_dti = 'dt_dti';                     % DTI diffusion tensor
fn_out_struc.dt_dti_eigval = 'dt_dti_eigval';       % DTI diffusion tensor eigenvalues
fn_out_struc.dt_dti_eigvec = 'dt_dti_eigvec';       % DTI diffusion tensor eigenvectors

% WMM outputs

fn_out_struc.awf = 'wmm_awf';                       % WMM axonal water fraction
fn_out_struc.da = 'wmm_da';                         % WMM axonal diffusivity
fn_out_struc.de_axial = 'wmm_de_ax';                % WMM axial extra-axonal diffusivity
fn_out_struc.de_radial = 'wmm_de_rad';              % WMM radial extra-axonal diffusivity
fn_out_struc.tortuosity = 'wmm_tort';               % WMM tortuosity

% -------------------------------------------
% Turn all warnings off
% -------------------------------------------

 warning('off','all');

% -------------------------------------------
% Process all subjects
% -------------------------------------------

for isubject = 1:length(subject_list)

    dir_subj = fullfile(studydir, subject_list{isubject});      % subject root folder

    diary off
    fn_diary = fullfile(dir_subj, 'dke.log');
    if exist(fn_diary, 'file')
        delete(fn_diary);
    end
    fid = fopen(fn_diary, 'w');
    if fid < 0
        error('Cannot open output file %s! Output directory does not exist or is write-protected.', fn_diary);
    end
    diary(fn_diary)

    fprintf('Start date and time: %s\n', datestr(now, 'mmmm dd, yyyy HH:MM:SS'))
    fprintf('%s\n',dkeVersion)% EM
    fnames = fieldnames(fn_out_struc);
    if isempty(subject_list{isubject})
        for ifield = 1:length(fnames)
            eval(['fn_out_struc_subject.' fnames{ifield} ' = fullfile(dir_subj, fn_out_struc.' fnames{ifield} ');']);
        end
    else
        for ifield = 1:length(fnames)
            eval(['fn_out_struc_subject.' fnames{ifield} ' = fullfile(dir_subj, [subject_list{isubject} ''_'' fn_out_struc.' fnames{ifield} ']);']);
        end
    end
    fn_out_struc_subject.dtype_out = fn_out_struc.dtype_out;
    
    if ~isempty(fn_noise)
        fn_subject_noise = fullfile(dir_subj, fn_noise);
    else
        fn_subject_noise = '';
    end
    
    if strcmpi(preprocess_options.format, 'dicom')
        eval(['!dke_preprocess_dicom "' dir_subj '" "' fn_params '"']);
        fn_subject_img_prefix = fullfile(dir_subj, 'intermediate_processing', 'combined', 'rdki');
        
    elseif strcmpi(preprocess_options.format, 'nifti')
        fn_subject_img_prefix = fullfile(dir_subj, preprocess_options.fn_nii);
    elseif strcmpi(preprocess_options.format, 'bruker')
        dke_preprocess_bruker(dir_subj, bval, preprocess_options);
        fn_subject_img_prefix = fullfile(dir_subj, 'intermediate_processing', 'combined', 'rdki');
    else
        error('Invalid input image format! Supported formats are DICOM and NIfTI.')
    end
    
    dke_estimate(fn_subject_img_prefix, idx_1st_img, bval, ndir, Kmin, NKmax, Kmin_final, Kmax_final, ...
        T, find_brain_mask_flag, dki_method, dti_method, fwhm_img, fn_subject_noise, fwhm_noise, ...
        fn_gradients, idx_gradients, fn_out_struc_subject, median_filter_method, map_interpolation_method, internal_flag)

end

diary off
