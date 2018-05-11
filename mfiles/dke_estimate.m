
function dke_estimate(fn_img_prefix, idx_1st_img, bval, ndir_orig, Kmin, NKmax, Kmin_final, Kmax_final, T, ...
    find_brain_mask_flag, dki_method, dti_method, fwhm_img, fn_noise, fwhm_noise, fn_gradients, ...
    idx_gradients, fn_out_struc, median_filter_param, map_interpolation_method, internal_flag)

% dke_estimate    Diffusional kurtosis estimator; estimates the diffusion and kurtosis
%       tensors and tensor-derived scalar measures of diffusivity and
%       diffusional kurtosis
%
% Syntax
%
%   dke(fn_img_prefix, idx_1st_img, bval, ndir_orig, Kmin, NKmax,
%       Kmin_final, Kmax_final, T, find_brain_mask_flag, dki_method,
%       dti_method, fn_noise, fwhm_img, fwhm_noise, fn_gradients,
%       idx_gradients, fn_out_struc, debug_flag)
%
% Inputs
%
%   fn_img_prefix       Prefix for diffuion-weighted image (DWI) file names (See read_img_set for more information)
%
%   idx_1st_img         Index of first DWI (typically 0 or 1)
%
%   bval                1-by-nbval vector of b values (in s/mm^2 units), where nbval is the number of b values,
%                       e.g., [0 1000 2000]
%   
%   ndir                1-by-(nbval-1) vector containing the number of gradient directions per b value, where nbval is the number 
%                       of b values; if ndir is a scalar, the number of directions will be set to that number for all b values
%
%   Kmin                Minimum acceptable directional kurtosis
%
%   NKmax               Parameter defining the maximum directional kurtosis as 
%                       Kmax = NKmax / (Di * bmax)
%                       where Di is directional diffusivity along direction i for the given voxel and bmax is the largest b value 
%
%   Kmin_final          Minimum allowable kurtosis value (mean, axial, radial, k1, k2)
%
%   Kmax_final          Maximum allowable kurtosis value (mean, axial, radial, k1, k2)
%
%   T                   Threshold for finding the brain mask
%
%   brain_mask_flag     Use connected component analysis to find brain mask
%
%   dki_method          Structure containing information for tensor fitting method
%                       linear_weighting    Whether to use unweighted (0) or weighted (1) linear least-squares (default: 0)
%                       linear_constrained  Whether to use unconstrained (0) or constrained (1) algorithm (default: 1)
%                       nonlinear           Whether (1) or not (0) to use nonlinear least-squares (only unconstrained 
%                                           solution) (default: 0)
%                       linear_violations   Whether (1) or not (0) to generate maps of constraint violations where 
%                                           unconstrained linear least squares solution violates directional diffusivities 
%                                           and kurtoses (default: 0)
%                       robust_option       Robust fitting option: (0) do not use robust fitting; (1) RESTORE-type algorithm 
%                                           (outlier detection and removal followed by tensor refitting) with a user-supplied 
%                                           'noise tolerance level' expressed as a fraction of the diffusion signal level; 
%                                           (2) RESTORE-type algorithm with a user-supplied 'significance level' used to 
%                                           determine the threshold for outlier detection based on fit residuals
%                       noise_tolerance     Noise tolerance level as fraction of signal level used for outlier detection 
%                                           (used with dki_method.robust_option = 1)
%                       significance_level  Significance level used to determine the threshold for outlier detection 
%                                           (used with dki_method.robust_option = 2)
%
%   dti_method          dti_flag            Whether or not to perform DTI computations
%                       dti_only            Whether or not to only perform DTI computations (no DKI)
%                       no_tensor           Estimate mean diffusivity using directional signal fits and 
%                                           no tensor estimation (default: 0)
%                       b_value             b-value for DTI computations
%                       directions          Indices of gradient directions for DTI computations these indices are relative 
%                                           to the indices specified in idx_gradients for dti_method.b_value
%
%   fn_noise            Noise map file name (in NIfTI format)
%
%   fwhm_img            FWHM (in mm units) for Gaussian smoothing applied to DWI's
%                       fwhm_img = 0 means no filtering
%
%   fwhm_noise          FWHM (in mm units) for Gaussian smoothing applied to the noise image
%                       fwhm_img = 0 means no filtering
%
%   fn_gradients        File name for gradient vectors; gradient vectors g
%                       will be normalized to ||g|| = 1
%
%   idx_gradients       1-by-(nbval-1) cell array containing indices of gradient directions used for estimation
%
%   fn_out_struc        Structure containing file names for output images and data files
%                       kmean           Mean kurtosis
%                       k2              Kurtosis along the eigenvector for 2nd largest DT eigenvalue
%                       k3              Kurtosis along the eigenvector for smallest DT eigenvalue
%                       kax             Axial kurtosis
%                       krad            Radial kurtosis
%                       dmean           Mean diffusivity as average of DT eigenvalues
%                       d2              Second largest DT eigenvalue
%                       d3              Smallest DT eigenvalue
%                       dax             Axial diffusivity
%                       drad            Radial diffusivity
%                       fa              Fractional anisotropy
%                       L               Eigenvalues of DT
%                       dtype_out       Data type for output files, e.g., 'float32'
%
%   debug_flag          If debug_flag is set, DTs, DT eigenvectors, DT eigenvalues, and kurtosis tensors are written to files 
%                       'diffusion_tensor.mat', 'diffusion_tensor_eigenvectors.mat', 'diffusion_tensor_eigenvalues.mat', 
%                       'kurtosis_tensor.mat', respectively. Also, the smoothed DWI's are saved with file names prefixed with 'sm_'.
%
% Outputs
%
%   Images and data file specified in fn_out_struc will be created. Note that the diffusivity images are saved in um^2/ms units. 

% Author: Ali Tabesh
% Version: 2.6.0
% Last modified: 11/25/14 by EM

%--------------------------------------------------------------------------
% Initialize parameters
%--------------------------------------------------------------------------

% Tolerance levels
err     = 1e-15;        % Tolerance in distinguishing pairs of eigenvalues of diffusion tensor (DT) (mm^2/s units)
errtol  = 1e-15;        % Tolerance in estimating R_D and R_F
Leps    = 1e-15;        % Smallest acceptable diffusion tensor eigenvalue (mm^2/s units)
Seps    = 1e-6;         % Smallest acceptable signal value for nonzero b-values

% parameters for white matter model
dki_method.ndir_Kmax = 10000;
dki_method.wm_model_flag = 1;
if internal_flag == 0
    dki_method.wm_model_flag = 0;
end

%--------------------------------------------------------------------------
% Check inputs
%--------------------------------------------------------------------------

% check if idx_gradient has the same number of elements as there are
% non-zero b-values

if length(idx_gradients) ~= length(bval) - 1
    error('Number of elements of idx_gradients must match the number of nonzero b-values in bval!')
end

% check if dti_method.directions has the same number of elements as there are
% non-zero b-values

if length(dti_method.directions) ~= length(dti_method.b_value)
    error('Number of elements of dti_method.directions must match the number of nonzero b-values in dti_method.b_value!')
end

% convert ndir_orig to a vector
if isscalar(ndir_orig)
    ndir_orig = ndir_orig * ones(length(bval) - 1, 1);
end

%--------------------------------------------------------------------------
% Read input files
%--------------------------------------------------------------------------

% read diffusion-weighted images
fprintf('Reading input images... ')

if exist(fn_img_prefix, 'file')  % assume 4d with complete filename
    [img_set hdr img_size voxel_size] = read_img_set_nii(fn_img_prefix);
elseif exist([fn_img_prefix '.nii'], 'file')  % assume 4d with missing '.nii'
    [img_set hdr img_size voxel_size] = read_img_set_nii([fn_img_prefix '.nii']);
elseif exist([fn_img_prefix '.nii.gz'], 'file')  % assume 4d with missing '.nii.gz'
    [img_set hdr img_size voxel_size] = read_img_set_nii([fn_img_prefix '.nii.gz']);
elseif exist([fn_img_prefix '_0.nii'], 'file') % assume prefix is for 3d series
    [img_set hdr img_size voxel_size] = read_img_set_nii(fn_img_prefix, ndir_orig, bval(2:end), idx_1st_img);
else
    error('Input NIfTI image %s does not exist!', fn_img_prefix)
end

fprintf('%s\n','complete')

% set ndir and ndir_total to correct values for the gradient subset 
ndir = zeros(1, length(idx_gradients));
for i = 1:length(idx_gradients)
    ndir(i) = length(idx_gradients{i});                % number of subsampled gradient directions
end

ndir_total = sum(ndir);

% take a subset of diffusion-weighted images
img_set2 = zeros([size(img_set, 1), size(img_set, 2), size(img_set, 3), 1 + ndir_total]);

img_set2(:, :, :, 1) = img_set(:, :, :, 1);
for i = 1:ndir(1)
    img_set2(:, :, :, i + 1) = img_set(:, :, :, 1 + idx_gradients{1}(i));
end

for ib = 2:(length(bval) - 1)
    for i = 1:ndir(ib)
        img_set2(:, :, :, i + sum(ndir(1:(ib - 1))) + 1) = img_set(:, :, :, 1 + sum(ndir_orig(1:(ib - 1))) + idx_gradients{ib}(i));
    end
end

img_set = img_set2;
clear img_set2;

% read and normalize gradient directions
g = read_gradients(fn_gradients, ndir, bval(2:end), idx_gradients);

% generate gradient set for white matter model

if isfield(dki_method, 'wm_model_flag')
    if dki_method.wm_model_flag == 1
        [dki_method.wm_model_DG, dki_method.wm_model_KG] = generate_wm_DG_KG(dki_method.ndir_Kmax);
    elseif dki_method.wm_model_flag ~= 0
        error('Invalid dki_method.wm_model_flag option! Parameter dki_method.wm_model_flag must be either 0 or 1.')
    end
else
    dki_method.wm_model_flag = 0;
end


%--------------------------------------------------------------------------
% Process images
%--------------------------------------------------------------------------

% find brain mask

if find_brain_mask_flag == 0                % find voxels with b = 0 signal > T
    mask = squeeze(img_set(:, :, :, 1)) > T;
elseif find_brain_mask_flag == 1
    mask = find_brain_mask(squeeze(img_set(:, :, :, 1)), T);  % find voxels with b = 0 signal > T and then refine results by connected component analysis
else
    error('Invalid find_brain_mask_flag option! Parameter find_brain_mask_flag can only be 0 or 1!')
end
mask = mask(:);

% Smooth gradient and noise images; remove noise

fprintf('Filtering input images... ')

x = filter_img_set(img_set, img_size, fn_noise, fwhm_img, fwhm_noise, voxel_size);

% if debug_flag   % if debug_flag is set, write smoothed images z{:} to file
%     [pathstr name ext] = fileparts(fn_img_prefix);
%     fn = fullfile(pathstr, ['sm_' name ext '_0']);
%     write_img(fn, z{1}, hdr, fn_out_struc.dtype_out)
%     for ibval = 2:length(bval)
%         for idir = 1:ndir(ibval-1)
%             fn = fullfile(pathstr, ['sm_' name ext '_' num2str(bval(ibval)) '_' num2str(idir-1)]);
%             write_img(fn, z{idir + sum(ndir(1:(ibval-2))) + 1}, hdr, fn_out_struc.dtype_out)
%         end
%     end
% end
% clear z

fprintf('complete.\n')

% form the coefficients matrix and its pseudoinverse, as well as the matrix
% needed for estimating the covariance matrix of estimated coefficients

[DG, KG, DKG] = coeff_mat(g, ndir, bval(2:end));
DKGpinv = pinv(DKG);
% if dki_method.robust_option == 2
%     DKGprodinv = inv(DKG' * DKG);
% else
DKGprodinv = [];
% end

if dki_method.wm_model_flag == 1
    DGpinv = pinv(DG);
end

% DKI coefficients matrix when no tensors are estimated
if dki_method.no_tensor == 1
    b = bval(2:end);
    A = [-b(:), (b(:) .^ 2) / 6];
    Apinv = pinv(A);
elseif dki_method.no_tensor ~= 0
    error('Invalid dki_method.no_tensor option! Parameter dki_method.no_tensor must be either 0 or 1.')
end

% form the coefficients matrix for DTI processing
if dti_method.dti_flag == 1

    DGdti = [];
    dti_method.DG_indices = [];
    
    for i = 1:length(dti_method.b_value)
        idx = find(bval(2:end) == dti_method.b_value(i));
        if isempty(idx)
            error('Cannot find b-value(s) specified for DTI computations in the vector of b-values bval!')
        end
        if length(dti_method.directions{i}) > ndir(idx)
            error('Number of directions in dti_method.directions cannot be larger than number of directions (ndir) for the b-value!')
        end
        if max(dti_method.directions{i}) > ndir(idx)
            error('Indices in dti_method.directions cannot be larger than number of directions (ndir) for the b-value!')
        end
        sdir = sum(ndir(1:(idx - 1)));     % starting row index for the b-value used for dte_core
        dti_method.DG_indices = [dti_method.DG_indices; sdir + dti_method.directions{i}(:)];     % indices of relevant rows of DG
        DGdti = [DGdti; -dti_method.b_value(i) * DG(sdir + dti_method.directions{i}, :)];
    end        

	% coefficient matrices
    if dti_method.no_tensor == 1
        b = dti_method.b_value;
        Adti = -b(:);
        Adtipinv = pinv(Adti);
    elseif dti_method.no_tensor == 0
        DGdtipinv = pinv(DGdti);
        DGdtiprodinv = inv(DGdti' * DGdti);
    else
        error('Invalid dti_method.no_tensor option! Parameter dti_method.no_tensor must be either 0 or 1.')
    end
end         % dti_method.dti_flag == 1

% inequality constraints matrix
C = [-DG, zeros(ndir_total, 15); zeros(ndir_total, 6), -KG; -DG * NKmax / max(bval) * bval(2), KG];

% equality constraints matrix and vector (minimum-norm constraint)
C_eq = null(DKG)';
c_eq = zeros(size(C_eq, 1), 1);
 
% if debug_flag is set, display condition numbers
 % if debug_flag
%     disp(['Condition number for diffusion coeff matrix: ' num2str(cond(DG))])
%     disp(['Condition number for kurtosis coeff matrix: ' num2str(cond(KG))])
%     disp(['Condition number for diffusion/kurtosis coeff matrix: ' num2str(cond(DKG))])
% end

% initialize output diffusion and kurtosis variables

nvoxel = size(x, 2);

if dti_method.dti_only == 0
    Dmean     = zeros(1, nvoxel);
    Dax       = zeros(1, nvoxel);
    Drad      = zeros(1, nvoxel);
    FA        = zeros(1, nvoxel);
    Kmean     = zeros(1, nvoxel);
    K2        = zeros(1, nvoxel);
    K3        = zeros(1, nvoxel);
    Kax       = zeros(1, nvoxel);
    Krad      = zeros(1, nvoxel);
    KFA       = zeros(1, nvoxel);
    DT        = zeros(6, nvoxel);
    L         = zeros(3, nvoxel);
    V         = zeros(9, nvoxel);
    KT        = zeros(15, nvoxel);
    nexcl     = zeros(1, nvoxel);
    D_viol    = zeros(1, nvoxel);
    Kmin_viol = zeros(1, nvoxel);
    Kmax_viol = zeros(1, nvoxel);
    fit_err   = zeros(1, nvoxel);
    MeanKT    = zeros(1, nvoxel);
    
    if dki_method.wm_model_flag == 1
        AWF        = zeros(1, nvoxel);
        Da         = zeros(1, nvoxel);
        De_axial   = zeros(1, nvoxel);
        De_radial  = zeros(1, nvoxel);
        tortuosity = zeros(1, nvoxel);
    end
    
end

if dti_method.dti_flag == 1
    Dax_dti     = zeros(1, nvoxel);
    Drad_dti    = zeros(1, nvoxel);
    Dmean_dti   = zeros(1, nvoxel);
    FA_dti      = zeros(1, nvoxel);
    DT_dti      = zeros(6, nvoxel);
    L_dti       = zeros(3, nvoxel);
    V_dti       = zeros(9, nvoxel);
    nexcl_dti   = zeros(1, nvoxel);
    fit_err_dti = zeros(1, nvoxel);
end

% variables for displaying progress 

I = (0:100) * nvoxel / 100;
J = I(1);
k = 1;

% set options for nlinfit (regular nonlinear least-squares)

if dki_method.nonlinear == 1
    nlinfit_options = statset('FunValCheck', 'off');
    % turn off warning messages
%    warning('off', 'stats:nlinfit:IterationLimitExceeded');
%    warning('off', 'stats:nlinfit:IllConditionedJacobian');
else
    nlinfit_options = [];
end

% set options for robust fitting

% if dki_method.robust_option ~= 0 || dti_method.robust_option ~= 0
%     warning('off', 'stats:statrobustfit:IterationLimit');
% end

% lower and upper limits from Carlson's codes for rd and rf
constants_struc.rd_lolim = 2 / (realmax('double') ^ (2 / 3));
constants_struc.rd_uplim = (0.1 * errtol / realmin('double')) ^ (2 / 3);
constants_struc.rf_lolim = realmin('double') * 5;
constants_struc.rf_uplim = realmax('double') / 5;

% loop over voxels
fprintf('Processing voxels... ')

if feature('numCores') > 1
    
    % Start a parallel pool, unless one has already been started
    % Before MATLAB 2013b (MATLAB version 8.2), 'matlabpool open' starts a parallel pool
    % For MATLAB 2013b and later, 'gcp' creates a parallel pool if one currently does not exist
    if verLessThan('matlab', '8.2')
        if matlabpool('size') <= 0
            evalc('matlabpool open');
        end
    else
        evalc('gcp');
    end
    
    if ~dti_method.dti_only
        if dki_method.no_tensor
            ndir_no_tensor = ndir(1);
            bval_no_tensor = bval(2:end);
        end
    end
    
    if dti_method.dti_flag
        xdti = x([1; dti_method.DG_indices(:) + 1], :);
        ndir_dti = length(dti_method.directions{1});
        b_value_dti = dti_method.b_value;
    end
    
    nbatch = 100;
    step = ceil(nvoxel / nbatch);
    
    for j = 1:nbatch
        
        start = (j - 1) * step + 1;
        stop = j * step;
        
        if ~dti_method.dti_only
            if dki_method.no_tensor
               parfor i = start:min(nvoxel, stop)
                    if mask(i) ~= 0
                        [Dax(i), Drad(i), FA(i), Dmean(i), Kmean(i), K2(i), K3(i), Kax(i), Krad(i), ...
                            DT(:, i), L(:, i), V(:, i), KT(:, i), nexcl(i), D_viol(i), Kmin_viol(i), Kmax_viol(i), fit_err(i)] = ...
                            dke_core_dirfit(x(:, i), ndir_no_tensor, bval_no_tensor, A, Apinv, Seps, Kmin, NKmax, Kmin_final, Kmax_final, dki_method);
                        [KFA(i),MeanKT(i)]= ComputeKFA(KT(:, i), x(:,i),Kmax_final,Kmin_final); %EM
                    end
                end
            else
                parfor i = start:min(nvoxel, stop)
                    if mask(i) ~= 0
                         warning('off','all');
                        [Dax(i), Drad(i), FA(i), Dmean(i), Kmean(i), K2(i), K3(i), Kax(i), Krad(i), ...
                            DT(:, i), L(:, i), V(:, i), KT(:, i), nexcl(i), D_viol(i), Kmin_viol(i), Kmax_viol(i), fit_err(i)] = ...
                                dke_core(x(:, i), ndir_total, ndir, bval, Kmin, NKmax, Kmin_final, Kmax_final, g, C, C_eq, c_eq, ...
                                DG, DKG, DKGpinv, DKGprodinv, dki_method, err, errtol, Leps, Seps, nlinfit_options, constants_struc);
                      [KFA(i),MeanKT(i)]= ComputeKFA(KT(:, i), x(:,i),Kmax_final,Kmin_final); %EM
                    end
                end
                if dki_method.wm_model_flag
                    wm_model_DG = dki_method.wm_model_DG;       % attempt to eliminate warning message
                    wm_model_KG = dki_method.wm_model_KG;       % attempt to eliminate warning message
                    parfor i = start:min(nvoxel, stop)
                        if mask(i) ~= 0
%                             [AWF(i), Da(i), De_axial(i), De_radial(i), tortuosity(i)] = ...
%                                 wm_model_estimate(DT(:, i), KT(:, i), Dax(i), Drad(i), Dmean(i), DG, DGpinv, KG, dki_method.wm_model_DG, dki_method.wm_model_KG);
                            [AWF(i), Da(i), De_axial(i), De_radial(i), tortuosity(i)] = ...
                                wm_model_estimate(DT(:, i), KT(:, i), Dax(i), Drad(i), Dmean(i), DG, DGpinv, KG, wm_model_DG, wm_model_KG);

                        end
                    end
                end
            end        
        end
        
        if dti_method.dti_flag
            if dti_method.no_tensor
                parfor i = start:min(nvoxel, stop)
                    if mask(i) ~= 0
                        [Dax_dti(i), Drad_dti(i), Dmean_dti(i), FA_dti(i), DT_dti(:, i), L_dti(:, i), V_dti(:, i), nexcl_dti(i), fit_err_dti(i)] = ...
                            dte_core_dirfit(xdti(:, i), ndir_dti, b_value_dti, Adti, Adtipinv, Seps, dti_method);
                    end
                end
            else
                parfor i = start:min(nvoxel, stop)
                    if mask(i) ~= 0
                        [Dax_dti(i), Drad_dti(i), Dmean_dti(i), FA_dti(i), DT_dti(:, i), L_dti(:, i), V_dti(:, i), nexcl_dti(i), fit_err_dti(i)] = ...
                            dte_core(xdti(:, i), DGdti, DGdtipinv, DGdtiprodinv, Leps, Seps, dti_method);
                    end
                end
            end
        end
        
        if j > 10
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b')
        elseif j > 1
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b')
        end
        fprintf('%d%% complete.', j * 100 / nbatch)
        
    end
    
    % Shut down the parallel pool
    if verLessThan('matlab', '8.2')
        evalc('matlabpool close');
    else
        poolobj = gcp('nocreate');
        evalc('delete(poolobj)');
    end

else
    
    for i = 1:nvoxel
        if mask(i) ~= 0
            if ~dti_method.dti_only
                if dki_method.no_tensor
                    [Dax(i), Drad(i), FA(i), Dmean(i), Kmean(i), K2(i), K3(i), Kax(i), Krad(i), ...
                        DT(:, i), L(:, i), V(:, i), KT(:, i), nexcl(i), D_viol(i), Kmin_viol(i), Kmax_viol(i), fit_err(i)] = ...
                        dke_core_dirfit(x(:, i), ndir(1), bval(2:end), A, Apinv, Seps, Kmin, NKmax, Kmin_final, Kmax_final, dki_method);
                else
                    [Dax(i), Drad(i), FA(i), Dmean(i), Kmean(i), K2(i), K3(i), Kax(i), Krad(i), ...
                        DT(:, i), L(:, i), V(:, i), KT(:, i), nexcl(i), D_viol(i), Kmin_viol(i), Kmax_viol(i), fit_err(i)] = ...
                        dke_core(x(:, i), ndir_total, ndir, bval, Kmin, NKmax, Kmin_final, Kmax_final, g, C, ...
                            C_eq, c_eq, DG, DKG, DKGpinv, DKGprodinv, dki_method, err, errtol, Leps, Seps, nlinfit_options, constants_struc);
                    if dki_method.wm_model_flag
                        [AWF(i), Da(i), De_axial(i), De_radial(i), tortuosity(i)] = ...
                            wm_model_estimate(DT(:, i), KT(:, i), Dax(i), Drad(i), Dmean(i), DG, DGpinv, KG, dki_method.wm_model_DG, dki_method.wm_model_KG);
                    end
                end
                [KFA(i),MeanKT(i)] = ComputeKFA(KT(:, i), x(:,i),Kmax_final,Kmin_final); %EM
            end
            
            if dti_method.dti_flag
                if dti_method.no_tensor
                    [Dax_dti(i), Drad_dti(i), Dmean_dti(i), FA_dti(i), DT_dti(:, i), L_dti(:, i), V_dti(:, i), nexcl_dti(i), fit_err_dti(i)] = ...
                        dte_core_dirfit(x([1; dti_method.DG_indices(:) + 1], i), length(dti_method.directions{1}), dti_method.b_value, Adti, Adtipinv, Seps, dti_method);
                else
                    [Dax_dti(i), Drad_dti(i), Dmean_dti(i), FA_dti(i), DT_dti(:, i), L_dti(:, i), V_dti(:, i), nexcl_dti(i), fit_err_dti(i)] = ...
                        dte_core(x([1; dti_method.DG_indices(:) + 1], i), DGdti, DGdtipinv, DGdtiprodinv, Leps, Seps, dti_method);
                end
            end
        end
        
        if i > J
            if k > 10
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b')
            elseif k > 1
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b')
            end
            fprintf('%d%% complete.', k)
            k = k + 1;
            J = I(k);
        end
    end
end

% convert diffusivities to units of um^2/ms

if ~dti_method.dti_only
    Dax   = 1000 * Dax;
    Drad  = 1000 * Drad;
    Dmean = 1000 * Dmean;
    DT    = 1000 * DT;
    L     = 1000 * L;
    if dki_method.wm_model_flag == 1
        Da        = 1000 * Da;
        De_axial  = 1000 * De_axial;
        De_radial = 1000 * De_radial;
    end  
end

if dti_method.dti_flag
    Dax_dti   = 1000 * Dax_dti;
    Drad_dti  = 1000 * Drad_dti;
    Dmean_dti = 1000 * Dmean_dti;
    DT_dti    = 1000 * DT_dti;
    L_dti     = 1000 * L_dti;
end

% median filter structure (ignored when generating diagnostic maps and dti maps)
median_filter_method.type = median_filter_param;


if ~dti_method.dti_only
    median_filter_method.img_ref  = reshape(Kmean, hdr.dim(2:4));
    median_filter_method.img_viol = reshape(Kmin_viol, hdr.dim(2:4));
    if median_filter_method.type == 1
        median_filter_method.T = 1 - (length(bval) - 1) * 15 / ndir_total;    % at least 15 good noncollinear directions
    elseif median_filter_method.type == 2
        median_filter_method.T = 0;
    end
end

% no median filtering for diagnostic maps (constraint violations, number of outliers in robust fitting)
median_filter_method_diagnostic_maps.type = 0;

% no median filtering for dti maps
median_filter_method_dti.type = 0;

% write images to file
fprintf('\n')
fprintf('\nWriting output images to files...\n')

if dti_method.dti_only == 0
    
    write_img(fn_out_struc.kmean, Kmean, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
    write_img(fn_out_struc.dmean, Dmean, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
    
    if dki_method.no_tensor == 0
        write_img(fn_out_struc.kax, Kax, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
        write_img(fn_out_struc.krad, Krad, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
        write_img(fn_out_struc.dax, Dax, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
        write_img(fn_out_struc.drad, Drad, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
        write_img(fn_out_struc.fa, FA, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
        write_img(fn_out_struc.kfa, KFA, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
        write_img(fn_out_struc.meankt, MeanKT, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
          
        %GRG------------------------------------------------------------------------------------------------------------
        write_img(fn_out_struc.kt, KT, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 2)
        write_img(fn_out_struc.dt, DT, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 2)
        %---------------------------------------------------------------------------------------------------------------
            
        if internal_flag == 1
            write_img(fn_out_struc.k2, K2, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
            write_img(fn_out_struc.k3, K3, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
            write_img(fn_out_struc.d2, L(2, :), hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
            write_img(fn_out_struc.d3, L(1, :), hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
            write_img(fn_out_struc.dt_eigval, L, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 1)
            write_img(fn_out_struc.dt_eigvec, V, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 1)
            make_fa_color_map(FA, V, img_size, fn_out_struc.fa_color_map)
            
        end

        if dki_method.linear_violations == 1
            write_img(fn_out_struc.d_viol, D_viol, hdr, fn_out_struc.dtype_out, median_filter_method_diagnostic_maps, map_interpolation_method, 0)
            write_img(fn_out_struc.kmin_viol, Kmin_viol, hdr, fn_out_struc.dtype_out, median_filter_method_diagnostic_maps, map_interpolation_method, 0)
            write_img(fn_out_struc.kmax_viol, Kmax_viol, hdr, fn_out_struc.dtype_out, median_filter_method_diagnostic_maps, map_interpolation_method, 0)
            write_img(fn_out_struc.fit_err, fit_err, hdr, fn_out_struc.dtype_out, median_filter_method_diagnostic_maps, map_interpolation_method, 1)
        end

        if dki_method.robust_option ~= 0
            write_img(fn_out_struc.nexcl, nexcl, hdr, fn_out_struc.dtype_out, median_filter_method_diagnostic_maps, map_interpolation_method, 0)
        end
        
        if dki_method.wm_model_flag == 1
            write_img(fn_out_struc.awf, AWF, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
            write_img(fn_out_struc.da, Da, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
            write_img(fn_out_struc.de_axial, De_axial, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
            write_img(fn_out_struc.de_radial, De_radial, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
            write_img(fn_out_struc.tortuosity, tortuosity, hdr, fn_out_struc.dtype_out, median_filter_method, map_interpolation_method, 0)
        end    
    end
end

if dti_method.dti_flag
    write_img(fn_out_struc.dmean_dti, Dmean_dti, hdr, fn_out_struc.dtype_out, median_filter_method_dti, map_interpolation_method, 0)

    if dti_method.no_tensor == 0
        write_img(fn_out_struc.dax_dti, Dax_dti, hdr, fn_out_struc.dtype_out, median_filter_method_dti, map_interpolation_method, 0)
        write_img(fn_out_struc.drad_dti, Drad_dti, hdr, fn_out_struc.dtype_out, median_filter_method_dti, map_interpolation_method, 0)

        if internal_flag == 1
            write_img(fn_out_struc.d2_dti, L_dti(2, :), hdr, fn_out_struc.dtype_out, median_filter_method_dti, map_interpolation_method, 0)
            write_img(fn_out_struc.d3_dti, L_dti(1, :), hdr, fn_out_struc.dtype_out, median_filter_method_dti, map_interpolation_method, 0)
            write_img(fn_out_struc.dt_dti, DT_dti, hdr, fn_out_struc.dtype_out, median_filter_method_dti, map_interpolation_method, 1)
            write_img(fn_out_struc.dt_dti_eigval, L_dti, hdr, fn_out_struc.dtype_out, median_filter_method_dti, map_interpolation_method, 1)
            write_img(fn_out_struc.dt_dti_eigvec, V_dti, hdr, fn_out_struc.dtype_out, median_filter_method_dti, map_interpolation_method, 1)
            write_img(fn_out_struc.fit_err_dti, fit_err_dti, hdr, fn_out_struc.dtype_out, median_filter_method_dti, map_interpolation_method, 1)
            make_fa_color_map(FA_dti, V_dti, img_size, fn_out_struc.fa_color_map_dti)
        end
        
        write_img(fn_out_struc.fa_dti, FA_dti, hdr, fn_out_struc.dtype_out, median_filter_method_dti, map_interpolation_method, 0)
        
        if dti_method.robust_option ~= 0
            write_img(fn_out_struc.nexcl_dti, nexcl_dti, hdr, fn_out_struc.dtype_out, median_filter_method_diagnostic_maps, map_interpolation_method, 0)
        end
    end
end

fprintf('Complete\n')

[pathstr, ~, ~] = fileparts(fn_out_struc.kmean);

fprintf('\n')
fprintf('Diffusional kurtosis maps are in folder %s\n', pathstr)
fprintf('DKE processing parameters file is %s\n', fullfile(pathstr, 'DKEParameters.dat'))
fprintf('DKE log file is %s\n\n', fullfile(pathstr, 'dke.log'))


%--------------------------------------------------------------------------
% read gradient vectors
%--------------------------------------------------------------------------

function g = read_gradients(fn_gradients, ndir, bval, idx_gradients)
g=cell(1,length(bval));
if iscell(fn_gradients)
    if length(fn_gradients) == 1
        fn_gradients = fn_gradients{1}; 
        g_orig = load(fn_gradients);%EM
        g_orig = g_orig ./ repmat(sqrt(sum(g_orig .^ 2, 2)), 1, 3);%EM
    
        for i = 1:length(bval)%EM
        g{i} = g_orig(idx_gradients{i}, :);  % take the subset of gradient set specified by idx_gradients
        end
        if size(g{1}, 1) ~= ndir(1)%EM
        error(['Number of gradient vectors in file ' fn_gradients ' does not match the user-specified number of directions!'])
        end
    
    else
        for i = 1:length(bval)
            g{i} = load(fn_gradients{i});
            g{i} = g{i}(idx_gradients{i}, :);  % take the subset of gradient set specified by idx_gradients
            if size(g{i}, 1) ~= ndir(i)
                error(['Number of gradient vectors in file ' fn_gradients{i} ' does not match the user-specified number of directions!'])
            end
            g{i} = g{i} ./ repmat(sqrt(sum(g{i} .^ 2, 2)), 1, 3);
        end
    end
else
    g_orig = load(fn_gradients);
    g_orig = g_orig ./ repmat(sqrt(sum(g_orig .^ 2, 2)), 1, 3);
    
    for i = 1:length(bval)
        g{i} = g_orig(idx_gradients{i}, :);  % take the subset of gradient set specified by idx_gradients
    end
    if size(g{1}, 1) ~= ndir(1)
        error(['Number of gradient vectors in file ' fn_gradients ' does not match the user-specified number of directions!'])
    end
end

g = cell2mat(g');



%--------------------------------------------------------------------------
% core function that estimates the tensors and tensor-derived measures
%--------------------------------------------------------------------------

 function [Dax, Drad, FA, Davg, Kmean, K2, K3, Kax, Krad, dt, L, V, kt, nexcl, D_viol, Kmin_viol, Kmax_viol, fit_err] = ...
    dke_core(x, n, ndir, b, Kmin, NKmax, Kmin_final, Kmax_final, g, C, C_eq, c_eq, DG, DKGu, DKGupinv, DKGprodinv, ...
    dki_method, err, errtol, Leps, Sbeps, nlinfit_options, const)

%dke_core   Core function to estimate diffusion and kurtosis parameters for a
%           single voxel
%
%Syntax
%
%   [Dax, Drad, FA, Davg, Kmean, K2, K3, Kax, Krad, dt, L, V, kt, nexcl, D_viol, Kmin_viol, Kmax_viol] = 
%   dke_core(x, n, ndir, b, Kmin, NKmax, Kmin_final, Kmax_final, g, C, C_eq, c_eq, DG, DKGu, DKGupinv, dki_method, err, errtol, Leps, Sbeps, nlinfit_options, const)
%
%Inputs
%
%   x           (ndir+1)-by-1 vector of signal values
%               The first row of x contains the b = 0 signal values, the next
%               ndir_1 rows contain the b = b_1 signal values, and so on, where
%               ndir_1 is the number of directions for b = b_1
%
%   ndir        Total number of gradient directions across b values
%
%   b           1-by-nbval vector of b values (in s/mm^2 units), where nbval is
%               the number of b values
%
%   Kmin        Minimum acceptable directional kurtosis
%
%   Kminf       Minimum allowable kurtosis value (mean, axial, radial)
%
%   Kmaxf       Maximum allowable kurtosis value (mean, axial, radial)
%
%   C           Matrix specifying the inequality constraints on directional diffusivities and kurtoses
%
%   C_eq        Matrix specifying the null space of DKG for imposing the minimum-norm constraint
%
%   c_eq        Vector of zeros for the minimum-norm constraint
%
%   DG          Matrix formed from gradient vectors to obtain diffusion tensor (DT)
%
%   DKG         Matrix formed from gradient vectors to obtain diffusion and kurtosis tensors
%
%   DKGpinv     Pseudoinverse of DKG
%
%   err         Tolerance in distinguishing pairs of eigenvalues of diffusion
%               tensor (DT) for identifying singularities (mm^2/s units) 
%
%   errtol      Tolerance in estimating R_D and R_F
%
%   Leps        Smallest acceptable DT eigenvalue (mm^2/s units)
%               Eigenvalues smaller than Leps are set to Leps.
%
%   Sbeps       Smallest acceptable signal value for nonzero b values
%               Signal values smaller than Sbeps are set
%               to Sbeps.
%
%   nlinfit_options options for unconstrained regular nonlinear least-squares
%
%   constants_struc structure containing lower and upper limits for inputs
%               to functions for computing rd and rf
%
%Outputs
%
%   Dax         Axial diffusivity (mm^2/s units)
%
%   Drad        Radial diffusivity (mm^2/s units)
%
%   FA          Fractional anisotropy
%
%   Davg        Mean diffusivity as average of DT eigenvalues (mm^2/s units)
%
%   Kmean       Mean kurtosis from analytical formula
%
%   K2          Kurtosis along the direction of the eigenvector corresponding
%               to the 2nd DT eigenvalue
%
%   K3          Kurtosis along the direction of minimum diffusivitiy
%
%   Kax         Axial kurtosis
%
%   Krad        Radial kurtosis
%
%   DT          6-by-1 vector containing elements of DT (mm^2/s units)
%               DT in matrix form denoted as DTfull can be obtained as 
%               DTfull = [DT(1) DT(4) DT(5); DT(4) DT(2) DT(6); DT(5) DT(6) DT(3)]
%
%   L           3-by-1 vector of DT eigenvalues (mm^2/s units) sorted in 
%               ascending order
%   
%   V           9-by-1 vector containing the eigenvectors of the DT
%               corresponding to eigenvalues in L
%
%   KT          15-by-1 vector containing elements of KT (mm^2/s units)
%               KT in array form denoted as KTfull can be obtained as 
%               KTfull(1,1,1,1) = KT(1);
%               KTfull(2,2,2,2) = KT(2);
%               KTfull(3,3,3,3) = KT(3);
%               KTfull(1,1,1,2) = KT(4);
%               KTfull(1,1,1,3) = KT(5);
%               KTfull(1,2,2,2) = KT(6);
%               KTfull(1,3,3,3) = KT(7);
%               KTfull(2,2,2,3) = KT(8);
%               KTfull(2,3,3,3) = KT(9);
%               KTfull(1,1,2,2) = KT(10);
%               KTfull(1,1,3,3) = KT(11);
%               KTfull(2,2,3,3) = KT(12);
%               KTfull(1,1,2,3) = KT(13);
%               KTfull(1,2,2,3) = KT(14);
%               KTfull(1,2,3,3) = KT(15);
%
%               The other elements of KTfull are obtained from the above 15
%               values noting that KTfull(i, j, k, l) is invariant to the order
%               of indices, e.g., KTfull(l, k, i, j) = KTfull(i, j, k, l)
%
%   nexcl       Number of excluded gradients directions (outliers) when fitting the
%               tensors 
%
%   D_viol      Proportion of violations of directional diffusivities
%
%   Kmin_viol   Proportion of violations of minimum directional kurtoses
%
%   Kmax_viol   Proportion of violations of maximum directional kurtoses

% Author: Ali Tabesh
% Last modified: 07/30/12

% Sb filter

if x(1) <= 0    % skip the voxel if b0 voxel is less than 0

    Dax   = 0;
    Drad  = 0;
    FA    = 0;
    Davg  = 0;
    Kmean = 0;
    K2    = 0;
    K3    = 0;
    Kax   = 0;
    Krad  = 0;
    dt = zeros(6, 1);
    L  = [0, 0, 0]';
    V  = zeros(9, 1);
    kt = zeros(15, 1);
    nexcl  = 0;
    D_viol = 0;
    Kmin_viol = 0;
    Kmax_viol = 0;
    fit_err   = 0;
    return
end

x(x < Sbeps) = Sbeps;    % set voxel values smaller than Sbeps to Sbeps

% form left-hand side of system of linear equations
% log (Sb / S0)

B = log(x(2:end)) - log(x(1));


%--------------------------------------------------------------------------
% step 1: compute diffusion and kurtosis tensors
%--------------------------------------------------------------------------

% initial unconstrained solution
if dki_method.linear_weighting == 0            % unweighted
    DKG = DKGu;
    t   = DKGupinv * B;
elseif dki_method.linear_weighting == 1        % weighting according to signal magnitude
    w   = x(2:end);
    DKG = diag(w) * DKGu;
    B   = w .* B;
    t   = mlsei(DKG, B);
    dki_method.w = w;
%     if dki_method.robust_option == 2
%         DKGprodinv = inv(DKG' * DKG);
%     end
else
    error('Invalid dki_method.linear_weighting option! Parameter dki_method.linear_weighting can only be 0 or 1!')
end

% outlier detection and removal

if dki_method.robust_option == 0
    nexcl = 0;

elseif dki_method.robust_option == 1 % || dki_method.robust_option == 2
    % determine prediction residuals and covariance matrix of estimated parameters (t)
    r = B - DKG * t;
    mse = sum(r .^ 2) / (size(DKG, 1) - size(DKG, 2));
    sigma = mse * DKGprodinv;
        
    % determine outlier detection threshold
    [log_signal_pred, outlier_threshold] = dki_predci(DKG, t, r, sigma, dki_method);
    
    % if the fitted signal does not match the actual signal, do a robust
    % fit with the regular nonlinear least-squares fit as the initial guess
    if any(abs(B - log_signal_pred) > outlier_threshold)

        % do a robust linear fit with the default Tukey's biweight function
        [t, r, sigma] = linrobustfit(DKG, B);

        % determine outlier detection threshold
        [log_signal_pred, outlier_threshold] = dki_predci(DKG, t, r, sigma, dki_method);

    end

    % identify 'inliers'
    idx_incl = find(abs(B - log_signal_pred) < outlier_threshold);
    nexcl = size(B, 1) - length(idx_incl);
    
    % update DKG and B by excluding the outliers and do a linear
    % least-squares fit using the updated DKG and B
    
    if nexcl > 0
        DKG = DKG(idx_incl, :);
        B = B(idx_incl);
        t = mlsei(DKG, B);
    end
else
    error('Invalid dki_method.robust_option option! Parameter dki_method.robust_option can only be 0 or 1!')
end
   
% count the number of constraint violations
% this is always done, but the maps are only written to file if dki_method.linear_violations == 1

if dki_method.linear_violations == 0 || dki_method.linear_violations == 1

    if Kmin ~= 0
        D = DG * t(1:6);    % directional diffusivities (D)
        D(D < 0) = 0;       % set negative D's to zero
    else
        D = zeros(n, 1);
    end

    % estimate c using D

    c = [zeros(n, 1); -D .^ 2 .* Kmin * b(2); zeros(n, 1)];
    
    % count the proportion of constraints that are violated

    viol = C * t > c;
    D_viol = sum(viol(1:n)) / n;
    Kmin_viol = sum(viol(n + 1:2 * n)) / n;
    Kmax_viol = sum(viol(2 * n + 1:3 * n)) / n;

else
    error('Invalid dki_method.linear_violations option! Parameter dki_method.linear_violations must be either 0 or 1!')
end

% run the constrained algorithm if desired (dki_method.linear_constrained == 1)
% and necessary (C*t > c)

if dki_method.linear_constrained == 1

    % estimate D using initial unconstrained solution

    if Kmin ~= 0
        D = DG * t(1:6);    % directional diffusivities (D)
        D(D < 0) = 0;       % set negative D's to zero
    else
        D = zeros(n, 1);
    end

    % estimate c using D

    c = [zeros(n, 1); -D .^ 2 .* Kmin * b(2); zeros(n, 1)];

    % if constraints are not satisfied, invoke iterative lsqlin

    if any(C * t > c)
        t = iterative_lsqlin(t, DKG, B, C, c, C_eq, c_eq, b, g, ndir, Kmin, NKmax);
    end

elseif dki_method.linear_constrained ~= 0
    
    error('Invalid dki_method.linear_constrained option! Parameter dki_method.linear_constrained must be either 0 or 1!')
    
end

% do a nonlinear least-squares fit if desired (dki_method.nonlinear == 1), with
% the unconstrained linear least-squares fit as the initial guess

if dki_method.nonlinear == 1

    if dki_method.linear_constrained == 1
        error('Invalid dki_method.linear_constrained option! When dki_method.nonlinear = 1, dki_method.linear_constrained must be 0.')
    end
    if dki_method.linear_weighting == 1
        error('Invalid dki_method.linear_weighting option! When dki_method.nonlinear = 1, dki_method.linear_weighting must be 0.')
    end
    t = nlinfit(DKG, exp(B), @dki_model, t, nlinfit_options);

elseif dki_method.nonlinear ~= 0

    error('Invalid dki_method.nonlinear option! Parameter dki_method.nonlinear must be either 0 or 1!')
    
end


% diffusion tensor
dt = t(1:6);

% form full diffusion tensor and compute its e'vals and e'vecs

DT = [dt(1), dt(4), dt(5); dt(4), dt(2), dt(6); dt(5), dt(6), dt(3)];
[V, L] = eig(DT);
L = diag(L);  % take the diagonal of L

% filter DT e'vals
L(L < Leps) = Leps;
    
% sort eigenvalues in ascending order
[L, idx] = sort(L);

% sort e'vecs according to e'vals
V = V(:, idx);      

% mean diffusivity
Davg = mean(L);             % average of DT e'vals

% kurtosis tensor
kt = t(7:21) / Davg ^ 2 / b(2);   % recall that columns of DKG corresponding to kt are scaled by 1/b(2)


%--------------------------------------------------------------------------
% step 2: rotate kurtosis tensor to reference frame
%--------------------------------------------------------------------------

% form full kurtosis tensor and rotate using e'vecs of DT

ktfull = zeros(3, 3, 3, 3);

% initialize 15 elements
% [W1111 W2222 W3333 W1112 W1113 W1222 W1333 W2223 W2333 W1122 W1133 W2233 W1123 W1223 W1233]
ktfull(1, 1, 1, 1) = kt(1);
ktfull(2, 2, 2, 2) = kt(2);
ktfull(3, 3, 3, 3) = kt(3);
ktfull(1, 1, 1, 2) = kt(4);
ktfull(1, 1, 1, 3) = kt(5);
ktfull(1, 2, 2, 2) = kt(6);
ktfull(1, 3, 3, 3) = kt(7);
ktfull(2, 2, 2, 3) = kt(8);
ktfull(2, 3, 3, 3) = kt(9);
ktfull(1, 1, 2, 2) = kt(10);
ktfull(1, 1, 3, 3) = kt(11);
ktfull(2, 2, 3, 3) = kt(12);
ktfull(1, 1, 2, 3) = kt(13);
ktfull(1, 2, 2, 3) = kt(14);
ktfull(1, 2, 3, 3) = kt(15);

% propagate values to other elements
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                I = sort([i, j, k, l,]);
                ktfull(i, j, k, l) = ktfull(I(1), I(2), I(3), I(4));
            end
        end
    end
end

% we only need Wt1111, Wt2222, Wt3333, Wt1122, Wt1133, Wt2233
kt_rot = zeros(6, 1);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                kt_rot(1) = kt_rot(1) + V(i, 1) * V(j, 1) * V(k, 1) * V(l, 1) * ktfull(i, j, k, l);  % Wt1111
                kt_rot(2) = kt_rot(2) + V(i, 2) * V(j, 2) * V(k, 2) * V(l, 2) * ktfull(i, j, k, l);  % Wt2222
                kt_rot(3) = kt_rot(3) + V(i, 3) * V(j, 3) * V(k, 3) * V(l, 3) * ktfull(i, j, k, l);  % Wt3333
                kt_rot(4) = kt_rot(4) + V(i, 1) * V(j, 1) * V(k, 2) * V(l, 2) * ktfull(i, j, k, l);  % Wt1122
                kt_rot(5) = kt_rot(5) + V(i, 1) * V(j, 1) * V(k, 3) * V(l, 3) * ktfull(i, j, k, l);  % Wt1133                
                kt_rot(6) = kt_rot(6) + V(i, 2) * V(j, 2) * V(k, 3) * V(l, 3) * ktfull(i, j, k, l);  % Wt2233                
            end
        end
    end
end


%--------------------------------------------------------------------------
% step 3: scalar parameters from diffusion and kurtosis tensors
%--------------------------------------------------------------------------

% scalars extracted from DT
Dax  = L(3);
Drad = 0.5 * (L(1) + L(2));
FA   = sqrt(((L(1) - L(2)) ^ 2 + (L(1) - L(3)) ^ 2 + (L(2) - L(3)) ^ 2) / (2 * sum(L .^ 2)));
V    = V(:);   % e'vecs of the DT

% scalars extracted from KT
% Kmean updated on 2009-03-25 based on the correction to the analytical formula

Kmean = 6 * A1122(L(1), L(2), L(3), err, errtol, const) * kt_rot(4) + ...
    6 * A1122(L(1), L(3), L(2), err, errtol, const) * kt_rot(5) + ...
    6 * A1122(L(2), L(3), L(1), err, errtol, const) * kt_rot(6) + ...
    A3333(L(2), L(3), L(1), err, errtol, const) * kt_rot(1) + ...
    A3333(L(1), L(3), L(2), err, errtol, const) * kt_rot(2) + ...
    A3333(L(1), L(2), L(3), err, errtol, const) * kt_rot(3);
K2   = sum(L) ^ 2 / (9 * L(2) ^ 2) * kt_rot(2);    
% K3 is defined according to the common convention (as opposed to the rest of this function)
% where lambda_3 is the smallest eigenvalue
K3   = sum(L) ^ 2 / (9 * L(1) ^ 2) * kt_rot(1);     
Kax  = sum(L) ^ 2 / (9 * L(3) ^ 2) * kt_rot(3);
Krad = 6 * C1122(L(1), L(2), L(3), err) * kt_rot(4) + ...
    C1111(L(1), L(2), L(3), err) * kt_rot(1) + ...
    C1111(L(2), L(1), L(3), err) * kt_rot(2);

% Kmean, Kax, Krad filter

Kmean(Kmean > Kmax_final) = Kmax_final;
Kax(Kax > Kmax_final)     = Kmax_final;
Krad(Krad > Kmax_final)   = Kmax_final;
K2(K2 > Kmax_final)       = Kmax_final;
K3(K3 > Kmax_final)       = Kmax_final;

Kmean(Kmean < Kmin_final) = Kmin_final;
Kax(Kax < Kmin_final)     = Kmin_final;
Krad(Krad < Kmin_final)   = Kmin_final;
K2(K2 < Kmin_final)       = Kmin_final;
K3(K3 < Kmin_final)       = Kmin_final;

% relative fit error (fraction of signal power unexplained by the model)

fit_err = sum((B - DKG * t) .^ 2) ./ sum(B .^ 2);
fit_err(isinf(fit_err)) = 0;



%--------------------------------------------------------------------------
% auxiliary functions for dke_core
%--------------------------------------------------------------------------

% compute A1122
function y = A1122(L1, L2, L3, err, errtol, const)
if abs(L2 - L1) > err
    y = ((L1 + L2 + L3) ^ 2) / (18 * (L1 - L2) ^ 2) * ((L1 + L2) / sqrt(L1 * L2) * ...
        compute_rf_mex(L3 / L1, L3 / L2, 1, errtol, const.rf_lolim, const.rf_uplim) + ...
        (2 * L3 - L1 - L2) / (3 * sqrt(L1 * L2)) * ...
        compute_rd_mex(L3 / L1, L3 / L2, 1, errtol, const.rd_lolim, const.rd_uplim) - 2);
elseif abs(L2 - L1) <= err &&  abs(L3 - L2) > err
    y = (2 * L1 + L3) ^ 2 / (144 * (L1 ^ 2) * (L3 - L1) ^ 2) * ...
        (L1 * (2 * L1 + L3) + L3 * (L3 - 4 * L1) * alpha(1 - L3 / L1));
else
    y = 1 / 15;
end


% compute alpha(x) for A1122 with singularity
function y = alpha(x)
if x > 0
    y = 1 / sqrt(x) * atanh(sqrt(x));
else
    y = 1 / sqrt(-x) * atan(sqrt(-x));
end


% compute A3333
function y = A3333(L1, L2, L3, err, errtol, const)
if abs(L3 - L1) > err && abs(L3 - L2) > err
    y = ((L1 + L2 + L3) ^ 2) / (18 * (L3 - L1) * (L3 - L2)) * ...
        (sqrt(L1 * L2) / L3 * compute_rf_mex(L3 / L1, L3 / L2, 1, errtol, const.rf_lolim, const.rf_uplim) + ...
        (3 * L3 ^ 2 - L1 * L2 - L1 * L3 - L2 * L3) / (3 * L3 * sqrt(L1 * L2)) * ...
        compute_rd_mex(L3 / L1, L3 / L2, 1, errtol, const.rd_lolim, const.rd_uplim) - 1);
elseif abs(L3 - L1) <= err && abs(L3 - L2) > err
    y = 3 * A1122(L1, L1, L2, err, errtol);
elseif abs(L3 - L1) > err && abs(L3 - L2) <= err
    y = 3 * A1122(L2, L2, L1, err, errtol);
else
    y = 0.2;
end


% compute C1111
function y = C1111(L1, L2, L3, err)
if abs(L1 - L2) > err
    y = (L1 + L2 + L3) ^ 2 / (18 * L1 * (L1 - L2) ^ 2) * (2 * L1 + (L2 ^ 2 - 3 * L1 * L2) / sqrt(L1 * L2));
else
    y = (2 * L1 + L3) ^ 2 / (24 * L1 ^ 2);
end


% compute C1122
function y = C1122(L1, L2, L3, err)
if abs(L1 - L2) > err
    y = (L1 + L2 + L3) ^ 2 / (18 * (L1 - L2) ^ 2) * ((L1 + L2) / sqrt(L1 * L2) - 2);
else
    y = (2 * L1 + L3) ^ 2 / (72 * L1 ^ 2);
end


%--------------------------------------------------------------------------
% white matter maps
%--------------------------------------------------------------------------

function [f, Da, De_axial, De_radial, tortuosity] = wm_model_estimate(dt, kt, Dax, Drad, Dmean, DG, DGpinv, KG, DG_maxK, KG_maxK)

f_lo = 0.0001;
f_hi = 0.9999;

% find maximum directional kurtosis
kmax = max_K(dt, kt, Dmean, DG_maxK, KG_maxK);

% axonal water fraction
f = kmax / (kmax + 3);

% constraints on f
if f < 0
    f = 0;
end

if isnan(f)
    f = 0;
end

if f < f_lo
    Da = 0;
    De_axial = Dax;
    De_radial = Drad;
    tortuosity = De_axial / De_radial;
elseif f > f_hi
    Da = 3 * Dmean;
    De_axial = 0;
    De_radial = 0;
    tortuosity = 0; 
else
    % directional diffusivities and kurtoses
    ddir = DG * dt;
    kdir = Dmean ^ 2 .* (KG * kt) ./ (ddir .^ 2);
    
    % extra-axonal diffusivities
    De = ddir .* (1 + sqrt(kdir * f / (3 * (1 - f))));
    Le = dt_estimate(De, DGpinv);
    
    Di = ddir .* (1 - sqrt(kdir * (1 - f) / (3 * f)));
    Li = dt_estimate(Di, DGpinv);
    
    Da = sum(Li);
    De_axial = Le(1);
    De_radial = (Le(2) + Le(3)) / 2;
    tortuosity = De_axial / De_radial;

end


% find minimum and maximum directional kurtoses
function maxK = max_K(dt, kt, Dmean, DG, KG)

K = Dmean ^ 2 * (KG * kt) ./ (DG * dt) .^ 2;
maxK = max(K, [], 1);


% solve for extra- and intra-axonal diffusion tensors
function L = dt_estimate(ddir, DGpinv)

dt = DGpinv * ddir;
DT = [dt(1), dt(4), dt(5); dt(4), dt(2), dt(6); dt(5), dt(6), dt(3)];
[~, L] = eig(DT);
% [V, L] = eig(DT);
L = diag(L);
L = sort(L, 'descend');


%--------------------------------------------------------------------------
% DKI signal model for nonlinear fitting
%--------------------------------------------------------------------------

function signal_pred = dki_model(t, DKG)
signal_pred = exp(DKG * t);


%--------------------------------------------------------------------------
% DKI log signal model for linear fitting
%--------------------------------------------------------------------------

function log_signal_pred = dki_log_signal(t, DKG)
log_signal_pred = DKG * t;


%--------------------------------------------------------------------------
% DTI log signal model for linear fitting
%--------------------------------------------------------------------------

function log_signal_pred = dti_log_signal(t, DG)
log_signal_pred = DG * t;


%--------------------------------------------------------------------------
% predicted DKI signal and confidence interval half-width
%--------------------------------------------------------------------------

function [log_signal_pred, outlier_threshold] = dki_predci(DKG, t, r, sigma, dki_method)
    
if dki_method.robust_option == 1
    log_signal_pred = DKG * t;
    if dki_method.linear_weighting == 1
        outlier_threshold = dki_method.noise_tolerance * dki_method.w;
    else
        outlier_threshold = dki_method.noise_tolerance;
    end
else
    [log_signal_pred, outlier_threshold] = nlpredci(@dki_log_signal, DKG, t, r, 'covar', sigma, 'alpha', dki_method.significance_level);
end


%--------------------------------------------------------------------------
% predicted DTI signal and confidence interval half-width
%--------------------------------------------------------------------------

function [log_signal_pred, outlier_threshold] = dti_predci(DG, t, r, sigma, dti_method)

if dti_method.robust_option == 1
    log_signal_pred = DG * t;
    if dti_method.linear_weighting == 1
        outlier_threshold = dti_method.noise_tolerance * dti_method.w;
    else
        outlier_threshold = dti_method.noise_tolerance;
    end
else
    [log_signal_pred, outlier_threshold] = nlpredci(@dti_log_signal, DG, t, r, 'covar', sigma, 'alpha', dti_method.significance_level);
end
    

%--------------------------------------------------------------------------
% robust linear fitt
%--------------------------------------------------------------------------

function [X, r, sigma] = linrobustfit(A, b)
[X, stats] = robustfit(A, b, 'bisquare', [], 'off');
r = b - A * X;
sigma = stats.coeffcorr .* (stats.se * stats.se');




%--------------------------------------------------------------------------
% core function to estimate diffusion tensor-derived parameters 
%--------------------------------------------------------------------------

function [Dax, Drad, Davg, FA, dt, L, V, nexcl, fit_err] = dte_core(x, DG, DGpinv, DGprodinv, Leps, Sbeps, dti_method)
% dte_core   Estimate diffusion tensor-derived parameters for a single voxel
%
% Syntax
%
%   [Dax, Drad, Davg, FA, DT, L, V] = dke_core(x, DG, DGpinv, DGprodinv, Leps, Sbeps, dti_method)
%
% Inputs
%
%   x           (ndir+1)-by-1 vector of signal values
%               The first row of x contains the b = 0 signal values, while the next
%               ndir rows contain the signal values for the single b used for
%               DTI computations
%
%   DG          Coefficients matrix
%
%   DGpinv      Pseudoinverse of coefficients matrix
%
%   DGprodinv   Inverse of (DG' * DG) used for estimating covariance matrix
%               of estimated DT
%
%   Leps        Smallest acceptable DT eigenvalue (mm^2/s units)
%               Eigenvalues smaller than Leps are set to Leps.
%
%   Sbeps       Smallest acceptable signal value for nonzero b values
%               Signal values smaller than Sbeps are set
%               to Sbeps.
%
%   dti_method  tensor estimation method
%
% Outputs
%
%   Dax     Axial diffusivity (mm^2/s units)
%
%   Drad    Radial diffusivity (mm^2/s units)
%
%   Davg    Mean diffusivity as average of DT eigenvalues (mm^2/s units)
%
%   FA      Fractional anisotropy
%
%   DT      6-by-1 vector containing elements of DT (mm^2/s units)
%           DT in matrix form denoted as DTfull can be obtained as 
%           DTfull = [DT(1) DT(4) DT(5); DT(4) DT(2) DT(6); DT(5) DT(6) DT(3)]
%
%   L       3-by-1 vector of DT eigenvalues (mm^2/s units) sorted in descending order
%   
%   V       9-by-1 vector containing the eigenvectors of the DT
%           corresponding to eigenvalues in L
%
%   nexcl   Number of excluded gradient directions

% Author: Ali Tabesh
% Last modified: 07/22/10

% Sb filter

if x(1) <= 0    % skip the voxel if b0 voxel is less than 0
    Dax     = 0;
    Drad    = 0;
    Davg    = 0;
    FA      = 0;
    dt      = zeros(6, 1);
    L       = zeros(3, 1);
    V       = zeros(9, 1);
    nexcl   = 0;
    fit_err = 0;
    return
end

% set voxel values smaller than Sbeps to Sbeps
x(x < Sbeps) = Sbeps;    

log_signal = log(x(2:end)) - log(x(1));

% compute diffusion tensor

if dti_method.linear_weighting == 1        % weighting according to signal magnitude
    w = x(2:end);
    DG = diag(w) * DG;
    DGpinv = pinv(DG);
    log_signal = w .* log_signal;
    dti_method.w = w;
elseif dti_method.linear_weighting ~= 0
    error('Invalid dti_method.linear_weighting option! Parameter dti_method.linear_weighting must be either 0 or 1.')
end

dt = DGpinv * log_signal;

% do a robust fit if desired (dti_method.robust_option == 1)

if dti_method.robust_option == 0
    nexcl = 0;
elseif dti_method.robust_option == 1 % || dti_method.robust_option == 2
    % determine covariance matrix of estimated tensor (dt)
    r = log_signal - DG * dt;
    mse = sum(r .^ 2) / (size(DG, 1) - size(DG, 2));
    sigma = mse * DGprodinv;
        
    % determine outlier detection threshold
    [log_signal_pred outlier_threshold] = dti_predci(DG, dt, r, sigma, dti_method);

    % if the fitted signal does not match the actual signal, do a robust
    % fit with the regular nonlinear least-squares fit as the initial guess
    if any(abs(log_signal - log_signal_pred) > outlier_threshold)

        % do a robust linear fit with the default Tukey's biweight function
        [dt r sigma] = linrobustfit(DG, log_signal);
        
        % determine outlier detection threshold
        [log_signal_pred outlier_threshold] = dti_predci(DG, dt, r, sigma, dti_method);
    end

    % identify 'inliers'
    idx_incl = find(abs(log_signal - log_signal_pred) < outlier_threshold);
    nexcl = size(log_signal, 1) - length(idx_incl);
    
    % update DG and signal array by excluding the outliers and do a linear
    % least-squares fit using the updated DG and signal array
    
    if nexcl > 0
        DG = DG(idx_incl, :);
        log_signal = log_signal(idx_incl);
        dt = mlsei(DG, log_signal);
    end     
else
    error('Invalid dti_method.robust_option option!')
end

% form full diffusion tensor and compute its e'vals and e'vecs
DT = [dt(1), dt(4), dt(5); dt(4), dt(2), dt(6); dt(5), dt(6), dt(3)];

[V, L] = eig(DT);
L = diag(L);  % take the diagonal of L

% filter DT e'vals

L(L < Leps) = Leps;
    
% sort eigenvalues in ascending order

[L idx] = sort(L);
V = V(:, idx);      % sort e'vecs according to e'vals
V = V(:);   % e'vecs of the DT

% scalar measures
Dax = L(3);
Drad = 0.5 * (L(1) + L(2));
Davg = mean(L);             % average of DT e'vals
FA = sqrt(((L(1) - L(2)) ^ 2 + (L(1) - L(3)) ^ 2 + (L(2) - L(3)) ^ 2) / (2 * sum(L .^ 2)));

% relative fit error (fraction of signal power unexplained by the model)

fit_err = sum((log_signal - DG * dt) .^ 2) ./ sum(log_signal .^ 2);
fit_err(isinf(fit_err)) = 0;



%--------------------------------------------------------------------------
% apply Gaussian smoothing to a set of input images
%--------------------------------------------------------------------------

function y = filter_img_set(x, siz_x, fn_noise, s_x, s_n, v)

% filter_img_set Filter noise from diffusion-weighted images using a noise map and
%       spatial smoothing
%
% Syntax
%
%   [y z] = filter_img_set(x, fn_noise, s_x, s_n, v)
%
% Inputs
%
%   x           Cell array containing diffusion-weighted images
%
%   siz_x       1-by-3 vector specifying the dimensions of the diffusion-weighted images
%
%   fn_noise    Noise file name
%
%   s_x         FWHM (in mm units) for Gaussian smoothing applied to diffusion-weighted images
%               s_x = 0 means no filtering
%
%   s_n         FWHM (in mm units) for Gaussian smoothing applied to noise image
%               s_n = 0 means no filtering
%
%   v           voxel dimensions (mm)
%
% Outputs
%
%   y           n-by-m matrix containing the filtered voxels
%               n = b * d + 1, where b is the number of b-values and d
%               is the number of gradient directions, and m is the
%               number of voxels in the input images. The first row of
%               x contains the b = 0 image, while the next d rows
%               contain the b = b_1 images, and so on.

% Author: Ali Tabesh
% Last modified: 05/25/12

% initialize

nimg = size(x, 4);

% read and optionally filter noise image

if ~isempty(fn_noise)
    if exist(fn_noise, 'file')
        [hdr, n] = read_nii(fn_noise);
        if any(hdr.dim(2:4) ~= siz_x)
            error('Noise image and diffusion-weighted images have different dimensions!')
        end
    else
        error('Noise image file %s does not exist!', fn_noise)
    end
else
    n = zeros(siz_x);
end

if s_n ~= 0
    n = gaussian_smooth(n, s_n, v);
end

% filter diffusion-weighted images
y = zeros(nimg, numel(x(:, :, :, 1)));

for i = 1:nimg
    if sum(s_x) ~= 0
        x(:, :, :, i) = gaussian_smooth(squeeze(x(:, :, :, i)), s_x, v);
    end
    tmp = x(:, :, :, i);
    y(i, :) = sqrt((tmp(:)) .^ 2 - (n(:)) .^ 2)';
end



%--------------------------------------------------------------------------
% find the brain in b = 0 image
%--------------------------------------------------------------------------

function y = find_brain_mask(x, T)

%find_brain_mask Create a binary brain mask from b = 0 image
%
%Syntax
%
%   y = find_brain_mask(x, T)
%
%Inputs/output
%
%   x   b = 0 image
%   T   Threshold for finding the brain mask 
%   y   Output binary mask
%
%Note
%
%   This function invokes MATLAB's Image Processing Toolbox functions
%   bwlabeln() and imfill().

% Author: Ali Tabesh
% Last modified: 04/23/09

% define the neighborhood matrix for bwlabeln
c = zeros(3, 3, 3);
c(:, :, 1) = [0 0 0; 0 1 0; 0 0 0];
c(:, :, 2) = [0 1 0; 1 1 1; 0 1 0];
c(:, :, 3) = [0 0 0; 0 1 0; 0 0 0];

% create initial brain mask by thresholding using T

y = x > T;                                      % apply threshold

% clean up the initial mask

[L, n] = bwlabeln(y, c);                         % label connected components
 A = zeros(1,n); %EM
for i = 1:n                                     % find areas of all connected components
    A(i) = length(find(L(:) == i));             
end
% [a, l] = max(A);                                 % find the largest connected component
[~, l] = max(A);                                 % find the largest connected component
y = L == l;                                     % keep only the largest connected component
for i = 1:size(y, 3)
    y(:, :, i) = imfill(y(:, :, i), 'holes');       % fill holes (in 2D) in the largest connected component
end



%--------------------------------------------------------------------------
% apply Gaussian smoothing to an input image
%--------------------------------------------------------------------------

function J = gaussian_smooth(I, s, v)
%gaussian_smooth Smooth image using a Gaussian filter
%
% Syntax
%
%   J = gaussian_smooth(I, s, v)
%
% Inputs/output
%
%   I   Input volume
%   s   Full width at half max 
%       s is a scalar for isotropic filtering or a 1-by-3 vector for
%       different degrees of smoothing along x, y and z
%   v   Voxel size (1-by-3 vector)
%   J   Output image
%
% Notes
%
% 1. Zero-padding along x and y is used for signal extension. Along z, the
%   filter support is truncated at the image boundary. The treatment is 
%   equivalent to that in SPM2.
%
% 2. Filter is always of odd length.
%
% 3. Filter kernel is truncated at 4*std dev of the Gaussian.

% Author: Ali Tabesh
% Last modified: 05/01/09

% initialize

I = double(I);                  % convert to double precision

c = 4;                       	% number of standard deviations to keep in the Gaussian kernel

if length(s) == 1
    s = [s, s, s]; 
end

s = s ./ v;                      % FWHM in voxel units
s = s / sqrt(8 * log(2));        % convert FWHM to std dev

% kernel support

t = ceil(c * s(1)); 
x = -t:t;                       % kernel support along x
t = ceil(c * s(2)); 
y = -t:t;                       % kernel support along y
t = ceil(c * s(3));
z = -t:t;                       % kernel support along z

% form filters
%
% sum of kernel coefficients = 1

if s(1) == 0
	x = 1;
else
	x = exp(-(x) .^ 2 / (2 * s(1) .^ 2));
end

if s(2) == 0
	y = 1;
else
    y = exp(-(y) .^ 2 / (2 * s(2) .^ 2));
end

if s(3) == 0
    z = 1;
else
    z = exp(-(z) .^ 2 / (2 * s(3) .^ 2));    
end

x = x / sum(x);
y = y / sum(y);                  
z = z / sum(z);

u = (length(x) - 1) / 2;
v = (length(y) - 1) / 2;
w = (length(z) - 1) / 2;

% separable 3d convolution 
%
% Extension is zero-padding along x and y; along z, the filter is truncated
% at the image boundary
%
% Extra zeros are added to enable keeping track of the center part of output

[nx ny nz] = size(I);

J = cat(1, I, zeros(u, ny, nz));        % zero-pad I along dimension 1 (x) to get J
J = filter(x, 1, J, [], 1);             % apply filter x to J
J = J((u+1):end, :, :);                 % discard the extra zeros

J = cat(2, J, zeros(nx, v, nz));
J = filter(y, 1, J, [], 2);
J = J(:, (v + 1):end, :);

% filter along z
%
% first compute the convolution with zero-padding 

I = cat(3, J, zeros(nx, ny, w));
J = filter(z, 1, I, [], 3);
J = J(:, :, (w + 1):end);

% then correct results for boundary voxels

for i = 1:min(w, size(J, 3))
    J(:, :, i) = J(:, :, i) ./ repmat(sum(z(((w+2-i):end))), nx, ny);
end

for i = max(nz - w + 1, 1):nz
    j = i - nz + w;
    J(:, :, i) = J(:, :, i) ./ repmat(sum(z((j + 1):end)), nx, ny);
end



%--------------------------------------------------------------------------
% read a set of diffusion-weighted images
%--------------------------------------------------------------------------

function [x hdr_b0 img_size voxel_size] = read_img_set_nii(fn_prefix, ndir, bval, start, numbering_pattern)

%read_img_set Read diffusion-weighted images from files
%
% Syntax
%   
%   [x hdr img_size voxel_size] = read_img_set(fn_prefix, ndir, bval, start, format)
%
% Inputs
%
%   fn_prefix       Prefix for input image ('img') files (e.g., 'rdki')
%
%                   Image file names are expected to follow the following
%                   naming convention: 
%                   "prefix underscore b-value underscore direction .img"
%                   Examples: rdki_500_1.img, rdki_1000_005.img
%
%   ndir            Number of gradient directions
%
%   bval            Vector of non-zero b values (e.g., [1000 2000])
%
%   start           (optional) gradient direction file name starting number
%                   (default = 0)
%
%                   This parameter is useful when the first gradient
%                   direction file number starts from 1. (rdki_500_1 rather
%                   than rdki_500_0)
%
%   format          (optional) format string for the file number (default =
%                   '1', '2', in 'rdki_500_1', 'rdki_500_2', ...)
%
%                   This option is useful when the file numbering has
%                   preceeding zeros, e.g., rdki_500_001, rdki_500_002, ...
%                   In the above example, format should be '%03d'.
%                   For more information on format, see help on function
%                   sprintf.
%
% Outputs
%   
%   x               Cell array containing the diffusion-weighted images                
%                   x{1} contains b = 0 image, x{2:(ndir+1)} contains b =
%                   b_1 image, and so on.
%
%   img_size        3-by-1 vector containing the dimensions of the input
%                   images
%   voxel_size      3-by-1 vector containing voxel dimensions

% Author: Ali Tabesh
% Last modified: 05/30/12

default_fn_numbering = 1;   % use default file numbering convention
startidx = 0;               % starting file number

if nargin == 1  % assume 4D nii
    [hdr_b0 x] = read_nii(fn_prefix);
    img_size = hdr_b0.dim(2:4);
    voxel_size = hdr_b0.pixdim(2:4);    
    return
elseif nargin < 3
    error('Insufficient number of input arguments!')
end
if nargin >= 4
    startidx = start;
end
if nargin == 5
    default_fn_numbering = 0;
end
if nargin > 5
    error('Too many input arguments!')
end

% number of non-zero b-values
nbval = length(bval);

% read b0 img
fn_img = [fn_prefix '_0.nii'];
hdr_b0 = read_nii(fn_img, 0);  % read header only

% determine image and voxel size
img_size = hdr_b0.dim(2:4);
voxel_size = hdr_b0.pixdim(2:4);    

% initialize and read x

x = zeros([img_size sum(ndir) + 1]);

[hdr_b0 img] = read_nii(fn_img);
x(:,:,:,1) = img;

for ibval = 1:nbval
    for idir = startidx:(startidx + ndir(ibval) - 1)
        if default_fn_numbering
            fn_img = [fn_prefix '_' num2str(bval(ibval)) '_' num2str(idir) '.nii'];
        else
            fn_img = [fn_prefix '_' num2str(bval(ibval)) '_' sprintf(numbering_pattern, idir) '.nii'];
        end
        [hdr, img] = read_nii(fn_img);
        if any(hdr.dim(2:4) ~= hdr_b0.dim(2:4))
            error('Diffusion-weighted images have inconsistent dimensions!')
        else
            x(:, :, :, idir + sum(ndir(1:(ibval-1))) - startidx + 2) = img;
        end
    end
end



%--------------------------------------------------------------------------
% create FA color maps
%--------------------------------------------------------------------------

function make_fa_color_map(FA, V, img_size, fn_out)

% make_fa_color_map Create an FA color map
%
% Syntax
%
%	make_fa_color_map(FA, V, img_size, fn_out)
%
% Inputs
%
%   FA          1-by-nvoxel FA map where nvoxel is the number of image voxels
%   V           9-by-nvoxel vector of diffusion tensor eigenvectors
%   img_size    1-by-3 vector specifying image dimensions
%   fn_out      Prefix for output FA color maps in TIFF format
%
% Output
%
%   A series of TIFF files will be created, where each TIFF file corresponds to a slice of the input volume

% Author: Ali Tabesh
% Last modified: 11/12/09

r = reshape(FA .* (V(7, :)) .^ 2, img_size);  % V(7, :) is 1st element of e'vec corresponding to largest e'val
g = reshape(FA .* (V(8, :)) .^ 2, img_size);  % V(8, :) is 2nd element of e'vec corresponding to largest e'val
b = reshape(FA .* (V(9, :)) .^ 2, img_size);  % V(9, :) is 3rd element of e'vec corresponding to largest e'val

for islice = 1:img_size(3)
    tif(:,:,1) = rot90(squeeze(r(:, :, islice)), 1);
    tif(:,:,2) = rot90(squeeze(g(:, :, islice)), 1);
    tif(:,:,3) = rot90(squeeze(b(:, :, islice)), 1);
    imwrite(tif, [fn_out sprintf('%03d', islice), '.tif']);
end



%--------------------------------------------------------------------------
% write image to NIfTI file
%--------------------------------------------------------------------------

function write_img(fn, img, hdr, datatype, filter_method, interpolation_method, tensor_flag)

if tensor_flag == 1 % 4D image
    % make sure hdr.dim corresponds to a 3D image (input image may be 4D whereas maps are 3D)
    hdr.dim(1) = 4;
    hdr.dim(5) = size(img, 1);
    % make sure img is a 4D array
    % NOTE: tensor data must be transposed.
    img = reshape(img', hdr.dim(2:5));

%GRG (changed else to else if: tensor_flag == 2 does not require dimension
%or image changes) --------------------------------------------------------
elseif tensor_flag == 0; %3D image 
%--------------------------------------------------------------------------
    % make sure hdr.dim corresponds to a 3D image (input image may be 4D whereas maps are 3D)
    hdr.dim(1) = 3;
    hdr.dim(5) = 1;
    % make sure img is a 3D array
    img = reshape(img, hdr.dim(2:4));
end

%GRG (added input: filter_method.img_res)----------------------------------
% filter image
if filter_method.type
    img = medianfilter(img, filter_method.img_ref, filter_method.img_viol, filter_method.T);
end
%--------------------------------------------------------------------------

% set output data type
if strcmpi(datatype, 'float32')
    hdr.datatype = 16;
    hdr.bitpix   = 32;
elseif strcmpi(datatype, 'int16')
    hdr.datatype = 4;
    hdr.bitpix   = 16;
else
    error('Invalid output data type!')
end

%GRG-----------------------------------------------------------------------
if tensor_flag == 2 %matlab
    if strfind(fn,'KT.mat'); KT = img; save(fn,'KT'); end
    if strfind(fn,'DT.mat'); DT = img; save(fn,'DT'); end
    fprintf('Writing %s\n', fn);
%--------------------------------------------------------------------------

else %nifti
% set output file name
hdr.fn = [fn, '.nii'];

% workaround for SPM bug; not sure if needed for in-house nifti read/write functions
if exist(fn, 'file')
    delete(fn);
end

if interpolation_method.flag && tensor_flag == 0
    fprintf('Writing and interpolating %s\n', hdr.fn);
else
    fprintf('Writing %s\n', hdr.fn);
end

write_nii(hdr, img);

% interpolate
if interpolation_method.flag && tensor_flag == 0
    %eval(['!map_interpolate "' hdr.fn '" ' num2str(interpolation_method.resolution) ' ' num2str(interpolation_method.order)])
    map_interpolate(hdr.fn,num2str(interpolation_method.resolution),num2str(interpolation_method.order))

end

end




%--------------------------------------------------------------------------
% core function to estimate mean diffusivity and mean kurtosis without
% estimating tensors
%--------------------------------------------------------------------------

function [Dax, Drad, FA, Davg, Kmean, K2, K3, Kax, Krad, dt, L, V, kt, nexcl, D_viol, Kmin_viol, Kmax_viol, fit_err] = ...
    dke_core_dirfit(x, n, b, A, Apinv, Sbeps, Kmin, NKmax, Kmin_final, Kmax_final, dki_method)

%dke_core_dirfit Core function to estimate mean diffusivity and mean kurtosis for a
%               single voxel without fitting tensors
%
%Syntax
%
%   [Dax, Drad, FA, Davg, Kmean, K2, K3, Kax, Krad, dt, L, V, kt, nexcl, D_viol, Kmin_viol, Kmax_viol] = dke_core_dirfit(x, n, b, A, Apinv, Sbeps, Kmin, NKmax, Kminf, Kmaxf, dki_method)
%
%Inputs
%
%   x           (ndir+1)-by-1 vector of signal values
%               The first row of x contains the b = 0 signal values, the next
%               ndir_1 rows contain the b = b_1 signal values, and so on, where
%               ndir_1 is the number of directions for b = b_1
%
%   ndir        Number of gradient directions for each b value
%
%   b           1-by-nbval vector of nonzero b values (in s/mm^2 units), where nbval is
%               the number of nonzero b values
%
%   A           coefficients matrix
%
%   Apinv       pseudoinverse of A
%
%   Sbeps       Smallest acceptable signal value for nonzero b values
%               Signal values smaller than Sbeps are set
%               to Sbeps.
%
%   Kmin        Minimum acceptable directional kurtosis
%
%   NKmax       Parameter defining the maximum directional kurtosis as
%               Kmax = NKmax / (Di * bmax)
%               where Di is directional diffusivity along direction i for the given voxel and bmax is the largest b value 
%
%   Kminf       Minimum allowable kurtosis value (mean, axial, radial)
%
%   Kmaxf       Maximum allowable kurtosis value (mean, axial, radial)
%
%   dki_method      Structure containing information for tensor fitting dki_method
%
%Outputs
%
%   Dax         Dummy output (set to zero)
%
%   Drad        Dummy output (set to zero)
%
%   FA          Dummy output (set to zero)
%
%   Davg        Mean diffusivity (mm^2/s units)
%
%   Kmean       Mean kurtosis
%
%   K2          Dummy output (set to zero)
%
%   K3          Dummy output (set to zero)
%
%   Kax         Dummy output (set to zero)
%
%   Krad        Dummy output (set to zero)
%
%   DT          Dummy output (set to zero)
%
%   L           Dummy output (set to zero)
%   
%   V           Dummy output (set to zero)
%
%   KT          Dummy output (set to zero)
%
%   nexcl       Number of excluded gradients directions (outliers) when fitting the
%               tensors 
%
%   D_viol      Dummy output (set to zero)
%
%   Kmin_viol   Dummy output (set to zero)
%
%   Kmax_viol   Dummy output (set to zero)

% Author: Ali Tabesh
% Last modified: 12/21/09

% Sb filter

Dax         = 0;
Drad        = 0;
FA          = 0;
K2          = 0;
K3          = 0;
Kax         = 0;
Krad        = 0;
dt          = zeros(6, 1);
L           = zeros(3, 1);
V           = zeros(9, 1);
kt          = zeros(15, 1);
nexcl       = 0;
D_viol      = 0;
Kmin_viol   = 0;
Kmax_viol   = 0;
fit_err     = 0;

if x(1) <= 0    % skip the voxel if b0 voxel is less than 0
    Davg  = 0;
    Kmean = 0;
    return
end

% set voxel values smaller than Sbeps to Sbeps
x(x < Sbeps) = Sbeps;    

% form left-hand side of system of linear equations
% log (Sb / S0)
B = log(x(2:end)) - log(x(1));
B = reshape(B, [n, length(b)])';

% find directional diffusivities and kurtoses (actually D^2 * K for each direction)

if dki_method.linear_weighting == 0
    X = Apinv * B;
elseif dki_method.linear_weighting == 1        % weighting according to signal magnitude
    w = reshape(x(2:end), [n length(b)])';
    B = w .* B;
   % X = zeros(2, n);%%% Check this , EM 
    for i = 1:n
        X(:, i) = mlsei(diag(w(:, i)) * A, B(:, i));
    end
else
    error('Invalid DKI weighting option!')
end

% impose constraints on directional diffusivities

X(1, X(1, :) < 0) = 0;

% estimate Dmean

Davg = mean(X(1, :));

% estimate directional kursoses
K = zeros(1, size(X, 2));           % initialize directional kurtoses to zero
idx = find(X(1, :) > 0);            % find directions with positive diffusivity
K(idx) = X(2, idx) ./ X(1, idx) .^ 2;      % calculate directional kurtoses for directions with positive diffusivity

% impose minimum constraint on directional kurtoses
K(K < Kmin) = Kmin;

% impose maximum constraint on directional kurtoses
Kmax = zeros(1, size(X, 2));                % initialize maximum directional kurtoses to zero
Kmax(idx) = NKmax ./ (max(b) * X(1, idx));  % set Kmax for directions with positive diffusivity
idx = find(K > Kmax);                       % find K's that exceed Kmax
K(idx) = Kmax(idx);                         % clip K's that exceed Kmax

% estimate Kmean
Kmean = mean(K);

% impose Kmin_final and Kmax_final constraints
Kmean(Kmean > Kmax_final) = Kmax_final;
Kmean(Kmean < Kmin_final) = Kmin_final;




%--------------------------------------------------------------------------
% core function to estimate mean diffusivity without estimating tensors
%--------------------------------------------------------------------------

function [Dax, Drad, Davg, FA, dt, L, V, nexcl, fit_err] = dte_core_dirfit(x, n, b, A, Apinv, Sbeps, dti_method)

% dte_core_dirfit   Estimate mean diffusivity for a single voxel
%
% Syntax
%
%   [Dax, Drad, Davg, FA, DT, L, V] = dte_core_dirfit(x, ndir, b, A, Apinv, Sbeps)
%
% Inputs
%
%   x           (ndir+1)-by-1 vector of signal values
%               The first row of x contains the b = 0 signal values, while the next
%               ndir rows contain the signal values for the single b used for
%               DTI computations
%
%   ndir        Number of gradient directions for each b value
%
%   b           Vector of b-values for DTI computations
%
%   A           coefficients matrix
%
%   Apinv       pseudoinverse of A
%
%   Sbeps       Smallest acceptable diffusion-weighted signal value
%               Signal values smaller than Sbeps are set to Sbeps
%
% Outputs
%
%   Dax     Dummy output (set to zero)
%
%   Drad    Dummy output (set to zero)
%
%   Davg    Mean diffusivity as average of DT eigenvalues (mm^2/s units)
%
%   FA      Dummy output (set to zero)
%
%   DT      Dummy output (set to zero)
%
%   L       Dummy output (set to zero)
%   
%   V       Dummy output (set to zero)
%
%   nexcl   Dummy output (set to zero)

% Author: Ali Tabesh
% Last modified: 12/29/10

% Sb filter

Dax     = 0;
Drad    = 0;
FA      = 0;
dt      = zeros(6, 1);
L       = [0, 0, 0]';
V       = zeros(9, 1);
nexcl   = 0;
fit_err = 0;

% skip the voxel if b0 voxel is less than 0

if x(1) <= 0
    Davg = 0;
    return
end

x(x < Sbeps) = Sbeps;    % set voxel values smaller than Sbeps to Sbeps

B = log(x(2:end)) - log(x(1));
B = reshape(B, [n, length(b)])';

% compute directional diffusivities

if dti_method.linear_weighting == 0
    D = Apinv * B;
elseif dti_method.linear_weighting == 1        % weighting according to signal magnitude
    w = reshape(x(2:end), [n length(b)])';
    B = w .* B;
    D = zeros(1, n); %EM
    for i = 1:n
        D(i) = mlsei(diag(w(:, i)) * A, B(:, i));
    end
else
    error('Invalid DTI weighting option!')
end

% impose non-negativity

D(D < 0) = 0;
    
% mean diffusivity

Davg = mean(D);             % average directional diffusivities



%--------------------------------------------------------------------------
% run constrained linear least-squares with iterative addition of
% constraints
%--------------------------------------------------------------------------

function t = iterative_lsqlin(t, DKG, B, C, c, C_eq, c_eq, b, g, ndir, Kmin, NKmax)

% t = iterative_lsqlin(t, DKG, B, C, c, C_eq, c_eq)
% Author: Ali Tabesh
% Last modified: 07/30/12
% maximum number of iterations

max_iter = 10;

% repeat addition of constraints until no constraint is violated or maximum
% number of iterations is exceeded

counter = 0;    % iteration counter

while any((C * t) > c) && counter < max_iter    % default tolerance in lsqlin is 100*eps
    t = mlsei(DKG, B, C, c, C_eq, c_eq);
    [g, ndir, C, c] = augment_constraints_mex(t, g, b, ndir, Kmin, NKmax);
    counter = counter + 1;
end



%--------------------------------------------------------------------------
% apply selective median filter to parametric maps
%--------------------------------------------------------------------------

function img_out = medianfilter(img_in, img_ref, img_viol, T)

w = 1;      % window size (w = 1 corresponds to 3x3x3 neighborhood)

dims = size(img_in);

% set tensor_flag if image is 4D
if length(dims) == 4
    tensor_flag = 1;
%GRG (added new case for tensors)------------------------------------------
elseif length(dims) == 2;   %Tensor input
    tensor_flag = 2; 
%--------------------------------------------------------------------------
else                        %3D image volume
    tensor_flag = 0;
end

img_out = img_in;

idx_viol = find(img_viol > T);

for i = 1:length(idx_viol)

    [I, J, K] = ind2sub(size(img_viol), idx_viol(i));    
    
    % take a 3D patch about the bad voxel
    patch_ref = img_ref(max(1, I - w):min(I + w, end), max(1, J - w):min(J + w, end), max(1, K - w):min(K + w, end));
    
    % take the corresponding violations patch
    patch_viol = img_viol(max(1, I - w):min(I + w, end), max(1, J - w):min(J + w, end), max(1, K - w):min(K + w, end));
    
    % identify good voxels (with fewer violations than T)
    idx = find(patch_viol <= T);
    [sorted, sortedidx] = sort(patch_ref(idx));
    
    if ~isempty(sorted)     % if empty (no good voxels in neighobrhood), leave voxel the way it is
        if length(sorted) == 2 * ceil(length(sorted) / 2)           % sorted has even length
           idx_median_voxel = idx(sortedidx((end / 2):(end/2 + 1)));
            
            %GRG-----------------------------------------------------------
            %Chose one closet to the mean
            ref_mean = mean(patch_ref(idx)); 
            ref_diff = abs(ref_mean-patch_ref(idx_median_voxel)); 
            idx_median_voxel = idx_median_voxel(find(ref_diff==min(ref_diff),1));
            %--------------------------------------------------------------
            
        else                                                        % sorted has odd length
            idx_median_voxel = idx(sortedidx((end + 1) / 2));
        end
        if tensor_flag == 1
            for j = 1:dims(4)
                patch_in = squeeze(img_in(max(1, I - w):min(I + w, end), max(1, J - w):min(J + w, end), max(1, K - w):min(K + w, end), j));
                img_out(I, J, K, j) = mean(patch_in(idx_median_voxel));
            end
            
        %GRG---------------------------------------------------------------
        elseif tensor_flag == 2
            r = max(1, I - w):min(I + w, size(img_viol,1)); 
            c = max(1, J - w):min(J + w, size(img_viol,2)); 
            p = max(1, K - w):min(K + w, size(img_viol,3));
            
            [C R P] = meshgrid(c,r,p); 
            
            patch_in_idx = (P(:)-1)*size(img_viol,1)*size(img_viol,2)+(C(:)-1)*size(img_viol,1)+R(:);       
            patch_in = img_in(:,patch_in_idx);
            img_out(:,idx_viol(i)) = patch_in(:,idx_median_voxel); 
        %------------------------------------------------------------------
        
        else
            patch_in = img_in(max(1, I - w):min(I + w, end), max(1, J - w):min(J + w, end), max(1, K - w):min(K + w, end));
            img_out(I, J, K) = mean(patch_in(idx_median_voxel));
        end
    end
end



%--------------------------------------------------------------------------
% generate a gradient set with n uniformly arranged vectors
%--------------------------------------------------------------------------

function [DG, KG] = generate_wm_DG_KG(n)

step = 1 / sqrt(n);

u = 0:step:(1 - step);

[theta, phi] = meshgrid(2 * pi * u, acos(1 - 2 * u));

theta = theta(:);
phi   = phi(:);

g = [cos(theta) .* sin(phi), sin(theta) .* sin(phi), cos(phi)];

DG = [g(:, 1) .^ 2 g(:, 2) .^ 2 g(:, 3) .^ 2, ...
    2 * g(:, 1) .* g(:, 2) 2*g(:, 1) .* g(:, 3), 2 * g(:, 2) .* g(:, 3)];

KG = [g(:, 1) .^ 4 g(:, 2) .^ 4 g(:, 3) .^ 4 ...
    4 * (g(:, 1) .^ 3) .* g(:, 2), ...
    4 * (g(:, 1) .^ 3) .* g(:, 3), ...
    4 * g(:, 1) .* (g(:, 2) .^ 3), ...
    4 * g(:, 1) .* (g(:, 3) .^ 3), ...
    4 * (g(:, 2) .^ 3) .* g(:, 3), ...
    4 * g(:, 2) .* (g(:, 3) .^ 3), ...
    6 * (g(:, 1) .^ 2) .* (g(:, 2) .^ 2), ...
    6 * (g(:, 1) .^ 2) .* (g(:, 3) .^ 2), ...
    6 * (g(:, 2) .^ 2) .* (g(:, 3) .^ 2) ...
    12 * (g(:, 1) .^ 2) .* g(:, 2) .* g(:, 3), ...
    12 * g(:, 1) .* (g(:, 2) .^ 2) .* g(:, 3), ...
    12 * g(:, 1) .* g(:, 2) .* (g(:, 3) .^ 2)];




%--------------------------------------------------------------------------
% form the coefficients matrices
%--------------------------------------------------------------------------

function [DG, KG, DKG] = coeff_mat(gmat, ndir, b)

% coeff_mat  Form coefficient matrices for dke_core
%
% Syntax
%
%   [DG KG DKG] = coeff_mat(gcell, ndir, b)
%
% Inputs
%
%   gcell   1-by-nb cell array containing gradient directions for each
%           b-value, where nb is the number of nonzero b-values
%
%   ndir    1-by-nb vector containing the numbers of gradient directions for
%           b-values
%
%   b       1-by-nb vector of nonzero b-values (in s/mm^2 units)
%
% Outputs
%
%   DG      diffusion matrix (for use in constraints matrix)
%
%   KG      kurtosis matrix (for use in constraints matrix)
%
%   DKG     complete coefficients matrix (actual coefficients matrix)

nb          = length(b);
ndir_total  = sum(ndir(1:nb));
DG          = zeros(ndir_total, 6);
DKG         = zeros(ndir_total, 21);
KG          = zeros(ndir_total, 15);

for ib = 1:nb
    sdir = sum(ndir(1:(ib - 1)));    % starting row index for the current b value
    g = gmat(sdir + (1:ndir(ib)), :);
    
    for idir = 1:ndir(ib)
        % [D11 D22 D33 D12 D13 D23]
        DG(sdir + idir, :) = [g(idir, 1) ^ 2, ...
            g(idir, 2) ^ 2, ...
            g(idir, 3) ^ 2, ...
            2 * g(idir, 1) * g(idir, 2), ...
            2 * g(idir, 1) * g(idir, 3), ...
            2 * g(idir, 2) * g(idir, 3)];
    
        % [W1111 W2222 W3333 W1112 W1113 W1222 W1333 W2223 W2333 W1122 W1133 W2233 W1123 W1223 W1233]
        KG(sdir + idir, :) = [g(idir, 1) ^ 4, ...
            g(idir, 2) ^ 4, g(idir, 3) ^ 4, ...
            4 * (g(idir, 1) ^ 3) * g(idir, 2), ...
            4 * (g(idir, 1) ^ 3) * g(idir, 3), ...
            4 * g(idir, 1) * (g(idir, 2) ^ 3), ...
            4 * g(idir, 1) * (g(idir, 3) ^ 3), ...
            4 * (g(idir, 2) ^ 3) * g(idir, 3), ...
            4 * g(idir, 2) * (g(idir, 3) ^ 3), ...
            6 * (g(idir, 1) ^ 2) * (g(idir, 2) ^ 2), ...
            6 * (g(idir, 1) ^ 2) * (g(idir, 3) ^ 2), ... 
            6 * (g(idir, 2) ^ 2) * (g(idir, 3) ^ 2), ...
            12 * (g(idir, 1) ^ 2) * g(idir, 2) * g(idir, 3), ...
            12 * g(idir, 1) * (g(idir, 2) ^ 2) * g(idir, 3), ...
            12 * g(idir, 1) * g(idir, 2) * (g(idir, 3) ^ 2)];
        
        DKG(sdir + idir, :) = [-b(ib) * DG(sdir + idir, :), (b(ib) ^ 2 / b(1) / 6) * KG(sdir + idir, :)];
    end         % for idir
end         % for ib


