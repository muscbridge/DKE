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

if options.gibbs_corr_flag == 0
    error('Not correcting for Gibbs artifact -- not supported yet')
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

b0_dir = fullfile(nifti_dir, 'B0');
mkdir(b0_dir);

eval(['!/Applications/MRIcroGL/dcm2niix -f %p ' basedir])

% move new NIfTI files to subdirectories of nifti_dir

% this part is needed when you have an additional separate b0 sequence
%     in=fullfile(basedir, '*B0*.nii');
%     copyfile(in, b0_dir);

dki1_dir = fullfile(nifti_dir, 'DKI1');
mkdir(dki1_dir);
file_in=dir(fullfile(basedir, '*DKI*.nii'));
in=fullfile(basedir, file_in(1).name);
out=fullfile(dki1_dir, 'dki.nii');
copyfile(in,out);

%--------------------------------------------------------------------------
% Split 4D nii into 3D nii
%--------------------------------------------------------------------------

% for separate b0
% Vdir = dir(fullfile(b0_dir, '*.nii'));
% V=fullfile(b0_dir, Vdir(1).name);
% Vo = spm_file_split(V, b0_dir);

Vdir = dir(fullfile(dki1_dir, '*.nii'));
V=fullfile(dki1_dir, Vdir(1).name);
Vo = spm_file_split(V,dki1_dir);

in=fullfile(dki1_dir, 'dki.nii');
out=fullfile(dki1_dir, '4D.nii');
movefile(in,out);

%--------------------------------------------------------------------------
% Denoise if denoise_flag = 1
%--------------------------------------------------------------------------

if options.denoise_flag == 1
    fprintf('Denoising data with dwidenoise (MRtrix)...  ')
    cd(dki1_dir);
    command=['/usr/local/mrtrix3/bin/dwidenoise 4D.nii 4D_DN.nii -noise noise.nii'];
    [status,cmdout] = system(command);
    if status == 0
        fprintf('done.\n')
    else
        fprintf('\nAn error occurred. Output of dwidenoise command:\n');
        cmdout
        error('Error running dwidenoise.')
    end
    fprintf('Calculating residuals with mrcalc (MRtrix)...  ')
    command=['/usr/local/mrtrix3/bin/mrcalc 4D.nii 4D_DN.nii -subtract res.nii'];
    [status,cmdout] = system(command);
    if status == 0
        fprintf('done.\n')
    else
        fprintf('\nAn error occurred. Output of mrcalc command:\n');
        cmdout
        error('Error running mrcalc.')
    end
else
    fprintf('Not denoising data\n')
end

%--------------------------------------------------------------------------
% Correct for Rician noise bias if rician_corr_flag = 1 and denoise_flag = 1
%--------------------------------------------------------------------------

if options.denoise_flag == 1
    if options.rician_corr_flag == 1
        fprintf('Correcting for Rician noise bias...  ')
        hdr_DN = spm_vol('4D_DN.nii');
        DN = spm_read_vols(hdr_DN);
        noise = spm_read_vols(spm_vol('noise.nii'));
        DN = sqrt(abs(DN.^2 - noise.^2));
        DN(isnan(DN)) = 0;
        make_4D_nii(hdr_DN, DN, '4D_DN_rician_corrected.nii');
        fprintf('done.\n')
    else
        fprintf('Not correcting for Rician noise bias\n')
    end
else
    if options.rician_corr_flag == 1
        fprintf('Not correcting for Rician noise bias (denoise_flag = 0)\n')
    end
end


%--------------------------------------------------------------------------
% Unring
%--------------------------------------------------------------------------

% If denoise_flag = 1 and rician_corr_flag = 1, use the denoised and Rician
% noise bias corrected volume 4D_DN_rician_corrected.nii.
% If denoise_flag = 1 and rician_corr_flag = 0, use the denoised
% volume 4D_DN.nii.
% Otherwise use 4D.nii.

if options.gibbs_corr_flag == 1
    fprintf('Correcting for Gibbs artifact...  ')
    if options.denoise_flag == 1
        if options.rician_corr_flag == 1
            DN=spm_read_vols(spm_vol(fullfile(dki1_dir,'4D_DN_rician_corrected.nii')));
        else
            DN=spm_read_vols(spm_vol(fullfile(dki1_dir,'4D_DN.nii')));
        end
    else
        DN=spm_read_vols(spm_vol(fullfile(dki1_dir,'4D.nii')));
    end

    list=dir(fullfile(dki1_dir,'*00*'));
    [dim1,dim2,dim3,dim4]=size(DN);
    parfor j=1:dim4
        img(:,:,:,j)=unring(DN(:,:,:,j));
        hdr=spm_vol(fullfile(dki1_dir,[list(j).name]));
        hdr.dt=[16 0];
        int=img(:,:,:,j);
        int(isnan(int))=0;
        spm_write_vol(hdr,int);
    end

    fprintf('done.\n')
else
    error('Not correcting for Gibbs artifact -- not supported yet')
end

%--------------------------------------------------------------------------
% Read DICOM headers from the DICOM files in basedir
%--------------------------------------------------------------------------

P = spm_select('fplist', basedir, '.*');  % fplist: list w/ full path
dicom_hdrs = spm_dicom_headers(P);

%--------------------------------------------------------------------------
% Get b values from the DICOM headers
% For every image,
%   If the sequence name from the header contains 'ep_b0', it's a b=0 image
%   Otherwise if the sequence name contains 'ep', get the b value
%     - b values are between 'b' and '#' (or 't') in the sequence name
%   Otherwise (if the sequence name does not contain 'ep'), show an error
%--------------------------------------------------------------------------

list=dir(fullfile(dki1_dir,'*00*.nii'));

for k=1:length(dicom_hdrs)
    hdr = dicom_hdrs{k};
    if ~isempty(strfind(hdr.SequenceName, 'ep_b0'))
        bval = '0';
    elseif ~isempty(strfind(hdr.SequenceName, 'ep'))
        bval_start = strfind(hdr.SequenceName, 'b') + 1;
        bval_end = strfind(hdr.SequenceName, '#') - 1;
        if isempty(bval_end)      % if '#' was not found, assume that DWI was acquired for only one direction (transversal)
            bval_end = strfind(hdr.SequenceName, 't') - 1;
            if isempty(bval_end)  % if 't' was not found, issue an error message
                error('Invalid gradient direction!')
            end
        end
        bval = hdr.SequenceName(bval_start:bval_end);
    else
        error('Invalid sequence name in DICOM header! Sequence name must contain ''ep''!')
    end
    list(k).bval = bval;
end

%--------------------------------------------------------------------------
% Rename NIfTI images (dki_vol#_bval.nii)
%--------------------------------------------------------------------------

for k=1:length(list)
    [list(k).newname] = strrep(list(k).name, '.nii', ['_b' list(k).bval '.nii']);
    movefile(fullfile(dki1_dir,list(k).name), fullfile(dki1_dir,list(k).newname));
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

