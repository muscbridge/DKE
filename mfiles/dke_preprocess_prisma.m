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

%--------------------------------------------------------------------------
% Convert DICOM images to NIfTI format
%--------------------------------------------------------------------------

b12root = [basedir '/intermediate_processing'];
mkdir(b12root);

eval(['!/Applications/MRIcroGL/dcm2niix -f %p ' basedir])

% move new NIfTI files to intermediate_processing/nifti

% this part is needed when you have an additional separate b0 sequence
%     mkdir([b12root '/nifti/B0'])
%     in=fullfile(basedir, '*B0*.nii');
%     out=fullfile(b12root, 'nifti/B0');
%     copyfile(in,out);

mkdir([b12root '/nifti/DKI1'])
file_in=dir(fullfile(basedir, '*DKI*.nii'));
in=fullfile(basedir, file_in(1).name);
out=fullfile(b12root, 'nifti/DKI1', 'dki.nii');
copyfile(in,out);

%--------------------------------------------------------------------------
% Split 4D nii into 3D nii
%--------------------------------------------------------------------------

% for separate b0
% Vdir = dir(fullfile(b12root, 'nifti/B0', '*.nii'));
% V=fullfile(b12root, 'nifti/B0', Vdir(1).name);
% Vo = spm_file_split(V,[b12root '/nifti/B0']);

Vdir = dir(fullfile(b12root, 'nifti/DKI1', '*.nii'));
V=fullfile(b12root, 'nifti/DKI1', Vdir(1).name);
Vo = spm_file_split(V,[b12root '/nifti/DKI1']);

in=fullfile(b12root, 'nifti/DKI1', 'dki.nii');
out=fullfile(b12root, 'nifti/DKI1', '4D.nii');
movefile(in,out);

mkdir([b12root '/dke']);

%--------------------------------------------------------------------------
% Denoise
%--------------------------------------------------------------------------

command=['/usr/local/mrtrix3/bin/dwidenoise ''' b12root '/nifti/DKI1/4D.nii'' ''' b12root '/nifti/DKI1/4D_DN.nii'' ' '-noise ''' b12root '/nifti/DKI1/noise.nii'''];
[status,cmdout] = system(command);

%--------------------------------------------------------------------------
% Unring
%--------------------------------------------------------------------------

DN=spm_read_vols(spm_vol(fullfile(b12root,'/nifti/DKI1','4D_DN.nii')));
list=dir(fullfile(b12root,'/nifti/DKI1','*00*'));
[dim1,dim2,dim3,dim4]=size(DN);
parfor j=1:dim4
    img(:,:,:,j)=unring(DN(:,:,:,j));
    hdr=spm_vol(fullfile(b12root,'nifti/DKI1',[list(j).name]));
    hdr.dt=[16 0];
    int=img(:,:,:,j);
    int(isnan(int))=0;
    spm_write_vol(hdr,int);
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

list=dir(fullfile(b12root,'/nifti/DKI1','*00*.nii'));

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
    movefile(fullfile(b12root,'/nifti/DKI1',list(k).name), fullfile(b12root,'/nifti/DKI1',list(k).newname));
end

% move b0s to folder 'intermediate_processing/nifti/b0'

clear list
list=dir(fullfile(b12root,'/nifti/DKI1','*b0*.nii'));
mkdir([b12root '/nifti/B0']);
for l=1:length(list)
    in=fullfile(b12root, 'nifti/DKI1',list(l).name);
    out=fullfile(b12root, '/nifti/B0', list(l).name);
    movefile(in,out);
end

%--------------------------------------------------------------------------
% Make gradient file
%--------------------------------------------------------------------------

cd(basedir);
name=dir('*.bvec');
A=importdata([name(1).name]);
B=A(:,any(A));
Gradient=B';
Gradient1=Gradient(1:(round(end/2)),:);
save(fullfile(b12root,'dke/gradient_dke.txt'),'Gradient1','-ASCII')

%--------------------------------------------------------------------------
% Coregister b0s to b0
%--------------------------------------------------------------------------

%ONLY USE WITH INTERLEAVED B0S
% dirb0=dir(fullfile(b12root,'nifti/B0/*.nii'));
% dirdki=dir(fullfile(b12root,'nifti/DKI1/*.nii'));
% 
% fprintf('Co-registering images...\n')
% fn_source = fullfile(b12root, 'nifti/B0', dirb0(1).name); % source file is the first b = 0 image in the series returned by the operating system
% 
%             fn_target = fullfile(b12root, 'nifti/DKI1', dirdki(1).name);
%             M=coregister(fn_target, fn_source, fullfile(b12root, 'nifti/B0'),'.nii');
% 
% delete(fullfile(b12root, 'nifti/DKI1', 'r*.nii'))
% movefile(fullfile(b12root, 'nifti/B0/r*.nii'),fullfile(b12root, 'nifti/B0_coreg')); 
% fprintf('Co-registration complete.\n')

%--------------------------------------------------------------------------
% Average b0's
%--------------------------------------------------------------------------
list = dir(fullfile(b12root,'nifti/B0/*.nii'));
hdr = spm_vol(fullfile(b12root,'nifti/B0',list(1).name));

imgavg = spm_read_vols(hdr);
for j = 2:length(list)
    hdr = spm_vol(fullfile(b12root,'nifti/B0',list(j).name));
    img = spm_read_vols(hdr);
    imgavg = imgavg + img;
end

imgavg = imgavg / (length(list));

mkdir(b12root, '/nifti/combined');
hdr.dt=[16 0];
hdr.fname = fullfile(b12root,'nifti/combined/b0_avg.nii');
imgavg(isnan(imgavg))=0;
spm_write_vol(hdr, imgavg);

%--------------------------------------------------------------------------
% Move DKI images and make 4D NIfTI images
%--------------------------------------------------------------------------
copyfile([b12root '/nifti/DKI1/*00*'], [b12root '/nifti/combined']);
files = dir([b12root '/nifti/combined/*.nii']);
make_4D_nii([b12root '/nifti/combined'],{files.name},'4D.nii');
movefile([b12root '/nifti/combined/4D.nii'],[b12root '/dke/4D.nii'])

img=spm_read_vols(spm_vol([b12root '/dke/4D.nii']));
img(isnan(img))=0;
make_4D_nii(spm_vol([b12root '/dke/4D.nii']),img,'4D.nii');

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

