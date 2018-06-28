function make_4D_nii(varargin)
%Create 4D image volume from 3D nifti images AND conserve header info, such
%as the hdr.mat vox_to_ras transformation
%
%EXAMPLES: 
%
%    MAKE_4D_NII with no input arguments will prompt user to select files
%
%    MAKE_4D_NII(D,FILES) where D is a directory and FILES is a cell
%       containing file names for nifti files within D. 
%
%    MAKE_4D_NII(hdr,img) where hdr is a nifti header (spm_vol) and img is 
%       a 4D image
%
%    MAKE_4D_NII(~,~,name) specifies output name
%
%This function was made to work with Jimmy Shen's make_nii and save_nii
%functions
%
%Author: Russell Glenn
%Medical University of South Carolina

%DECIDE WHAT THE USER WANTS TO DO----------------------------------------

read_images = 0; 
if nargin == 0||nargin==1
    [files d] = uigetfile('*','Select Nifti Files to Convert','MultiSelect','on'); 
    if d~=0
        read_images = 1; 
    else
        error('Nothing to do')
    end
    
    name = []; 
    if nargin==1; name = varargin{1}; end
    if ~ischar(name);
        [d f e] = fileparts(fullfile(d,files{1})); 
        name = [f '_4D.nii']; 
    end
    
    if length(files)<2
	error('Nothing to do')
    end
    
else
    v1 = varargin{1}; 
    v2 = varargin{2};  
    
    if nargin>2
        name = varargin{3}; 
        if ~ischar(name); name = 'img_4D.nii'; end
        if ~strcmp(name(end-3:end),'.nii'); name = [name '.nii']; end
    end
        
    if ischar(v1)&&isdir(v1)
        read_images = 1; 
        d = v1; 
        if iscell(v2)
            files = v2; 
        else
            error('files must be a cell array of strings')
        end
    else
        if size(v2,4)>1
            hdr = v1(1); img = v2; 
            if size(hdr,1)==1; hdr = repmat(hdr,size(img,4),1); end
%             if size(hdr,1)~=size(img,4); error('hdr and img size are not compatible'); end
        else
            error('img must be a 4D image volume')
        end 
    end        
end

%FORMAT NEW HEADER---------------------------------------------------------

if read_images
    hdr = spm_vol(fullfile(d,files{1}));  
    nii0 = read_nii(fullfile(d,files{1})); 

%     vox = sqrt(sum(hdr.mat(1:3,1:3).^2)); 
%     A = bsxfun(@rdivide,hdr.mat(1:3,1:3),vox); 
%     vo = hdr.mat(1:3,4); 
    
    img = zeros(hdr.dim(1),hdr.dim(2),hdr.dim(3),length(files));
    for i = 1:length(files)
        img(:,:,:,i) = spm_read_vols(spm_vol(fullfile(d,files{i})));
    end
    
else
    d = fileparts(hdr(1).fname); 
    nii0 = read_nii(hdr(1).fname); 
end
    vox = sqrt(sum(hdr(1).mat(1:3,1:3).^2)); 
    nii = make_nii(img); 
    
%     srow = hdr.mat(1:3,4)+A*vox'; 
% 
%     nii.hdr.dime.pixdim = [-1 vox 0 0 0 0];
%     nii.hdr.dime.vox_offset = hdr.private.dat.offset;
%     nii.hdr.dime.scl_slope = hdr.private.dat.scl_slope;
%     nii.hdr.dime.scl_inter = hdr.private.dat.scl_inter; 
% 
%     nii.hdr.hist.descrip = hdr.descrip;
%     nii.hdr.hist.qform_code = 2; 
%     nii.hdr.hist.sform_code = 2; 
%     nii.hdr.hist.qoffset_x = srow(1); 
%     nii.hdr.hist.qoffset_y = srow(2); 
%     nii.hdr.hist.qoffset_z = srow(3); 
%     nii.hdr.hist.srow_x = [hdr.mat(1,1:3) srow(1)]; 
%     nii.hdr.hist.srow_y = [hdr.mat(2,1:3) srow(2)]; 
%     nii.hdr.hist.srow_z = [hdr.mat(3,1:3) srow(3)]; 
%     nii.hdr.hist.magic = 'n+1 '; 
    
    nii.hdr.hk.sizeof_hdr = nii0.sizeof_hdr;
    nii.hdr.hk.db_name = nii0.db_name;
    nii.hdr.hk.extents = nii0.extents;
    nii.hdr.hk.session_error = nii0.session_error;
    nii.hdr.hk.regular = nii0.regular;
    nii.hdr.hk.dim_info = nii0.dim_info;

%     nii.hdr.dime.dim = [4 nii0.dim(2:end)]; 
    nii.hdr.dime.intent_p1 = nii0.intent_p1;
    nii.hdr.dime.intent_p2 = nii0.intent_p2; 
    nii.hdr.dime.intent_p3 = nii0.intent_p3;
    nii.hdr.dime.intent_code = nii0.intent_code;
    nii.hdr.dime.slice_start = nii0.slice_start;
%     nii.hdr.dime.pixdim = [1 vox 0 0 0 0];
    nii.hdr.dime.pixdim = nii0.pixdim; 
    nii.hdr.dime.vox_offset = nii0.vox_offset;
    nii.hdr.dime.scl_slope = nii0.scl_slope;
    nii.hdr.dime.scl_inter = nii0.scl_inter;
    nii.hdr.dime.slice_end = nii0.slice_end;
    nii.hdr.dime.slice_code = nii0.slice_code;
    nii.hdr.dime.xyzt_units = nii0.xyzt_units;
    nii.hdr.dime.cal_max = nii0.cal_max;
    nii.hdr.dime.cal_min = nii0.cal_min;
    nii.hdr.dime.slice_duration = nii0.slice_duration;
    nii.hdr.dime.toffset = nii0.toffset;
    nii.hdr.dime.glmax = nii0.glmax;
    nii.hdr.dime.glmin = nii0.glmin;
    
    nii.hdr.hist.descrip = nii0.descrip;
    nii.hdr.hist.aux_file = nii0.aux_file;
    nii.hdr.hist.qform_code = nii0.qform_code;
    nii.hdr.hist.sform_code = nii0.sform_code;
    nii.hdr.hist.quatern_b = nii0.quatern_b;
    nii.hdr.hist.quatern_c = nii0.quatern_c;
    nii.hdr.hist.quatern_d = nii0.quatern_d;
    nii.hdr.hist.qoffset_x = nii0.quatern_x;
    nii.hdr.hist.qoffset_y = nii0.quatern_y;
    nii.hdr.hist.qoffset_z = nii0.quatern_z;
    nii.hdr.hist.srow_x = nii0.srow_x;
    nii.hdr.hist.srow_y = nii0.srow_y;
    nii.hdr.hist.srow_z = nii0.srow_z;
    nii.hdr.hist.intent_name = nii0.intent_name;
    nii.hdr.hist.magic = nii0.magic;
    
    

    save_nii(nii,fullfile(d,name))

end