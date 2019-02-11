function SEED = random_seed_FT(FT_struct)
%RANDOM_SEED_FT(FT_STRUCT) generates random seed points indicated by the
%tracking mask indicated by the fields in FT_struct, including the fa threshold 
%and additional tracking mask used. 
%
%The DKI FT scripts expect LPS orientation whereas nifti outputs LAS. Thus
%the second dimension should typically be inverted. 
%
%FT_Struct is defined in the optimize_DKI_dODF.m script and includes the
%following fields: 
%     tractography_flg = 1;   %Perform tractography
%     fa_threshold = 0.1;     %FA threshold
%     angle_threshold = 35;   %Angle threshold in degrees
%     trk_length = 20;        %Minimum tract length in mm
%     step_size = 1;          %Step size in mm (0 defaults to half of the voxel length)
%     seednum = 1E5;          %Number of random seed points in the tracking mask; 
%     trk_mask = '';          %Tracking mask to apply in addition to other criteria defined above
%     shift = 0.5;            %Shift applied to voxel coordinates in .trk file
%     permute_odf = [1 2 3];  %Permutation of [x,y,z] dimensions in case orientation of images is off
%     invert_odf = [1 -1 1];  %+/- 1 to flip the given dimension. The tractography script expects LPS orientaion
%     hdr                     %SPM header for fa.nii image

%GET FA IMAGE
if isempty(FT_struct.seed_mask)
    fa = permute(spm_read_vols(FT_struct.hdr),FT_struct.permute_img); 
    inv_dim = find(FT_struct.invert_img==-1); 
    for i = inv_dim; fa = flipdim(fa,i); end

    %GET BRAIN MASK
    if ~isempty(FT_struct.trk_mask)
        mask = permute(spm_read_vols(spm_vol(FT_struct.trk_mask)),FT_struct.permute_img); 
        for i = inv_dim; mask = flipdim(mask,i); end
        mask = logical(mask)&fa>FT_struct.fa_threshold; 
    else
        mask = fa>FT_struct.fa_threshold;
    end
else
    mask = logical(permute(spm_read_vols(spm_vol(FT_struct.seed_mask)),FT_struct.permute_img)); 
    inv_dim = find(FT_struct.invert_img==-1); 
    for i = inv_dim; mask = flipdim(mask,i); end
end

%GENERATE RANDOM SEED POINTS IN BRAIN MASK---------------------------------
idx = find(mask);
[ROW COL PAGE] = ind2sub(FT_struct.hdr.dim(FT_struct.permute_img),idx); %RCP representation of brain mask

%Randomly choose voxels within the brain mask
rint = round(rand(FT_struct.seednum,1)*length(idx)+0.5); 

%Randomly distribute each random integer within the voxel (ie between i-0.5
%and i+0.5 for each dimension)
Ri = rand(FT_struct.seednum,1)-0.5;
Ci = rand(FT_struct.seednum,1)-0.5;
Pi = rand(FT_struct.seednum,1)-0.5;

SEED = [ROW(rint)+Ri, COL(rint)+Ci, PAGE(rint)+Pi]; %Seed points