function EULER_DKE(FT_struc,C,name)
%This script uses DKE output files to perfom WM fiber tractography with the
%EULER algorithm, using either the DKI or DTI odf reconstructions
%
%INPUTS: Designed to work with DKE 2.6.0 and ODF_DKE.m
%   FT_STRUCT - FT_STRUCT defined with optimize_DKI_dODF.m containes the 
%   following fields:
%       tractography_flg  %Perform tractography
%       fa_threshold      %FA threshold
%       angle_threshold   %Angle threshold in degrees
%       trk_length        %Minimum tract length in mm
%       step_size         %Step size in mm (0 defaults to half of the voxel length)
%       seednum           %Number of random seed points in the tracking mask; 
%       trk_mask          %Tracking mask to apply in addition to other criteria defined above
%       shift             %Shift applied to voxel coordinates in .trk file
%       permute_odf       %Permutation of [x,y,z] dimensions in case orientation of images is off
%       invert_odf        %+/- 1 to flip the given dimension. The tractography script expects LPS orientaion
%       hdr               %SPM header for fa.nii image ~ expected to be in LAS orientaitno from SPM
%       %outdir           %Path to directory to output tracts (subj_dir)
%       %reset_memory     %Option if parfor causes memory issues

%
%DEFAULT VARIABLES: Change internally (for now)
%   ANGLE THRESHOLD - 35. Stop tracking tract turns by more than this
%   FA THRESHOLD - 0.1. Stop tracking if fa falls below
%   TRACK LENGTH - 20 (mm). Only consider something a track if it is longer than
%   20 mm in length
%
%OUPUTS: 
%   FT_NAME.trk - tracks in TrackVis format
%   FT_NAME.mat - tracks in matlab format - each element of the TRK cell
%   contains an Nx3 matrix with the xyz components of the track in voxel
%   space, ie [1 1 1] refers to the middle of the first voxel, [2 1 1]
%   refers to the middle of the second voxel, and so on. 
%
%ORIENATION NOTES: 
%   This program expects both image files and orientation estimates to be
%   in LPS coordinate system. This means dimension 1,2,and 3 of the image
%   volume corresponds to L,P,and S, respectively and the first,second, and
%   third elements of each orientation vector correspond to L,P,and S,
%   respectively. 

%Author: Russell Glenn
%Medical University of South Carolina
%04/08/2015


%Initialize things---------------------------------------------------------
%     aT = 35;                %Angle threshold in degrees
%     trk_length = 20;        %Length in mm required to count as track (affects trackDensity)
%     faT = 0.1;              %Fa thrshold
%     step_size = 1;           %Step Size in mm: Default 0 will re-calculate as 1/2 voxel dimensions
%     shift = 0.5;            %Shift required to take the coordinates from my 'voxel space' where
%                             %the middle of each voxel is located at it's index (with the bottom 
%                             %of the image volume at (0.5,0.5,0.5) to the space used by TrackVis,  
%                             %with the bottom of the image volume at (0,0,0). 
% 
%     %Change orientation vectors from whatever coordinate system they were
%     %calculated in to LPS-------
%     permute_odf = [1 2 3];  %permute x,y,z
%     invert_odf = [1 1 1];   %invert x,y,z
%---------------------------
    
%NOTE permute_odf and invert_odf are added because I don't know how other
%people will orient their gradient tables relative to how spm flips
%images around in the nifti format, so added these options. To see their
%effects, see viewODF_DKE.m

%NOTE: These are currently set in the script, but can be changed to be
%supplied as function arguments if you want or of you want to convert to
%some type of GUI... Eventually, I think it would be nice if these were
%contained within the dke_params file
%--------------------------------------------------------------------------

%GET FA IMAGE ~ THIS PROGRAM EXPECTS LPS
hdr = FT_struc.hdr; %Use as a template later
fa = permute(spm_read_vols(FT_struc.hdr),FT_struc.permute_img); 
C = permute(C,FT_struc.permute_img); 
inv_dim = find(FT_struc.invert_img==-1); 
for i = inv_dim; fa = flipdim(fa,i); C = flipdim(C,i); end

dim = hdr.dim(FT_struc.permute_img)'; %Image matrix size
vox = sqrt(sum(hdr.mat(FT_struc.permute_img,1:3).^2)); %Voxel Dimensions
vox = round(vox'*1E3)/1E3;  


%COMBINE WM TERMINATION CRITERIA TO MAKE WM TRACKING MASK
if ~isempty(FT_struc.trk_mask)
     mask = permute(spm_read_vols(spm_vol(FT_struc.trk_mask)),FT_struc.permute_img); 
     for i = inv_dim; mask = flipdim(mask,i); end
else mask = ones(size(fa)); 
end

mask = logical(mask)&fa>FT_struc.fa_threshold; 

%CONVERT C
for i = 1:numel(C); 
    ci = C{i}; 
    if~isempty(ci);             
        C{i} = bsxfun(@times,ci(FT_struc.permute_odf,:)',FT_struc.invert_odf)';
    end; 
end

TRK_C = cell(size(FT_struc.SEED,1),1);

numtrks = 0;
maxIter = 1000;

nbatch = 20; 
% if FT_struc.reset_memory; nbatch = 10; 
% else nbatch = 100;
% end


step = ceil(size(FT_struc.SEED,1) / nbatch);
fpb = fprintf('Processing %s Tractography... ',name); 

%NOTE: Coordinates are treated in 'voxel space' meaning if you round the
%current point it will tell you which voxel, however, to apply the correct
%step size, an appropriate scaling factor is applied, ie it's scaled to mm
%space, step is applied and then the voxel is determined by rescaling
%Also weird way for handling TRK / TRK_C is due to parfor loop restrictions

for k = 1:nbatch
% for k = 1; 
start = (k - 1) * step + 1;
stop = k * step;

    parfor si = start:min(size(FT_struc.SEED,1),stop)
%     for si = 1:size(FT_struct.SEED,1)
%     for si = 19
        %Initialize Things  
        r0 = FT_struc.SEED(si,:)'.*vox; %initial point
        ci = round(FT_struc.SEED(si,:));  %current index
        idx = (ci(3)-1)*prod(dim(1:2))+(ci(2)-1)*dim(1)+ci(1); %sub2ind
        V = C{idx}; %Vector directions; 
        trk_c = {}; %Hold tracts
        for i = 1:size(V,2); %Check all vectors
            v0 = V(:,i); 
            R = []; %Clear Points
            for dir = [-1 1]; %Check both directions 
                vi = dir*v0; %initialize vector direction
                ri = r0+FT_struc.step_size*vi; %Initial Point
                ci = round(ri./vox); 
                if ci(1)>=1&&ci(1)<=dim(1)&&ci(2)>=1&&ci(2)<=dim(2)&&ci(3)>=1&&ci(3)<=dim(3) %Inside image boundary
                    Ri = [r0 ri]; %Points for vector i along this direction
                    iter = 0; %Number of iterations completed
                    run = 1; %Run flag

                    while run&&iter<maxIter
                        ci = round(ri./vox); 
                        if ci(1)>=1&&ci(1)<=dim(1)&&ci(2)>=1&&ci(2)<=dim(2)&&ci(3)>=1&&ci(3)<=dim(3) %Inside image boundary
                            idx = (ci(3)-1)*prod(dim(1:2))+(ci(2)-1)*dim(1)+ci(1); %sub2ind
                            if fa(idx)>FT_struc.fa_threshold&&mask(idx)&&~isempty(C{idx}) %FA and extra mask criteria  

                                Vi = C{idx};
                                Vi = [Vi -Vi]; %Check all directions
                                deg = real(acosd(vi'*Vi)); %Round off error can make imaginary if vi=Vi
                                if min(deg)<FT_struc.angle_threshold
                                   vi = Vi(:,find(deg==min(deg),1)); 
                                   ri = ri+FT_struc.step_size*vi; 
                                    ci = round(ri./vox); 
                                    if ci(1)>=1&&ci(1)<=dim(1)&&ci(2)>=1&&ci(2)<=dim(2)&&ci(3)>=1&&ci(3)<=dim(3) %Inside image boundary
                                        Ri = [Ri ri]; 
                                    else run = 0; 
                                    end
                                else run = 0;                                                    
                                end
                            else run = 0;  
                            end
                        else run = 0;
                        end                    
                        iter = iter+1; 
                    end

                    if dir==-1&&size(Ri,2)>1; 
                        R = fliplr(Ri(:,2:end)); 
                    elseif size(Ri,2)>1
                        R = [R Ri]; 
                    end     
                end
            end
            
            len = sum(sqrt(sum((R(:,1:end-1)-R(:,2:end)).^2)));
            if len > FT_struc.trk_length
                trk_c{1,end+1} = bsxfun(@rdivide,R,vox)-FT_struc.shift; 
                numtrks = numtrks+1; 
            end            
        end
        TRK_C{si} = trk_c; 
    end

    fprintf(repmat('\b',1,fpb)); 
    fpb = fprintf('Processing %s Tractography: %d%% complete.',name,k * 100 / nbatch); 
        
end

fprintf(repmat('\b',1,fpb)); 
fpb = fprintf('Saving Data...');  

%BUILD OUTPUT STRUCTURES---------------------------------------------------

TRK = cell(numtrks,1); 
trknum = zeros(dim'); 
trk_i = 1; 
for i = 1:size(TRK_C,1)
        TRK_i = TRK_C{i}; 
        for j = 1:size(TRK_i,2)
            trk = TRK_i{j}; 
            TRK{trk_i} = trk';
            ci = ceil(trk);
            idx = (ci(3,:)-1).*prod(dim(1:2))+(ci(2,:)-1).*dim(1)+ci(1,:); %sub2ind
            trknum(idx) = trknum(idx)+1; 
            trk_i = trk_i+1; 
        end
end

%--------------------------------------------------------------------------

%SAVE AS TRK FILE: View with TrackVis

%SAVE TRK IN TRACK VIS FORMAT
    hdr_trk = struct('dim',[1 dim'],'pixdim',[1 vox']);
    
    %Get vox_to_ras (vox origin [0 0 0] and vox order LPS) and IOP (orientation of image volume)
    vox_to_ras = FT_struc.hdr.mat([FT_struc.permute_img 4],:)*[diag(FT_struc.invert_img) (FT_struc.invert_img==-1)'.*dim+FT_struc.shift*[1;1;1];0 0 0 1];  
%     iop = [FT_struc.invert_img(1) 0 0 0 FT_struc.invert_img(2) 0]; 
    
    matlab2trk(fullfile(FT_struc.outdir, [FT_struc.pre_name 'FT_' name FT_struc.post_name '.trk']),TRK,vox_to_ras,hdr_trk); 
%     matlab2trk(fullfile(FT_struc.outdir, [FT_struc.pre_name 'FT_' name FT_struc.post_name '.trk']),TRK,vox_to_ras,hdr_trk,iop); 
%     matlab2trk(fullfile(FT_struc.outdir, [FT_struc.pre_name 'FT_' name FT_struc.post_name '.trk']),TRK,hdr.mat,hdr_trk); 
    

%SAVE TRK IN MATLAB FORMAT    
    save(fullfile(FT_struc.outdir, [FT_struc.pre_name 'FT_' name FT_struc.post_name '.mat']),'TRK'); 

%SAVE TRACK DENSITY IN NIFTI FORMAT
    td = trknum./prod(vox); 
    
    %Undo inversions and permutations
    inv_dim = find(FT_struc.invert_img==-1); 
    for i = inv_dim; td = flipdim(td,i); end
    td = permute(td,FT_struc.permute_img); 

    hdr.fname = fullfile(FT_struc.outdir, [FT_struc.pre_name 'TrackDensity_' name FT_struc.post_name '.nii']);
    spm_write_vol(hdr,td); 

clear hdr fa C inv_dim dim vox mask ci TRK_C numtrks maxIter nbatch step fbp 
clear TRK start stop si ro ci idx V trk_c v0 R vi ri ci Ri run deg dir R len hdr_trk td

fprintf(repmat('\b',1,fpb));