function fg = trk2afq(FT_struc)
%Convert from .trk file to the fiber group structure used by AFQ using the
%information defined in the FT structure defined in the dke_ft modeule 
%
%Author: Russell Glenn
%Medical University of South Carolina

fn_trk = fullfile(FT_struc.outdir,FT_struc.name); 

[hdr trks] = trk_read(fn_trk); 

fibers = cell(length(trks),1); 
for i = 1:length(trks); 
    fi = hdr.vox_to_ras*[bsxfun(@rdivide,double(trks(i).matrix),hdr.voxel_size)';ones(1,trks(i).nPoints)]; 
    fibers{i} = fi(1:3,:); 
end

% params.faThresh = FT_struc.fa_threshold;
% params.lengthThreshMm = [FT_struc.trk_length 1000];
% params.stepSizeMm = FT_struc.step_size; 

params = {'faThresh',FT_struc.fa_threshold,'lengthThreshMm',[FT_struc.trk_length 500],...
    'stepSizeMm',FT_struc.step_size}; 

seeds = hdr.vox_to_ras*[FT_struc.SEED';ones(1,size(FT_struc.SEED,1))];

fg.name = 'wholeBrain'; 
fg.colorRgb = [20 90 200]; 
fg.thickness = -0.5;
fg.visible = 1; 
fg.seeds = seeds(1:3,:)';   
fg.seedRadius = 0; 
fg.seedVoxelOffsets = 'random'; 
fg.params = params; 
fg.fibers = fibers; 
fg.query_id = -1; 

% %%
% delete(gca)
% hold on 
% 
% idx = randperm(20,20); 
% 
% for i = [15 16 19 20]
%     
%     xi = x(i).fibers; 
%     
%     for j = 1:length(xi)        
%         f = xi{j}; 
%         plot3(f(1,:),f(2,:),f(3,:),'color',C(idx(i),:))
%     end
% end
%     