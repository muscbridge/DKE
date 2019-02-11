function matlab2trk(fn_trk, tracts,A,hdr)
% function matlab2trk(fn_trk, tracts,A,hdr,iop)
% matlab2trk(fn_trk, tracts, hdr)matlab2trk
%
% fn_trk    string containing the output track file name
% tracts    cell array where each cell is an n-by-3 matrix containing the coordinates of a track
% A         vox_to_ras transformation (LPS > RAS)
% iop       image orientation patient (from invert_img)

%Adapted from trk_write----------------------------------------------

%TRK_WRITE - Write TrackVis .trk files
%
% Syntax: trk_write(header,tracks,savePath)
%
% Inputs:
%    header   - Header information for .trk file [struc]
%    tracks   - Track data struc array [1 x nTracks]
%      nPoints  - # of points in each track
%      matrix   - XYZ coordinates (in mm) and associated scalars [nPoints x 3+nScalars]
%      props    - Properties of the whole tract
%    savePath - Path where .trk file will be saved [char]
%
% Output files:
%    Saves .trk file to disk at location given by 'savePath'.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Apr 2010


tracts_out = cell2struct(tracts, {'matrix'}, 2);

for i = 1:length(tracts_out)

    tracts_out(i).matrix(:,1) = tracts_out(i).matrix(:, 1) * hdr.pixdim(2);
    tracts_out(i).matrix(:,2) = tracts_out(i).matrix(:, 2) * hdr.pixdim(3);
    tracts_out(i).matrix(:,3) = tracts_out(i).matrix(:, 3) * hdr.pixdim(4);

    tracts_out(i).nPoints = size(tracts_out(i).matrix, 1);

end

%BUILD HEADER


header.dim = hdr.dim(2:4); 
header.hdr_size = 1000;
header.id_string = ['TRACK' char(0)];
header.invert_x = 0;
header.invert_y = 0;
header.invert_z = 0;
header.n_count = length(tracts_out);
header.n_properties = 0;
header.n_scalars = 0;
header.origin = [0 0 0];
header.pad1 = char([0 0]);
header.pad2 = ['LPS' char(0)];
header.property_name = char(zeros(10,20));
header.reserved = char(zeros(444,1));
header.scalar_name = char(zeros(10,20));
header.swap_xy = 0;
header.swap_yz = 0;
header.swap_zx = 0;
header.version = 2;
header.voxel_order = ['LPS' char(0)];
header.voxel_size = hdr.pixdim(2:4);
header.vox_to_ras = A;
% header.image_orientation_patient = iop; 
% header.image_orientation_patient = [1 0 0 0 1 0]; %Tracts were calculated from image volume in LPS
header.image_orientation_patient = [[-header.vox_to_ras(1:2, 1); header.vox_to_ras(3, 1)]/hdr.pixdim(2); [-header.vox_to_ras(1:2, 2); header.vox_to_ras(3, 2)]/hdr.pixdim(3)]';
% header.image_orientation_patient = [A(1:3,1)'./hdr.pixdim(2) A(1:3,2)'./hdr.pixdim(3)]; 

%WRITE TRK FILE

fid = fopen(fn_trk, 'w');

fwrite(fid, header.id_string, '*char');
fwrite(fid, header.dim, 'short');
fwrite(fid, header.voxel_size, 'float');
fwrite(fid, header.origin, 'float');
fwrite(fid, header.n_scalars , 'short');
fwrite(fid, header.scalar_name', '*char');
fwrite(fid, header.n_properties, 'short');
fwrite(fid, header.property_name', '*char');
fwrite(fid, header.vox_to_ras', 'float');
% fwrite(fid, header.vox_to_ras2', 'float');
fwrite(fid, header.reserved, '*char');
fwrite(fid, header.voxel_order, '*char');
fwrite(fid, header.pad2, '*char');
fwrite(fid, header.image_orientation_patient, 'float');
fwrite(fid, header.pad1, '*char');
fwrite(fid, header.invert_x, 'uchar');
fwrite(fid, header.invert_y, 'uchar');
fwrite(fid, header.invert_z, 'uchar');
fwrite(fid, header.swap_xy, 'uchar');
fwrite(fid, header.swap_yz, 'uchar');
fwrite(fid, header.swap_zx, 'uchar');
fwrite(fid, header.n_count, 'int');
fwrite(fid, header.version, 'int');
fwrite(fid, header.hdr_size, 'int');

% Check orientation
[tmp ix] = max(abs(header.image_orientation_patient(1:3)));
[tmp iy] = max(abs(header.image_orientation_patient(4:6)));
iz = 1:3;
iz([ix iy]) = [];

% Write body
for iTrk = 1:header.n_count
    % Modify orientation back to LPS for display in TrackVis
    header.dim        = header.dim([ix iy iz]);
    header.voxel_size = header.voxel_size([ix iy iz]);
    coords = tracts_out(iTrk).matrix(:,1:3);
    coords = coords(:,[ix iy iz]);
    if header.image_orientation_patient(ix) < 0
        coords(:,ix) = header.dim(ix)*header.voxel_size(ix) - coords(:,ix);
    end
    if header.image_orientation_patient(3+iy) < 0
        coords(:,iy) = header.dim(iy)*header.voxel_size(iy) - coords(:,iy);
    end
    tracts_out(iTrk).matrix(:,1:3) = coords;
    
    fwrite(fid, tracts_out(iTrk).nPoints, 'int');
    fwrite(fid, tracts_out(iTrk).matrix', 'float');
    if header.n_properties
        fwrite(fid, tracts_out(iTrk).props, 'float');
    end
end

fclose(fid);
