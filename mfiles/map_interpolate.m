function map_interpolate(fn, vox, order)

% map_interpolate interpolate parametric maps to specified voxel size

% Author: Ali Tabesh
% Version: 2.6.0
% Last modified: 08/18/14 by mvh

% Portions of this code are based on John's Gems:
%
% http://www-personal.umich.edu/~nichols/JohnsGems.html

if nargin ~= 3
    fprintf('\n')
    fprintf('Usage: map_interpolate fn vox order\n')
    fprintf('fn  input filename\n')
    fprintf('vox  target voxel size in mm\n')
    fprintf('order  interpolation polynomial order; 0 is nearest neighbor and 1 is trilinear interpolation\n\n')
    return
end

warning('off')

vox   = str2double(vox);
order = str2double(order);

% interpolate parametric maps

hdr_in = spm_vol(fn);

d = hdr_in.dim(1:3);

% corners in voxel-space
c = [1,    1,    1,    1;
     1,    1,    d(3), 1;
     1,    d(2), 1,    1;
     1,    d(2), d(3), 1;
     d(1), 1,    1,    1;
     d(1), 1,    d(3), 1;
     d(1), d(2), 1,    1;
     d(1), d(2), d(3), 1]';

% corners in world-space
tc = hdr_in.mat(1:3, 1:4) * c;

% bounding box (world) min and max
mn = min(tc, [], 2)';
mx = max(tc, [], 2)';

mat = spm_matrix([mn, 0, 0, 0, vox, vox, vox]) * spm_matrix([-1, -1, -1]);

dims = ceil(mat \ [mx, 1]' - 0.1)';

% output image
hdr_out = hdr_in;

% [path, name, ext]  = fileparts(hdr_in.fname);
[path, ~, ~]     = fileparts(hdr_in.fname);
hdr_out.fname    = fullfile(path, 'tmp.nii');
hdr_out.dim(1:3) = dims(1:3);
hdr_out.mat      = mat;
hdr_out = spm_create_vol(hdr_out);

for i = 1:dims(3)
    M = inv(spm_matrix([0, 0, -i]) * inv(hdr_out.mat) * hdr_in.mat);
    img = spm_slice_vol(hdr_in, M, dims(1:2), order);
    spm_write_plane(hdr_out, img, i);
end

movefile(hdr_out.fname, hdr_in.fname, 'f')


