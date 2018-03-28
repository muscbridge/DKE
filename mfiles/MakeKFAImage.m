
% KFA-----------------------------------------------------------------------
[hdr, img] = read_nii('D:\Projects\DKE\DKE_MATLAB\testBed\4Drdki_0.nii');
hdr.fn = 'kfatest1.nii';
hdr.dim(1) = 3;
hdr.dim(5) = 1;

% imageSize = [82, 82, 45];
imageSize = hdr.dim(2:4);

nvoxels = prod(imageSize);

kfa = zeros(nvoxels, 1);

load('KT.mat');


% I4 = zeros(3, 3, 3, 3);
% for i = 1:3
%     for j = 1:3
%         for k = 1:3
%             for l = 1:3
%                 I4(i, j, k, l) = (((i == j) && (k == l)) + ((i == k) && (j == l)) + ((i == l) && (j == k))) / 3;
%             end
%         end
%     end
% end

% initialize 15 elements
% [W1111 W2222 W3333 W1112 W1113 W1222 W1333 W2223 W2333 W1122 W1133 W2233 W1123 W1223 W1233]

% kurtosis tensor
% KTfull = zeros(3, 3, 3, 3);

if feature('numCores') > 1,
    if matlabpool('size') <= 0
        evalc('matlabpool open');
    end
    tic
    parfor n = 1:nvoxels,
        % form full kurtosis tensor
%         KTfull = MakeFullKT(KT(:, n));
%         normKT = norm(KTfull(:));
%         if normKT > 1e-6,
%             % mean of the kurtosis tensor
%             meanKT = (KTfull(1, 1, 1, 1) + KTfull(2, 2, 2, 2) + KTfull(3, 3, 3, 3) + 2 * KTfull(1, 1, 2, 2) + ...
%                 2 * KTfull(1, 1, 3, 3) + 2 * KTfull(2, 2, 3, 3)) / 5;
%             
%             % Isotropic component of the kurtosis tensor
%             KTI    = KTfull - meanKT * I4;
%             kfa(n) = norm(KTI(:)) / normKT;
             %kfa(n) = ComputeKFAX(KT(:, n));
             kfa(n) = ComputeKFA(KT(:, n));
%         end
        
    end
    toc
    matlabpool('close');
end

kfa = reshape(kfa, imageSize(1), imageSize(2), imageSize(3));
write_nii(hdr, round(1000*kfa));




