function [KFA,MeanKT] = ComputeKFA(KT,x,Kmax_final,Kmin_final)
%
% computes the KFA given a 15-vector of kurtosis tensor values
% Ali ordering
%
% Author: Mark Van Horn
% Version: 2.6.0
% Last modified: 11/25/14 by EM 

if x(1) <= 0    % skip the voxel if b0 voxel is less than 0
  KFA = 0;
  MeanKT= 0;
    return
end

W1111 = KT(1);
W2222 = KT(2);
W3333 = KT(3);
W1112 = KT(4);
W1113 = KT(5);
W1222 = KT(6);
W1333 = KT(7);
W2223 = KT(8);
W2333 = KT(9);
W1122 = KT(10);
W1133 = KT(11);
W2233 = KT(12);
W1123 = KT(13);
W1223 = KT(14);
W1233 = KT(15);

W_F = sqrt(W1111 ^ 2 + W2222 ^ 2 + W3333 ^ 2 + 6 * W1122 ^ 2 + 6 * W1133 ^ 2 + 6 * W2233 ^ 2 + ...
    4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + 4 * W1333 ^ 2 + ...
    4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);

Wbar=ComputeMeanKT(KT);

if W_F < 1e-3,
    KFA = 0;
else
    W_diff_F = sqrt((W1111 - Wbar) ^ 2 + (W2222 - Wbar) ^ 2 + (W3333 - Wbar) ^ 2 + ...
        6 * (W1122 - Wbar / 3) ^ 2 + 6 * (W1133 - Wbar / 3) ^ 2 + 6 * (W2233 - Wbar / 3) ^ 2 + ...
        4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + ...
        4 * W1333 ^ 2 + 4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
    
    KFA = W_diff_F / W_F;
end

MeanKT=Wbar;
MeanKT(MeanKT > Kmax_final) = Kmax_final;
MeanKT(MeanKT < Kmin_final) = Kmin_final;

end


% function KFA = ComputeKFA(KT)
% %
% % computes the KFA given a 15-vector of kurtosis tensor values
% %
% % Author: Mark Van Horn
% % Version: 2.6.0
% % Last modified: 08/18/14
% 
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
% 
% KTfull = MakeFullKT(KT);
% normKT = norm(KTfull(:));
% if normKT > 1e-6,
%     % mean of the kurtosis tensor
%     meanKT = (KTfull(1, 1, 1, 1) + KTfull(2, 2, 2, 2) + KTfull(3, 3, 3, 3) + ...
%         2 * (KTfull(1, 1, 2, 2) + KTfull(1, 1, 3, 3) + KTfull(2, 2, 3, 3))) / 5;
%     
%     % isotropic component of the kurtosis tensor
%     KTI = KTfull - meanKT * I4;
%     KFA = norm(KTI(:)) / normKT;
% else
%     KFA = 0;
% end
% 
