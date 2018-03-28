
function [KFA2, KFA3, KFA4] = latestTest

KT = [ 0.192007426269198
   0.541929820430468
   0.643516593494627
  -0.159663710495804
   0.141490229284791
  -0.060816989967150
   0.239696609520885
  -0.026389619555339
   0.224699964963515
   0.100528956963001
   0.178067538569462
   0.241285405654159
  -0.001819606484770
  -0.092319698903758
   0.000683635970370];

% permute DKE ordering to Jens' ordering
indx2 = [1, 2, 3, 10, 11, 12, 4, 5, 6, 8, 7, 9, 13, 14, 15];

k = KT(indx2);

KFA2 = ComputeKFA2(k);  % Jens code and indexing
KFA3 = ComputeKFA3(k); % Russell code and DKE indexing
KFA4 = ComputeKFA4(KT); % Jens code and DKE indexing


function KFA = ComputeKFA2(KT)
%
% computes the KFA given a 15-vector of kurtosis tensor values
% Jens ordering
%
% Author: Mark Van Horn
% Version: 2.6.0
% Last modified: 08/18/14
W1111 = KT(1);
W2222 = KT(2);
W3333 = KT(3);
W1122 = KT(4);
W1133 = KT(5);
W2233 = KT(6);
W1112 = KT(7);
W1113 = KT(8);
W1222 = KT(9);
W2223 = KT(10);
W1333 = KT(11);
W2333 = KT(12);
W1123 = KT(13);
W1223 = KT(14);
W1233 = KT(15);

W_F = sqrt(W1111 ^ 2 + W2222 ^ 2 + W3333 ^ 2 + 6 * W1122 ^ 2 + 6 * W1133 ^ 2 + 6 * W2233 ^ 2 + ...
    4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + 4 * W1333 ^ 2 + ...
    4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
if W_F < 1e-3,
    KFA = 0;
else
    Wbar = (W1111 + W2222 + W3333 + 2 * (W1122 + W1133 + W2233)) / 5;
    
    W_diff_F = sqrt((W1111 - Wbar) ^ 2 + (W2222 - Wbar) ^ 2 + (W3333 - Wbar) ^ 2 + ...
        6 * (W1122 - Wbar / 3) ^ 2 + 6 * (W1133 - Wbar / 3) ^ 2 + 6 * (W2233 - Wbar / 3) ^ 2 + ...
        4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + ...
        4 * W1333 ^ 2 + 4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
    
    KFA = W_diff_F / W_F;
end


function KFA = ComputeKFA3(KT)
%
% Russell version
%

I4 = zeros(3, 3, 3, 3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                I4(i, j, k, l) = (((i == j) && (k == l)) + ((i == k) && (j == l)) + ((i == l) && (j == k))) / 3;
            end
        end
    end
end


KTfull(1, 1, 1, 1) = KT(1);
KTfull(2, 2, 2, 2) = KT(2);
KTfull(3, 3, 3, 3) = KT(3);
KTfull(1, 1, 2, 2) = KT(4);
KTfull(1, 1, 3, 3) = KT(5);
KTfull(2, 2, 3, 3) = KT(6);
KTfull(1, 1, 1, 2) = KT(7);
KTfull(1, 1, 1, 3) = KT(8);
KTfull(1, 2, 2, 2) = KT(9);
KTfull(2, 2, 2, 3) = KT(10);
KTfull(1, 3, 3, 3) = KT(11);
KTfull(2, 3, 3, 3) = KT(12);
KTfull(1, 1, 2, 3) = KT(13);
KTfull(1, 2, 2, 3) = KT(14);
KTfull(1, 2, 3, 3) = KT(15);

% propagate values to other elements
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                I = sort([i, j, k, l]);
                KTfull(i, j, k, l) = KTfull(I(1), I(2), I(3), I(4));
            end
        end
    end
end

normKT = norm(KTfull(:));

if normKT > 1e-3,
    % mean of the kurtosis tensor
    meanKT = (KTfull(1, 1, 1, 1) + KTfull(2, 2, 2, 2) + KTfull(3, 3, 3, 3) + ...
        2 * (KTfull(1, 1, 2, 2) + KTfull(1, 1, 3, 3) + KTfull(2, 2, 3, 3))) / 5;
    
    % isotropic component of the kurtosis tensor
    KTI = KTfull - meanKT * I4;
    KFA = norm(KTI(:)) / normKT;
else
    KFA = 0;
end


function KFA = ComputeKFA4(KT)
%
% computes the KFA given a 15-vector of kurtosis tensor values
% Ali ordering
%
% Author: Mark Van Horn
% Version: 2.6.0
% Last modified: 08/19/14

W1111 = KT(1);
W2222 = KT(2);
W3333 = KT(3);

W1122 = KT(10);
W1133 = KT(11);
W2233 = KT(12);

W1112 = KT(4);
W1113 = KT(5);
W1222 = KT(6);

W2223 = KT(8);
W1333 = KT(7);
W2333 = KT(9);

W1123 = KT(13);
W1223 = KT(14);
W1233 = KT(15);

W_F = sqrt(W1111 ^ 2 + W2222 ^ 2 + W3333 ^ 2 + 6 * W1122 ^ 2 + 6 * W1133 ^ 2 + 6 * W2233 ^ 2 + ...
    4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + 4 * W1333 ^ 2 + ...
    4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
if W_F < 1e-3,
    KFA = 0;
else
    Wbar = (W1111 + W2222 + W3333 + 2 * (W1122 + W1133 + W2233)) / 5;
    
    W_diff_F = sqrt((W1111 - Wbar) ^ 2 + (W2222 - Wbar) ^ 2 + (W3333 - Wbar) ^ 2 + ...
        6 * (W1122 - Wbar / 3) ^ 2 + 6 * (W1133 - Wbar / 3) ^ 2 + 6 * (W2233 - Wbar / 3) ^ 2 + ...
        4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + ...
        4 * W1333 ^ 2 + 4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
    
    KFA = W_diff_F / W_F;
end
