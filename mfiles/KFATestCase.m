


function KFATestCase

N = 1000;

KFATest = zeros(7, 3);

KFATest(1, 1) = 0;
KFATest(2, 1) = .6324555322;
KFATest(3, 1) = .6831300509;
KFATest(4, 1) = .5773502692;
KFATest(5, 1) = .9979777532;
KFATest(6, 1) = .9866759405;
KFATest(7, 1) = .6491753010;

% W1111, W2222, W3333, W1122, W1133, W2233, W1112, W1113, W1222, W2223, W1333, W2333, W1123, W1223, W1233
t1 = [1, 1, 1, 1 / 3, 1 / 3, 1 / 3, 0, 0, 0, 0, 0, 0, 0, 0, 0];
t2 = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
t3 = [2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
t4 = [1, 2, 1, 0, 0, 1 / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0];
t5 = [1, 0, 1, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0];
t6 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
t7 = [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1];

%  indx converts the kurtosis vector from Jens' ordering to the existing
%  ordering in the DKE software
indx = [1 2 3 7 8 9 11 10 12 4 5 6 13 14 15];

T1 = t1(indx);
T2 = t2(indx);
T3 = t3(indx);
T4 = t4(indx);
T5 = t5(indx);
T6 = t6(indx);
T7 = t7(indx);

tstart = tic;
for n=1:N
    KFATest(1, 2) = ComputeKFA1(1, 1, 1, 1 / 3, 1 / 3, 1 / 3, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    KFATest(2, 2) = ComputeKFA1(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    KFATest(3, 2) = ComputeKFA1(2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    KFATest(4, 2) = ComputeKFA1(1, 2, 1, 0, 0, 1 / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    KFATest(5, 2) = ComputeKFA1(1, 0, 1, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0);
    KFATest(6, 2) = ComputeKFA1(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
    KFATest(7, 2) = ComputeKFA1(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1);
end
et1 = toc(tstart);


tstart = tic;
for n=1:N
    KFATest(1, 3) = ComputeKFA2(t1);
    KFATest(2, 3) = ComputeKFA2(t2);
    KFATest(3, 3) = ComputeKFA2(t3);
    KFATest(4, 3) = ComputeKFA2(t4);
    KFATest(5, 3) = ComputeKFA2(t5);
    KFATest(6, 3) = ComputeKFA2(t6);
    KFATest(7, 3) = ComputeKFA2(t7);
end
et2 = toc(tstart);


tstart = tic;
for n=1:N
    KFATest(1, 4) = ComputeKFA3(t1);
    KFATest(2, 4) = ComputeKFA3(t2);
    KFATest(3, 4) = ComputeKFA3(t3);
    KFATest(4, 4) = ComputeKFA3(t4);
    KFATest(5, 4) = ComputeKFA3(t5);
    KFATest(6, 4) = ComputeKFA3(t6);
    KFATest(7, 4) = ComputeKFA3(t7);
end
et3 = toc(tstart);

% Jens ordering
% W1111, W2222, W3333,
% W1122, W1133, W2233,
% W1112, W1113, W1222,
% W2223, W1333, W2333,
% W1123, W1223, W1233

% Ali ordering
% W1111, W2222, W3333, 
% W1112, W1113, W1222, 
% W1333, W2223, W2333,
% W1122, W1133, W2233,
% W1123, W1223, W1233

tstart = tic;
for n=1:N
    KFATest(1, 5) = ComputeKFA4(T1);
    KFATest(2, 5) = ComputeKFA4(T2);
    KFATest(3, 5) = ComputeKFA4(T3);
    KFATest(4, 5) = ComputeKFA4(T4);
    KFATest(5, 5) = ComputeKFA4(T5);
    KFATest(6, 5) = ComputeKFA4(T6);
    KFATest(7, 5) = ComputeKFA4(T7);
end
et4 = toc(tstart);


fprintf('\n')
for n=1:7
    disp([KFATest(n, 1) KFATest(n, 2) KFATest(n, 3) KFATest(n, 4) KFATest(n, 5)])
end

fprintf('method 1 - %f sec.\n', et1);
fprintf('method 2 - %f sec.\n', et2);
fprintf('method 3 - %f sec.\n', et3);
fprintf('method 4 - %f sec.\n', et4);

end     % function


% DKE ordering
% KTfull(1, 1, 1, 1) = KT(1);
% KTfull(2, 2, 2, 2) = KT(2);
% KTfull(3, 3, 3, 3) = KT(3);
% KTfull(1, 1, 1, 2) = KT(4);
% KTfull(1, 1, 1, 3) = KT(5);
% KTfull(1, 2, 2, 2) = KT(6);
% KTfull(1, 3, 3, 3) = KT(7);
% KTfull(2, 2, 2, 3) = KT(8);
% KTfull(2, 3, 3, 3) = KT(9);
% KTfull(1, 1, 2, 2) = KT(10);
% KTfull(1, 1, 3, 3) = KT(11);
% KTfull(2, 2, 3, 3) = KT(12);
% KTfull(1, 1, 2, 3) = KT(13);
% KTfull(1, 2, 2, 3) = KT(14);
% KTfull(1, 2, 3, 3) = KT(15);





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

% use with t1 (unpermuted for testing)
% W1111 = KT(1);
% W2222 = KT(2);
% W3333 = KT(3);
% W1122 = KT(4);
% W1133 = KT(5);
% W2233 = KT(6);
% W1112 = KT(7);
% W1113 = KT(8);
% W1222 = KT(9);
% W2223 = KT(10);
% W1333 = KT(11);
% W2333 = KT(12);
% W1123 = KT(13);
% W1223 = KT(14);
% W1233 = KT(15);

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
end




function KFA = ComputeKFA1(W1111, W2222, W3333, W1122, W1133, W2233, W1112, W1113, W1222, W2223, ...
    W1333, W2333, W1123, W1223, W1233)

W_F = sqrt(W1111 ^ 2 + W2222 ^ 2 + W3333 ^ 2 + 6 * W1122 ^ 2 + 6 * W1133 ^ 2 + 6 * W2233 ^ 2 + ...
    4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + 4 * W1333 ^ 2 + ...
    4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);

if W_F < 1e-3
    KFA = 0;
else
    Wbar = (W1111 + W2222 + W3333 + 2 * (W1122 + W1133 + W2233)) / 5;
    
    W_diff_F = sqrt((W1111 - Wbar) ^ 2 + (W2222 - Wbar) ^ 2 + (W3333 - Wbar) ^ 2 + ...
        6 * (W1122 - Wbar / 3) ^ 2 + 6 * (W1133 - Wbar / 3) ^ 2 + 6 * (W2233 - Wbar / 3) ^ 2 + ...
        4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + ...
        4 * W1333 ^ 2 + 4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
    
    KFA = W_diff_F / W_F;
end
end


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

end





%
% function KFATestCase
%
% N = 1000;
%
% KFATest = zeros(7, 5);
%
% KFATest(1, 1) = 0;
% KFATest(2, 1) = .6324555322;
% KFATest(3, 1) = .6831300509;
% KFATest(4, 1) = .5773502692;
% KFATest(5, 1) = .9979777532;
% KFATest(6, 1) = .9866759405;
% KFATest(7, 1) = .6491753010;
%
% tic
% for n=1:N
%     KFATest(1, 2) = ComputeKFA(1, 1, 1, 1 / 3, 1 / 3, 1 / 3, 0, 0, 0, 0, 0, 0, 0, 0, 0);
%     KFATest(2, 2) = ComputeKFA(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
%     KFATest(3, 2) = ComputeKFA(2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
%     KFATest(4, 2) = ComputeKFA(1, 2, 1, 0, 0, 1 / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0);
%     KFATest(5, 2) = ComputeKFA(1, 0, 1, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0);
%     KFATest(6, 2) = ComputeKFA(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
%     KFATest(7, 2) = ComputeKFA(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1);
% end
% toc
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
% t1 = [1, 1, 1, 1 / 3, 1 / 3, 1 / 3, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% t2 = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% t3 = [2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% t4 = [1, 2, 1, 0, 0, 1 / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% t5 = [1, 0, 1, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0];
% t6 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
% t7 = [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1];
%
% T1 = t1([1, 2, 3, 10, 11, 12, 4, 5, 6, 8, 7, 9, 13, 14, 15]);
% T2 = t2([1, 2, 3, 10, 11, 12, 4, 5, 6, 8, 7, 9, 13, 14, 15]);
% T3 = t3([1, 2, 3, 10, 11, 12, 4, 5, 6, 8, 7, 9, 13, 14, 15]);
% T4 = t4([1, 2, 3, 10, 11, 12, 4, 5, 6, 8, 7, 9, 13, 14, 15]);
% T5 = t5([1, 2, 3, 10, 11, 12, 4, 5, 6, 8, 7, 9, 13, 14, 15]);
% T6 = t6([1, 2, 3, 10, 11, 12, 4, 5, 6, 8, 7, 9, 13, 14, 15]);
% T7 = t7([1, 2, 3, 10, 11, 12, 4, 5, 6, 8, 7, 9, 13, 14, 15]);
%
%
% tic
% for n=1:N
%     KFATest(1, 3) = ComputeKFA2(T1);
%     KFATest(2, 3) = ComputeKFA2(T2);
%     KFATest(3, 3) = ComputeKFA2(T3);
%     KFATest(4, 3) = ComputeKFA2(T4);
%     KFATest(5, 3) = ComputeKFA2(T5);
%     KFATest(6, 3) = ComputeKFA2(T6);
%     KFATest(7, 3) = ComputeKFA2(T7);
% end
% toc
%
% tic
% for n=1:N
%     KFATest(1, 4) = ComputeKFA3(t1);
%     KFATest(2, 4) = ComputeKFA3(t2);
%     KFATest(3, 4) = ComputeKFA3(t3);
%     KFATest(4, 4) = ComputeKFA3(t4);
%     KFATest(5, 4) = ComputeKFA3(t5);
%     KFATest(6, 4) = ComputeKFA3(t6);
%     KFATest(7, 4) = ComputeKFA3(t7);
% end
% toc
%
% tic
% for n=1:N
%     KFATest(1, 5) = ComputeKFA4(T1);
%     KFATest(2, 5) = ComputeKFA4(T2);
%     KFATest(3, 5) = ComputeKFA4(T3);
%     KFATest(4, 5) = ComputeKFA4(T4);
%     KFATest(5, 5) = ComputeKFA4(T5);
%     KFATest(6, 5) = ComputeKFA4(T6);
%     KFATest(7, 5) = ComputeKFA4(T7);
% end
% toc
%
% fprintf('\n')
% for n=1:7
%     disp([KFATest(n, 1) KFATest(n, 2) KFATest(n, 3) KFATest(n, 4) KFATest(n, 5)])
% end
%
% end
%
%
% % original
% function KFA = ComputeKFA(W1111, W2222, W3333, W1122, W1133, W2233, W1112, W1113, W1222, W2223, ...
%     W1333, W2333, W1123, W1223, W1233)
%
% W_F = sqrt(W1111 ^ 2 + W2222 ^ 2 + W3333 ^ 2 + 6 * W1122 ^ 2 + 6 * W1133 ^ 2 + 6 * W2233 ^ 2 + ...
%     4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + 4 * W1333 ^ 2 + ...
%     4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
%
% Wbar = (W1111 + W2222 + W3333 + 2 * (W1122 + W1133 + W2233)) / 5;
%
% W_diff_F = sqrt((W1111 - Wbar) ^ 2 + (W2222 - Wbar) ^ 2 + (W3333 - Wbar) ^ 2 + ...
%     6 * (W1122 - Wbar / 3) ^ 2 + 6 * (W1133 - Wbar / 3) ^ 2 + 6 * (W2233 - Wbar / 3) ^ 2 + ...
%     4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + ...
%     4 * W1333 ^ 2 + 4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
%
% KFA = W_diff_F / W_F;
%
% end
%
%
%
% function KFA = ComputeKFA2(KT)
% %
% % computes the KFA given a 15-vector of kurtosis tensor values
% %
% % Author: Mark Van Horn
% % Version: 2.6.0
% % Last modified: 08/18/14
%
% % makes almost no diffference if I4 is calculated here or passed to the
% % function
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
% % original ordering
% % KTfull(1, 1, 1, 1) = KT(1);
% % KTfull(2, 2, 2, 2) = KT(2);
% % KTfull(3, 3, 3, 3) = KT(3);
% % KTfull(1, 1, 1, 2) = KT(4);
% % KTfull(1, 1, 1, 3) = KT(5);
% % KTfull(1, 2, 2, 2) = KT(6);
% % KTfull(1, 3, 3, 3) = KT(7);
% % KTfull(2, 2, 2, 3) = KT(8);
% % KTfull(2, 3, 3, 3) = KT(9);
% % KTfull(1, 1, 2, 2) = KT(10);
% % KTfull(1, 1, 3, 3) = KT(11);
% % KTfull(2, 2, 3, 3) = KT(12);
% % KTfull(1, 1, 2, 3) = KT(13);
% % KTfull(1, 2, 2, 3) = KT(14);
% % KTfull(1, 2, 3, 3) = KT(15);
%
% % ordering from Jens' code
% KTfull(1, 1, 1, 1) = KT(1);
% KTfull(2, 2, 2, 2) = KT(2);
% KTfull(3, 3, 3, 3) = KT(3);
% KTfull(1, 1, 2, 2) = KT(4);
% KTfull(1, 1, 3, 3) = KT(5);
% KTfull(2, 2, 3, 3) = KT(6);
% KTfull(1, 1, 1, 2) = KT(7);
% KTfull(1, 1, 1, 3) = KT(8);
% KTfull(1, 2, 2, 2) = KT(9);
% KTfull(2, 2, 2, 3) = KT(10);
% KTfull(1, 3, 3, 3) = KT(11);
% KTfull(2, 3, 3, 3) = KT(12);
% KTfull(1, 1, 2, 3) = KT(13);
% KTfull(1, 2, 2, 3) = KT(14);
% KTfull(1, 2, 3, 3) = KT(15);
%
% % propagate values to other elements
% for i = 1:3
%     for j = 1:3
%         for k = 1:3
%             for l = 1:3
%                 I = sort([i, j, k, l]);
%                 KTfull(i, j, k, l) = KTfull(I(1), I(2), I(3), I(4));
%             end
%         end
%     end
% end
%
% normKT = norm(KTfull(:));
%
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
% end
%
% % W1111, W2222, W3333, W1122, W1133, W2233, W1112, W1113, W1222, W2223,
% %     W1333, W2333, W1123, W1223, W1233
%
% % fastest method
% function KFA = ComputeKFA3(KT)
%
% W1111 = KT(1);
% W2222 = KT(2);
% W3333 = KT(3);
% W1122 = KT(4);
% W1133 = KT(5);
% W2233 = KT(6);
% W1112 = KT(7);
% W1113 = KT(8);
% W1222 = KT(9);
% W2223 = KT(10);
% W1333 = KT(11);
% W2333 = KT(12);
% W1123 = KT(13);
% W1223 = KT(14);
% W1233 = KT(15);
%
% W_F = sqrt(W1111 ^ 2 + W2222 ^ 2 + W3333 ^ 2 + 6 * W1122 ^ 2 + 6 * W1133 ^ 2 + 6 * W2233 ^ 2 + ...
%     4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + 4 * W1333 ^ 2 + ...
%     4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
%
% Wbar = (W1111 + W2222 + W3333 + 2 * (W1122 + W1133 + W2233)) / 5;
%
% W_diff_F = sqrt((W1111 - Wbar) ^ 2 + (W2222 - Wbar) ^ 2 + (W3333 - Wbar) ^ 2 + ...
%     6 * (W1122 - Wbar / 3) ^ 2 + 6 * (W1133 - Wbar / 3) ^ 2 + 6 * (W2233 - Wbar / 3) ^ 2 + ...
%     4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + ...
%     4 * W1333 ^ 2 + 4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
%
% KFA = W_diff_F / W_F;
%
% end
%
%
% % fastest method
% function KFA = ComputeKFA4(KT)
%
% W1111 = KT(1);
% W2222 = KT(2);
% W3333 = KT(3);
% W1112 = KT(4);
% W1113 = KT(5);
% W1222 = KT(6);
% W1333 = KT(7);
% W2223 = KT(8);
% W2333 = KT(9);
% W1122 = KT(10);
% W1133 = KT(11);
% W2233 = KT(12);
% W1123 = KT(13);
% W1223 = KT(14);
% W1233 = KT(15);
%
% W_F = sqrt(W1111 ^ 2 + W2222 ^ 2 + W3333 ^ 2 + 6 * W1122 ^ 2 + 6 * W1133 ^ 2 + 6 * W2233 ^ 2 + ...
%     4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + 4 * W1333 ^ 2 + ...
%     4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
%
% Wbar = (W1111 + W2222 + W3333 + 2 * (W1122 + W1133 + W2233)) / 5;
%
% W_diff_F = sqrt((W1111 - Wbar) ^ 2 + (W2222 - Wbar) ^ 2 + (W3333 - Wbar) ^ 2 + ...
%     6 * (W1122 - Wbar / 3) ^ 2 + 6 * (W1133 - Wbar / 3) ^ 2 + 6 * (W2233 - Wbar / 3) ^ 2 + ...
%     4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + ...
%     4 * W1333 ^ 2 + 4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
%
% KFA = W_diff_F / W_F;
%
% end
