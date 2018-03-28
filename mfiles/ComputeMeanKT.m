function [MeanKT]=ComputeMeanKT(KT)
% computes the MeanKT given a 15-vector of kurtosis tensor values
% Ali ordering
%
% Author: Emilie McKinnon
% Version: 2.6.0
% Last modified: 11/25/14

%   KT 15-by-1 vector containing elements of KT (mm^2/s units)
%              
W1111 = KT(1);
W2222 = KT(2);
W3333 = KT(3);
W1122 = KT(10);
W1133 = KT(11);
W2233 = KT(12);

% mean kurtosis tensor: W=1/5 * (W1111 + W2222 + W3333 + 2*W1122+ 2*W1133+ 2*W2233) 

MeanKT= 1/5 * (W1111 + W2222+ W3333 + 2 * (W1122 + W1133 + W2233) );

