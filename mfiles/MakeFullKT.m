
function ktfull = MakeFullKT(KT)
% 
% generate the full kurtosis tensor given the 15-vector of tensor values
% 
% Author: Mark Van Horn
% Version: 2.6.0
% Last modified: 08/18/14

ktfull(1, 1, 1, 1) = KT(1);
ktfull(2, 2, 2, 2) = KT(2);
ktfull(3, 3, 3, 3) = KT(3);
ktfull(1, 1, 1, 2) = KT(4);
ktfull(1, 1, 1, 3) = KT(5);
ktfull(1, 2, 2, 2) = KT(6);
ktfull(1, 3, 3, 3) = KT(7);
ktfull(2, 2, 2, 3) = KT(8);
ktfull(2, 3, 3, 3) = KT(9);
ktfull(1, 1, 2, 2) = KT(10);
ktfull(1, 1, 3, 3) = KT(11);
ktfull(2, 2, 3, 3) = KT(12);
ktfull(1, 1, 2, 3) = KT(13);
ktfull(1, 2, 2, 3) = KT(14);
ktfull(1, 2, 3, 3) = KT(15);

% propagate values to other elements
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                I = sort([i, j, k, l]);
                ktfull(i, j, k, l) = ktfull(I(1), I(2), I(3), I(4));
            end
        end
    end
end