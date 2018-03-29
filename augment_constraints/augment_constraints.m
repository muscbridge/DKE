function [gmat_aug ndir_aug C_aug c_aug] = augment_constraints(t, gmat, b, ndir, Kmin, NKmax)

% calculate diffusion tensor

dt = t(1:6);
DT = [dt(1), dt(4), dt(5); dt(4), dt(2), dt(6); dt(5), dt(6), dt(3)];
[V L] = eig(DT);
L = diag(L);        % take the diagonal of L
[L idx] = sort(L);
V = V(:, idx);      % sort e'vecs according to e'vals

% augment gradient set

gmat_aug = gmat;
ndir_aug = ndir;
for ibval = 1:length(b(2:end))
    lastdir = sum(ndir_aug(1:ibval));    % last row index of gradient set for the current b value
    gmat_aug = [gmat_aug(1:lastdir, :); real(V'); gmat_aug((lastdir + 1):end, :)];
    ndir_aug(ibval) = ndir(ibval) + 3;
end

% for ibval = 1:length(b(2:end))
%     g_aug{ibval} = [g{ibval}; V'];
%     ndir_aug(ibval) = ndir(ibval) + 3;
% end

% calculate coefficient matrices for augmented gradient set

[DG_aug KG_aug DKG_aug] = coeff_mat(gmat_aug, ndir_aug, b(2:end));
n_aug = sum(ndir_aug);

% calculate inequality constraint matrix C_aug and vector c_aug

C_aug = [-DG_aug,                    zeros(n_aug, 15);
        zeros(n_aug, 6),             -KG_aug; ...
        -DG_aug * NKmax / max(b) * b(2),   KG_aug];

if Kmin ~= 0
    D = DG_aug * t(1:6);    % directional diffusivities (D)
    D(D < 0) = 0;           % set negative D's to zero
else
    D = zeros(n_aug, 1);
end

c_aug = [zeros(n_aug, 1); -D.^2 .* Kmin * b(2); zeros(n_aug, 1)];



%--------------------------------------------------------------------------
% form the coefficients matrices
%--------------------------------------------------------------------------

function [DG KG DKG] = coeff_mat(gmat, ndir, b)

%coeff_mat  Form coefficient matrices for dke_core
%
%Syntax
%
%   [DG KG DKG] = coeff_mat(gcell, ndir, b)
%
%Inputs
%
%   gcell   1-by-nb cell array containing gradient directions for each
%           b-value, where nb is the number of nonzero b-values
%
%   ndir    1-by-nb vector containing the numbers of gradient directions for
%           b-values
%
%   b       1-by-nb vector of nonzero b-values (in s/mm^2 units)
%
%Outputs
%
%   DG      diffusion matrix (for use in constraints matrix)
%
%   KG      kurtosis matrix (for use in constraints matrix)
%
%   DKG     complete coefficients matrix (actual coefficients matrix)

nb = length(b);
ndir_total = sum(ndir(1:nb));
DG = zeros(ndir_total, 6);
DKG = zeros(ndir_total, 21);
KG = zeros(ndir_total, 15);

for ib = 1:nb

    sdir = sum(ndir(1:(ib-1)));    % starting row index for the current b value
    g = gmat(sdir + (1:ndir(ib)), :);
    for idir = 1:ndir(ib)

        % [D11 D22 D33 D12 D13 D23]
        DG(sdir + idir, :) = [g(idir, 1)^2 g(idir, 2)^2 g(idir, 3)^2 2*g(idir, 1)*g(idir, 2) 2*g(idir, 1)*g(idir, 3) 2*g(idir, 2)*g(idir, 3)];
    
        % [W1111 W2222 W3333 W1112 W1113 W1222 W1333 W2223 W2333 W1122 W1133 W2233 W1123 W1223 W1233]
        KG(sdir + idir, :) = [g(idir, 1)^4 g(idir, 2)^4 g(idir, 3)^4 ...
            4*(g(idir, 1)^3)*g(idir, 2) 4*(g(idir, 1)^3)*g(idir, 3) 4*g(idir, 1)*(g(idir, 2)^3) 4*g(idir, 1)*(g(idir, 3)^3) 4*(g(idir, 2)^3)*g(idir, 3) 4*g(idir, 2)*(g(idir, 3)^3)...
            6*(g(idir, 1)^2)*(g(idir, 2)^2) 6*(g(idir, 1)^2)*(g(idir, 3)^2) 6*(g(idir, 2)^2)*(g(idir, 3)^2) ...
            12*(g(idir, 1)^2)*g(idir, 2)*g(idir, 3) 12*g(idir, 1)*(g(idir, 2)^2)*g(idir, 3) 12*g(idir, 1)*g(idir, 2)*(g(idir, 3)^2)];
        
        DKG(sdir + idir, :) = [-b(ib) * DG(sdir + idir, :) (b(ib)^2 / b(1) / 6) * KG(sdir + idir, :)];

    end

end
