%--------------------------------------------------------------------------
% compute the incomplete elliptic integral of the second kind (RF)
%--------------------------------------------------------------------------

function rf = compute_rf(x, y, z, errtol, lolim, uplim)

%RF Computes the incomplete elliptic integral of the second kind
%
%Usage
%
%   r = compute_rf(x, y, z, errtol, lolim, uplim)
%
%Inputs/outputs
%
%   x, y, z      Input arguments
%   errtol       Error tolerance
%   lolim, uplim Lower and upper limits on inputs
%   r            Output
%
%The following description and the algorithm are from Carlson's Algorithm 577 
%available at
%
%http://portal.acm.org/citation.cfm?doid=355958.355970 (Supplement)
%
%RF(x,y,z) = integral from zero to infinity of
%
%                      -1/2     -1/2     -1/2
%            (1/2)(t+x)    (t+y)    (t+z)    dt,
%
%where x, y, and z are nonnegative and at most one of them
%is zero.  if one of them is zero, the integral is complete.
%the duplication theorem is iterated until the variables are
%nearly equal, and the function is then expanded in taylor
%series to fifth order.  
%
%Reference: B. C. Carlson, Computing elliptic integrals by duplication, 
%Numer. Math. 33 (1979), 1-16.
%
%Coded by B. C. Carlson and Elaine M. Notis, Ames Laboratory-DOE, 
%Iowa State University, Ames, Iowa 50011. March 1, 1980.
%
%Check by addition theorem: rf(x,x+z,x+w) + rf(y,y+z,y+w)
%= rf(0,z,w), where x,y,z,w are positive and x * y = z * w.
%
%lolim and uplim determine the range of valid arguments.
%lolim is not less than the machine minimum multiplied by 5.
%uplim is not greater than the machine maximum divided by 5.
%
%errtol is set to the desired error tolerance.
%relative error due to truncation is less than
%errtol ** 6 / (4 * (1 - errtol)).

% Author: Ali Tabesh
% Last modified: 03/25/10

% ensure x, y, z are column vectors
x = x(:);
y = y(:);
z = z(:);

if any(min([x y z], [], 2) < 0 | min([x+y x+z y+z], [], 2) < lolim | max([x y z], [], 2) > uplim)
    error('Invalid input arguments to function compute_rf!')
end

xn = x;
yn = y;
zn = z;

% the following loop is a bit different than Carlson's code; though
% essentially doing the same thing
mu = (xn + yn + zn) / 3;
xndev = 2 - (mu + xn) ./ mu;
yndev = 2 - (mu + yn) ./ mu;
zndev = 2 - (mu + zn) ./ mu;
epslon = max(abs([xndev; yndev; zndev]));
while epslon > errtol
    xnroot = sqrt(xn);
    ynroot = sqrt(yn);
    znroot = sqrt(zn);
    lamda = xnroot .* (ynroot + znroot) + ynroot .* znroot;
    xn = (xn + lamda) * 0.25;
    yn = (yn + lamda) * 0.25;
    zn = (zn + lamda) * 0.25;
    mu = (xn + yn + zn) / 3;
    xndev = 2 - (mu + xn) ./ mu;
    yndev = 2 - (mu + yn) ./ mu;
    zndev = 2 - (mu + zn) ./ mu;
    epslon = max(abs([xndev; yndev; zndev]));
end

c1 = 1 / 24;
c2 = 3 / 44;
c3 = 1 / 14;
e2 = xndev .* yndev - zndev .* zndev;
e3 = xndev .* yndev .* zndev;
s = 1 + (c1 .* e2 - 0.1 - c2 .* e3) .* e2 + c3 .* e3;
rf = s ./ sqrt(mu);
