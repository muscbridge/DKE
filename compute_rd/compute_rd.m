%--------------------------------------------------------------------------
% compute the incomplete elliptic integral of the second kind (RD)
%--------------------------------------------------------------------------

function rd = compute_rd(x, y, z, errtol, lolim, uplim)

%RD Computes the incomplete elliptic integral of the second kind
%
%Usage
%
%   r = compute_rd(x, y, z, errtol, lolim, uplim)
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
%RD(x,y,z) = Integral from zero to infinity of
%
%                                -1/2     -1/2     -3/2
%                      (3/2)(t+x)    (t+y)    (t+z)    dt,
%
%where x and y are nonnegative, x + y is positive, and z is
%positive.  if x or y is zero, the integral is complete.
%the duplication theorem is iterated until the variables are
%nearly equal, and the function is then expanded in taylor
%series to fifth order.  
%
%Reference: B. B. Carlson, Computing elliptic integrals by duplication, Numer. 
%Math. 33 (1979), 1-16.  
%Coded by B. C. Carlson and Elaine M. Notis, Ames Laboratory-DOE, Iowa State 
%University, Ames, Iowa 50011. March 1, 1980.
%
%Check: rd(x,y,z) + rd(y,z,x) + rd(z,x,y) = 3 / dsqrt(x * y * z), where x, y, 
%and z are positive.
%
%lolim and uplim determine the range of valid arguments.
%lolim is not less than 2 / (machine maximum) ** (2/3).
%uplim is not greater than (0.1 * errtol / machine
%minimum) ** (2/3), where errtol is described below.
%
%errtol is set to the desired error tolerance.
%relative error due to truncation is less than
%3 * errtol ** 6 / (1-errtol) ** 3/2.

% Author: Ali Tabesh
% Last modified: 03/25/10

% ensure x, y, z are column vectors
x = x(:);
y = y(:);
z = z(:);

if any(min([x y], [], 2) < 0 | min([x+y z], [], 2) < lolim | max([x y z], [], 2) > uplim)
    error('Invalid input arguments to function compute_rd!')
end

xn = x;
yn = y;
zn = z;
sigma = 0;
power4 = 1;

% the following loop is a bit different than Carlson's code; though
% essentially doing the same thing
mu = (xn + yn + 3 * zn) * 0.2;
xndev = (mu - xn) ./ mu;
yndev = (mu - yn) ./ mu;
zndev = (mu - zn) ./ mu;
epslon = max(abs([xndev; yndev; zndev]));
while epslon > errtol
    xnroot = sqrt(xn);
    ynroot = sqrt(yn);
    znroot = sqrt(zn);
    lamda = xnroot .* (ynroot + znroot) + ynroot .* znroot;
    sigma = sigma + power4 ./ (znroot .* (zn + lamda));
    power4 = power4 * 0.25;
    xn = (xn + lamda) * 0.25;
    yn = (yn + lamda) * 0.25;
    zn = (zn + lamda) * 0.25;
    mu = (xn + yn + 3 * zn) * 0.2;
    xndev = (mu - xn) ./ mu;
    yndev = (mu - yn) ./ mu;
    zndev = (mu - zn) ./ mu;
    epslon = max(abs([xndev; yndev; zndev]));
end

c1 = 3 / 14;
c2 = 1 / 6;
c3 = 9 / 22;
c4 = 3 / 26;
ea = xndev .* yndev;
eb = zndev .* zndev;
ec = ea - eb;
ed = ea - 6 * eb;
ef = ed + ec + ec;
s1 = ed .* (-c1 + 0.25 * c3 * ed - 1.5 * c4 * zndev .* ef);
s2 = zndev .* (c2 * ef + zndev .* (-c3 * ec + c4 * zndev .* ea));
rd = 3 * sigma + power4 * (1 + s1 + s2) ./ (mu .* sqrt(mu));
