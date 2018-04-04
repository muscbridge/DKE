function [x, rnorme, rnorml, exitflag] = mlsei(C, d, A, b, Aeq, beq, lb, ub)
% mlsei: Constrained linear least squares based on the lsei algoithm from
% the slatec fortran library.
%
% X = mlsei(C,d,A,b) attempts to solve the least-squares problem
% 
%         min  0.5*(NORM(C*x-d)).^2       subject to    A*x <= b
%          x
% 
% where C is m-by-n.
% 
% X = mlsei(C,d,A,b,Aeq,beq) solves the least-squares
% (with equality constraints) problem:
% 
%         min  0.5*(NORM(C*x-d)).^2    subject to 
%          x                               A*x <= b and Aeq*x = beq
% 
% X = mlsei(C,d,A,b,Aeq,beq,LB,UB) defines a set of lower and upper
% bounds on the design variables, X, so that the solution
% is in the range LB <= X <= UB. Set LB(i) = -Inf if X(i) is unbounded 
% below; set UB(i) = Inf if X(i) is unbounded above.
%
%
% mlsei is a gateway function to the fortran routine dlsei, the
% double-precision version of the lsei algorithm from the slatec
% fortran mathemetical library available here:
%
%   http://www.netlib.org/slatec/
%
% For a detailed description of the algorithm used see:
%
% Richard J. Hanson and Karen H. Haskell. 1982. Algorithm 587: Two
% Algorithms for the Linearly Constrained Least Squares Problem. ACM Trans.
% Math. Softw. 8, 3 (September 1982), 323-333. DOI=10.1145/356004.356010
% http://doi.acm.org/10.1145/356004.356010 
%
% The required functions have been compiled into mex functions for several
% platforms and are accessed through the mexlsei routine. Unfortunately
% this means mlsei currently only works on the following architectures:
% win32, win64, glx, glxa64.
%

% Copyright Richard Crozier 2011

    % get the number of variables
    n = size(C, 2);
    
    if nargin == 2
        
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        
    elseif nargin == 4
        
        Aeq = [];
        beq = [];
    
    elseif nargin > 6
        
        % Bounds on the variables can be handled as additional 
        % linear constaints
        
        if length(lb) ~= n
            error('lb vector must be of same size as number of variables, set unbounded variables to -InF.')
        elseif length(ub) ~= n
            error('ub vector must be of same size as number of variables, set unbounded variables to InF')
        end
        
        % get the members of lb that are not equal to -inf (or inf)
        noninflb = ~isinf(lb);

        if nnz(noninflb) > 0

            % create an n x n matrix with diagonal line
            % of negative ones and zeros elswhere
            lbmatrix = -eye(n);

            % add non-Inf bounds to C matrix
            A = [A; lbmatrix(noninflb, :)];

            % add the constraints to the b vector
            b = [b; -lb(noninflb)];

        end

        % upper bounds
        noninfub = ~eq(ub,inf);
        
        if nnz(noninfub) > 0
            
            % create an n x n matrix with diagonal line
            % of one and zeros elswhere
            ubmatrix = eye(n);
            
            % add non-Inf bounds to C matrix
            A = [A; ubmatrix(noninfub,:)];
            
            b = [b; ub(noninfub)];
            
        end

    end

    % Build the w matrix for dlsei
    w = [Aeq, beq;
         C,   d;
        -A,  -b];
    
    mdw = size(w,1);
    me = size(Aeq,1);
    ma = size(C,1);
    mg = size(A,1);
    
    % the options for lsei
    prgopt = 1;
    
    % dlsei requires working space to be provided in
    % ws
    k = max(ma+mg, n);
    wsize = 2*(me+n)+k+(mg+2)*(n+7);
    ws = zeros(1,wsize);
    
    % dlsei requires working space to be provided in ip 
    %
    % From lsei doc:
    %
    % Before computation
    %
    %     IP(1),       The amounts of working storage actually 
    %     IP(2)        allocated for the working arrays WS(*) and 
    %                  IP(*), respectively.  These quantities are 
    %                  compared with the actual amounts of storage 
    %                  needed by DLSEI( ).  Insufficient storage allocated for either WS(*) or IP(*) is an */
    %                  error.  This feature was included in DLSEI( ) */
    %                  because miscalculating the storage formulas */
    %                  for WS(*) and IP(*) might very well lead to */
    %                  subtle and hard-to-find execution errors. */
    %
    %                  The length of WS(*) must be at least */
    %
    %                  LW = 2*(ME+N)+K+(MG+2)*(N+7) */
    %
    %                  where K = max(MA+MG,N) */
    %                  This test will not be made if IP(1).LE.0. */
    %
    %                  The length of IP(*) must be at least */
    %
    %                  LIP = MG+2*N+2 */
    %                  This test will not be made if IP(2).LE.0. */
    %
    % IP   The integer working array has three entries that provide rank 
    %      and working array length information after completion.
    %
    %      IP(1) = rank of equality constraint matrix.              
    %
    %      IP(2) = rank of reduced least squares problem.              
    %
    %      IP(3) = the amount of storage in the working array WS(*) that
    %              was actually used by the subprogram. The formula given
    %              above for the length of WS(*) is a necessary
    %              overestimate. If exactly the same problem matrices are
    %              used in subsequent executions, the declared dimension of
    %              WS(*) can be reduced to this output value.
    
    ipsize = mg+2*n+2;
    ip = zeros(1, ipsize);
    
    % We pass in the size of the 
    ip(1) = wsize;
    ip(2) = ipsize;
    
    % construct the input vector W
    w = reshape(w, 1, []);

    % call dlsei via the mex gateway function mexlsei
    [x, rnorme, rnorml, exitflag] = mexlsei(w, mdw, me, ma, mg, n, prgopt, ws, ip);
    
    if exitflag
        exitflag = -exitflag;
    end

end
