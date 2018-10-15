function [ODF_maxima, v1, gfa, gfa_RGB, nfd, A] = ODF_DKE_Optimization(dt,kt,alpha,S,IDX,area,BFGS_options)
%ODF_DKE_OPTIMIZATION evaluates the kurtosis ODF from DKE output
%
%Author: Russell Glenn
%Medical University of South Carolina
%Sept. 23, 2014
    
    D = [dt(1) dt(4) dt(5); dt(4) dt(2) dt(6); dt(5) dt(6) dt(3)]; %Diffusion Tensor

    W = zeros(3,3,3,3); %Kurtosis Tensor
    W(1,1,1,1) = kt(1);
    W(2,2,2,2) = kt(2);
    W(3,3,3,3) = kt(3);
    W(1,1,1,2) = kt(4);  W(1,1,2,1) = W(1,1,1,2); W(1,2,1,1) = W(1,1,1,2); W(2,1,1,1) = W(1,1,1,2);
    W(1,1,1,3) = kt(5);  W(1,1,3,1) = W(1,1,1,3); W(1,3,1,1) = W(1,1,1,3); W(3,1,1,1) = W(1,1,1,3);
    W(1,2,2,2) = kt(6);  W(2,1,2,2) = W(1,2,2,2); W(2,2,1,2) = W(1,2,2,2); W(2,2,2,1) = W(1,2,2,2);
    W(1,3,3,3) = kt(7);  W(3,1,3,3) = W(1,3,3,3); W(3,3,1,3) = W(1,3,3,3); W(3,3,3,1) = W(1,3,3,3);
    W(2,2,2,3) = kt(8);  W(2,2,3,2) = W(2,2,2,3); W(2,3,2,2) = W(2,2,2,3); W(3,2,2,2) = W(2,2,2,3);
    W(2,3,3,3) = kt(9);  W(3,2,3,3) = W(2,3,3,3); W(3,3,2,3) = W(2,3,3,3); W(3,3,3,2) = W(2,3,3,3);
    W(1,1,2,2) = kt(10); W(1,2,1,2) = W(1,1,2,2); W(1,2,2,1) = W(1,1,2,2); W(2,1,1,2) = W(1,1,2,2); W(2,1,2,1) = W(1,1,2,2); W(2,2,1,1) = W(1,1,2,2);
    W(1,1,3,3) = kt(11); W(1,3,1,3) = W(1,1,3,3); W(1,3,3,1) = W(1,1,3,3); W(3,1,1,3) = W(1,1,3,3); W(3,1,3,1) = W(1,1,3,3); W(3,3,1,1) = W(1,1,3,3);
    W(2,2,3,3) = kt(12); W(2,3,2,3) = W(2,2,3,3); W(2,3,3,2) = W(2,2,3,3); W(3,2,2,3) = W(2,2,3,3); W(3,2,3,2) = W(2,2,3,3); W(3,3,2,2) = W(2,2,3,3);
    W(1,1,2,3) = kt(13); W(1,1,3,2) = W(1,1,2,3); W(1,2,1,3) = W(1,1,2,3); W(1,2,3,1) = W(1,1,2,3); W(1,3,1,2) = W(1,1,2,3); W(1,3,2,1) = W(1,1,2,3); W(2,1,1,3) = W(1,1,2,3); W(2,1,3,1) = W(1,1,2,3); W(2,3,1,1) = W(1,1,2,3); W(3,1,1,2) = W(1,1,2,3); W(3,1,2,1) = W(1,1,2,3); W(3,2,1,1) = W(1,1,2,3);
    W(1,2,2,3) = kt(14); W(1,2,3,2) = W(1,2,2,3); W(1,3,2,2) = W(1,2,2,3); W(2,1,2,3) = W(1,2,2,3); W(2,1,3,2) = W(1,2,2,3); W(2,2,1,3) = W(1,2,2,3); W(2,2,3,1) = W(1,2,2,3); W(2,3,1,2) = W(1,2,2,3); W(2,3,2,1) = W(1,2,2,3); W(3,1,2,2) = W(1,2,2,3); W(3,2,1,2) = W(1,2,2,3); W(3,2,2,1) = W(1,2,2,3);
    W(1,2,3,3) = kt(15); W(1,3,2,3) = W(1,2,3,3); W(1,3,3,2) = W(1,2,3,3); W(2,1,3,3) = W(1,2,3,3); W(2,3,1,3) = W(1,2,3,3); W(2,3,3,1) = W(1,2,3,3); W(3,1,2,3) = W(1,2,3,3); W(3,1,3,2) = W(1,2,3,3); W(3,2,1,3) = W(1,2,3,3); W(3,2,3,1) = W(1,2,3,3); W(3,3,1,2) = W(1,2,3,3); W(3,3,2,1) = W(1,2,3,3);
    
    Davg = trace(D)/3;
    U = Davg*D^-1;

    A = ODF_coeff(U,W,alpha); %Get ODF coefficients
    
    %GET ODF FUNCTION FROM COEFFICIENTS IN SPHERICAL FORM
    %NOTE: minus sign has been added for optimization
    eval(sprintf(['fK = @(x)-(1/((sin(x(1))*cos(x(2)))^2*%d+(sin(x(1))*sin(x(2)))^2*%d+cos(x(1))^2*%d+2*(sin(x(1))*cos(x(2)))*(sin(x(1))*sin(x(2)))*%d+',...
    '2*(sin(x(1))*cos(x(2)))*cos(x(1))*%d+2*(sin(x(1))*sin(x(2)))*cos(x(1))*%d))^((%d+1)/2)*(1+(%d+',...
    '(%d*(sin(x(1))*cos(x(2)))^2+%d*(sin(x(1))*cos(x(2)))*(sin(x(1))*sin(x(2)))+%d*(sin(x(1))*cos(x(2)))*cos(x(1))+%d*(sin(x(1))*sin(x(2)))^2+',...
    '%d*(sin(x(1))*sin(x(2)))*cos(x(1))+%d*cos(x(1))^2)/((sin(x(1))*cos(x(2)))^2*%d+(sin(x(1))*sin(x(2)))^2*%d+cos(x(1))^2*%d+',...
    '2*(sin(x(1))*cos(x(2)))*(sin(x(1))*sin(x(2)))*%d+2*(sin(x(1))*cos(x(2)))*cos(x(1))*%d+2*(sin(x(1))*sin(x(2)))*cos(x(1))*%d)+',...
    '(%d*(sin(x(1))*cos(x(2)))^4+%d*(sin(x(1))*cos(x(2)))^3*(sin(x(1))*sin(x(2)))+%d*(sin(x(1))*cos(x(2)))^3*cos(x(1))+%d*(sin(x(1))*cos(x(2)))^2*(sin(x(1))*sin(x(2)))^2+',...
    '%d*(sin(x(1))*cos(x(2)))^2*(sin(x(1))*sin(x(2)))*cos(x(1))+%d*(sin(x(1))*cos(x(2)))^2*cos(x(1))^2+%d*(sin(x(1))*cos(x(2)))*(sin(x(1))*sin(x(2)))^3+',...
    '%d*(sin(x(1))*cos(x(2)))*(sin(x(1))*sin(x(2)))^2*cos(x(1))+%d*(sin(x(1))*cos(x(2)))*(sin(x(1))*sin(x(2)))*cos(x(1))^2+%d*(sin(x(1))*cos(x(2)))*cos(x(1))^3+%d*(sin(x(1))*sin(x(2)))^4+',...
    '%d*(sin(x(1))*sin(x(2)))^3*cos(x(1))+%d*(sin(x(1))*sin(x(2)))^2*cos(x(1))^2+%d*(sin(x(1))*sin(x(2)))*cos(x(1))^3+%d*cos(x(1))^4)/((sin(x(1))*cos(x(2)))^2*%d+(sin(x(1))*sin(x(2)))^2*%d+',...
    'cos(x(1))^2*%d+2*(sin(x(1))*cos(x(2)))*(sin(x(1))*sin(x(2)))*%d+2*(sin(x(1))*cos(x(2)))*cos(x(1))*%d+2*(sin(x(1))*sin(x(2)))*cos(x(1))*%d)^2)/24);'],...
    A(23),A(24),A(25),A(26),A(27),A(28),A(29),A(1),A(2),A(3),A(4),A(5),A(6),A(7),A(23),A(24),A(25),A(26),A(27),A(28),A(8),A(9),A(10),A(11),A(12),A(13),A(14),A(15),A(16),A(17),A(18),A(19),A(20),A(21),A(22),A(23),A(24),A(25),A(26),A(27),A(28)))

    
%--------------------------------------------------------------------------
%FIND PEAKS OVER SPHERICAL GRID
%--------------------------------------------------------------------------
    
    ODF = zeros(size(S,1),1);   %ODF VALUES
    ODF_maxima = [];            %Local maxima pair orientations (normalized)
    
    % Evaluate ODFs at each coordinate
    for i = 1:size(S,1); 
        ODF(i) = fK(S(i,:)); 
    end
    
    %Find local maxima pair (min here because you inverted fK to work with
    %fminunc)
    for i = 1:size(IDX,1)
        if ODF(IDX(i,1))==min(ODF(IDX(i,:)));
            ODF_maxima = [ODF_maxima; S(IDX(i,1),:)]; %At this point, they're still in spherical coordinates
        end
    end

%--------------------------------------------------------------------------
%REFINE PEAK ESTIMATE USING QUASI-NEWTON METHOD
%--------------------------------------------------------------------------
%Start with the ODF_maxima you found, and replace with refined estimates
    seed = ODF_maxima; ODF_maxima = []; peaks = []; %ODF value at local maxima (peak)
    for i = 1:size(seed,1); 
        try
            [p,fval,ef] = fminunc(fK,seed(i,:),BFGS_options);
            if ef==1||ef==3||ef==5
                peaks = [peaks abs(fval)];
                ODF_maxima = [ODF_maxima [sin(p(1))*cos(p(2));sin(p(1))*sin(p(2));cos(p(1))]]; %Cartesian coordinates [x;y;z];
            end
        catch 
        end
    end

%SORT PEAKS BY MAGNITUDE: This will always be a permutation of the indices to avoid
%the scenario of accidentally removing a peak that has exactly the same
%magnitude as another (as may occur in sumulations)
    peaks_sorted = sort(peaks,'descend');
    idx = [];
    for i = 1:size(peaks_sorted,2)
        idxi = find(peaks==peaks_sorted(i)); 
        for j = 1:length(idxi);
            if isempty(idx)||~sum(ismember(idx,idxi(j)))
                idx = [idx idxi(j)];
            end
        end
    end
    ODF_maxima = ODF_maxima(:,idx);
    
% ODF_maxima
%CHECK FOR DUPLICATES: This can occur for example due to the relative
%geometry of the sampling distribution used in estimating peaks to the ODF,
%but the BFGS algorithm will converge on the same peak, which will affect
%inflate fiber number
rFlag = zeros(1,size(ODF_maxima,2));    %remove flag
T=1;                            %Angle threshold in degrees (Quasi-Newton will converge)
for i = 1:size(ODF_maxima,2)-1
    for j = i+1:size(ODF_maxima,2)
        d1 = acosd(dot(ODF_maxima(:,i),ODF_maxima(:,j)));
        d2 = acosd(dot(ODF_maxima(:,i),-ODF_maxima(:,j)));
        if d1<T||d2<T; rFlag(j)=1; end
    end
end
ODF_maxima = ODF_maxima(:,~rFlag);

%Get principal eigenvector of the diffusion tensor to estimate maxima of
%the diffusion ODF
[evecs evals] = eig(D); evals = diag(evals); 
v1 = evecs(:,find(evals==max(evals),1));

%Update fractional anisotropy measures (n = 1281 for sampling distribution)

%For GFA
    smk = (sum(area(IDX(:,1)).*ODF(IDX(:,1))))^2;    %squared mean of kurtosis ODF
    msk = sum(area(IDX(:,1)).*ODF(IDX(:,1)).^2);     %mean squared of kurtosis ODF
    gfa = sqrt(1-smk/msk); 

    gfa_RGB = gfa*permute(ODF_maxima(:,1),[2 1 3])';

%Number of fiber directions
    nfd = size(ODF_maxima,2); 
