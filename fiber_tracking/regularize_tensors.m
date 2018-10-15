function [dt kt] = regularize_tensors(dt,kt,fa_t)
%In rare cases the diffusion tensor may be extremely isotropic with very
%small eigenvalues, causing the kurtosis dODF to have erratic behavior with
%very large values, as the kurtosis dODF evaluates the inverse of D. 
%
%Since the voxel is extremely isotropic, preserving the the principal
%orientation of the diffusion tensor is preserved for tractography. 
%
%This issue may arrise from constrained fitting, for example, if one
%eigenvalue is less than zero. Since the eigenvalues affect computation of
%the kurtosis tensor, the kurtosis tensor may be unreliable in these voxels
%as well. 
%
%Consequently we only utilize the diffusion tensor orientations. 
%
%This issue is very rare. Currently this is only implemented when FA >0.98,
%in which case the fiber bundles are estimated to be highly aligned. 
%
%Author: Russell Glenn
%Medical University of South Carolina

D = [dt(1) dt(4) dt(5); dt(4) dt(2) dt(6); dt(5) dt(6) dt(3)]; 
[V L] = eig(D); 
[L idx] = sort(diag(L),'descend'); 
V = V(:,idx); 

x = roots([2*(1-2*fa_t^2)/3,-4*L(1)/3,2*(1-fa_t^2)/3*L(1)^2]);


if ~isempty(x(x>0&x<L(1))); L(2:3) = x(x>0&x<L(1)); 
else Davg = trace(D)/3; L(L<0.1*Davg)=0.1*Davg; 
end

D = V*diag(L)*V^-1; 

dt = [D(1);D(5);D(9);D(2);D(3);D(6)]; 
kt = zeros(size(kt)); 