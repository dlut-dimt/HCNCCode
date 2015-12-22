function [BH,mean_dist]=hc_compute(Bsamp,Tsamp,mean_dist,nbins_theta,nbins_r,r_inner,r_outer,out_vec)
% [BH,mean_dist]=sc_compute(Bsamp,Tsamp,mean_dist,nbins_theta,nbins_r,r_inner,r_outer,out_vec);
%
% compute (r,theta) histograms for points along boundary 
%
% Bsamp is 2 x nsamp (x and y coords.)
% Tsamp is 1 x nsamp (tangent theta)
% out_vec is 1 x nsamp (0 for inlier, 1 for outlier)
%
% mean_dist is the mean distance, used for length normalization
% if it is not supplied, then it is computed from the data
%
% outliers are not counted in the histograms, but they do get
% assigned a histogram
%

nsamp=size(Bsamp,2);
in_vec=out_vec==0;

% compute height functions
r_array=zeros(nsamp,nsamp);
r_array(3:nsamp-1,:) = compu_contour_HF(Bsamp');
r_array=abs(r_array)';

% compute r,theta arrays
% r_array=real(sqrt(dist2(Bsamp',Bsamp'))); % real is needed to
                                          % prevent bug in Unix version
theta_array_abs=atan2(Bsamp(2,:)'*ones(1,nsamp)-ones(nsamp,1)*Bsamp(2,:),Bsamp(1,:)'*ones(1,nsamp)-ones(nsamp,1)*Bsamp(1,:))';
theta_array=theta_array_abs-Tsamp'*ones(1,nsamp);

% create joint (r,theta) histogram by binning r_array and
% theta_array

% normalize distance by mean, ignoring outliers
if isempty(mean_dist)
   tmp=r_array(in_vec,:);
   tmp=tmp(:,in_vec);
   mean_dist=mean(tmp(:));
end
r_array_n=r_array/mean_dist;

% use a log. scale for binning the distances
r_bin_edges=logspace(log10(r_inner),log10(r_outer),nbins_r);
r_array_q=zeros(nsamp,nsamp);
% r_array_height=zeros(nsamp,nsamp);
for m=1:nbins_r
    r_array_q=r_array_q+(r_array_n<r_bin_edges(m));%%%统计个数
%     if r_array_n<r_bin_edges(m)
%         r_array_height = r_array_height + r_array_n;  
%     end
end
fz=r_array_q>0; % flag all points inside outer boundary

% put all angles in [0,2pi) range
theta_array_2 = rem(rem(theta_array,2*pi)+2*pi,2*pi);
% quantize to a fixed set of angles (bin edges lie on 0,(2*pi)/k,...2*pi
theta_array_q = 1+floor(theta_array_2/(2*pi/nbins_theta));

nbins=nbins_theta*nbins_r;
BH=zeros(nsamp,nbins);
for n=1:nsamp
   fzn=fz(n,:)&in_vec;
   Sn=sparse(theta_array_q(n,fzn),r_array_q(n,fzn),1,nbins_theta,nbins_r);
%    BH(n,:)= Sn(:)';
%    for i=1:nsamp
%       Dn(theta_array_q(n,i),r_array_q(n,i)) = Dn(theta_array_q(n,i),r_array_q(n,i)) + r_array_n();
%    end
   Dn = sparse(theta_array_q(n,fzn),r_array_q(n,fzn),r_array_n(n,fzn),nbins_theta,nbins_r);
   BH(n,:)= Dn(:)';
%    temp1(n,:)=Sn(:)';
%    temp2(n,:)=Dn(:)';
%    for i=1:nbins
%        if (temp1(n,i) > 0)
%            BH(n,i)= temp2(n,i) / temp1(n,i);
%        end
%    end
%    scTest(Bsamp,r_inner,mean_dist,r_outer,nbins_theta,nbins_r,n,Sn);%画图
end
