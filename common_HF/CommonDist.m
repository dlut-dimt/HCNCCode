function d = CommonDist(X1,X2,dist_type)
% CommonDist calculates dist or similarity between two histograms. This
% function implements the distances (esp. the BRD-dissimilarities) used 
% in the following paper
%   Use Bin-Ratio Information for Category and Scene Classification
%   N. Xie, H. Ling, W. Hu, and X. Zhang
%   IEEE Conf. on Computer Vision and Pattern Recognition (CVPR), San Francisco, 2010. %   X1 and X2 are row vectors
%   
% dist_type can be 'L1','HI','Chi-square','L1-BRD','BRD','Chi-square-BRD'
%   Note that this function does NOT include normalization
%
% Nianhua Xie and Haibin Ling, July 2010
n=size(X1,2);
switch dist_type
    case 'L1'
        d = sum(abs(X1-X2));
%     d=0;
%     for t=1:n
% %       w = 1/min(t,n-t);
%         w=1;
%         d = d + w*abs(X1(1,t)-X2(1,t));        
%     end
    case 'HI'
        d = sum(min(X1,X2));
    case 'Chi-square'
        d = 0.5*sum(((X1-X2).^2)./((X1+X2)+1e-10));
    case 'BRD'
        tmp1 = abs(1-X1*X2');
        tmp2 = ((X1-X2).^2 + 2*tmp1*(X1.*X2))./((X1+X2).^2+1e-10);
        d = abs(sum(tmp2));
    case 'L1-BRD'
        tmp1 = abs(1-X1*X2');
        tmp2 = ((X1-X2).^2 + 2*tmp1*(X1.*X2)).*abs(X1-X2)./((X1+X2).^2+1e-10);
        d = abs(sum(tmp2));
    case 'Chi-square-BRD'
        tmp1 = abs(1-X1*X2');
        tmp2 = ((X1-X2).^2 + 2*tmp1*(X1.*X2)).*(X1-X2).^2./((X1+X2).^3+1e-10);
        d = abs(sum(tmp2));
    otherwise
        d = -1;
        error('unknown dist type!');
end