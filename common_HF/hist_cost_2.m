function HC=hist_cost_2(BH1,BH2)
% HC=hist_cost_2(BH1,BH2);
%
% same as hist_cost.m but BH1 and BH2 can be of different lengths

%labels	= ceil((nbins_theta*nbins_r:-1:1)/nbins_theta);
BH1=BH1';
BH2=BH2';
[nsamp1,nbins]=size(BH1);
[nsamp2,nbins]=size(BH2);

BH1n=BH1./repmat(sum(BH1,2)+eps,[1 nbins]);
BH2n=BH2./repmat(sum(BH2,2)+eps,[1 nbins]);
% 
 tmp1=repmat(permute(BH1n,[1 3 2]),[1 nsamp2 1]);
 tmp2=repmat(permute(BH2n',[3 2 1]),[nsamp1 1 1]);
%  HC=0.5*sum(((tmp1-tmp2).^2)./(tmp1+tmp2+eps),3);
 HC = sum(abs(tmp1-tmp2),3);
% for i=1:nsamp1
%     BH1n(i,:)=BH1(i,:)/norm(BH1(i,:));
%     BH2n(i,:)=BH2(i,:)/norm(BH2(i,:));
% %     BH1n_l1(i,:)=BH1(i,:)/sum(BH1(i,:));
% %     BH2n_l1(i,:)=BH2(i,:)/sum(BH2(i,:));
% end
% % %%%ÃÌº”»®÷ÿ%%%%%
% % % nbins_theta=12;
% % % labels=zeros(nsamp1,nsamp1,nbins);
% % % for i=1:nsamp1
% % %     for j=1:nsamp1
% % %        labels(i,j,:)= 0.1*ceil((nbins:-1:1)/nbins_theta);
% % %     end
% % % end
% % %%%
% HC=zeros(nsamp1,nsamp1);
% for i=1:nsamp1
%    for j=1:nsamp1
%        HC(i,j)=CommonDist(BH1n(i,:),BH2n(j,:),'L1');
%        %HC(j,i)=HC(i,j);
%    end
% end

% %%%%%%%%%%BRD (each row calculates BRD seperately)%%%%%%
% n_dist	= 5;
% n_theta	= 12;
% temp1=zeros(1,n_theta);
% temp2=zeros(1,n_theta);
% for i=1:nsamp1
%    for j=1:nsamp1
%       for k=1:n_dist
%           temp1=BH1(i,(k-1)*n_theta+1:k*n_theta);
%           temp2=BH2(i,(k-1)*n_theta+1:k*n_theta);
%           norm1=norm(temp1(1,:));
%           norm2=norm(temp2(1,:));
%           if norm1~=0
%             temp1=temp1/norm1;
%           end
%           if norm2~=0
%             temp2=temp2/norm2;
%           end
%           HC(i,j)=HC(i,j)+CommonDist(temp1,temp2,'BRD');
%       end
%       HC(i,j)=HC(i,j)/n_dist;
%    end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%% X2 (each row calculates BRD seperately)%%%%%%
% n_dist	= 5;
% n_theta	= 12;
% temp1=zeros(1,n_theta);
% temp2=zeros(1,n_theta);
% for i=1:nsamp1
%    for j=1:nsamp1
%       for k=1:n_dist
%           temp1=BH1(i,(k-1)*n_theta+1:k*n_theta);
%           temp2=BH2(i,(k-1)*n_theta+1:k*n_theta);
%           sum1=sum(temp1(1,:));
%           sum2=sum(temp2(1,:));
%           if sum1~=0
%             temp1=temp1/sum1;
%           end
%           if sum2~=0
%             temp2=temp2/sum2;
%           end
%           HC(i,j)=HC(i,j)+CommonDist(temp1,temp2,'Chi-square');
%       end
%       HC(i,j)=HC(i,j)/n_dist;
%    end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%L1-BRD(bin ratio dissimilarity)%%%%%%%%%%%%%%%
% for i=1:nsamp1
%    for j=1:nsamp1
%        temp1=sum(abs(BH1n(i,:)-BH2n(j,:)));
%        temp2=sum(abs(BH1n(i,:)+BH2n(j,:)).^2);
%        temp3=0;
%        for k=1:nbins
%            if BH1n(i,k)==0&&BH2n(j,k)==0
%                temp=0;         
%            else
%                temp=(abs(BH1n(i,k)-BH2n(j,k))*BH1n(i,k)*BH2n(j,k))/((BH1n(i,k)+BH2n(j,k))^2);
%           end
%            temp3=temp3+temp;
%        end
%        HC(i,j)=temp1-temp2*temp3;
%    end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%X2-BRD(bin ratio dissimilarity)%%%%%%%%%%%%%%%
% for i=1:nsamp1
%    for j=1:nsamp1
%        temp2=sum(abs(BH1n(i,:)+BH2n(j,:)).^2);
%        temp3=0;
%        temp1=0;
%        temp=0;  
%        temp4=0;
%        for k=1:nbins
%            if BH1n(i,k)==0&&BH2n(j,k)==0
%                temp=0;
%                temp4=0;
%            else
%                temp=((BH1n(i,k)-BH2n(j,k))^2*BH1n(i,k)*BH2n(j,k))/((BH1n(i,k)+BH2n(j,k))^3);
%                temp4=0.5*(BH1n(i,k)-BH2n(j,k))^2/(BH1n(i,k)+BH2n(j,k));
%           end
%            temp3=temp3+temp;
%            temp1=temp1+temp4;
%        end
%        HC(i,j)=temp1-2*temp2*temp3;
%    end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%X2-BRD with different normlization (bin ratio dissimilarity)%%%%%%%%%%%%%%%
% for i=1:nsamp1
%    for j=1:nsamp1
%        temp2=sum(abs(BH1n(i,:)+BH2n(j,:)).^2);
%        temp3=0;
%        for k=1:nbins
%            if BH1n(i,k)==0&&BH2n(j,k)==0
%                temp=0;
%            else
%                temp=((BH1n_l1(i,k)-BH2n_l1(j,k))^2*BH1n(i,k)*BH2n(j,k))/((BH1n(i,k)+BH2n(j,k))^2*(BH1n_l1(i,k)+BH2n_l1(j,k)));
%           end
%            temp3=temp3+temp;
%        end
%        x2_cost=CommonDist(BH1n_l1(i,:),BH2n_l1(j,:),'Chi-square');
%        HC(i,j)=x2_cost-2*temp2*temp3;
%    end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%L1-BRD with different normlization (bin ratio dissimilarity)%%%%%%%%%%%%%%%
% for i=1:nsamp1
%    for j=1:nsamp1
%        ll_cost=sum(abs(BH1n_l1(i,:)-BH2n_l1(j,:)));
%        temp2=sum(abs(BH1n(i,:)+BH2n(j,:)).^2);
%        temp3=0;
%        for k=1:nbins
%            if BH1n(i,k)==0&&BH2n(j,k)==0
%                temp=0;         
%            else
%                temp=(abs(BH1n_l1(i,k)-BH2n_l1(j,k))*BH1n(i,k)*BH2n(j,k))/((BH1n(i,k)+BH2n(j,k))^2);
%           end
%            temp3=temp3+temp;
%        end
%        HC(i,j)=ll_cost-temp2*temp3;
%    end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%X2-L1(bin ratio dissimilarity)%%%%%%%%%%%%%%%
% for i=1:nsamp1
%    for j=1:nsamp1  
%        temp1=0;
%        temp=0;  
%        for k=1:nbins
%            if BH1n_l1(i,k)==0&&BH2n_l1(j,k)==0
%                temp1=0;
%            else
%                temp1=0.5*abs(BH1n_l1(i,k)-BH2n_l1(j,k))*(BH1n_l1(i,k)-BH2n_l1(j,k))^2/(BH1n_l1(i,k)+BH2n_l1(j,k));
%            end
%            temp=temp+temp1;
%        end
%        HC(i,j)=temp;
%    end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
