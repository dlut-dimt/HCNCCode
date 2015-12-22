function testMultiIntevals_CHN(mirror,order,dirname)
%%Calculating FCNs in shapes
addpath common_HF;

sImage  = strcat('./',dirname,'/');

f_structure = dir(sImage);
m = length(f_structure);

ifname = cell(1, m-2);
for i = 3 : m  % batch process
      ifname{i-2} = strcat(sImage, f_structure(i).name);
end

n_contour  = 100;
times=25;
resultdata = cell(1,m-2);

for num = 1:m-2
        im =  double(imread(ifname{num}));
%         im = double(imread('./bell.gif'));
        % im = imresize(im,[1200 1200]);
        % imshow(im);
        im = im(:,:,1);          % image preprocess 1: for binary images, im(:,:,1) is simply the image itself. here is to avoid colour image data, which can't be processed by the 'extract_longest_cont' function
       % imshow(im);
        % im(im < 100) = 0;          
        im(im > 0) = 255;        % image preprocess 2: omit noise pixel values
        im=255-im;                % additional process: for 255 background and 0 foreground (check the boundary, or may cause fatal error!)
       
        
        [len, wid] = size(im);
        im2 = ones(len+2, wid+2).*im(len,1);
        im2(2:len+1,2:wid+1) = im;  % image preprocess 3: add empty lines around the object, avoiding contour error (which will make all later work useless!)
        if(mirror ==1)
             im2 = im2';%mirror image
        end

        %%% ????????
        edgedata = extract_longest_cont(im2, n_contour);
    %    disp(edgedata); 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        edgedata = edgedata';

        %%%%%%%%%%%%%%%%%%%
        edgedata = repmat(edgedata,1,4);
       
        featureVec = zeros(times,n_contour);
        %%% Calculate FCNs
        for k = (n_contour+1) : 2*n_contour
            p3 = edgedata(:,k);
            for  intervals=1:times              
                p2 = edgedata(:,k-intervals);
                p4= edgedata(:,k+intervals);
                
                
                if(strcmp(order,'order1'))
                     p1 = edgedata(:,k-2*intervals);
                     p5 = edgedata(:,k+2*intervals);
                else
                    % exchange p1 p5
                     p5 = edgedata(:,k-2*intervals);
                     p1 = edgedata(:,k+2*intervals);
                end
               
  %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [A(1,:), A(2,:)] = linecrosspoint(p1,p3,p2,p4);
                if isnan(A(1,1))
                       featureVec(times+1-intervals,k-n_contour) = 0;
                       continue;
                end     
                [Inner(1,:), Inner(2,:)] = linecrosspoint(p3,A,p1,p5);
                if isnan(Inner(1,1))
                       featureVec(times+1-intervals,k-n_contour) = 0;
                       continue;
                end     
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                crosspoints = [Inner,p2,p3,p4,p5];
                p1 = Inner;
                r=threePointsColinear(crosspoints(1,:),crosspoints(2,:));         
                if r==1
                    featureVec(times+1-intervals,k-n_contour) = 0; 
%                     fprintf('i=%d,there is collinear condition.\n',num);
                else 
                    [P(1,:), P(2,:)] = linecrosspoint(p1,p2,p4,p5);
                    if isnan(P(1,1))
                       featureVec(times+1-intervals,k-n_contour) = 0;
                       continue;
                    end
                    [Q(1,:), Q(2,:)] = linecrosspoint(p1,p2,p3,p4);
                    if isnan(Q(1,1))
                       featureVec(times+1-intervals,k-n_contour) = 0;
                       continue;
                    end                    
                    [R(1,:), R(2,:)] = linecrosspoint(p2,p3,p4,p5);
                    if isnan(R(1,1))
                       featureVec(times+1-intervals,k-n_contour) = 0;
                       continue;
                    end                    
                    [K(1,:), K(2,:)] = linecrosspoint(p1,R,p5,Q);
                    if isnan(K(1,1))
                       featureVec(times+1-intervals,k-n_contour) = 0;
                       continue;
                    end                    
                    [M(1,:), M(2,:)] = linecrosspoint(P,K,p1,p5);
                     if isnan(M(1,1))
                       featureVec(times+1-intervals,k-n_contour) = 0;
                       continue;
                    end                   
                    [N(1,:), N(2,:)] = linecrosspoint(P,p3,p1,p5);
                    if isnan(N(1,1))
                       featureVec(times+1-intervals,k-n_contour) = 0;
                       continue;
                    end                    
%                     charatio = 1;
                    chra1 = characterRatio(P,p1,p2,Q);
                    chra2 = characterRatio(p1,p5,M,N);
                    chra3 = characterRatio(p5,P,p4,R);
                    charatio = chra1 * chra2 * chra3;
%                     featureVec(times+1-intervals,k-n_contour) = charatio;
%                     if  isnan(charatio) || abs(charatio) > 1 
%                         featureVec(times+1-intervals,k-n_contour) = 0;% ???????????
%                     else
                    featureVec(times+1-intervals,k-n_contour) = charatio;  
%                     end
                end        
                
            end
        end
     
        resultdata{1,num} = featureVec;
%         save('bell dot pi feature.mat','featureVec');
     fprintf('%d\n',num);
end  
% close all;
if(mirror==1)
    save(strcat(dirname,order,' img1.mat'),'resultdata');
else
    save(strcat(dirname,order,' img0.mat'),'resultdata');
end
end

%%% 
function r=threePointsColinear(x,y)
ii=nchoosek(1:length(x),3);
xx=x(ii);
yy=y(ii);
cc=((yy(:,2)-yy(:,1)).*(xx(:,3)-xx(:,1))-(xx(:,2)-xx(:,1)).*(yy(:,3)-yy(:,1)));
r=any(cc==0);
end

function [X Y]= linecrosspoint(X1,Y1,X2,Y2)
    if X1(1)==Y1(1)
        X=X1(1);
        k2=(Y2(2)-X2(2))/(Y2(1)-X2(1));
        b2=X2(2)-k2*X2(1); 
        Y=k2*X+b2;
        return;
    end
    if X2(1)==Y2(1)
        X=X2(1);
        k1=(Y1(2)-X1(2))/(Y1(1)-X1(1));
        b1=X1(2)-k1*X1(1);
        Y=k1*X+b1;
        return;
    end
    if X1(1)~=Y1(1) && X2(1)~=Y2(1)
        k1=(Y1(2)-X1(2))/(Y1(1)-X1(1));
        k2=(Y2(2)-X2(2))/(Y2(1)-X2(1));
        b1=X1(2)-k1*X1(1);
        b2=X2(2)-k2*X2(1);
        if k1==k2
           X=NaN;
           Y=NaN;
        else
        X=(b2-b1)/(k1-k2);
        Y=k1*X+b1;
        end
    end
end

