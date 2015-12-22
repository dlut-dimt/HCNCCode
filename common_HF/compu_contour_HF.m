% [sshr] = compu_contour_SSHR( cont )
% 
% Compute the global segment height representation of the given contour points.
%
% Output	
%	sshr		: global segment height representation, N columns, the i-th column is a N-3 
%				vector indicating the shr feature at the i-th point
%
% Input:	
%	cont	: Nx2 matrix, input N contour points
%
% how to get the height of every point related to the chord:
% triangle area = 0.5*det(triple) = 0.5*chord * height
%
%	Junwei Wang, MC lab, EI, hust.edu.cn  2010.08.27
%
function [hf] = compu_contour_HF( cont )
											
%------ Parameters ----------------------------------------------
n_pt= size(cont,1);
X	= cont(:,1);
Y	= cont(:,2);  % X,Y are both column vectors
hf = zeros(n_pt-3, n_pt);
% hf = zeros(n_pt, n_pt);
%-- Geodesic distances between all landmark points ---------
Xs			= repmat(X,1,n_pt);  % Xs = [X X ... X];
dXs			= Xs-Xs';
Ys			= repmat(Y,1,n_pt);  % Xs,Ys are both square matrix 
dYs			= Ys-Ys';
dis_mat		= sqrt(dXs.^2+dYs.^2);    % this data representation and algorithm to get all distances is great! use matrix operations naturally.
diameter    = max(dis_mat(:));  % max or mean both try

%-- SAR for every landmark point ---------
X3 = repmat(X,3,1);
Y3 = repmat(Y,3,1);

for p_index = 1+n_pt : n_pt+n_pt
    scale_index = n_pt/2-1; 
    left = p_index-scale_index;
    right = p_index+scale_index;
%     chord = pdist([X3(left) Y3(left); X3(right) Y3(right)]); % chord length
    height_vector = zeros(2*scale_index-1, 1);
%  height_vector = zeros(100, 1);
    %%%%%%%%%%%%%%%
%          clf;
%         figure(1);
%         hold on;
%     tangentVec = [X3(right)-X3(left)  Y3(right)-Y3(left)];
    for i = left+1 : right-1
%      for i = left-1 : right
%         tempVec = [X3(i)-X3(left)  Y3(i)-Y3(left)];
%         angle = acos(dot(tangentVec,tempVec)/(norm(tangentVec)*norm(tempVec)))*180/pi;
        eudist = pdist([X3(i) Y3(i); X3(left-1) Y3(left-1)]);
% %         
        temp = det([X3(left) Y3(left) 1; X3(i) Y3(i) 1; X3(right) Y3(right) 1]);
        if temp < 0
            eudist = -eudist;
        end
        height_vector(i-left) = eudist;
        
        %%%%%%%%%%% paint %%%%%%%%%%
%       if p_index == 89+n_pt && i == left+49 %|| p_index == 28+n_pt || p_index == 89+n_pt
%    
% %       plot(Y,X,'-g','linewidth',1);
%         plot(cont(:,1),cont(:,2),'r.');
%         plot(X3(left-1),Y3(left-1),'ms','linewidth',2);      
% %         line([X3(left-1) X3(i)],[Y3(left-1) Y3(i)]);
%         plot([X3(left) X3(i)],[Y3(left) Y3(i)],'-','linewidth',0.5);
%         plot([X3(left) X3(right)],[Y3(left) Y3(right)],'-','linewidth',0.5);
%         plot([X3(right) X3(i)],[Y3(right) Y3(i)],'-','linewidth',0.5);
%       end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
%         height_vector(i-left) = height_vector(i-left) /  pdist([X3(i) Y3(i); X3(left-1) Y3(left-1)]);
%           temp = det([X3(left) Y3(left) 1; X3(i) Y3(i) 1; X3(right) Y3(right) 1]);
%         acos(dot(A,B)/(norm(A)*norm(B)))*180/pi 
       %%%%%%%%%%%%%%%%%%%%%% The Triangle Descriptor%%%%%%%%%%%
%        dist1=pdist([X3(i) Y3(i); X3(left) Y3(left)]);
%        dist2=pdist([X3(i) Y3(i); X3(right) Y3(right)]);
%        s=(dist1+dist2+chord)/2;
%        delta=sqrt(s*(s-dist1)*(s-dist2)*(s-chord));
%        height_vector(i-left)=delta/s;
%        if temp < 0
%            height_vector(i-left) = -height_vector(i-left);
%        end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 计算曲率 %%%%%
%         tangentLeft = [X3(i)-X3(i-1) Y3(i)-Y3(i-1)];
%         tangentRight = [X3(i+1)-X3(i) Y3(i+1)-Y3(i)];
%         chordLeft = pdist([X3(i) Y3(i); X3(i-1) Y3(i-1)]);
%         chordRight = pdist([X3(i) Y3(i); X3(i+1) Y3(i+1)]);
%         angle = acos(dot(tangentLeft,tangentRight)/(norm(tangentLeft)*norm(tangentRight)))*180/pi;
%         ratio = angle / (chordLeft + chordRight);
%         
%         height_vector(i-left) = det([X3(left) Y3(left) 1; X3(i) Y3(i) 1; X3(right) Y3(right) 1]); % signed area
        %%%%%%%%%%%%%%%
         %% 角度
%          lineVec = [X3(i)-X3(left-1) Y3(i)-Y3(left-1)];             
         
%           if temp >= 0
% %              height_vector(i-left) =  pdist([X3(i) Y3(i); X3(left-1) Y3(left-1)]);
%                height_vector(i-left) = dot(tangentVec,lineVec) / (norm(tangentVec) * norm(lineVec));
%           else
% %              height_vector(i-left) =  - pdist([X3(i) Y3(i); X3(left-1) Y3(left-1)]);
%                height_vector(i-left) = - dot(tangentVec,lineVec) / (norm(tangentVec) * norm(lineVec));
%           end
% %            height_vector(i-left) =  pdist([X3(i) Y3(i); X3(left-1) Y3(left-1)]);
            
    end
%         hold off;
%         axis off;
%     height_vector = height_vector / chord; % signed height 每个点的特征对应一列
    hf(:, p_index-n_pt) = height_vector;
    
%     if  p_index == 89+n_pt
%         figure(4);
%         plot(1:100, [0 -15 height_vector' -18]); 
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
     figure(1);
        hold on;

   
%       plot(Y,X,'-g','linewidth',1);
        plot(cont(:,1),cont(:,2),'r.');
        plot(X3(left-1),Y3(left-1),'ms','linewidth',2);      
%         line([X3(left-1) X3(i)],[Y3(left-1) Y3(i)]);
        hold off;
        axis off;
%       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Normalize tar with the shape diameter --------------------
diameter  = max(hf(:)); 
hf = hf / diameter;

return;
