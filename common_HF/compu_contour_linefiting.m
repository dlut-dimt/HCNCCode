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
function [hf] = compu_contour_linefiting( cont )
											
%------ Parameters ----------------------------------------------
n_pt= size(cont,1);
X	= cont(:,1);
Y	= cont(:,2);  % X,Y are both column vectors
hf = zeros(40, n_pt);
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
%     height_vector = zeros(2*scale_index-1, 1);

    % 旋转形状，使得该点的切线与X轴重合，保证旋转不变
    xvector = [1 0];
    tangentVec = [X3(right)-X3(left)  Y3(right)-Y3(left)];
    angle = acos(dot(tangentVec,xvector)/(norm(tangentVec)*norm(xvector)))*180/pi;
    Xinv=X'; Yinv=Y'; 
    if tangentVec(1,2) > 0 
       rotateMat=[cosd(-angle) -sind(-angle);sind(-angle) cosd(-angle)]*[Xinv;Yinv];
    else
       rotateMat=[cosd(angle) -sind(angle);sind(angle) cosd(angle)]*[Xinv;Yinv];
    end

    rotateX3 = repmat(rotateMat(1,:)',3,1);
    rotateY3 = repmat(rotateMat(2,:)',3,1);
    % 平移坐标，以当前点为中心点，保证平移不变
    rotateX3 = rotateX3 - rotateX3(left-1);
    rotateY3 = rotateY3 - rotateY3(left-1);
    
    % 显示旋转结果
%     hold on;
%     plot(X,Y,'r',rotateMat(1,:),rotateMat(2,:));    
%     plot(X3(left-1),Y3(left-1),'ms','linewidth',2);
%     plot(rotateMat(1,left-1),rotateMat(1,left-1),'ms','linewidth',2);
%     axis equal;
%     legend('原图像','旋转后的图像');

    % 计算拟合直线方程的系数
    lineVector = zeros(40, 1);
    for i=1:20      
        inx1=left-1+(i-1)*5;
        if inx1>300
            inx1=mod(inx1,300);
        end
        inx2=right+(i-1)*5;
        if inx2>300
            inx2=mod(inx2,300);
        end
        inx3=right-1+(i-1)*5;
        if inx3>300
            inx3=mod(inx3,300);
        end
        inx4=right-2+(i-1)*5;
        if inx4>300
            inx4=mod(inx4,300);
        end
        inx5=right-3+(i-1)*5;
        if inx5>300
            inx5=mod(inx5,300);
        end
        
        point1=[rotateX3(inx1) rotateY3(inx1)];
        point2=[rotateX3(inx2) rotateY3(inx2)];
        point3=[rotateX3(inx3) rotateY3(inx3)];
        point4=[rotateX3(inx4) rotateY3(inx4)];
        point5=[rotateX3(inx5) rotateY3(inx5)];
        [a,b]=linearfit([point1(1,1) point2(1,1) point3(1,1) point4(1,1) point5(1,1)],[point1(1,2) point2(1,2) point3(1,2) point4(1,2) point5(1,2)]);
        lineVector(2*i-1,1)=a;
        lineVector(2*i,1)=b;
    end
    hf(:, p_index-n_pt) = lineVector;
end
return;

% 最小二乘法拟合直线
function [a,b]=linearfit(x,y)
    xy=x.*y;
    x2=x.^2;
    x_mean=mean(x);
    y_mean=mean(y);
    xy_mean=mean(xy);
    x2_mean=mean(x2);
    b = (xy_mean-x_mean*y_mean)/(x2_mean-x_mean^2);
    
    if x2_mean~=x_mean^2        
        a=y_mean-b*x_mean;
    else
        a = x_mean;
    end
    
return;
