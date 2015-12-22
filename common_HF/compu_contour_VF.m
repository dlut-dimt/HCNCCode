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
function [hf] = compu_contour_VF( cont )
											
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
    left = p_index-scale_index;  %◊Û¡⁄”Ú
    right = p_index+scale_index; %”“¡⁄”Ú
%     chord = pdist([X3(left) Y3(left); X3(right) Y3(right)]); % chord length
    height_vector = zeros(2*scale_index-1, 1);
    if Y3(left)~= Y3(right)
       k = (Y3(left) - Y3(right))/(X3(left) - X3(right));
       k=-1/k;
%        b= k*X3(p_index)-Y3(p_index);
       b=k*X3(left-1)-Y3(left-1);
    else
       k=0;
    end
    
    for i = left+1 : right-1
        if k ~= 0
            height_vector(i-left) = (k*X3(i)-Y3(i)-b)/sqrt(k^2+1);
        else
            height_vector(i-left) = X3(i) - X3(left-1);
        end
    end
    hf(:, p_index-n_pt) = height_vector;
end

%-- Normalize tar with the shape diameter --------------------
hf = hf / diameter;

return;
