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
function [hf] = compu_contour_ID( cont, im)
											
%------ Parameters ----------------------------------------------
n_pt= size(cont,1);
X	= cont(:,1);
Y	= cont(:,2);  % X,Y are both column vectors
hf = zeros(n_pt, n_pt);
% hf = zeros(n_pt, n_pt);
%-- Geodesic distances between all landmark points ---------
Xs			= repmat(X,1,n_pt);  % Xs = [X X ... X];
dXs			= Xs-Xs';
Ys			= repmat(Y,1,n_pt);  % Xs,Ys are both square matrix 
dYs			= Ys-Ys';
% dis_mat		= sqrt(dXs.^2+dYs.^2);    % this data representation and algorithm to get all distances is great! use matrix operations naturally.
% diameter    = max(dis_mat(:));  % max or mean both try

%-- SAR for every landmark point ---------
X3 = repmat(X,3,1);
Y3 = repmat(Y,3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------ ¼ÆËã inner distance ----------------------------------------------
%- the following build the graph, in that each pair of points has an edge
%between them if they can see each other inside the shape boundary
im	= double(im);
fg_mask	= double(im>.5);
E	= build_graph_contour_C(X,Y,fg_mask,1);
E	= E';
% disp_graph(V,E);		keyboard;
[dis_mat,ang_mat] = bellman_ford_ex_C(X,Y,E);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p_index = 1 : n_pt

    left = p_index-1;
    right = p_index+1;
    
    if p_index==1
        left=n_pt;
    end
    
    if p_index==n_pt
        right=1;
    end        
%     tangentVec = [X3(right)-X3(left)  Y3(right)-Y3(left)];
    for i = p_index : p_index+n_pt-1
          if i>n_pt
              i=mod(i,n_pt);
          end
          
          temp = det([X(left) Y(left) 1; X(i) Y(i) 1; X(right) Y(right) 1]);
          if temp >= 0
             height_vector(i) =  ang_mat(p_index,i);
          else
             height_vector(i) =  - ang_mat(p_index,i);
          end           
    end

    height_vector = circshift(height_vector, [0,-(p_index-1)]);
    hf(:, p_index) = height_vector';
end

%-- Normalize tar with the shape diameter --------------------
diameter  = max(hf(:)); 
hf = hf / diameter;

return;
