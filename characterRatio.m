function [ ratio ] = characterRatio( u ,v, p1, p2, p3)
%CHARACTER calculate character ratio between a serial points in a line of spacial
%object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   *Input*
%   u:  2*N matrix ,N is  the number of points to calculate.or just the
%   firstpoint
%   v,p1,p2,p3: the other points.
%
%   *Output*
%   ratio:  character ratio value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 2
        u(:,end+1)=v;
    case 3
        u(:,end+1)=v;
        u(:,end+1)=p1;
    case 4
        u(:,end+1)=v;
        u(:,end+1)=p1;
        u(:,end+1)=p2;
    case 5
        u(:,end+1)=v;
        u(:,end+1)=p1;
        u(:,end+1)=p2;
        u(:,end+1)=p3;
end

%remove duplicate columns from u 
i=1;
while(i<size(u,2))
    j=i+1;
    while(j<=size(u,2))
        if(abs(sum(u(:,i)-u(:,j)))<1.0e-8)
            u(:,j)=[];
        else
            j=j+1;
        end
    end
    i=i+1;
end

n=size(u,2);
if(n<1)
    error('Arguments not enough.character.');
elseif(n<3) % 
    ratio=0;
    return;
end

ratio=1;
for i=3:n
    C=cross(u(1,[1 2 i]),u(2,[1 2 i]));
    if(any(C==0) && size(u,1)==3)
            C=cross(u(1,[1 2 i]),u(3,[1 2 i]));
    end
    if(any(C==0))
            C=cross(u(1,[1 2 i]),[1 1 1]);
    end
    if(any(C==0))
        C=cross(u(2,[1 2 i]),[1 1 1]);
    end
    if(any(C==0))
        error('Character:zero in C.');
    end
    %ratio=ratio * (C(2)+0.001) ./ (C(1)+0.001);
    c=(C(2)) ./ (C(1));
    ratio=ratio * c;
end