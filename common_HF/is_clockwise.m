function bClk = is_clockwise(cont,fg_mask)

if 1	% new robust approach
	X	= cont(:,1);
	Y	= cont(:,2);
	Xs	= [X; X(1)];
	Ys	= [Y; Y(1)];
	dX	= diff(Xs);
	dY	= diff(Ys);
	
	angs	= atan2(dY,dX);
	angs	= [angs; angs(1)];
	dA		= diff(angs);
	id		= find(dA>pi);
	dA(id)	= dA(id)-2*pi;
	id		= find(dA<-pi);
	dA(id)	= dA(id)+2*pi;
	bClk	= sum(dA)<0;
% 	keyboard;
	
else	
	X	= cont(:,1);
	Y	= cont(:,2);
	Xs	= [X; X(1)];
	Ys	= [Y; Y(1)];
	dX	= diff(Xs);
	dY	= diff(Ys);
	d	= sqrt(dX.^2+dY.^2);
	coss 	= dX./d;
	sins	= dY./d;
	%thetas	= atan2(dY,dX);
	
	Xs	= round(X+1.5*sins);
	Ys	= round(Y-1.5*coss);
	[nRow,nCol]	= size(fg_mask);
	id_gd		= find(Xs>0 & Xs<=nCol & Ys>0 & Ys<=nRow);
	
	Xs	= Xs(id_gd);
	Ys	= Ys(id_gd);	
	idx	= round(nRow*(Xs-1) + Ys);
	
 	n	= length(X);
	f	= fg_mask(idx);
	bClk= sum(f)>(n/2);
end

return;
