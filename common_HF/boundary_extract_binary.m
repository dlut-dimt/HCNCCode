function [Cs]=boundary_extract_binary(im)

c		= contourc(im,[.8 .8]);
len_c	= size(c,2);

if len_c>0
	n_max	= length(find(c(1,:)==.8));
	Cs		= cell(n_max,1);
	n_C		= 0;
	iHead	= 1;
	while iHead<len_c
		n_len	= c(2,iHead);
		iTail	= iHead+n_len;
		n_C		= n_C+1;
		Cs{n_C}	= c(:,iHead+1:iTail);
		iHead	= iTail+1;
	
		if 0	% test
			figure;	clf;
			imshow(im);	hold on;
			plot(Cs{n_C}(1,:), Cs{n_C}(2,:),'b');
			pause;
		end
	end
	
	if n_C<n_max
		Cs	= Cs(1:n_C);
	end
	
else
	Cs	= [];
end

