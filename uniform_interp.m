function [XIs,YIs] = uniform_interp(Xs,Ys,n_samp)

% 	% test
% 	Xs		= (1:10)';
% 	Ys		= [3 4 12 4 3 5 6 12 2 9]';
% 	n_samp	= 11;
	
	XIs		= zeros(n_samp,1);
	YIs		= zeros(n_samp,1);
	
	len_all	= 0;
	difXs	= diff(Xs);
	difYs	= diff(Ys);
	seg_lens	= [0; sqrt(difXs.^2+difYs.^2)];
	len_all	= sum(seg_lens);
	d_len	= len_all/(n_samp+1);
	
	% interpolate to get each point
	n_pt	= length(Xs);
	for(ii=2:n_pt)
		seg_lens(ii)	= seg_lens(ii)+seg_lens(ii-1);
	end
	
	i_fill	= 0;
	cur_len	= d_len;
	for(ii=2:n_pt)
		while cur_len <= seg_lens(ii)  & i_fill<n_samp
			% interpolate a point
			r	= (cur_len-seg_lens(ii-1)) / (seg_lens(ii)-seg_lens(ii-1));
			x	= Xs(ii-1)+r*(Xs(ii)-Xs(ii-1));
			y	= Ys(ii-1)+r*(Ys(ii)-Ys(ii-1));
			i_fill	= i_fill+1;
			XIs(i_fill)	= x;
			YIs(i_fill)	= y;
			cur_len	= cur_len+d_len;
		end
	end
	
	% test
	if 0
		figure(97);	clf; hold on;
		plot(Xs,Ys,'r.-');
		plot(XIs,YIs,'b+');
		title(['n_samp=' i2s(n_samp) ', len=' i2s(length(XIs))]);
		keyboard;
	end