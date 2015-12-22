function [C,T]=DPMatching(A, bCircular, thre)
%
%[C,T]=DPMatching(A, thre)
%	A 	- a square cost matrix.
%	thre- use average*thre as the occlusion
%	C 	- the optimal assignment.
%	T 	- the cost of the optimal assignment.
%	
%	bCircular

if ~exist('thre')		thre		= .2;		end			
if ~exist('bCircular')	bCircular	= 1;		end			

[M,N] = size(A);
A2 	= fliplr(A);		% two direction

if bCircular
	[C,T]=DPMatchingCircular(A, thre);
else
	[C,T]=DPMatchingSimple(A, thre);
end
return;

% if bCircular
% 	[C1,T1]=DPMatchingCircular(A, thre);
% 	[C2,T2]=DPMatchingCircular(A2, thre);
% else
% 	[C1,T1]=DPMatchingSimple(A, thre);
% 	[C2,T2]=DPMatchingSimple(A2, thre);
% end
% 
% if T1<=T2
% 	T = T1;
% 	C = C1;
% else
% 	T = T2;
% 	C = C2;
% 	good = find(C<=N);
% 	C(good) = N-C(good)+1;
% end		
 
%-----------------------------------------------------------------------------
% Circular mapping !!
function [C,T]=DPMatchingCircular(A, thre)
	[M,N]	= size(A);
	Anew	= A;
	T		= 100000;
	for iStart=1:N
		[CC,TT]	= DPMatchingSimple(Anew, thre, T);
		if T-0.00001 > TT
			T 	= TT;
			C 	= CC;
			ind	= iStart;
			disp(['start=' int2str(iStart) '  T=' num2str(T) '   # of match= ' int2str(length(find(C<1.5*N)))]);
		end
		Anew	= [Anew(:,2:N),Anew(:,1)];
	end
	
	for ii=1:M
		if C(ii) <= N		% not occluded
			C(ii) = C(ii)+ind-1;
			if C(ii)>N
				C(ii) = C(ii)-N;
			end
		end
	end


%-----------------------------------------------------------------------------
function [C,T]=DPMatchingSimple(A, thre, uplimit)
	
	if ~exist('uplimit')		uplimit = 1000000;	end
	
	avg 	= mean(A(:));	
	pen1	= avg*thre;				% ocllusion penalty

	D 		= zeros(size(A));		% DP matrix
	links	= zeros(size(D));
	[M,N]	= size(D);
		
	% initialize
	[D(1,1),links(1,1)] = min([A(1,1),pen1]);
	for jj=2:N
		tmp = [A(1,jj)+(jj-1)*pen1, pen1+D(1,jj-1)];
		[D(1,jj), id] = min(tmp);
		if id==1		links(1,jj) = 1;
		else 			links(1,jj) = 2;
		end
	end
	for ii=2:M
		tmp = [A(ii,1)+(ii-1)*pen1, pen1+D(ii-1,1)];
		[D(ii,1), id] = min(tmp);
		if id==1		links(ii,1) = 1;
		else 			links(ii,1) = 2;
		end
	end
	
	for ii=2:M
		for jj=2:N
			tmp = [A(ii,jj)+D(ii-1,jj-1), pen1+D(ii-1,jj), pen1+D(ii,jj-1)];
			[D(ii,jj), links(ii,jj)] = min(tmp);
		end
	end
	
	% Get the mapping result
	OCL	= 4*N;		% ocllusion using this index
	C	= OCL*ones(1,M);
	T 	= D(M,N);
	if T > uplimit		% terminate before final result
		C = [];
		%disp(['Terminate without computing C,  T=', num2str(T)]);
		return;
	end

	ii 	= M;
	jj 	= N;
	while ii>0 & jj>0
		switch links(ii,jj)
		case 1
			C(ii)	= jj;
			ii 		= ii-1;
			jj		= jj-1;
			
		case 2
			ii 		= ii-1;
			
		case 3
			jj 		= jj-1;
		end
	end
	
% 	keyboard;
% 	T,C(1:10)