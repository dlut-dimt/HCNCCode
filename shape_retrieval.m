function shape_retrieval(dirname)
%%Calculate similarities between any two shapes

addpath common_HF;
%clear;
%clc;

normMat = load(strcat(dirname,' img0.mat'));
normMat = normMat.resultdata;
normMat2 = load(strcat(dirname,' img1.mat'));
normMat2 = normMat2.resultdata;

delete(strcat(dirname,' img0.mat'));
delete(strcat(dirname,' img1.mat'));


m=length(normMat);

%- matching parameters
num_start = 100;
search_step	= 1;
thre = 0.1;

Score=Inf(m,m);
% TA1 = fix(clock); % start time
for k1 = 1:m%1:m-2  % the diagonal also need matching!

     f1=normMat{k1};

    for k2 = k1:m  % matching the two shapes with their feature by DP
        f2=normMat{k2};
        f3=normMat2{k2};
        [costmat] = weighted0_tar_cost(f1, f2);
        [costmats] = weighted0_tar_cost(f1, f3);
        
        %- MATCHING
        %- in current order
        [cvec1, match_cost1] = mixDPMatching_C(costmat, thre, num_start, search_step);
%         %- in reverse order
        [cvec2, match_cost2] = mixDPMatching_C(costmats, thre, num_start, search_step);

        %- get the best result
        Score(k1, k2) = min(match_cost1, match_cost2);
        Score(k2, k1) = Score(k1, k2);
%          Score(k1, k2) = match_cost1;
    end
    
    fprintf('Process line %i end \n', k1);
  
end
% TA2 = fix(clock); % end time
  save(strcat(dirname,' score.mat'), 'Score');   
end
