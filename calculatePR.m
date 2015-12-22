
function calculatePR(dirname,classnum,viewnum)
%%% calculate precision and recall 
clc;

simmatrix = load(strcat(dirname,' score.mat'));

simmatrix = simmatrix.Score;

imgnum = classnum * viewnum;
precisionmat = zeros(1,imgnum);
recallmat = zeros(1,imgnum);
%ccount=zeros(1,11);
all=0;
n=1;

for T =2:2:classnum*2%[2, 5, 8, 10, 15, 20, 25, 50]; 
    
for class = 1:classnum
    for view = 1:viewnum
        num = (class-1)*viewnum+view;
        querysim = simmatrix(num,:);
        [~,indexsim]= sort(querysim);  
        counters = 0;
        for i = 1:T
            left = (class-1)*viewnum + 1;
            right = class*viewnum;
            if indexsim(1,i)>=left && indexsim(1,i)<=right
                counters = counters + 1;             
            end
            
        end
        all=all+counters;
        precisionmat(1,num) = counters / T;
        recallmat(1,num) = counters / viewnum;
     %   fprintf('%d,%d,%d,%d\n',T,class,view,counters);
%         disp(precisionmat(1,num));
    end
end
    
% 
% save('precision.mat','precisionmat');
% save('recall.mat','recallmat');
%%% ????????
% disp('precision');
precisionAvg = mean(precisionmat);
% disp(precisionAvg);
% disp('recall');
recallAvg = mean(recallmat);
% disp('--------------');
disp([T,precisionAvg,recallAvg]);
PRv(n,:)=[T,precisionAvg,recallAvg];
n=n+1;
end
save(strcat(dirname,'PR.mat'),'PRv');
figure;
plot(PRv(:,3),PRv(:,2));

end %%% file end

