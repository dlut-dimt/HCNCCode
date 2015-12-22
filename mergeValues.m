function mergeValues(morror,dirname)
%%Generat HCNC descriptors with FCN
clc;
data1 = load(strcat(dirname,'order1 img',num2str(morror),'.mat'));
data1 = data1.resultdata;
data2 = load(strcat(dirname,'order2 img',num2str(morror),'.mat'));
data2 = data2.resultdata;

delete(strcat(dirname,'order1 img',num2str(morror),'.mat'));
delete(strcat(dirname,'order2 img',num2str(morror),'.mat'));

num = length(data1); %%% 

resultdata = cell(1,num);
for i=1:num
    feature1=data1{1,i};
    feature2=data2{1,i};
    tempfeature = feature1 ./ feature2;
%     tempfeature = feature1;
    rows=size(tempfeature,1);
    cols=size(tempfeature,2);
    
    for m=1:rows
       for n=1:cols   
          if isinf(abs(tempfeature(m,n))) || isnan(abs(tempfeature(m,n))) 
              tempfeature(m,n)=0;   
              continue;
          elseif  abs(tempfeature(m,n))>1
                tempfeature(m,n) = sign(tempfeature(m,n))*1;   
          end 
       end
    end
    
    resultdata{1,i} = tempfeature;
    
    disp(i);   
end

save(strcat(dirname,' img',num2str(morror),'.mat'),'resultdata');

end %%% file end