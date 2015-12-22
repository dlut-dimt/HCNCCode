 
 function startHCNC
 
     dirname='data';%database directory
     
     disp('Calculating FCNs of original shapes');
     testMultiIntevals_CHN(0,'order1',dirname);
     testMultiIntevals_CHN(0,'order2',dirname);
     
     disp('Calculating FCNs of mirror shapes');
     testMultiIntevals_CHN(1,'order1',dirname);
     testMultiIntevals_CHN(1,'order2',dirname);
     
     disp('Generating descriptors');
     mergeValues(0,dirname);
     mergeValues(1,dirname);

     disp('Calculating similarities');
     shape_retrieval(dirname);
      
     %precision and recall 
%      classnum = 8;%%number of shape classes
%      viewnum = 5; %%instances in each class
%      calculatePR(dirname,classnum,viewnum);
 end