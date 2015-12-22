
This package includes all necessary source codes for shape matching with Hierarchical Projective Invariant Contexts. 

Some of these functions are implemented by directly using or slightly changing the source codes of Inner-Distance Shape Context (IDSC) provided online by Ling and Jacobs, and the source codes of height functions(HF) provided online by Xiang Bai.

The source codes are permitted for non-profit research use only.

To run the program, please put shape images into directory 'data' and use the command line below in Matlab: 

startHCNC

This program will output the similarity matrix 'Score' into file 'data score.mat', where Score(i,j) is the similarity between the ith shape and the jth shape read from 'data'.


For algorithm details, please refer to 
Qi Jia, Xin Fan, Yu Liu, Haojie Li, Zhongxuan Luo, and He Guo, "Hierarchical Projective Invariant Contexts for Shape Recognition", Pattern Recognition. 
(doi:10.1016/j.patcog.2015.11.003)
http://www.sciencedirect.com/science/article/pii/S0031320315004239
