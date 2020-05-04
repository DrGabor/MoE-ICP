# MoE-ICP

The backbone MATLAB code of MoE-ICP [1] is illustrated in this repository, and the registration can be seen by runing the file of main_demo.m. The MoE-ICP employs a variational inference framework and includes widely-used metric distances such as point-to-point, point-to-plane and plane-to-plane, which can be viewed as an extension of MiNoM [2]. If you cannot run the code, some missing functions like Loc2Glo() and SkewFun() are in the repository of DrGabor/LiDAR/CommonFunctions. The implementation of VBMoEFun() and IRLSFun() will be released in the near future. Albeit, the readers can inspect the pipeline and details of MoE-ICP as follows:

1. The residual errors and information matrices of plane-to-plane distance are illustrated in compute_statisticsFun() and computeCovFun().

2. The vector-to-vector operations with the help of MATLAB symbolic toolbox are illustrated in CalHbFun(). 

For background of variational inference, the manuscript titled "A Brief Note on Variational Inference" is provided. 

[1] Wang, D., MoE-ICP, under review. 

[2] Wang, D., Xue, J., Tao, Z., Zhong, Y., Cui, D., Du, S. and Zheng, N., 2018, October. Accurate Mix-Norm-Based Scan Matching. In 2018 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) (pp. 1665-1671). IEEE.
