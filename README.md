# MoE-ICP

The backbone MATLAB code of MoE-ICP [1] is illustrated in this repository. The MoE-ICP employs a variational inference framework and includes widely-used metric distances such as point-to-point, point-to-plane and plane-to-plane, which can be viewed as an extension of MiNoM [2]. Currently, the variational inference (VI) in VBMoEFun.m is missing owing to the review process of MoE-ICP, and we will release it as soon as the change of status. Albeit, the readers can inspect the pipeline and details of MoE-ICP as follows:

1. The residual errors and information matrices of plane-to-plane distance are illustrated in compute_statisticsFun.m.

2. The iteratively re-weighted least-squares (IRLS) is illustrated in IRLSFun.m.

3. The vector-to-vector operations with the help of MATLAB symbolic toolbox are illustrated in CalHbFun.m. 

For background of variational inference, the manuscript titled "A Brief Note on Variational Inference" is provided. 

[1] Wang, D., MoE-ICP, under review. 

[2] Wang, D., Xue, J., Tao, Z., Zhong, Y., Cui, D., Du, S. and Zheng, N., 2018, October. Accurate Mix-Norm-Based Scan Matching. In 2018 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) (pp. 1665-1671). IEEE.
