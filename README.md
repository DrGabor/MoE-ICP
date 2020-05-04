# MoE-ICP

The backbone MATLAB code of MoE-ICP is illustrated in this repository. Currently, the variational inference (VI) in VBMoEFun.m is missing owing to the paper of MoE-ICP is still under review, and we will release it as soon as the paper is accepted. Albeit, the readers can inspect the pipeline and details of MoE-ICP as follows:

1. The residual errors and information matrices of plane-to-plane distance are illustrated in compute_statisticsFun.m.

2. The iteratively re-weighted least-squares (IRLS) is illustrated in IRLSFun.m.

3. The vector-to-vector operations with the help of MATLAB symbolic toolbox are illustrated in CalHbFun.m. 

For background of variational inference, the manuscript titled "A Brief Note on Variational Inference" is provided. 
