# Multivariate-Intrinsic-Local-Polynomial-Regression-on-Isometric-Riemannian-Manifolds


![relation](https://user-images.githubusercontent.com/102549261/234109560-67556ca0-12ef-43db-9e90-044f2345ac68.png)


This is a toolkit for calculating multivariate intrinsic local polynomial regression on matrix Riemannian manifolds with Euclidean pullback metrics. It can be used to perform nonparametric regression analysis on data from symmetric positive definite (SPD) manifolds, such as diffusion tensor imaging (DTI) data.

# Installation
To use this toolkit, you need to have MATLAB R2022a or later installed on your computer. You also need to download and install Manopt, a MATLAB toolbox for optimization on manifolds https://github.com/NicolasBoumal/manopt. In particular, you need the function positive_definite_karcher_mean.m from Manopt to compute the Karcher mean of SPD matrices.

# Contents

This toolkit contains the following folders and files:

CholeskyManifold: This folder contains the program CholeskyManifold.m that creates the structure of the SPD manifold with the log-Cholesky Riemannian metric. This structure allows the calculation of Riemannian operations such as distance, geodesic, exponential map, and logarithmic map on the SPD manifold.

DataGeneration: This folder contains the program spd_high_dimensional.m that generates true simulated data for Monte Carlo simulations. The data are generated from a multivariate intrinsic local polynomial model on the SPD manifold with a specified signal-to-noise ratio (SNR).

Figures: This folder contains codes that serve as examples of how to use the regression method and the code used to generate the figures of the paper and to conduct Monte Carlo simulations. The figures illustrate the performance of the proposed method on simulated and real DTI data sets.

MILPR: This folder contains the implementation of the multivariate intrinsic local polynomial regression (MILPR) estimator for a Euclidean pullback metric (EPM). The main function is MILPR.m, which takes as input a sample of SPD matrices and their corresponding covariates, and returns the estimated regression function and its derivatives at any given point on the manifold. The function also allows the user to specify the bandwidth, the degree of the polynomial, and the type of EPM to use.

Noise_Model: This folder contains the implementation of a Riemannian log-normal model to add noise to a sample of SPD data. The function High_Dimension_Riem_Noise.m takes as input a sample of SPD matrices and a noise level parameter, and returns a noisy sample of SPD matrices with the same mean and covariance structure as the original sample.

Riem t-SNE: This folder contains the implementation of log-Cholesky Riemannian t-distributed stochastic neighbor embedding (Rie-SNE) on the SPD manifold. The function riem_tsne.m takes as input a sample of SPD matrices and returns a low-dimensional embedding that preserves the local structure of the data.

# Examples

To illustrate how to use this toolkit, we provide two examples: one on simulated DTI data and one on high-dimensional data.

DTI simulated data 

![DTI](https://user-images.githubusercontent.com/102549261/234109220-5f8a778d-4d34-452a-966a-84e4c860a677.png)

High Dimensional preformance 

![trade_off](https://user-images.githubusercontent.com/102549261/234109349-70bee79b-b862-4305-b51b-e676c6e32630.png)

Riem t-SNE

![Rie_SNE](https://user-images.githubusercontent.com/102549261/234109427-8aacbbd5-bd47-4ab9-9920-bfd7c28b2756.png)
