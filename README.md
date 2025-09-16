# A computational framework and implementation for implicit priors in Bayesian inverse problems

#### [arXiv](https://arxiv.org/abs/2509.11781)

This repository contains supplementary materials to our manuscript.

## Notebooks
### Pedagogical examples
- [Simplest linear Bayesian inverse problem (Section 4.2-4.4)](simplest_linear/simplest_linear.ipynb), RLRTO
- [Simplest nonlinear Bayesian inverse problem (Section 4.5)](simplest_nonlinear/simplest_nonlinear.ipynb), Langevin

### Application case studies
- [Deblurring with constraints (Section 5.1)](deblurring/staircase.ipynb), RLRTO
- Various inverse problems with the Poisson equation (Section 5.2)
  - [Source term 1](pde_source/source_1d_03.ipynb) [2](pde_source/source_1d_1.ipynb) [3](pde_source/source_1d_3.ipynb), RLRTO
  - [Boundary value](pde_boundary_value/boundary_value_2d.ipynb), RLRTO
  - [Conductivity](pde_conductivity/Poisson_2D_MYULA_Part1.ipynb) together with [comparing different TV regularization strength](pde_conductivity/Poisson_2D_MYULA_Part2.ipynb), Langevin
- Image inpainting (Section 5.3)
  - [Wavelet denoiser](inpainting/inpainting_wavelet.ipynb), Langevin
  - [DnCNN denoiser](inpainting/inpainting_DnCNN.ipynb), Langevin

## Environment

Core libraries being used and their respective versions:
- [CUQIpy](https://github.com/CUQI-DTU/CUQIpy) 1.3.0
- [CUQIpy-FEniCS](https://github.com/CUQI-DTU/CUQIpy-FEniCS) 0.8.0
- [Scikit-image](https://github.com/scikit-image/scikit-image) 0.23.2
- [Deepinverse](https://github.com/deepinv/deepinv) 0.3.0
- [FEniCS](https://anaconda.org/conda-forge/fenics) 2019