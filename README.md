# A computational framework and implementation for implicit priors in Bayesian inverse problems

#### [arXiv to-be-filled](https://arxiv.org/abs/xxxx)

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
  - [Conductivity](pde_myula/Poisson_2D_MYULA_Part1.ipynb) together with [comparing different TV regularization strength](pde_myula/Poisson_2D_MYULA_Part2.ipynb), Langevin
- Image inpainting (Section 5.3)
  - [Wavelet denoiser](inpainting/inpainting_wavelet.ipynb), Langevin
  - [DnCNN denoiser](inpainting/inpainting_DnCNN.ipynb), Langevin
