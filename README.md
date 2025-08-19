# A computational framework and implementation for implicit priors in Bayesian inverse problems

#### [arXiv to-be-filled](https://arxiv.org/abs/xxxx)

This repository contains supplementary materials to our manuscript.

## Notebooks
### Pedagogical examples
- [Simplest linear Bayesian inverse problem](simplest_linear/simplest_linear.ipynb), linear, RLRTO
- [Simplest nonlinear Bayesian inverse problem](simplest_nonlinear/simplest_nonlinear.ipynb), nonlinear, Langevin

### Application case studies
- [Deblurring with constraints](showcase_regularizedGaussian/showcase_regGauss.ipynb), linear, RLRTO
- Various inverse problems with the Poisson equation
  - [Source term](pde_source/right_hand_side_1d_demo.ipynb), linear, RLRTO
  - [Boundary value](pde_boundary_value/boundary_value_demo.ipynb), linear, RLRTO
  - [Conductivity](pde_myula/Poisson_2D_MYULA_short.ipynb), nonlinear, Langevin
- [Image inpainting](inpainting/inpainting.ipynb), linear, Langevin
