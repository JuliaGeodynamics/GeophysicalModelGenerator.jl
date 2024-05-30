# pTatin

We provide a few routines that allow exporting GMG input data sets to a format that can be read by the 3D geodynamic modelling software [pTatin](https://bitbucket.org/dmay/ptatin-total-dev.git/src), which is an open-source Cartesian code to perform crustal and lithospheric-scale simulations. 

The `Q1Data` format can be saved directly to pTatin.

```@docs
GeophysicalModelGenerator.write_pTatin_mesh
GeophysicalModelGenerator.swap_yz_dims
```