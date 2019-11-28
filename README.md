# REKINDLE
REKINDLE algorithm as published in [Tax et al., MRM 2015](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.25165)

How to use:

1) First set some parameters:
```
par.iter1 = 5;
par.iter2 = 20;
par.con = 1.0000e-03;
kappa = 6;
DKI = false;
```
This implementation has the option to fit DTI and DKI, by setting the `DKI` parameter to false or true, respectively.

2) Load some data, for example using [this NIFTI reader](https://uk.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image):
```
nii = load_untouch_nii('DWI.nii');
mask = load_untouch_nii('mask.nii');
bval = load('DWI.bval'); bval = bval'; % optional transpose to make N x 1
g = load('DWI.bvec'); bvec = bvec'; % optional transpose to make N x 3
```

3) Run REKINDLE iterative reweighting:
```
res = REKINDLE(nii.img, bval, g, logical(mask.img), DKI, par);
```

4) Run REKINDLE outlier rejection and final fit, adapt kappa if necessary:
```
[outlier, DWIB0, DT, varargout] = REKINDLE_Outlier_Removal_Fit(nii.img, bval, g, logical(mask.img), DKI, res, kappa);
```

5) Save residuals for future adaptation of kappa:
```
nii.img = res; save_untouch_nii(nii, 'res.nii');
```
