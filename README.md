# jfreg2: Joint fMRI Registration (Version 2)

`jfreg2` is a Python script that is intended to facilitate analysis of
functional magnetic resonance imaging (fMRI) datasets. `jfreg2` provides a
pipeline for registering at least one fMRI time series to standard space. The
pipeline includes motion correction, distortion correction using field maps,
functional-structural registration, and registration of the structural dataset
to standard space. Transformations and warps from the individual steps are
combined into a single transformation and warp that brings the fMRI time series
into standard space.

The `jfreg2` pipeline uses the FMRI Software Library (FSL).

# License

Copyright (c) 2021, Jeffrey M. Engelmann

`jfreg2` is released under the revised (3-clause) BSD license.
For details, see [LICENSE.txt](LICENSE.txt).

The [FMRI Software Library](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) (FSL) is a
comprehensive library of analysis tools for brain imaging data that is
developed by the Analysis Group at the Wellcome Centre for Integrative
Neuroimaging at the University of Oxford. FSL is released for non-commercial
use under an open-source
[license](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence).
