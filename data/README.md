# Folder `data`
The folder `data` contains the following items:

* `ref`: folder containing reference meshes of regular prisms.
  * `sidesX_diam1um.faces`: faces of a reference regular prism with X-sided bases (X = 4, 5, 6) and unitary size, according to the Polygon File Format (PLY) mesh format.
  * `sidesX_diam1um.header`: header of a reference PLY file for a regular prism with X-sided bases (X = 4, 5, 6) and unitary size.
  * `sidesX_diam1um.ply`: reference PLY file for a regular prism with X-sided bases (X = 4, 5, 6) and unitary size.
  * `sidesX_diam1um.vertices`: vertices of a regular prism with X-sided bases (X = 4, 5, 6) and unitary size, according to the PLY mesh format.
  * `sidesX_diamLum.vertices`: vertices of a regular prism with X-sided bases (X = 4, 5, 6) and size of L (L = 11, ..., 60 um), according to the PLY mesh format.
* `perturbed`: folder containing the perturbations of the regular prisms contained in `ref`, which are used as models of hepatocytes.
  * `sidesX_diamLum_fpert0.1_npertN.ply`: PLY file containing a mesh obtained by perturbing a reference regular prism. The generic `sidesX_diamLum_fpert0.1_npertN.ply` file corresponds to the N-th perturbation (N = 1, ..., 5) of the regular prism with X-sided bases (X = 4, 5, 6) and size L (L = 11, ..., 60 um). Perturbations are drawn independently for each mesh vertex from the normal distribution of standard deviation P = 0.1L.
  * `sidesX_diamLum_fpert0.1_npertN.vertices`: vertices contained in `sidesX_diamLum_fpert0.1_npertN.ply`.
* `perturbed_randomwalks`: folder that would contain the diffusion random walks. The total size of this folder would be on the order of 1TByte and is not included as it exceeds the maximum repository space in GitHub (10 GByte).

Note: PLY files can be rendered in 3D easily with free software such as [MeshLab](https://www.meshlab.net).


