# First-Rank Voronoi Tessellation in MATLAB

This repository contains a MATLAB implementation of the **First-Rank Voronoi Tessellation** method, as described by Levet et al. (2015) 
in their publication:

📄 Levet, F. et al. (2015)
SR-Tesseler: a method to segment and quantify localization-based super-resolution microscopy data 
Nature Methods, 12(11), 1065–1071. [https://doi.org/10.1038/nmeth.3579]

🧠 About

This code computes local point densities from a 2D localization dataset using Voronoi tessellation, focusing on **first-rank neighbors** 
(i.e., immediate neighboring regions). It statistically compares observed densities with those from a random distribution to define a threshold, 
allowing identification of significantly high-density regions.

This approach is particularly useful in **single-molecule localization microscopy (SMLM)**, **spatial biology**, and **point-pattern analysis**.

📌 Features

- Computes Voronoi tessellation and area-based density
- Uses first-rank neighborhood statistics for robust density estimation
- Compares with random spatial distribution to identify significant clusters
- MATLAB-compatible and easy to integrate into analysis pipelines

🚀 Getting Started

### Requirements

- MATLAB R2018b or newer (may work with older versions)
- Delaunay/Voronoi functions (`delaunayTriangulation`, `voronoiDiagram`)
- The helper function: `AreaVoronoiRegions_1stRank.m` 

### Usage

```matlab
[Data_voronoi, Data_voronoi_high, threshold] = VoronoiTessel_1st_rank_density_th(XY, thL);
```
INPUTS:

  XY   - Nx2 or Nx3 matrix with 2D coordinates (and optional intensity)
  thL  - Threshold level (quantile) for identifying high-density points
OUTPUTS:

  Data_voronoi       - Voronoi-based metrics for input points:
                      [X Y (optional Z) density norm. rank score]
  Data_voronoi_high  - Subset of Data_voronoi above the density threshold
  threshold          - Density threshold (in log10 scale) for filtering
