# Wildcard LiDAR Pipeline

High-performance LiDAR processing pipeline designed for voxel-based forest structure analysis.

The pipeline runs in an HPC (PBS) environment, stages data to scratch storage, performs cloud homogenization, preprocessing, and computes structural metrics using R.

---

## Overview

The pipeline processes a single `.laz` file clipped by a polygon `.shp` and computes voxel-based forest metrics.

Main components:

- Bash orchestration (`WildcardPipeline.sh`)
- PBS job wrapper (`Run.sh`)
- Point cloud homogenization (shell scripts)
- Preprocessing and metric computation (R scripts)

---

## Execution Environment

- PBS scheduler
- Local scratch storage (`SCRATCHDIR`)
- R environment
- `s5cmd` for S3 transfers
- AWS-compatible S3 endpoint

---

## Input Parameters

`WildcardPipeline.sh` expects:
```
<path_to_laz> <path_to_polygon_shp> <nOfElements> <voxelSize> <site_id> <coordinates>
```

Key parametrs:
- **LAZ file** – input point cloud
- **Polygon SHP** – clipping geometry
- **nOfElements** – voxel subdivision parameter
- **voxelSize** – voxel resolution (e.g. 0.2)
- **coordinates** – CRS (e.g. EPSG:5514)

---

## Processing Steps

1. **Scratch Staging**
   - Input data copied to node-local scratch
   - Temporary directories isolated per job

2. **Cloud Homogenization**
   - `CloudHomogenization.sh`
   - Optional geometric nearest-neighbor variant
   - Produces density-normalized cloud

3. **Preprocessing (R)**
   - `Preprocessing.R`
   - `Preprocessing_long.R`
   - Voxelization and structural preparation

4. **Metric Computation (R)**

   Implemented metrics:

   - `MaxHeight98`
   - `HeightVariabilityMetrics`
   - `VerticalStructureMetrics`
   - `CanopyDensity`
   - `CanopyVolumeMetrics`
   - `GappinesMetrics`

5. **Results Export**
   - Processed outputs copied back from scratch
   - Scratch cleaned on successful completion
   - Logs written per run

---

## HPC Job Wrapper

`Run.sh`:

- Pulls LAZ data from S3 (`wp2` bucket)
- Derives `SITE_ID` from filename
- Submits PBS job
- Configurable:
  - Voxel size
  - Number of elements
  - CRS

---

## Logging & Fault Handling

- Timestamped log per run
- Execution timing tracking
- Automatic copy-back on exit
- Scratch preserved on failure for debugging

---

## Output

Results are written to:
```
Processed/<fileName>/
Results/<fileName>/
```


Outputs include voxel-level structural forest metrics derived from homogenized point clouds.

---

## Purpose

Designed for large-scale forest structural analysis using voxel-based LiDAR metrics in distributed HPC environments.