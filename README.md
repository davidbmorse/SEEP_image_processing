# SEEP_image_processing
code to reproduce image processing and analysis employed by SEEP method

## Data Structure
#### \Data
##### c1 : all c1 (raw imaging data in plate form)
##### c1_working_plane : c1 at working plane
##### c1_working_plane_mask : spheroid mask of c1 at working plane c1_working_plane_normalized : normalized c1
##### c1_working_plane_shell : shells of c1

## Step 1: Compute bottom/working plane
#### \Code
S0_bottom_plane_computation.m
All c1 at working plane are saved in \c1_working_plane

## Step 2: Spheroid segmentation
Open source code (Python) https://github.com/matterport/Mask_RCNN
All spheroid mask of c1 are saved in \c1_working_plane_mask.

## Step 3: Normalization and shell extraction
#### \Code
S1_shell_extraction.m
