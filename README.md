# VesselAnalysis

`vesselanalysis.m` is a MATLAB pipeline for semi-automatic blood vessel analysis from 2D images.  

Pipeline:

- Batches through all images in a folder  
- Uses a single, user-defined global threshold 
- Segments vessels, fills holes, and bridges gaps  
- Skeletonizes the vessel network and splits it into individual segments  
- Samples multiple cross-sections perpendicular to each vessel segment  
- Outputs per-segment vessel diameter, length, and approximate area, plus QC overlays

---

## Overview

The script is designed for analysis of 2D vessel images (e.g., fluorescence). It provides:

- Semi-automatic control over segmentation via interactive threshold selection  
- Consistent scaling into physical units (µm and µm²)  
- Per-segment metrics: mean diameter, variance, length, and approximate 2D area  

---

## Requirements

- MATLAB  
- Image Processing Toolbox  
- Computer Vision Toolbox (for `insertText`, used in QC overlays)

---

## Workflow

1. **Run the script in MATLAB**

   Make sure `vesselanalysis.m` is on your MATLAB path, then run:
   ```matlab
   vesselanalysis

2. **Select the input folder**

   A dialog asks you to choose the folder containing your vessel images.
   Supported file types:
   *.jpg, *.jpeg, *.png, *.tif, *.tiff, *.bmp

3. **Choose an example image for thresholding**

   A list dialog shows all image names in the folder.
   Pick one representative image to use for setting the global threshold.

4. **Set threshold using the preview GUI**
   A GUI opens with:
    Left panel: original grayscale image
    Right panel: binary mask at the current threshold
    Adjust the slider until the vessels are well captured with minimal noise.
    Click “Use This Threshold” to apply it to all images.
    Click Cancel to abort the analysis.

5. **Enter the scale (pixels per micron)**
   A dialog prompts for px / µm (pixels per micron).

6. **Automatic batch processing**
   For each image, the script performs:
   A. Preprocessing & Segmentation
     a. Convert to grayscale (im2double)
     b. Threshold using your chosen value
     c. Ensure vessels are white (invert if necessary)
     d. Remove tiny specks (bwareaopen)
     e. Fill internal holes (imfill)
     f. Bridge small gaps (bwmorph with bridge)
     g. Apply morphological closing with a disk structuring element to fill weak gaps

   B. Skeletonization & Segmentation
     a. Skeletonize vessels (bwskel or bwmorph('skel'))
     b. Remove small spurs (bwmorph('spur'))
     c. Identify branch points and remove them
     d. Label individual skeleton segments (bwlabel)

   C. Cross-section Diameter Sampling
     a. For each skeleton segment (above a minimum length):
         Sample skeleton points at intervals (sampleStep)
         For each sampled point:
          Find local skeleton neighborhood (neighborRadiusPixels)
          Use PCA to determine the local vessel direction (tangent)
          Compute the normal vector (perpendicular direction)
          March along the normal in both directions until leaving the vessel mask
          Sum distances to get a diameter in pixels
         Convert all diameters to µm
         Compute segment length from skeleton pixel count (in µm)
         Estimate segment area using a rectangular approximation.
   
    D. QC Overlay Creation
       Generate an RGB overlay:
        Grayscale background
        Red: vessel edges (bwperim)
        Green: vessel skeleton
       Overlay segment IDs (yellow labels) using insertText at each segment centroid.

  7. **Results are saved to folder**
    A CSV file containing all measurement data
    A folder of QC images for visual inspection

