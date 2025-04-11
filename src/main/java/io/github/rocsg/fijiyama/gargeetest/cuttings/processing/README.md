# Cuttings T1 T2 Data Registration Pipeline

This repository contains the image registration pipeline for aligning and preparing T1 and T2 MRI image sequences for downstream analysis using **Fijirelax**. The pipeline is implemented within the **Fijiyama** framework.

---

## 📂 Overview

This pipeline processes multi-sequence MRI data through a combination of manual, automatic, and dense registration techniques, followed by normalization. It ensures spatial and intensity consistency across all images and timepoints before feeding them to the **Fijirelax** tool.

---

## 🧬 Pipeline Steps

### 1. **Input Check**
- **Input:** Raw image (can be a hypermap or regular sequence)
- If hypermap:
  - Select sequence where `TR = 2400`, `TE = 12`
- If not:
  - Locate folder `TR002400` and open file `TE000012`

---


### **TODO: Need to Automate Step 2-5, As of now the registration is done using the fijiyama gui (serial registration)**


### 2. **Manual Registration**
- Manually register T1, T2, T3, and T4 of the selected sequence
- **Output:** `TR matrix (1)`

---

### 3. **Auto Registration**
- **Input:** Raw sequences + `TR matrix (1)` × 4 (timepoints)
- Automatically aligns sequences
- **Output:** `TR matrix (2)`

---

### 4. **Dense Registration**
- **Input:** Raw sequences + `TR matrix (1)` + `TR matrix (2)` × 4
- Refines alignment
- **Output:** `TR matrix (3)`

---

### 5. **Inoc Alignment**
- **Input:** Raw sequences + all TR matrices (1–2) × 4 + coordinate CSV file (Dense not applied yet)
- Applies additional alignment based on inoculation points
- **Output:** `TR matrix (4)`

---

### 6. **Apply Transformations**
- **Input:** 19 channels × 4 timepoints

This step applies a transformation matrix to each image in a temporal sequence (e.g. MRI scans across time points). The matrix is typically computed in a previous registration step. Here's what the class does:

- Loads the raw image for each timestamp.

- Reads the corresponding transformation matrix file.

- Applies the transformation to the sequence.

- Saves the registered (aligned) sequence to disk.

- This ensures that all images across time are spatially aligned for further processing.


- **Output:** Registered hyperstack `(c = 19, t = 4)`

---

### 7. **Normalization**
This step normalizes the intensity values of T1/T2 MRI sequences to ensure consistency across channels and time points.

- Loads registered T1/T2 image sequences for a given specimen.

- Iterates through each channel and normalizes slices based on a reference “capillary” intensity. (Mean cap val file available on drive)

- Reconstructs normalized channels into a combined hyperstack.

- Applies visualization enhancements like LUTs and intensity scaling.

- Adjusts spatial calibration units (e.g., voxel size in mm).

- Performs **sigma correction** by updating image metadata using predefined normalization factors.

- Saves the normalized and corrected multi-channel image as a TIFF file.

- **Output:** Hyperstack with pixel intensities ranging from `0 to 1`

---

### 8. **Output**

This step builds a final "HyperMap" that combines all previously normalized and registered T1/T2 sequences over multiple timepoints.

- Iterates through each timepoint in the experiment.

- Loads the normalized T1/T2 MRI data corresponding to each timepoint.

- Concatenates these into a single 5D hyperstack (channels × slices × frames).

- Converts the concatenated stack into a proper hyperstack format using ImageJ utilities.

- Saves the final HyperMap as a `.tif` file for downstream analysis or visualization.

- Feed the final normalized hyperstack into **Fijirelax** for analysis.

---


