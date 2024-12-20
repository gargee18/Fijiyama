/*
 * This Java class, CropWithRoi, is designed to crop specific regions from a series of images based on predefined coordinates.
 *  
 * Methods:
 * - readCropCoordsAndDimensionsFromCSVFile: Reads crop coordinates and dimensions from CSV files, applies them to crop the images, and saves the cropped images.
 * - crop: Crops an image based on the specified region of interest (ROI) and returns the cropped image.
 * 
 * Input Parameters:
 * - year: The year of the data, used in directory paths and file names.
 * - mainDir: The main directory where the input and output directories are located.
 * - specimen: An array of specimen identifiers used to locate the specific images.
 * - inputDir: The directory containing the raw images to be cropped.
 * - outputDir: The directory where the cropped images will be saved.
 * - pathToCSVFileforCropCoords: The path to the CSV file containing crop coordinates.
 * - pathToCSVFileforCropSize: The path to the CSV file containing crop dimensions.
 * 
 * Output Parameters:
 * - imgCropped: The cropped image with the specified dimensions.
 * 
 */

package io.github.rocsg.fijiyama.gargeetest.ceps;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Roi;
import io.github.rocsg.fijiyama.common.VitimageUtils;

public class CropImagesWithRoi {
    public static String year = "2024";

    public static void main(String[] args) {
        // Open an image using ImageJ
        ImageJ ij = new ImageJ();

        // Define the specimen identifiers and directories
        String[] specimen = GeneralUtils.specimenListXR;
        String mainDir = "/home/phukon/Desktop/";
        String outDir = mainDir + "cropdata" + year + "/";
        String inputDir = mainDir + "XR_" + year + "_raw/";
        String pathToCSVFileforCropCoords = "/home/phukon/Desktop/" + year + "_cropcoords.csv";
        String pathToCSVFileforCropSize = "/home/phukon/Desktop/croppedImageDimensions" + year + ".csv";

        // Perform cropping based on ROI
        readCropCoordsAndDimensionsFromCSVFile(specimen, inputDir, outDir, pathToCSVFileforCropCoords, pathToCSVFileforCropSize);
    }

    // Reads crop coordinates and dimensions from CSV files, applies them to crop the images, and saves the cropped images
    public static ImagePlus readCropCoordsAndDimensionsFromCSVFile(String[] specimen, String inputDir, String outputDir, String pathToCSVFileforCropCoords, String pathToCSVFileforCropSize ) {
        ImagePlus imgCropped = null;
        // Read CSV data for crop coordinates and crop dimensions
        String[][] cropCoordsData = VitimageUtils.readStringTabFromCsv(pathToCSVFileforCropCoords);
        String[][] sizeCropData = VitimageUtils.readStringTabFromCsv(pathToCSVFileforCropSize);

        int numberOfSpecimens = specimen.length;
        int numberOfCoordinates = 3; // x,y,z

        int[][] trCoords = new int[numberOfSpecimens + 1][numberOfCoordinates + 1];
        int[][] cropImgdim = new int[numberOfSpecimens + 1][numberOfCoordinates + 1];

        for (int i = 1; i < numberOfSpecimens + 1; i++) {
            trCoords[i][1] = (int) Double.parseDouble(cropCoordsData[i][1]); // x
            trCoords[i][2] = (int) Double.parseDouble(cropCoordsData[i][2]); // y
            cropImgdim[i][1] = (int) Double.parseDouble(sizeCropData[i][1]); // width
            cropImgdim[i][2] = (int) Double.parseDouble(sizeCropData[i][2]); // height
            Roi roi = new Roi(trCoords[i][1], trCoords[i][2], cropImgdim[i][1], cropImgdim[i][2]);
            ImagePlus img = new ImagePlus(inputDir + "CEP_" + specimen[i - 1] + "_" + year + "_XR.tif");
            imgCropped = crop(img, roi, "stack");
            IJ.saveAsTiff(imgCropped, outputDir + "CEP_" + specimen[i - 1] + "_" + year + "_XR_crop.tif");
        }
        // Print crop coordinates
        System.out.println("Crop coordinates:");
        for (int i = 1; i < numberOfSpecimens + 1; i++) {
            System.out.println("Specimen " + specimen[i - 1] + ": x = " + trCoords[i][1] + ", y = " + trCoords[i][2]
                    + ", z = " + trCoords[i][3]);
        }
        // Print crop dimensions    
        System.out.println("Crop Dimensions:");
        for (int i = 1; i < numberOfSpecimens + 1; i++) {
            System.out.println("Specimen " + specimen[i - 1] + ": x = " + cropImgdim[i][1] + ", y = " + cropImgdim[i][2]
                    + ", z = " + cropImgdim[i][3]);
        }
        return imgCropped;
    }
    // Crops an image based on the specified region of interest (ROI) and returns the cropped image
    public static ImagePlus crop(ImagePlus img, Roi roi, String options) {
        ImagePlus cropImps = new ImagePlus();
        Roi cropRoi = roi;
        String name = cropRoi.getName();  
        // Set the slice if options contain "slice" and image has multiple slices
        if (options.contains("slice") && img.getStackSize() > 1) {
            int position = cropRoi.getPosition();
            img.setSlice(position); // no effect if roi position is undefined (=0), ok
        }
        // Set ROI and crop the image
        img.setRoi(cropRoi);
        ImagePlus cropped = img.crop(options);
        // Reads crop coordinates and dimensions from CSV files, applies them to crop the images, and saves the cropped images
        if (cropRoi.getType() != Roi.RECTANGLE) {
            Roi cropRoi2 = (Roi) cropRoi.clone();
            cropRoi2.setLocation(0, 0);
            cropped.setRoi(cropRoi2);
        }
        // Set the title for the cropped image
        String name2 = IJ.pad(1, 3) + "_" + img.getTitle();
        cropped.setTitle(name != null ? name : name2);
        cropped.setOverlay(null);
        cropImps = cropped;
        return cropImps;
    }

}