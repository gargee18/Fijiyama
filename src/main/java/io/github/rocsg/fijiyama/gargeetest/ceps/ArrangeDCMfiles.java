/*
 * This Java class, ArrangeDCMfiles, is designed to organize and manage XR/MRI DICOM files.
 * The class reads DICOM files from a specified input directory, extracts slice location parameters from the metadata of each file,
 * and arranges them in the correct order based on these parameters. The ordered files are then saved into a new output directory,
 * with their filenames appended by slice numbers for clarity.
 * 
 * Methods:
 * - filesInRightOrder: Reads DICOM files from the input directory, extracts and sorts them based on slice location parameters.
 * - getParamsFromDCM: Extracts specified parameters from the metadata of a DICOM image.
 * - saveFiles: Saves the ordered DICOM files into the output directory, renaming them to include slice numbers.
 * 
 * This program utilizes the ImageJ library for handling DICOM image processing and metadata extraction.
 * 
 * Usage:
 * 1. Set the `year`, `specimen`, `mod`, `dir`, and `dirOutput` variables to appropriate values.
 * 2. Run the program to process and organize the DICOM files.
 */

package io.github.rocsg.fijiyama.gargeetest.ceps;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ij.IJ;
import ij.ImagePlus;

public class ArrangeDCMfiles {
     // Define the year of the data
    public static String year = "";

     // Define the specimen and modality from the GeneralUtils class
    public static String specimen = GeneralUtils.specimenListMRI[7];
    public static String mod = GeneralUtils.modalities[0];

    // Define the input and output directories
    public static String dir = "";
    public static String dirOutput = "";

    public static void main(String[] args) {
        // Get the list of files arranged in the correct order
        List<FileSliceLocationPair> filesArranged = filesInRightOrder(dir, mod, GeneralUtils.paramsWanted[5]);
        // Save the ordered files to the output directory
        saveFiles(filesArranged, dirOutput, specimen);
    }

    // Extracts specified parameters (e.g., slice location) from the DICOM image metadata
    public static String[] getParamsFromDCM(ImagePlus imgPar, String paramsWanted) {
        if (imgPar == null)
            System.out.println(
                    "Warning in TransformUsingMetadata: image used for parameters detection is detected as null image. Computation will fail");
        String paramGlob = imgPar.getInfoProperty();
        String[] paramLines = paramGlob.split("\n");
        String[] ret = new String[paramsWanted.length()];
        for (int i = 0; i < paramsWanted.length(); i++) {
            ret[i] = "NOT DETECTED";
            for (int j = 0; j < paramLines.length; j++) {
                if (paramLines[j].split(": ")[0].indexOf(paramsWanted) >= 0)
                    ret[i] = paramLines[j].split(": ")[1];
            }
        }
        return ret;
    }
    // Reads DICOM files from the input directory, extracts and sorts them based on slice location parameters
    public static List<FileSliceLocationPair> filesInRightOrder(String inputDirectoryOfDCMFiles, String modality, String parameters) {
        File[] dcmFiles = new File(inputDirectoryOfDCMFiles).listFiles();
        if (dcmFiles != null && dcmFiles.length > 0) {
            List<FileSliceLocationPair> fileList = new ArrayList<>();
            for (File dcmFile : dcmFiles) {
                ImagePlus img = IJ.openImage(dcmFile.getPath());
                if (img != null) {
                    String[] sliceLocation = getParamsFromDCM(img, parameters);
                    if (sliceLocation != null && sliceLocation.length > 0) {
                        fileList.add(new FileSliceLocationPair(dcmFile, Double.parseDouble(sliceLocation[0])));
                    } else {
                        System.out.println("No slice location found for " + dcmFile.getName());
                    }
                } else {
                    System.out.println("Failed to open DICOM file: " + dcmFile.getName());
                }
            }
            // Sort the list based on slice location values
            Collections.sort(fileList);
            // Display the sorted files
            for (FileSliceLocationPair pair : fileList) {
            System.out.println("Slice Location for " + pair.getFile().getName() + ": " + pair.getSliceLocation());
            }
            return fileList;
        } else {
            System.out.println("No DICOM files found in the directory: " + inputDirectoryOfDCMFiles);
        }
        return null;
    }
    // Saves the ordered DICOM files into the output directory, renaming them to include slice numbers
    public static void saveFiles(List<FileSliceLocationPair> fileList, String dirOutput, String specimen) {
        for (FileSliceLocationPair pair : fileList) {
            File sourceFile = pair.getFile();
            File destFile = new File(dirOutput, specimen + "_slicenumber_" + (fileList.indexOf(pair) + 1) + "_"
                    + pair.getFile().getName() + ".dcm");
            try {
                Files.copy(sourceFile.toPath(), destFile.toPath());
                System.out.println("Copied: " + sourceFile.getName() + " to " + destFile.getName());
            } catch (IOException e) {
                System.out.println("Failed to copy file: " + sourceFile.getName());
                e.printStackTrace();
            }
        }
    }

    // Helper class to hold file and slice location pair
    static class FileSliceLocationPair implements Comparable<FileSliceLocationPair> {
        private File file;
        private double sliceLocation;
        public FileSliceLocationPair(File file, double sliceLocation) {
            this.file = file;
            this.sliceLocation = sliceLocation;
        }
        public File getFile() {
            return file;
        }
        public double getSliceLocation() {
            return sliceLocation;
        }
        @Override
        public int compareTo(FileSliceLocationPair other) {
            return Double.compare(this.sliceLocation, other.sliceLocation);
        }
    }
}