package io.github.rocsg.fijiyama.gargeetest;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ij.IJ;
import ij.ImagePlus;

public class ArrangeDCMfiles {
    public static String year = "2023";
    public static String specimen = GeneralUtils.specimenListMRI[7];
    public static String mod = GeneralUtils.modalities[0];
    public static String dir = "/home/phukon/Desktop/MRI/02_ceps/2023_CEPS_suivi_clinique/0_RAW_data_MRI_2023/"+specimen+"/"+mod+"/";
    public static String dirOutput = "/home/phukon/Desktop/MRI/02_ceps/2023_CEPS_suivi_clinique/0_RAW_data_MRI_2023/"+specimen+"/"+mod+"_ordered";

    public static void main(String[] args) { 
        List<FileSliceLocationPair> filesArranged = filesInRightOrder(dir,mod,GeneralUtils.paramsWanted[5]);
        saveFiles(filesArranged, dirOutput, specimen);
        // String[] arrFiles = filesArranged.stream()
        //                        .map(pair -> "Slice Location for " + pair.getFile().getName() + ": " + pair.getSliceLocation())
        //                        .toArray(String[]::new);
  
    }


    public static String[] getParamsFromDCM(ImagePlus imgPar,String paramsWanted){
        if(imgPar==null) System.out.println("Warning in TransformUsingMetadata: image used for parameters detection is detected as null image. Computation will fail");
        String paramGlob=imgPar.getInfoProperty();
        String[]paramLines=paramGlob.split("\n");
        String[]ret=new String[paramsWanted.length()];
        for(int i=0;i<paramsWanted.length();i++) {
            ret [i]="NOT DETECTED";
            for(int j=0;j<paramLines.length ;j++) {
                if(paramLines[j].split(": ")[0].indexOf(paramsWanted)>=0)ret[i]=paramLines[j].split(": ")[1];
            }
        }
        return ret; 
	}

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
            // for (FileSliceLocationPair pair : fileList) {
            //     System.out.println("Slice Location for " + pair.getFile().getName() + ": " + pair.getSliceLocation());
            // }
            return fileList;   
        } else {
            System.out.println("No DICOM files found in the directory: " + inputDirectoryOfDCMFiles);
        }
        return null;
    }

    public static void saveFiles(List<FileSliceLocationPair> fileList, String dirOutput, String specimen){

        for (FileSliceLocationPair pair : fileList) {


            File sourceFile = pair.getFile();
            File destFile = new File(dirOutput,specimen +"_slicenumber_" + (fileList.indexOf(pair) + 1)+ "_"+pair.getFile().getName() + ".dcm");
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


