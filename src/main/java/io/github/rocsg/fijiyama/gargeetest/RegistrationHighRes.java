/*
 * 
 * This Java program implements a registration pipeline for aligning high-resolution images acquired in different years
 * of the same biological specimens. The registration process involves computing transformation matrices to align images
 * captured in consecutive years and generating a high-resolution registered image.
 * 
 * The registration pipeline consists of the following steps:
 * 1. Load the reference and moving images corresponding to consecutive years.
 * 2. Compute transformation matrices to align the moving image to the reference image.
 * 3. Combine the transformation matrices to produce a single transformation that accounts for changes in both translation and rotation.
 * 4. Apply the combined transformation to the moving image to generate a high-resolution registered image.
 * 5. Save the registered image and the transformation matrix for further analysis and visualization.
 * 
 * This code is specifically designed to process high-resolution images of biological specimens captured in consecutive years.
 * It provides a comprehensive registration solution to ensure accurate alignment and analysis of temporal changes in the specimens.
 * 
 * Input:
 * - Data of biological specimens captured in different years.
 * 
 * Output:
 * - High-resolution moving image after registration.
 * - Transformation matrix representing the registration transformation.
 * - The computed matrices and images are saved in a specific directory (that is hard-coded)
 * 
 * Dependencies:
 * - ImageJ library for image processing and manipulation.
 * - Custom library of Fijiyama (io.github.rocsg.fijiyama.common.VitimageUtils) for common utility functions.
 * 
 * Usage:
 * - Modify the 'mainDir', 'year1', and 'year2' variables to specify the directory and years of image data.
 * - Compile and execute the program to perform registration and generate high-resolution registered images.
 * 
 * Author: Gargee PHUKON, Romain FERNANDEZ
 * 
 */



package io.github.rocsg.fijiyama.gargeetest;

//java import
import java.util.ArrayList;

//specific libraries
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;

//my libraries
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.registration.ItkTransform;



public class RegistrationHighRes {
    public static String mainDir="/home/phukon/Desktop/registration/";// Path to the directory consisting the data
    public static String year1="2022";
    public static String year2="2023";


    //Main input points
    public static void main(String[] args) {
            
        // Specify the single specimen here or set it to null to get the entire list
        String singleSpecimen = "1191"; 
        
        // Build a list of all specimen names (with respect with the folder s name)
        ArrayList<String> specimenList = getSpecimen(singleSpecimen);

        //streamLineTheVerificationOfThatEverythingProbablyWentGoodAndThRegistrationAreAccurate(singleSpecimen);
        //if(true)return;

        ImageJ ij=new ImageJ();
        System.out.println("Now running the test...");

        // Iterate through the list
        for (String specimen : specimenList) {
            if(!specimen.contains("1191"))continue;
            /*---------Step 1 and 2 : compute the transformation matrix corresponding to the inverse of the crop of year 2022
             Load reference image and the moving image (the image that will be transformed to match the other one)*/
            ImagePlus imgHighRes22 = IJ.openImage(mainDir+"cep_"+specimen+"/"+"raw"+"/"+specimen+"_"+year1+".tif");
            double[] tr22Inv = giveCropInfoAsTranslationVectorDependingOnImageNameAndYear(specimen,year1);
            ItkTransform tCrop22Inv=ItkTransform.array16ElementsToItkTransform(new double[]{1,0,0,tr22Inv[0],   0,1,0,tr22Inv[1],   0,0,1,tr22Inv[2],   0,0,0,1 });
          
            /*---------Step 3 and 4 : compute the transformation matrix corresponding to the crop of year 2023
            // Load reference image and image that will be transformed to match the other one for the second year*/
            ImagePlus imgHighRes23 = IJ.openImage(mainDir+"cep_"+specimen+"/"+"raw"+"/"+specimen+"_"+year2+".tif");
            double[] tr23 = giveCropInfoAsTranslationVectorDependingOnImageNameAndYear(specimen,year2);
            ItkTransform tCrop23=ItkTransform.array16ElementsToItkTransform(new double[]{1,0,0,tr23[0],   0,1,0,tr23[1],   0,0,1,tr23[2],   0,0,0,1 });

            /*---------Step 5 and 6 : compute the transformation matrix corresponding to registering 23 high res (as moving) over 22 high res (as fixed)
            Produce the combined transform, combining the cropping steps and the registration matrix that was computed at low resolution*/
            //Import the transform global matrix 
            String pathToGlobalTransformAtLowRes=mainDir+"cep_"+specimen+"/"+"result/Exported_data/transform_global_img_moving.txt";
            ItkTransform tGlobalLowRes= ItkTransform.readTransformFromFile(pathToGlobalTransformAtLowRes);

            // Transformation matrix for high res
            ItkTransform trTotal=tCrop23;
            trTotal.addTransform(tGlobalLowRes);
            trTotal.addTransform(tCrop22Inv);

            // Compute the result (high res registered image), and save it
            ImagePlus finalResult=trTotal.transformImage(imgHighRes22,imgHighRes23);
            IJ.saveAsTiff(finalResult,mainDir + "registrationHighRes_results/"+specimen+"_img_moving_after_registration.tif");
            System.out.println("High resolution moving image saved to " + mainDir + "registrationHighRes_results/ as "+specimen+"_img_moving_after_registration.tif");
            
            // Save the high res transform matrix 
            trTotal.writeMatrixTransformToFile(mainDir+"cep_"+specimen+"/"+"result/testMatrices/tRigHigh.txt");
            System.out.println("Transformation Matrix saved to " + mainDir+"cep_"+specimen+"/"+"result/testMatrices/tRigHigh.txt");

            // Show single specimen (testing)
            if(singleSpecimen != null){
                 VitimageUtils.compositeNoAdjustOf(imgHighRes22, finalResult).show();
            }
                   
        }
        System.out.println("Test complete!");
    }
       

    //Big entries for computing
    public static double[] giveCropInfoAsTranslationVectorDependingOnImageNameAndYear(String specimenName,String year){
        double[]vect=new double[3];
        if(specimenName.contains("1181")){
            if(year.contains("2022"))vect=new double[]{-84.96,     -95.51,   -0.22};
            if(year.contains("2023"))vect=new double[]{204.92,     178.24,    1.58};
        }
        else if(specimenName.contains("313B")){
            if(year.contains("2022"))vect=new double[]{-73.43,    -68.8,   -0.59};
            if(year.contains("2023"))vect=new double[]{155.63,	  124.69,   0.6};
        }
        else if(specimenName.contains("318")){
            if(year.contains("2022"))vect=new double[]{-64.69,   -90,       -0.94};
            if(year.contains("2023"))vect=new double[]{127,      164.38,    -0.6};
        }
        else if(specimenName.contains("322")){
            if(year.contains("2022"))vect=new double[]{-121.7,   -140.51,   -0.78};
            if(year.contains("2023"))vect=new double[]{213.12,    178.25,    1.48};
        }
        else if(specimenName.contains("323")){
            if(year.contains("2022"))vect=new double[]{-79.1,   -65.62,    -0.15};
            if(year.contains("2023"))vect=new double[]{146.49,  192.38,     0.94};
        }
        else if(specimenName.contains("330")){
            if(year.contains("2022"))vect=new double[]{-126.9,   -105.99,   -0.64};
            if(year.contains("2023"))vect=new double[]{139.65,    211.91,    0.94};
        }
        else if(specimenName.contains("335")){
            if(year.contains("2022"))vect=new double[]{-86.13,    -113.2,   -0.8};
            if(year.contains("2023"))vect=new double[]{144.14,     189.31,  -0.27};
        }
        else if(specimenName.contains("368B")){
            if(year.contains("2022"))vect=new double[]{-26.42,    -86.16,   -0.61};
            if(year.contains("2023"))vect=new double[]{102.69,    142.41,    1.04};
        }
        else if(specimenName.contains("378A")){
            if(year.contains("2022"))vect=new double[]{-51.23,    -42.58,   -0.73};
            if(year.contains("2023"))vect=new double[]{136.72,    151.37,    1.8};
        }
        else if(specimenName.contains("378B")){
            if(year.contains("2022"))vect=new double[]{-98.16,    -48.4,    -0.58};
            if(year.contains("2023"))vect=new double[]{156.78,    126.87,    0.86};
        }
        else if(specimenName.contains("380A")){
            if(year.contains("2022"))vect=new double[]{-105.18,   -156.51,   0};
            if(year.contains("2023"))vect=new double[]{176.1,      144.51,   0};
        }
        else if(specimenName.contains("764B")){
            if(year.contains("2022"))vect=new double[]{-93.24,    -72.94,   -0.24};
            if(year.contains("2023"))vect=new double[]{178.71,    163.08,    0.6};
        }
        else if(specimenName.contains("988B")){
            if(year.contains("2022"))vect=new double[]{-42.99,   -117.76,   -0.58};
            if(year.contains("2023"))vect=new double[]{87.79,     122.37,   -0.18};
        }
        else if(specimenName.contains("1186A")){
            if(year.contains("2022"))vect=new double[]{-51.87,    -23.61,   -0.12};
            if(year.contains("2023"))vect=new double[]{173.83,    120.12,    0.6};
        }
        else if(specimenName.contains("1189")){
            if(year.contains("2022"))vect=new double[]{-131.25,  -129.88,   -0.8};
            if(year.contains("2023"))vect=new double[]{145.46,    200.98,    0.49};
        }
        else if(specimenName.contains("1191")){
            if(year.contains("2022"))vect=new double[]{-83.82,    -79.74,    0};
            if(year.contains("2023"))vect=new double[]{159.51,    131.31,    0.6};
        }
        else if(specimenName.contains("1193")){
            if(year.contains("2022"))vect=new double[]{-51.13,    -62.09,    0.26};
            if(year.contains("2023"))vect=new double[]{91.79,     143.56,    0.75};
        }
        else if(specimenName.contains("1195")){
            if(year.contains("2022"))vect=new double[]{-81.52,    -93.66,   -0.65};
            if(year.contains("2023"))vect=new double[]{74.22,      80.08,    1.2};
        }
        else if(specimenName.contains("1266A")){
            if(year.contains("2022"))vect=new double[]{-62.14,    -81.21,    0.9};
            if(year.contains("2023"))vect=new double[]{92.7,       89.41,   -0.12};
        }
        else if(specimenName.contains("2184A")){
            if(year.contains("2022"))vect=new double[]{-82.64,    -68.13,   -0.2};
            if(year.contains("2023"))vect=new double[]{148.44,    126.96,    0.26};
        }
        

        else {
            IJ.log("Warning : specimen name or year not found : "+specimenName+" "+year);
        }
        return vect;

    }

    
    //Tiny helpers

    public static ArrayList<String> getSpecimen(String singleSpecimen) {
        ArrayList<String> specimenList = new ArrayList<>();
        if (singleSpecimen != null && !singleSpecimen.isEmpty()) {
            // If a single specimen is specified, add only that specimen to the list
            specimenList.add(singleSpecimen);
        } else {
            // Otherwise, add all specimens to the list
            specimenList = getSpecimenList();
        }
        return specimenList;
    }

    public static ArrayList<String> getSpecimenList(){
        ArrayList<String> specimenList=new ArrayList<>();
        specimenList.add("313B");
        specimenList.add("318");
        specimenList.add("322");
        specimenList.add("323");
        specimenList.add("330");
        specimenList.add("335");
        specimenList.add("368B");
        specimenList.add("378A");
        specimenList.add("378B");
        specimenList.add("380A");
        specimenList.add("764B");
        specimenList.add("988B");
        specimenList.add("1181");
        specimenList.add("1186A");
        specimenList.add("1189");
        specimenList.add("1191");
        specimenList.add("1193");
        specimenList.add("1195");
        specimenList.add("1266A");
        specimenList.add("2184A");
        return  specimenList;
    }

    public static void streamLineTheVerificationOfThatEverythingProbablyWentGoodAndThRegistrationAreAccurate(String singleSpecimen){
        ArrayList<String> specimenList=getSpecimen(singleSpecimen);
        for(String specimen:specimenList){
            String year1="2022";
            ImagePlus imgRef = IJ.openImage(mainDir+"cep_"+specimen+"/"+"raw"+"/"+specimen+"_"+year1+".tif");
            ImagePlus imgMov = IJ.openImage(mainDir + "registrationHighRes_results/"+specimen+"_img_moving_after_registration.tif");
            ImagePlus imgFus = VitimageUtils.compositeNoAdjustOf(imgRef, imgMov);
            imgFus.setTitle(specimen);
            VitimageUtils.imageChecking(imgFus,specimen,5); 
            imgFus.close();
            imgRef.close();
            imgMov.close();
        }            
    }

}



