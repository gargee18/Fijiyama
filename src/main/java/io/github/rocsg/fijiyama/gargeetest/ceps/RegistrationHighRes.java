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



package io.github.rocsg.fijiyama.gargeetest.ceps;

//specific libraries
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
//my libraries
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.registration.ItkTransform;



public class RegistrationHighRes {
    
    public static String year1=GeneralUtils.years[1];
    public static String year2=GeneralUtils.years[2];
    //public static String mainDir="/home/phukon/Desktop/registration/CEP_2022_XR_2023_XR_2024_XR/";// Path to the directory consisting the data
    public static String mainDir="/home/phukon/Desktop/";// Path to the directory consisting the data


    //Main input points
    public static void main(String[] args) {
            
        // Specify the single specimen here or set it to null to get the entire list
        String singleSpecimen ="330"; 
        
        // Build a list of all specimen names (with respect with the folder s name)
        String[] specimenList = getSpecimen(singleSpecimen);

        //if(true)return;

        ImageJ ij=new ImageJ();
        
        System.out.println("Now running the test...");

        //Iterate through the list
        for (String specimen : specimenList) {
            getHighResTransform(specimen);
           // streamLineTheVerificationOfThatEverythingProbablyWentGoodAndThRegistrationAreAccurate(specimen);

        }



       
    }   
        public static ItkTransform getHighResTransform(String specimen){
            //if(!specimen.contains("1191"))continue;
            /*---------Step 1 and 2 : compute the transformation matrix corresponding to the inverse of the crop of year 2022
             Load reference image and the moving image (the image that will be transformed to match the other one)*/
            ImagePlus imgHighRes22 = IJ.openImage("/home/phukon/Desktop/Results_2024/Cep_Monitoring/XR/4_Registered_Images/16bit_Registered_2022_2023_2024/CEP_"+specimen+"_2023_XR.tif");
            //ImagePlus imgHighRes22 = IJ.openImage("/home/phukon/Desktop/Results_2024/Cep_Monitoring/MRI/XR_MRI_2023_16bit/"+specimen+"_TransformedImage_T1.tif");
            double[] tr22Inv = giveCropInfoAsTranslationVectorDependingOnImageNameAndYear23and24(specimen,year1);
            // double[] tr22Inv = giveCropInfoAsTranslationVectorForXRandMRI24(specimen,"XR");
            ItkTransform tCrop22Inv=ItkTransform.array16ElementsToItkTransform(new double[]{1,0,0,-tr22Inv[0],   0,1,0,-tr22Inv[1],   0,0,1,-tr22Inv[2],   0,0,0,1 });
          
            /*---------Step 3 and 4 : compute the transformation matrix corresponding to the crop of year 2023
            // Load reference image and image that will be transformed to match the other one for the second year*/
            ImagePlus imgHighRes23 = IJ.openImage("/home/phukon/Desktop/Results_2024/Cep_Monitoring/XR/4_Registered_Images/16bit_Registered_2022_2023_2024/CEP_"+specimen+"_2024_XR.tif");
            //ImagePlus imgHighRes23 = IJ.openImage("/home/phukon/Desktop/MRI/02_ceps/2022_CEPS_suivi_clinique/1_STACK_16-bits/cep_"+specimen+"_"+year1+"_MRI_T1.tif");
            // double[] tr23 = giveCropInfoAsTranslationVectorForXRandMRI24(specimen,"MRI");
            double[] tr23 = giveCropInfoAsTranslationVectorDependingOnImageNameAndYear23and24(specimen,year2);
            ItkTransform tCrop23=ItkTransform.array16ElementsToItkTransform(new double[]{1,0,0,tr23[0],   0,1,0,tr23[1],   0,0,1,tr23[2],   0,0,0,1 });

            /*---------Step 5 and 6 : compute the transformation matrix corresponding to registering 23 high res (as moving) over 22 high res (as fixed)
            Produce the combined transform, combining the cropping steps and the registration matrix that was computed at low resolution*/
            //Import the transform global matrix 
            //String pathToGlobalTransformAtLowRes="/home/phukon/Desktop/master_repo/mrit1_2023_mrit2_2023/cep_"+specimen+"/res/Exported_data/"+specimen+"_trMatrix_2023_2024.txt"; 
            String pathToGlobalTransformAtLowRes="/home/phukon/Desktop/master_repo/xray_2023_2024/cep_"+specimen+"/res/Exported_data/"+specimen+"_trMatrix_2023_2024.txt";
            ItkTransform tGlobalLowRes= ItkTransform.readTransformFromFile(pathToGlobalTransformAtLowRes);

            // Transformation matrix for high res
            // ItkTransform trTotal=tCrop23;
            // trTotal.addTransform(tGlobalLowRes);
            // trTotal.addTransform(tCrop22Inv);
            
            ItkTransform trTotal=tCrop23;
            trTotal.addTransform(tGlobalLowRes);
            trTotal.addTransform(tCrop22Inv);
            
            



            // Compute the result (high res registered image), and save it
            ImagePlus finalResult=trTotal.transformImage(imgHighRes22,imgHighRes23);
            // IJ.saveAsTiff(finalResult,mainDir+"cep_"+specimen+"/result_"+year1+"_"+year2+"_highRes/"+specimen+"_img_moving_after_registration_test.tif");
            IJ.saveAsTiff(finalResult,"/home/phukon/Desktop/Results_2024/Cep_Monitoring/MRI/XR_MRI_2024_16bit/"+specimen+"_TransformedImage_T1.tif");
            //System.out.println("High resolution moving image saved to : /home/phukon/Desktop/XR_2022_raw/"+specimen+"_TransformedImage.tif");
            
            // Save the high res transform matrix 
            trTotal.writeMatrixTransformToFile("/home/phukon/Desktop/Results_2024/Cep_Monitoring/MRI/XR_MRI_2024_16bit/"+specimen+"_img_moving_after_registration_highres_T1.txt");
            //System.out.println("Transformation Matrix saved to :  "+ "/home/phukon/Desktop/XR_2022_raw/"+specimen+"_img_moving_after_registration_highres.txt");


            ImagePlus imgFus = VitimageUtils.compositeNoAdjustOf(imgHighRes22, finalResult);
            // IJ.saveAsTiff(imgFus, "/home/phukon/Desktop/xray_2022_2023_hires/"+specimen+"_Composite.tif");
            IJ.saveAsTiff(imgFus, "/home/phukon/Desktop/Results_2024/Cep_Monitoring/MRI/XR_MRI_2024_16bit/"+specimen+"_Composite_T1.tif");

            return trTotal;
        }


        public static double[] giveCropInfoAsTranslationVectorForMRIT1andT222(String specimenName,String mod){
            double[]vect=new double[3];
            if(specimenName.contains("318")){
                if(mod.contains("T1"))vect=new double[]{64.69,	90,	-0.24};
                if(mod.contains("T2"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("322")){
                if(mod.contains("T1"))vect=new double[]{121.77,	140.51,	0.34};
                if(mod.contains("T2"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("323")){
                if(mod.contains("T1"))vect=new double[]{79.1,	65.62,	0.16};
                if(mod.contains("T2"))vect=new double[]{0,	0,	0};
            }
            else if(specimenName.contains("330")){
                if(mod.contains("T1"))vect=new double[]{126.9,	105.98,	-0.12};
                if(mod.contains("T2"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("335")){
                if(mod.contains("T1"))vect=new double[]{85.31,	113.2,	0.8};
                if(mod.contains("T2"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("764B")){
                if(mod.contains("T1"))vect=new double[]{93.24,	72.94,	-0.18};
                if(mod.contains("T2"))vect=new double[]{0,0,0};
 
            }
        
            else if(specimenName.contains("1181")){
                if(mod.contains("T1"))vect=new double[]{85.54,	95.5,	-0.41};
                if(mod.contains("T2"))vect=new double[]{97.53,	0,	0};

            }
           
            else if(specimenName.contains("1193")){
                if(mod.contains("T1"))vect=new double[]{51.13,	61.36,	0.17};
                if(mod.contains("T2"))vect=new double[]{0,0,0};

            }
           
            else {
                IJ.log("Warning : specimen name or mod not found : "+specimenName+" "+mod);
            }
            return vect;
        }

        public static double[] giveCropInfoAsTranslationVectorForMRIT1andT223(String specimenName,String mod){
            double[]vect=new double[3];
            if(specimenName.contains("318")){
                if(mod.contains("T1"))vect=new double[]{64.12,	90,	2.12};
                if(mod.contains("T2"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("322")){
                if(mod.contains("T1"))vect=new double[]{122.63,	140.51,	-0.84};
                if(mod.contains("T2"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("323")){
                if(mod.contains("T1"))vect=new double[]{79.1,	65.62,	-0.64};
                if(mod.contains("T2"))vect=new double[]{0,	0,	0};
            }
            else if(specimenName.contains("330")){
                if(mod.contains("T1"))vect=new double[]{126.9,	105.98,	0.28};
                if(mod.contains("T2"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("335")){
                if(mod.contains("T1"))vect=new double[]{86.14,	113.2,	2};
                if(mod.contains("T2"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("764B")){
                if(mod.contains("T1"))vect=new double[]{93.24,	72.94,	1.42};
                if(mod.contains("T2"))vect=new double[]{0,0,0};
 
            }
        
            else if(specimenName.contains("1181")){
                if(mod.contains("T1"))vect=new double[]{84.96,	96.09,	-0.39};
                if(mod.contains("T2"))vect=new double[]{103.9,	0	,0};

            }
           
            else if(specimenName.contains("1193")){
                if(mod.contains("T1"))vect=new double[]{50.4,	62.09,	-0.22};
                if(mod.contains("T2"))vect=new double[]{0,0,0};

            }
           
            else {
                IJ.log("Warning : specimen name or mod not found : "+specimenName+" "+mod);
            }
            return vect;
        }
        
        
        public static double[] giveCropInfoAsTranslationVectorForXRandMRI22(String specimenName,String mod){
            double[]vect=new double[3];
            if(specimenName.contains("318")){
                if(mod.contains("XR"))vect=new double[]{64.69,89.44,0};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("322")){
                if(mod.contains("XR"))vect=new double[]{121.78,	140.51,	0};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("323")){
                if(mod.contains("XR"))vect=new double[]{79.1,65.63,0};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("330")){
                if(mod.contains("XR"))vect=new double[]{126.9,105.93,0};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("335")){
                if(mod.contains("XR"))vect=new double[]{86.13,113.21,0};
                if(mod.contains("v"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("764B")){
                if(mod.contains("XR"))vect=new double[]{93.24,72.94,0};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
 
            }

            else if(specimenName.contains("1181")){
                if(mod.contains("XR"))vect=new double[]{84.96,	95.51,	0.99};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};

            }
           
            else if(specimenName.contains("1193")){
                if(mod.contains("XR"))vect=new double[]{51.14,61.36,0};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};

            }
           
            else {
                IJ.log("Warning : specimen name or mod not found : "+specimenName+" "+mod);
            }
            return vect;
        }
            
    
        public static double[] giveCropInfoAsTranslationVectorForXRandMRI23(String specimenName,String mod){
            double[]vect=new double[3];
            if(specimenName.contains("318")){
                if(mod.contains("XR"))vect=new double[]{64.13,	90,	2.54};
                if(mod.contains("MRI"))vect=new double[]{67.68,0,0};
            }
            else if(specimenName.contains("322")){
                if(mod.contains("XR"))vect=new double[]{121.77,	140.51,	0.65};
                if(mod.contains("MRI"))vect=new double[]{79.98,	3.42,	0};
            }
            else if(specimenName.contains("323")){
                if(mod.contains("XR"))vect=new double[]{79.11,	65.62,	-0.47};
                if(mod.contains("MRI"))vect=new double[]{45.81,	0,	0};
            }
            else if(specimenName.contains("330")){
                if(mod.contains("XR"))vect=new double[]{126.9,	105.98,	0.68};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("335")){
                if(mod.contains("XR"))vect=new double[]{86.13,	113.21,	2};
                if(mod.contains("v"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("764B")){
                if(mod.contains("XR"))vect=new double[]{93.24,	72.94,	1.04};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
 
            }
            else if(specimenName.contains("988B")){
                if(mod.contains("XR"))vect=new double[]{87.79,	121.48,	0};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
            }

            else if(specimenName.contains("1181")){
                if(mod.contains("XR"))vect=new double[]{84.96,	95.96,	-0.44};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};

            }
           
            else if(specimenName.contains("1193")){
                if(mod.contains("XR"))vect=new double[]{51.13,	62.09,	0.54};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};

            }
           
            else {
                IJ.log("Warning : specimen name or mod not found : "+specimenName+" "+mod);
            }
            return vect;
        }


        public static double[] giveCropInfoAsTranslationVectorForXRandMRI24(String specimenName,String mod){
            double[]vect=new double[3];
            if(specimenName.contains("318")){
                if(mod.contains("XR"))vect=new double[]{64.13,	90,	0.97};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("322")){
                if(mod.contains("XR"))vect=new double[]{121.78,	140.51,	0.73};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("323")){
                if(mod.contains("XR"))vect=new double[]{84.96,	73.24,	0.35};
                if(mod.contains("MRI"))vect=new double[]{0,	0,	0};
            }
            else if(specimenName.contains("330")){
                if(mod.contains("XR"))vect=new double[]{126.9,	105.98,	-0.52};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("335")){
                if(mod.contains("XR"))vect=new double[]{86.13,	111.56,	2.4};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
            }
            else if(specimenName.contains("764B")){
                if(mod.contains("XR"))vect=new double[]{90.23,	75.94,	0.26};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};
 
            }
            else if(specimenName.contains("1181")){
                if(mod.contains("XR"))vect=new double[]{83.79,	94.34,	-0.37};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};

            }
           
            else if(specimenName.contains("1193")){
                if(mod.contains("XR"))vect=new double[]{51.13,	62.09,	0.54};
                if(mod.contains("MRI"))vect=new double[]{0,0,0};

            }
           
            else {
                IJ.log("Warning : specimen name or mod not found : "+specimenName+" "+mod);
            }
            return vect;
        }
            
    
       

    //Big entries for computing
    public static double[] giveCropInfoAsTranslationVectorDependingOnImageNameAndYear22and23(String specimenName,String year){
        double[]vect=new double[3];
        if(specimenName.contains("1181")){
            if(year.contains("2022"))vect=new double[]{84.96,       95.51,    0.22};
            if(year.contains("2023"))vect=new double[]{204.92,     178.24,    1.58};
        }

        else if(specimenName.contains("313B")){
            if(year.contains("2022"))vect=new double[]{73.43,      68.8,    0.59};
            if(year.contains("2023"))vect=new double[]{155.63,    124.69,	0.6};
        }
        else if(specimenName.contains("318")){
            if(year.contains("2022"))vect=new double[]{64.69,     90,        0.94};
            if(year.contains("2023"))vect=new double[]{127,      164.38,    -0.6};
        }
        else if(specimenName.contains("322")){
            if(year.contains("2022"))vect=new double[]{121.7,     140.51,    0.78};
            if(year.contains("2023"))vect=new double[]{213.12,    178.25,    1.48};
        }
        else if(specimenName.contains("323")){
            if(year.contains("2022"))vect=new double[]{79.1,      65.62,     0.15};
            if(year.contains("2023"))vect=new double[]{146.49,   192.38,     0.94};
        }
        else if(specimenName.contains("330")){
            if(year.contains("2022"))vect=new double[]{126.9,     105.99,    0.64};
            if(year.contains("2023"))vect=new double[]{139.65,    211.91,    0.94};
        }
        else if(specimenName.contains("335")){
            if(year.contains("2022"))vect=new double[]{86.13,      113.2,    0.8};
            if(year.contains("2023"))vect=new double[]{144.14,     189.31,  -0.27};
        }
        else if(specimenName.contains("368B")){
            if(year.contains("2022"))vect=new double[]{26.42,      86.16,    0.61};
            if(year.contains("2023"))vect=new double[]{102.69,    142.41,    1.04};
        }
        else if(specimenName.contains("378A")){
            if(year.contains("2022"))vect=new double[]{51.23,      42.58,    0.73};
            if(year.contains("2023"))vect=new double[]{136.72,    151.37,    1.8};
        }
        else if(specimenName.contains("378B")){
            if(year.contains("2022"))vect=new double[]{98.16,      48.4,     0.58};
            if(year.contains("2023"))vect=new double[]{156.78,    126.87,    0.86};
        }
        else if(specimenName.contains("380A")){
            if(year.contains("2022"))vect=new double[]{105.18,     156.51,   0};
            if(year.contains("2023"))vect=new double[]{176.1,      144.51,   0};
        }
        else if(specimenName.contains("764B")){
            if(year.contains("2022"))vect=new double[]{93.24,      72.94,    0.24};
            if(year.contains("2023"))vect=new double[]{178.71,    163.08,    0.6};
        }
        else if(specimenName.contains("988B")){
            if(year.contains("2022"))vect=new double[]{42.99,     117.76,    0.58};
            if(year.contains("2023"))vect=new double[]{87.79,     122.37,   -0.18};
        }
        else if(specimenName.contains("1186A")){
            if(year.contains("2022"))vect=new double[]{51.87,      23.61,    0.12};
            if(year.contains("2023"))vect=new double[]{173.83,    120.12,    0.6};
        }
        else if(specimenName.contains("1189")){
            if(year.contains("2022"))vect=new double[]{131.25,    129.88,    0.8};
            if(year.contains("2023"))vect=new double[]{145.46,    200.98,    0.49};
        }
        else if(specimenName.contains("1191")){
            if(year.contains("2022"))vect=new double[]{83.82,      79.74,    0};
            if(year.contains("2023"))vect=new double[]{159.51,    131.31,    0.6};
        }
        else if(specimenName.contains("1193")){
            if(year.contains("2022"))vect=new double[]{51.13,      62.09,   -0.26};
            if(year.contains("2023"))vect=new double[]{91.79,     143.56,    0.75};
        }
        else if(specimenName.contains("1195")){
            if(year.contains("2022"))vect=new double[]{81.52,      93.66,    0.65};
            if(year.contains("2023"))vect=new double[]{74.22,      80.08,    1.2};
        }
        else if(specimenName.contains("1266A")){
            if(year.contains("2022"))vect=new double[]{62.14,      81.21,   -0.9};
            if(year.contains("2023"))vect=new double[]{92.7,       89.41,   -0.12};
        }
        else if(specimenName.contains("2184A")){
            if(year.contains("2022"))vect=new double[]{82.64,      68.13,    0.2};
            if(year.contains("2023"))vect=new double[]{148.44,    126.96,    0.26};
        }
        

        else {
            IJ.log("Warning : specimen name or year not found : "+specimenName+" "+year);
        }
        return vect;

    }

    public static double[] giveCropInfoAsTranslationVectorDependingOnImageNameAndYear23and24(String specimenName,String year){
        double[]vect=new double[3];
        if(specimenName.contains("1181")){
            if(year.contains("2023"))vect=new double[]{84.96,	95.51,	-0.43};
            if(year.contains("2024"))vect=new double[]{125.51,     107.46,   -0.13};
        }
        else if(specimenName.contains("313B")){
            if(year.contains("2023"))vect=new double[]{72.85,	  69.37,	1.8};
            if(year.contains("2024"))vect=new double[]{93.52, 	  91.05,	0.41};
        }
        else if(specimenName.contains("318")){
            if(year.contains("2023"))vect=new double[]{64.13,     90,        1.32};
            if(year.contains("2024"))vect=new double[]{100.08,	 102.54,	 0.6};
        }
        else if(specimenName.contains("322")){
            if(year.contains("2023"))vect=new double[]{121.71,    140.5,    0.56};
            if(year.contains("2024"))vect=new double[]{167.34,	  140.27,	 0.55};
        }
        else if(specimenName.contains("323")){
            if(year.contains("2023"))vect=new double[]{79.1,      65.63,     0.34};
            if(year.contains("2024"))vect=new double[]{106.64,	 118.12,    -0.24};
        }
        else if(specimenName.contains("330")){
            if(year.contains("2023"))vect=new double[]{127.06,    105.98,    0.31};
            if(year.contains("2024"))vect=new double[]{105.82,	   92.7, 	 0.9};
        }
        else if(specimenName.contains("335")){
            if(year.contains("2023"))vect=new double[]{86.14,      113.2,    1.6};
            if(year.contains("2024"))vect=new double[]{116.9,	   107.23,	-0.57};
        }
        else if(specimenName.contains("368B")){
            if(year.contains("2023"))vect=new double[]{26.41,      86.16,    0.57};
            if(year.contains("2024"))vect=new double[]{143.19,	  143.18, 	 0.6};
        }
        else if(specimenName.contains("378A")){
            if(year.contains("2023"))vect=new double[]{51.22,      42.59,   -0.49};
            if(year.contains("2024"))vect=new double[]{73.82,	   55.78,	 0.3};
        }
        else if(specimenName.contains("378B")){
            if(year.contains("2023"))vect=new double[]{98.16,	48.4,	0.27};
            if(year.contains("2024"))vect=new double[]{110.75,	  107.46,	-0.3};
        }
        else if(specimenName.contains("380A")){
            if(year.contains("2023"))vect=new double[]{105.17,	156.51,	0.03};
            if(year.contains("2024"))vect=new double[]{125.07,	    43.5,	 0.5};
        }
        else if(specimenName.contains("764B")){
            if(year.contains("2023"))vect=new double[]{93.25,	72.94,	-0.47};
            if(year.contains("2024"))vect=new double[]{100.2,	   43.89,	 0.7};
        }
        else if(specimenName.contains("988B")){
            if(year.contains("2023"))vect=new double[]{42.36,	118.38,	  1};
            if(year.contains("2024"))vect=new double[]{94.34,	  109.1,	 0.04};
        }
        else if(specimenName.contains("1186A")){
            if(year.contains("2023"))vect=new double[]{52.07,	24.21,	0.51};
            if(year.contains("2024"))vect=new double[]{112.38,	   77.11,	 0.6};
        }
        else if(specimenName.contains("1189")){
            if(year.contains("2023"))vect=new double[]{130.57,	129.88,	1.2};
            if(year.contains("2024"))vect=new double[]{114.62,	   87.5,	 0.26};
        }
        else if(specimenName.contains("1191")){
            if(year.contains("2023"))vect=new double[]{84.4,	83.23,	2.13};
            if(year.contains("2024"))vect=new double[]{141.37,	  100.59,	 0.3};
        }
        else if(specimenName.contains("1193")){
            if(year.contains("2023"))vect=new double[]{51.13,	62.09,	-0.27};
            if(year.contains("2024"))vect=new double[]{51.2,  	   10.59,	 1};
        }
        else if(specimenName.contains("1195")){
            if(year.contains("2023"))vect=new double[]{81.52,	93.66,	0.69};
            if(year.contains("2024"))vect=new double[]{95.14,	    0,	     0.6};
        }
        else if(specimenName.contains("1266A")){
            if(year.contains("2023"))vect=new double[]{62.14,	81.21,	-0.1};
            if(year.contains("2024"))vect=new double[]{114.84,	   82.85,	 0.3};
        }
        else if(specimenName.contains("2184A")){
            if(year.contains("2023"))vect=new double[]{82.65,	68.13,	0.99};
            if(year.contains("2024"))vect=new double[]{81.22,	  100.07,	 0.56};
        }
        

        else {
            IJ.log("Warning : specimen name or year not found : "+specimenName+" "+year);
        }
        return vect;

    }
    
    //Tiny helpers
    public static String[] getSpecimen(String singleSpecimen) {
        String[] specimenList;
        if (singleSpecimen != null && !singleSpecimen.isEmpty()) {
            // Select one specimen
            specimenList = new String[]{singleSpecimen};
        } else {
            // Select the entire list
            specimenList = GeneralUtils.specimenListMRI;
        }
        return specimenList;
    }

    
    public static void streamLineTheVerificationOfThatEverythingProbablyWentGoodAndThRegistrationAreAccurate(String singleSpecimen){
        String[] specimenList=getSpecimen(singleSpecimen);
        for(String specimen:specimenList){
            ImagePlus imgRef = IJ.openImage(mainDir+"cep_"+specimen+"/raw/"+specimen+"_"+year1+".tif");
            ImagePlus imgMov = IJ.openImage(mainDir+"cep_"+specimen+"/res/"+specimen+"_TransformedImage.tif");
            ImagePlus imgFus = VitimageUtils.compositeNoAdjustOf(imgRef, imgMov);
            imgFus.setTitle(specimen);
            VitimageUtils.imageChecking(imgFus,specimen,5); 
            imgFus.close();
            imgRef.close();
            imgMov.close();
        }            
    }

}



