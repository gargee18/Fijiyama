package io.github.rocsg.fijiyama.gargeetest;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.fijiyamaplugin.RegistrationAction;
import io.github.rocsg.fijiyama.registration.BlockMatchingRegistration;
import io.github.rocsg.fijiyama.registration.ItkTransform;

public class TestGargee {


    // public static void testRegistration(){
    //     ImagePlus imgRef=null;
    //     ImagePlus imgMov=null;
    //     RegistrationAction regAct=new RegistrationAction();
    //     regAct.defineSettingsFromTwoImages(imgRef, imgMov, null, false);
    //     regAct.bhsX=3;
    //     regAct.bhsY=3;
        
    //     BlockMatchingRegistration bm = BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef, imgMov,regAct);
    //     ItkTransform tr=bm.runBlockMatching(null, false);

    //     //Apply the transform 
    //     ImagePlus imgRegistered=tr.transformImage(imgRef, imgMov);
        
    // }
    static String mainDir="/home/phukon/Desktop/registration/";
       
    public static void main(String[] args) {
        ImageJ ij=new ImageJ();
        System.out.println("Now running the test...");


        //makeDemo("902");return;
        testMri();if(true)return;

        
        // Path to the directory consisting the data
        String specimen="318";

        //---------Step 1 and 2 : compute the transformation matrix corresponding to the inverse of the crop of year 2022
        // Load reference image and image that will be transformed to match the other one for the first year
        String year1="2022";
        ImagePlus imgHighRes22 = IJ.openImage(mainDir+"cep_"+specimen+"/"+"raw"+"/"+specimen+"_"+year1+".tif");
        ImagePlus imgCrop22 = IJ.openImage(mainDir+"cep_"+specimen+"/"+"raw"+"/"+specimen+"_"+year1+"_crop_sub_z.tif");
        double[] tr22Inv = giveCropInfoAsTranslationVectorDependingOnImageNameAndYear(specimen,year1);
        ItkTransform tCrop22Inv=ItkTransform.array16ElementsToItkTransform(new double[]{1,0,0,tr22Inv[0],   0,1,0,tr22Inv[1],   0,0,1,tr22Inv[2],   0,0,0,1 });


        //---------Step 3 and 4 : compute the transformation matrix corresponding to the crop of year 2023
        // Load reference image and image that will be transformed to match the other one for the second year
        String year2="2023";
        ImagePlus imgHighRes23 = IJ.openImage(mainDir+"cep_"+specimen+"/"+"raw"+"/"+specimen+"_"+year2+".tif");
        ImagePlus imgCrop23 = IJ.openImage(mainDir+"cep_"+specimen+"/"+"raw"+"/"+specimen+"_"+year2+"_crop_sub_z.tif");
        double[] tr23 = giveCropInfoAsTranslationVectorDependingOnImageNameAndYear(specimen,year2);
        ItkTransform tCrop23=ItkTransform.array16ElementsToItkTransform(new double[]{1,0,0,tr23[0],   0,1,0,tr23[1],   0,0,1,tr23[2],   0,0,0,1 });


        //---------Step 5 and 6 : compute the transformation matrix corresponding to registering 23 high res (as moving) over 22 high res (as fixed)

        //Produce the combined transform, combining the cropping steps and the registration matrix that was computed at low resolution
        //Import the transform global matrix 
        String pathToGlobalTransformAtLowRes=mainDir+"cep_"+specimen+"/"+"result/Exported_data/transform_global_img_moving.txt";
        ItkTransform tGlobalLowRes= ItkTransform.readTransformFromFile(pathToGlobalTransformAtLowRes);
    
        // Display Low res transform
        // ImagePlus testLowRes=tGlobalLowRes.transformImage(imgCrop22,imgCrop23);
        // testLowRes.show(); testLowRes.setTitle("LowResTransform");
        // imgCrop23.show();
        // imgCrop22.show(); 
        // VitimageUtils.compositeNoAdjustOf(imgCrop22, testLowRes).show();
        
         
        // Transformation matrix for high res
        ItkTransform trTotal=tCrop23;
        trTotal.addTransform(tGlobalLowRes);
        trTotal.addTransform(tCrop22Inv);
        
        // Display High res transform
        // ImagePlus finalResult=trTotal.transformImage(imgHighRes22,imgHighRes23);
        // finalResult.show(); finalResult.setTitle("HighResTransform");
        // imgHighRes23.show();
        // imgHighRes22.show(); 
        // VitimageUtils.compositeNoAdjustOf(imgHighRes22, finalResult).show();

        // Save the high res transform matix 
        trTotal.writeMatrixTransformToFile(mainDir+"cep_"+specimen+"/"+"result/testMatrices/tRigHigh.txt");
        System.out.println("Transformation Matrix saved to " + mainDir+"cep_"+specimen+"/"+"result/testMatrices/tRigHigh.txt");
        System.out.println("Test complete!");

 
    }
    
    
    public static void testMri(){
        // MRI data
        String specimen="318";
        String year1="2022";
        ImagePlus imgHighRes22 = IJ.openImage(mainDir+"cep_"+specimen+"/"+"raw"+"/"+specimen+"_"+year1+".tif");
        ImagePlus imgMRI = IJ.openImage("/home/phukon/Desktop/MRI/02_ceps/2022-03_CEPS_suivi_clinique/trial_XR_MRI/Reslice of CEP_318_2022_MRI_T1_mirrored_x.tif");
        RegistrationAction regAct=new RegistrationAction();
        System.out.println("The registration action at very start have these parameters :\n"+regAct);
        regAct.defineSettingsFromTwoImages(imgHighRes22, imgMRI, null, false);
        System.out.println("The registration action after define settings have these parameters :\n"+regAct);

        regAct.bhsX=3;
        regAct.bhsY=3;
        regAct.bhsZ=2;
        regAct.neighX=2;
        regAct.neighY=2;
        regAct.neighZ=2;
        regAct.strideX=8;
        regAct.strideY=8;
        regAct.strideZ=11;
        regAct.iterationsBMLin=8;
        regAct.selectScore=95;
        regAct.selectLTS=80;
        System.out.println("The registration action after custom settings have these parameters :\n"+regAct);

        BlockMatchingRegistration bm = BlockMatchingRegistration.setupBlockMatchingRegistration(imgHighRes22, imgMRI,regAct);
        ItkTransform tr=bm.runBlockMatching(null, false);

        //Apply the transform 
        ImagePlus imgRegistered=tr.transformImage(imgHighRes22, imgMRI);
        imgRegistered.show(); imgRegistered.setTitle("HighResTransform");
        VitimageUtils.compositeNoAdjustOf(imgHighRes22, imgRegistered).show();
        System.out.println("Toto end");
        // VitimageUtils.waitFor(100000);
    }


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
            if(year.contains("2022"))vect=new double[]{-64.69,   -90,      -0.94};
            if(year.contains("2023"))vect=new double[]{127,      164.38,   -0.6};
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
            if(year.contains("2022"))vect=new double[]{-86.13,    -113.2,     -0.8};
            if(year.contains("2023"))vect=new double[]{144.14,     189.31,    -0.27};
        }
        else if(specimenName.contains("368B")){
            if(year.contains("2022"))vect=new double[]{-26.42,    -86.16,     -0.61};
            if(year.contains("2023"))vect=new double[]{102.69,    142.41,      1.04};
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
            if(year.contains("2022"))vect=new double[]{-83.81,    -80.32,    0};
            if(year.contains("2023"))vect=new double[]{159.51,    131.31,   -0.6};
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
            IJ.showMessage("Warning : specimen name or year not found : "+specimenName+" "+year);
        }
        return vect;
    }












































    // 
    public static double[][] cropMatrixRef(double translationX, double translationY, double translationZ) {
        double[][] Tcrop22 = {
            {1.0, 0.0, 0.0, -translationX},
            {0.0, 1.0, 0.0, -translationY},
            {0.0, 0.0, 1.0, -translationZ},
            {0.0, 0.0, 0.0, 1.0}
        };
        return Tcrop22;
    }
    
    // Print Transfomation Matrix
    public static void printMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    // Method definition
    /*
     * This method takes two parameters: inputImage, which is the original image to be transformed, and transformationMatrix, matrix representing the transformation to be applied.
       It returns an ImagePlus object representing the transformed image.
     */
    public static ImagePlus applyTransformation(ImagePlus inputImage, double[][] transformationMatrix) {
        // Dimensions of the input image
        int width = inputImage.getWidth();;  
        int height =  inputImage.getHeight();  
        int depth = inputImage.getStackSize();
        
        // Blank image to store the transformed image
        ImagePlus transformedImage = IJ.createHyperStack("Transformed Image", width, height, 1, depth, 1, 8);
        


        // Get the pixel data from the input image
        for (int z = 1; z <= depth; z++) { // iterates over each slice of the input image
            inputImage.setSlice(z); // sets the current slice of the input image.
            ImageProcessor inputProcessor = inputImage.getProcessor();
            ImageProcessor transformedProcessor = inputProcessor.duplicate();
            
            /*
            * Nested loops iterate over each pixel in the current slice. 
            * For each pixel, the transformation matrix is applied to determine the new coordinates after transformation.
            */
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    double[] result = applyTransformationMatrix(transformationMatrix, x, y, z); // Applies the transformation matrix to the current pixel coordinates.
                    int newX = (int) result[0]; // Extracts the new x and y coordinates after transformation.
                    int newY = (int) result[1];
                    int pixelValue = inputProcessor.getPixel(newX, newY); // Retrieves the pixel value at the transformed coordinates from the input image.
                    transformedProcessor.putPixel(x, y, pixelValue);//  Sets the pixel value in the transformed image processor at the original pixel coordinates.
                }
            }
            
            transformedImage.getStack().setProcessor(transformedProcessor, z);
        }
        
        return transformedImage;
    }


    /*public static void testForWritingStuff(){
        ItkTransform t=new ItkTransform();
        //t= receive siomlet
        t=ItkTransform.array16ElementsToItkTransform(transformationMatrix);
        double[]blabla=t.transformPoint(new double[]{x,y,z});

    }*/
    

    private static double[] applyTransformationMatrix(double[][] transformationMatrix, int x, int y, int z) {
        double[] result = new double[3];
        result[0] = transformationMatrix[0][0] * x + transformationMatrix[0][1] * y + transformationMatrix[0][2] * z + transformationMatrix[0][3];
        result[1] = transformationMatrix[1][0] * x + transformationMatrix[1][1] * y + transformationMatrix[1][2] * z + transformationMatrix[1][3];
        result[2] = transformationMatrix[2][0] * x + transformationMatrix[2][1] * y + transformationMatrix[2][2] * z + transformationMatrix[2][3];
        return result;
    }


















    public static void test() {
        //ItkTransform myRigidBodyTransform=ItkTransform.readTransformFromFile("Please write here your path");
        //or
        //ItkTransform myRigidBodyTransform=ItkTransform.array16ElementsToItkTransform(new double[]{1,0,0,1,  3,1,4,2,  6,8,1,3,  0,0,0,1});
        ItkTransform cropMatrixRef=ItkTransform.array16ElementsToItkTransform(new double[]{1,0,0,1,  3,1,4,2,  6,8,1,3,  0,0,0,1});
        //doSomething(myRigidBodyTransform);
        System.out.println("Test passed!");

        double[] aMatrixtoVector = cropMatrixRef.from2dMatrixto1dVector();

        //System.out.println(aMatrixtoVector);   
       // for(int i=0;i<6 ; i++)System.out.println(aMatrixtoVector[i]);
        
        for(int i=0;i<6 ; i++)System.out.print(aMatrixtoVector[i]+", ");
        System.out.println();
    }



    public static void doSomething(ItkTransform myRigidBodyTransform) {
        System.out.println("I am doing something with the transform");
    }

    

    

}
