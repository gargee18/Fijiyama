package io.github.rocsg.fijiyama.gargeetest.ceps;


import java.util.ArrayList;

//specific libraries
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.fijiyamaplugin.RegistrationAction;
import io.github.rocsg.fijiyama.registration.BlockMatchingRegistration;
import io.github.rocsg.fijiyama.registration.ItkTransform;
import io.github.rocsg.fijiyama.registration.Transform3DType;

import java.time.Duration;
import java.time.Instant;


public class RegistrationAutomate {
 
 
    /* The intend of this class is to automate actions of registration (rigid and dense) over a set of specimens
    * corresponding to 20 grapevine trunk, observed in 2022 in both MRI(7) and Xray(20)
    * The registration will be computed with Xray as fixed image an MRI as moving image or Xray of year 1 as 
    * fixed and Xray of year 2 as moving or MRI T1 as fixed and MRI T2 as moving for the same year 
    * autoLinearRegistration performs Rigid registration and autoNonLinearRegistration performs Dense field registration. Both return the transformation matrix. 
    * The result (the ItkTransform) is saved in a designated directory.
    * runRigidRegistration and runDenseRegistration calls its respective registration methods and saves the transformed moving image
    * and composite image to a designated directory. 
    */

    //public static  String mainDir="/home/phukon/Desktop/registration/CEP_2022_XR_2023_XR_2024_XR/";// Path to the directory consisting the data
    public static String mainDir = "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/Cuttings_MRI_registration/";
    public static String year1 = GeneralUtils.years[0];
    public static String year2 = GeneralUtils.years[1];
    static String[] timestamps = new String[]{"J001", "J029", "J077", "J141"};
   
   
    
    // public RegistrationAutomate(RegistrationAutomate otherRegToBeCopied){    this.year1=otherRegToBeCopied.year1; }

    public static void main(String[] args) {
        ImageJ ij=new ImageJ();
        //String year1="2023";
        // String chooseSpecimen = "B_201";
    
        String[] specimenList = GeneralUtils.getSpecimenListMRI();
        // for (String specimen : specimenList) {
        //      runRigidRegistration(specimen);
        // }
        for (int i = 201; i <= 201; i++) {
            if (i == 204) {
                continue;
            }
            String chooseSpecimen = "B_" + i;
            Instant start = Instant.now();
            runRigidRegistration(chooseSpecimen);
            Instant end = Instant.now();
            Duration duration = Duration.between(start, end);
            long seconds = duration.getSeconds();
            System.out.println("Elapsed time in seconds: " + seconds);
        }
        //seeResultsOfRigidRegistration(chooseSpecimen);
        //System.exit(0);
    }

    public static String getPathToImageLowRes(String specimen,int step){
              return mainDir+specimen+"/raw_subsampled/"+specimen+"_"+timestamps[step]+"_sub222.tif";
              
          }
    
    public static String getPathToImageHighRes(String specimen,int step){
                return mainDir+specimen+"/raw/"+specimen+"_"+timestamps[step]+".tif";
                
    }
    
    public static void seeResultsOfRigidRegistration(String specimen){
        ImagePlus imgRef = IJ.openImage(getPathToImageHighRes(specimen,0));
        ImagePlus imgMov = IJ.openImage(getPathToImageHighRes(specimen,1));
        ItkTransform tr=ItkTransform.readTransformFromFile(mainDir+specimen+"/transforms_corrected/Automatic_Transform_2to0.txt");
        ImagePlus movRegistered=tr.transformImage(imgRef, imgMov);
        VitimageUtils.compositeNoAdjustOf(imgRef, movRegistered, "composite").show();
        

    }

    public RegistrationAutomate(){ }
    
    public RegistrationAutomate(String year1){    
        this.year1=year1;
    }
    
    public static String getPathToReferenceImage(String specimen){
        return mainDir+specimen+"/raw/"+specimen+"_J001_aligned.tif";
    }

    public static String getPathToMovingImage(String specimen){
        return mainDir+specimen+"/raw/"+specimen+"_J029_normalized.tif";
    }

    public static String getPathToTrInit(String specimen){
        return mainDir+specimen+"/transforms_corrected/Final_transform_001_029.txt";
    }

    public static String getPathToTrFinal(String specimen){
        return mainDir+specimen+"/transforms_corrected/testMaskedCambiumRegistration.txt";
    }

    public static String getPathToMask(String specimen){
        return mainDir+specimen+"/raw/"+specimen+"_mask_cambium.tif";
    }


    public static void runRigidRegistration(String specimen){

        System.out.println("Now running Rigid Registration for specimen: " + specimen);
        System.out.println("-----------------------------------------------");
    
        // Open Reference and Moving Images
        ImagePlus imgRef = IJ.openImage(getPathToReferenceImage(specimen));
        ImagePlus imgMov = IJ.openImage(getPathToMovingImage(specimen));
        // ImagePlus imgRef = IJ.openImage(getPathToImageLowRes(specimen,0));
        // System.out.println(getPathToImageLowRes(specimen,0));
        // ImagePlus imgMov = IJ.openImage(getPathToImageLowRes(specimen,1));
        // System.out.println(getPathToImageLowRes(specimen,1));
        imgRef.show();
        imgMov.show();

        // Get the global matrix from manual registration and run automatic registration step (rigid body)
        ItkTransform trInit=ItkTransform.readTransformFromFile(getPathToTrInit(specimen));
        ItkTransform trFinal = autoLinearRegistration(imgRef,imgMov, trInit,specimen);
        ImagePlus finalResult=trFinal.transformImage(imgRef,imgMov);

        // IJ.saveAsTiff(finalResult, mainDir+specimen+"/outputs/"+specimen+"_001_029_TransformedImage.tif");
        // System.out.println(mainDir+specimen+"/outputs/"+specimen+"_001_029_TransformedImage.tif");

        //IJ.saveAsTiff(finalResult,mainDir+specimen+"_XR_MRI/result_automate/"+specimen+"_TransformedImageXRMRI.tif");
        ImagePlus imgComposite=VitimageUtils.compositeNoAdjustOf(imgRef, finalResult);
        // IJ.saveAsTiff(imgComposite,mainDir+specimen+"/outputs/"+specimen+"_001_029_Composite.tif");
        // System.out.println(mainDir+specimen+"/outputs/"+specimen+"_001_029_Composite.tif");
        //IJ.saveAsTiff(imgComposite,mainDir+specimen+"_XR_MRI/result_automate/"+specimen+"_CompositeImageXRMRI.tif");g        imgComposite.show();
        imgComposite.close();
        imgRef.close();
        imgMov.close();
        
    }
    
    // Test automatic linear rigid-body registration
    public static ItkTransform autoLinearRegistration(ImagePlus imgRef, ImagePlus imgMov, ItkTransform trInit, String specimen){
        RegistrationAction regAct = new RegistrationAction();
        regAct.defineSettingsFromTwoImages(imgRef, imgMov, null, false);
        regAct.typeTrans=Transform3DType.RIGID;
        regAct.typeAutoDisplay=0;   
        regAct.higherAcc=0; 
        regAct.levelMaxLinear=3;
        regAct.levelMinLinear=3;
        regAct.bhsX=3;
        regAct.bhsY=3;
        regAct.bhsZ=3;
        regAct.strideX=8;
        regAct.strideY=8;
        regAct.strideZ=8;
        regAct.iterationsBMLin=8;
        regAct.neighX=2;
        regAct.neighX=2;
        regAct.neighZ=2;
        VitimageUtils.showWithParams(imgMov, specimen, 0, 0, 0);
        imgRef.show();
        imgMov.show();
        VitimageUtils.waitFor(5000);
        imgRef.hide();
        imgMov.hide();
        BlockMatchingRegistration bmRegistration = BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef, imgMov, regAct);
        bmRegistration.mask=IJ.openImage(getPathToMask(specimen));
        ItkTransform trFinal=bmRegistration.runBlockMatching(trInit, false);
        trFinal.writeMatrixTransformToFile(getPathToTrFinal(specimen));
        bmRegistration.closeLastImages();
        bmRegistration.freeMemory();
        
        return trFinal;
    }

    public static void runDenseRegistration(String specimen){
     
        System.out.println("Now running Dense Field registration for specimen: " + specimen);
        System.out.println("-----------------------------------------------");
    

        // Open Reference and Moving Images
        ImagePlus imgRef = IJ.openImage(getPathToReferenceImage(specimen));
        ImagePlus imgMov = IJ.openImage(getPathToMovingImage(specimen));
        

        // Get the global matrix from manual registration and run automatic registration step (dense field)
        ItkTransform trInit=ItkTransform.readTransformFromFile(getPathToTrInit(specimen));
        
        ItkTransform trFinal = autoNonLinearRegistration(imgRef,imgMov,trInit,specimen);
        ImagePlus finalResult=trFinal.transformImage(imgRef,imgMov);
        IJ.saveAsTiff(finalResult,mainDir+"cep_"+specimen+"/result_"+year1+"_"+year2+"_auto_reg/TransformedImage.tif");
        ImagePlus imgComposite=VitimageUtils.compositeNoAdjustOf(imgRef, finalResult);
        IJ.saveAsTiff(imgComposite,mainDir+"cep_"+specimen+"/result_"+year1+"_"+year2+"_auto_reg/CompositeImage.tif");
        imgComposite.show();
        imgComposite.close();
        imgRef.close();
        imgMov.close();
   
        
    }

    public static ItkTransform autoNonLinearRegistration(ImagePlus imgRef, ImagePlus imgMov, ItkTransform trInit, String specimen){
        RegistrationAction regAct = new RegistrationAction();
        regAct.defineSettingsFromTwoImages(imgRef, imgMov, null, false);
        regAct.typeTrans=Transform3DType.DENSE;
        regAct.typeAutoDisplay=2;   
        regAct.higherAcc=0; 
        regAct.levelMaxDense=3;
        regAct.levelMinDense=1;
        regAct.bhsX=3;
        regAct.bhsY=3;
        regAct.bhsZ=1;
        regAct.strideX=8;
        regAct.strideY=8;
        regAct.strideZ=8;
        regAct.sigmaDense=50;
        regAct.iterationsBMDen=5;
        regAct.neighX=2;
        regAct.neighX=2;
        regAct.neighZ=2;
        VitimageUtils.showWithParams(imgMov, specimen, 0, 0, 0);
        imgRef.show();
        imgMov.show();
        VitimageUtils.waitFor(5000);
        imgRef.hide();
        imgMov.hide();
        BlockMatchingRegistration bmRegistration = BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef, imgMov, regAct);
        bmRegistration.mask=IJ.openImage(getPathToMask(specimen));
        ItkTransform trFinal=bmRegistration.runBlockMatching(trInit, false);
        trFinal.writeAsDenseField(getPathToTrFinal(specimen),imgRef);
        bmRegistration.closeLastImages();
        bmRegistration.freeMemory();
        
        return trFinal;
    }


    public static String[] getSpecimen(String singleSpecimen) {
        String[] specimenList;
        if (singleSpecimen != null && !singleSpecimen.isEmpty()) {
            // Select one specimen
            specimenList = new String[]{singleSpecimen};
        } else {
            // Select the entire list
            specimenList = GeneralUtils.getSpecimenListXR();
        }
        return specimenList;
    }

    
}



        