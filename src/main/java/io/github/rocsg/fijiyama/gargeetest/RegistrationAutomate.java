package io.github.rocsg.fijiyama.gargeetest;


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
    public static String mainDir = "/home/phukon/Desktop/mrit1_2022_mrit2_2022/";
    public static String year1 = GeneralUtils.years[0];
    public static String year2 = GeneralUtils.years[1];
   
   
    
    // public RegistrationAutomate(RegistrationAutomate otherRegToBeCopied){    this.year1=otherRegToBeCopied.year1; }

    public static void main(String[] args) {
        ImageJ ij=new ImageJ();
        //String year1="2023";
        //String chooseSpecimen = "1181";

    
        String[] specimenList = GeneralUtils.getSpecimenListMRI();
        for (String specimen : specimenList) {
             runRigidRegistration(specimen);
        }
        // runRigidRegistration(chooseSpecimen);
        //System.exit(0);
    }

    public RegistrationAutomate(){ }
    
    public RegistrationAutomate(String year1){    
        this.year1=year1;
    }
    
    public static String getPathToReferenceImage(String specimen){
        System.out.println(mainDir+"cep_"+specimen+"/raw/"+specimen+"_imgRef.tif");
        //System.out.println(mainDir+"cep_"+specimen+"/raw/"+specimen+"_"+year1+"_crop_sub_z.tif");
        //return mainDir+"cep_"+specimen+"/raw/"+specimen+"_"+year1+"_crop_sub_z.tif"; //
        return mainDir+"cep_"+specimen+"/raw/"+specimen+"_imgRef.tif";
    }

    public static String getPathToMovingImage(String specimen){
        //System.out.println(mainDir+"cep_"+specimen+"/raw/"+specimen+"_"+year2+"_crop_sub_z.tif");
        System.out.println(mainDir+"cep_"+specimen+"/raw/"+specimen+"_imgMov.tif");
        //return mainDir+"cep_"+specimen+"/raw/"+specimen+"_"+year2+"_crop_sub_z.tif"; 
        return mainDir+"cep_"+specimen+"/raw/"+specimen+"_imgMov.tif";
    }

    public static String getPathToTrInit(String specimen){
        System.out.println(mainDir+"cep_"+specimen+"/res/Exported_data/transform_global_img_moving.txt");
        //return mainDir+"cep_"+specimen+"/result_"+year1+"_"+year2+"_manual_reg/Exported_data/transform_global_img_moving.txt";

        
       
        return mainDir+"cep_"+specimen+"/res/Exported_data/transform_global_img_moving.txt";
    }

    public static String getPathToTrFinal(String specimen){
       
        System.out.println(mainDir+"cep_"+specimen+"/res/Exported_data/"+specimen+"_trMatrix_2023_2024.txt");
        //return mainDir+specimen+"_XR_MRI/result_automate/"+specimen+"_trMatrix.txt";
        return mainDir+"cep_"+specimen+"/res/Exported_data/"+specimen+"_trMatrix_2023_2024.txt";
    }

 
    public static String getPathToMask(String specimen){
        System.out.println(mainDir+"cep_"+specimen+"/raw/"+specimen+"_mask.tif");
        return mainDir+"cep_"+specimen+"/raw/"+specimen+"_mask.tif";
    }

    

    public static void runRigidRegistration(String specimen){
     
        System.out.println("Now running Rigid Registration for specimen: " + specimen);
        System.out.println("-----------------------------------------------");
    

        // Open Reference and Moving Images
        ImagePlus imgRef = IJ.openImage(getPathToReferenceImage(specimen));
        ImagePlus imgMov = IJ.openImage(getPathToMovingImage(specimen));
        

        // Get the global matrix from manual registration and run automatic registration step (rigid body)
        ItkTransform trInit=ItkTransform.readTransformFromFile(getPathToTrInit(specimen));
        
        ItkTransform trFinal = autoLinearRegistration(imgRef,imgMov, trInit,specimen);
        ImagePlus finalResult=trFinal.transformImage(imgRef,imgMov);
        IJ.saveAsTiff(finalResult, mainDir+"cep_"+specimen+"/res/"+specimen+"_TransformedImage.tif");
        System.out.println(mainDir+"cep_"+specimen+"/res/"+specimen+"_TransformedImage.tif");
        //IJ.saveAsTiff(finalResult,mainDir+specimen+"_XR_MRI/result_automate/"+specimen+"_TransformedImageXRMRI.tif");
        ImagePlus imgComposite=VitimageUtils.compositeNoAdjustOf(imgRef, finalResult);
        IJ.saveAsTiff(imgComposite,mainDir+"cep_"+specimen+"/res/"+specimen+"_Composite.tif");
        System.out.println(mainDir+"cep_"+specimen+"/res/"+specimen+"_Composite.tif");
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
        regAct.typeAutoDisplay=2;   
        regAct.higherAcc=0; 
        regAct.levelMaxLinear=3;
        regAct.levelMinLinear=1;
        regAct.bhsX=3;
        regAct.bhsY=3;
        regAct.bhsZ=1;
        regAct.strideX=8;
        regAct.strideY=8;
        regAct.strideZ=8;
        regAct.iterationsBMLin=12;
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
        regAct.iterationsBMDen=8;
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



        