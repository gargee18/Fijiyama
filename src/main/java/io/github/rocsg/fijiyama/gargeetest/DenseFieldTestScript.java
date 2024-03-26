/*
 * Test script for understanding the functioning of dense field and the causes for singulairities and to be able to use mask appropriately
 * 
 */


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




public class DenseFieldTestScript{


    public static String mainDir="/home/phukon/Desktop/MRI/02_ceps/2022-03_CEPS_suivi_clinique/Registered_XR_MRI/";
   
    public static void main(String[]args) {

		ImageJ ij=new ImageJ();
        String singleSpecimen = "1193"; 
        ArrayList<String> specimenList = getSpecimen(singleSpecimen);
        for (String specimen : specimenList) {
            System.out.println("Now running the test for specimen: " + specimen);
            registerXRayandMRI(specimen);
        }
    }


    public static void registerXRayandMRI(String specimenName){
       
        // Open XRay image(ref) and MRI image(mov)
        ImagePlus imgRef = IJ.openImage(mainDir+specimenName+"/"+specimenName+"_imgRef.tif");
        ImagePlus imgMov = IJ.openImage(mainDir+specimenName+"/"+specimenName+"_imgMov.tif");
        ImagePlus mask = IJ.openImage(mainDir+specimenName+"/"+specimenName+"_maskTrunk.tif");
        ItkTransform trInit=ItkTransform.readTransformFromFile(mainDir+specimenName+"/"+specimenName+"_trRigid.txt");
        
        //Create a new registration action
        RegistrationAction regAct = new RegistrationAction();
        regAct.typeTrans=Transform3DType.DENSE;
        regAct.typeAutoDisplay=2;   
        regAct.higherAcc=0; 
        regAct.levelMaxDense=2;
        regAct.levelMinDense=2;
        regAct.bhsX=3;
        regAct.bhsY=3;
        regAct.bhsZ=1;
        regAct.strideX=8;
        regAct.strideY=8;
        regAct.strideZ=3;
        regAct.sigmaDense=50;
        regAct.iterationsBMDen=10;
        regAct.neighX=2;
        regAct.neighX=2;
        regAct.neighZ=1;
        

        ImagePlus imgMovAfterRigidBody=trInit.transformImage(imgRef, imgMov);
        ImagePlus imgFusAfterRigidBodyRegistration=VitimageUtils.compositeNoAdjustOf(imgRef, imgMovAfterRigidBody);
        imgFusAfterRigidBodyRegistration.show();
        imgFusAfterRigidBodyRegistration.setTitle("After rigid "+specimenName);


        // Start blockmatching
        BlockMatchingRegistration bmRegistration = BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef, imgMov, regAct);
        bmRegistration.mask=mask;
        ItkTransform trFinal=bmRegistration.runBlockMatching(trInit, false);
        bmRegistration.closeLastImages();
        trFinal.writeAsDenseField(mainDir + specimenName+"/"+specimenName+"_trDense.tif", imgRef);
        
        ImagePlus imgMovAfterRigidBodyAndDense=trFinal.transformImage(imgRef, imgMov);
        ImagePlus imgFusAfterRigidBodyAndDenseRegistration=VitimageUtils.compositeNoAdjustOf(imgRef, imgMovAfterRigidBodyAndDense);
        imgFusAfterRigidBodyAndDenseRegistration.show();
        imgFusAfterRigidBodyAndDenseRegistration.setTitle("After dense "+specimenName);


        // Apply final ITK transform and save as tiff
        IJ.saveAsTiff(imgMovAfterRigidBodyAndDense,mainDir + specimenName+"/"+specimenName+"_regDense_sigma50.tif");
        System.out.println("Finished");
        

    }

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
        specimenList.add("318");
        specimenList.add("322");
        specimenList.add("323");
        specimenList.add("330");
        specimenList.add("335");
        specimenList.add("764B");
        specimenList.add("1181");
        specimenList.add("1193");
        return  specimenList;
    }
}



        
