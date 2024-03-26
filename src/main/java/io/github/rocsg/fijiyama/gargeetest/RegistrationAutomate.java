package io.github.rocsg.fijiyama.gargeetest;


import java.util.ArrayList;

//specific libraries
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.fijiyamaplugin.RegistrationAction;
import io.github.rocsg.fijiyama.registration.BlockMatchingRegistration;
import io.github.rocsg.fijiyama.registration.ItkTransform;
import io.github.rocsg.fijiyama.registration.Transform3DType;




public class RegistrationAutomate {
 
 
    /* The intend of this class is to automate actions of registration over a set of specimens
    * corresponding to 7 grapevine trunk, observed in 2022 in both MRI and Xray
    * The registration will be computed with Xray as fixed image an MRI as moving image.
    * The result (the ItkTransform) will be saved in some $COMPUTED_DATA
    * 
    */

    public  String mainDir="/home/phukon/Desktop/MRI/02_ceps/2022-03_CEPS_suivi_clinique/";// Path to the directory consisting the data
    public  String year1;
    public  String year2="2023";
    public  String specimen = "318";
    

    public static void main(String[] args) {
        ImageJ ij=new ImageJ();
        String year1="2022";
        RegistrationAutomate regAuto=new RegistrationAutomate(year1);
        regAuto.runStuff();
    }

    
    public String getPathToReferenceImage(){
        return mainDir+specimen+"_XR_MRI/"+specimen+"_imgRefcut.tif";
    }

    public String getPathToMovingImage(){
        return mainDir+specimen+"_XR_MRI/"+specimen+"_imgMov.tif";
    }

    public String getPathToTrInit(){
        return mainDir+specimen+"_XR_MRI/result_manual/Exported_data/transform_global_img_moving.txt";
    }

    public String getPathToTrFinal(){
        return mainDir+specimen+"_XR_MRI/result_automate/transform_automatic_rigid-body_step2.txt";
    }

    public RegistrationAutomate (String year1){
        this.year1=year1;
    }

    public String getPathToMask(){
        return mainDir+specimen+"_XR_MRI/mask.tif";
    }

    public void runStuff(){
     
        System.out.println("Now running the test for specimen: " + specimen);
    

        // Open XRay image and MRI image
        ImagePlus imgRef = IJ.openImage(getPathToReferenceImage());
        ImagePlus imgMov = IJ.openImage(getPathToMovingImage());

        // Get the global matrix from manual registration and run automatic registration step (rigid body auto)
        ItkTransform trInit=ItkTransform.readTransformFromFile(getPathToTrInit());
        this.autoLinearRegistration(imgRef,imgMov,trInit);
        
    }
    
    // Test automatic linear rigid-body registration
    public void autoLinearRegistration(ImagePlus imgRef, ImagePlus imgMov, ItkTransform trInit){
        RegistrationAction regAct = new RegistrationAction();
        regAct.defineSettingsFromTwoImages(imgRef, imgMov, null, false);
        regAct.typeTrans=Transform3DType.RIGID;
        regAct.typeAutoDisplay=2;   // Display at the end
        regAct.higherAcc=0; 
        regAct.levelMaxLinear=3;
        regAct.levelMinLinear=1;
        regAct.bhsX=3;
        regAct.bhsY=3;
        regAct.bhsZ=1;
        regAct.strideX=8;
        regAct.strideY=8;
        regAct.strideZ=3;
        regAct.iterationsBMLin=8;
        regAct.neighX=2;
        regAct.neighX=2;
        regAct.neighZ=1;
        BlockMatchingRegistration bmRegistration = BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef, imgMov, regAct);
        bmRegistration.mask=IJ.openImage(getPathToMask());
        ItkTransform trFinal=bmRegistration.runBlockMatching(trInit, false);
        bmRegistration.closeLastImages();

        trFinal.writeMatrixTransformToFile(getPathToTrFinal());
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



        