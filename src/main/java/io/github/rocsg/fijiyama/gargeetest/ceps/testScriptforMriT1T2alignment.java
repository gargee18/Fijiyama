/*
 * Test script
 * 
 */


package io.github.rocsg.fijiyama.gargeetest.ceps;


import java.io.File;
//specific libraries
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.fijiyamaplugin.RegistrationAction;
import io.github.rocsg.fijiyama.registration.BlockMatchingRegistration;
import io.github.rocsg.fijiyama.registration.ItkTransform;
import io.github.rocsg.fijiyama.registration.Transform3DType;
import io.github.rocsg.fijiyama.registration.TransformUtils;




public class testScriptforMriT1T2alignment{

    
    public static String year = GeneralUtils.years[1];
    public static String mainDir="/home/phukon/Desktop/MRI/02_ceps/"+year+"_CEPS_suivi_clinique/";
    public static void main(String[]args) {

		ImageJ ij=new ImageJ();
        
        String singleSpecimen = null; 
        int[] indexSpec = {0,1,2,3,4,5,6,7};
        String[] specimenList = getSpecimen(singleSpecimen);
        for(int i = 0; i<specimenList.length; i++){
            System.out.println("Now running the test for specimen: " + specimenList[i] +" with index for translation: "+ i);
            registerXRayandMRI(specimenList[i],indexSpec[i]);
        }

    }

    static double[]requiredAdditionalTranslationForWeirdOneX=new double[]{ 
        0 /* 318*/,
        0 /* 322*/,
        0 /* 323*/,
        23 /* 330*/,
        41 /* 335*/,
        0 /* 764B*/,
        -5 /* 1181*/,
        0 /* 1193*/,
    };

    static double[]requiredAdditionalTranslationForWeirdOneY=new double[]{
        28 /* 318*/,
        0 /* 322*/,
        0 /* 323*/,
        17 /* 330*/,
        33 /* 335*/,
        0 /* 764B*/,
        -5 /* 1181*/,
        23 /* 1193*/,
    };
    // coords for 2023
    static double[]requiredAdditionalTranslationForWeirdOneZ=new double[]{
        -18 /* 318*/,
        -30 /* 322*/,
        -25 /* 323*/,
        -20 /* 330*/,
        -23 /* 335*/,
        -23 /* 764B*/,
        -23 /* 1181*/,
        -23.5 /* 1193*/,
    };

    // coords for 2024
    // static double[]requiredAdditionalTranslationForWeirdOneZ=new double[]{
    //     -0 /* 318*/,
    //     0 /* 322*/,
    //     0 /* 323*/,
    //     -10 /* 330*/,
    //     0 /* 335*/,
    //     -18 /* 764B*/,
    //     -5 /* 1181*/,
    //     0 /* 1193*/,
    // };

    public static void registerXRayandMRI(String specimenName, int indexSpec){

        
        // Open XRay image(ref) and MRI image(mov) or MRI T1(ref) and MRI T2(mov)
        String pathtoRef = mainDir+"raw_T1_T2/CEP_"+specimenName+"_"+year+"_MRI_T1.tif";
        System.out.println(pathtoRef);
        ImagePlus imgRef = IJ.openImage(mainDir+"raw_T1_T2/CEP_"+specimenName+"_"+year+"_MRI_T1.tif");
        ImagePlus imgMov = IJ.openImage(mainDir+"raw_T1_T2/CEP_"+specimenName+"_"+year+"_MRI_T2.tif");


        //Create a new registration action
        RegistrationAction regAct = new RegistrationAction();
        regAct.typeTrans=Transform3DType.TRANSLATION;
        regAct.typeAutoDisplay=2;   
        regAct.higherAcc=0; 
        regAct.levelMaxLinear=3;
        regAct.levelMinLinear=1;
        regAct.bhsX=2;
        regAct.bhsY=2;
        regAct.bhsZ=1;
        regAct.strideX=8;
        regAct.strideY=8;
        regAct.strideZ=8;
        regAct.neighX=2;
        regAct.neighX=2;
        regAct.neighZ=2;
        regAct.iterationsBMLin = 8;
        //regAct.typeAction=RegistrationAction.TYPEACTION_EVALUATE;


        // Start blockmatching
        BlockMatchingRegistration bmRegistration = BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef, imgMov, regAct);
        //bmRegistration.mask=mask;
        ItkTransform trAddAlignementZ=ItkTransform.array16ElementsToItkTransform(new double[]{1,0,0,requiredAdditionalTranslationForWeirdOneX[indexSpec],0,1,0,requiredAdditionalTranslationForWeirdOneY[indexSpec],0,0,1,requiredAdditionalTranslationForWeirdOneZ[indexSpec]});
        ItkTransform trFinal=bmRegistration.runBlockMatching(trAddAlignementZ, false);
        bmRegistration.closeLastImages();
        trFinal.writeMatrixTransformToFile(mainDir + "MRI_T1_T2_reg/"+specimenName+"_T1_T2_"+year+".txt");
        trFinal.writeMatrixTransformToFile(mainDir + "_T1_T2_"+year+".txt");
        ImagePlus imgMovAfterRigidBody=trFinal.transformImage(imgRef, imgMov);
        //imgMovAfterRigidBody.show();
        imgMovAfterRigidBody.setTitle("Rigid registration of "+specimenName);
        

        IJ.saveAsTiff(imgMovAfterRigidBody,mainDir + "MRI_T1_T2_reg/"+specimenName+"_TransformedImage.tif");
        ImagePlus imgComposite=VitimageUtils.compositeNoAdjustOf(imgRef, imgMovAfterRigidBody);
        IJ.saveAsTiff(imgComposite,mainDir + "MRI_T1_T2_reg/"+specimenName+"_CompositeImage.tif");

        bmRegistration.freeMemory();
        

    }

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

    // public static String[] getSpecimenList(){
    //     String[] specimenList=new String[]{
    //         "318",
    //         "322",
    //         "323",
    //         "330",
    //         "335",
    //         "764B",
    //         "1181",
    //         "1193",
        
    //     };
    //     return  specimenList;
    // }

    public static double getSlicesToPerformVoxelCalculation(String mod, String specimen)
    {
        String inputDir="/home/phukon/Desktop/MRI/02_ceps/2024_CEPS_suivi_clinique/0_RAW_data_MRI_"+year+"/"+specimen+"/"+mod;
        File f=new File(inputDir);
        String[]tabNames=f.list();
        for (String tab: tabNames) {
            System.out.println(tab);
            
        }
        System.out.println("First slice = "+tabNames[0]);
        System.out.println("First slice = "+tabNames[tabNames.length - 1]);


        int nbSlices=Integer.parseInt(tabNames[0])-Integer.parseInt(tabNames[tabNames.length - 1]);
        System.out.println(nbSlices);
        String pathToLastSlice=tabNames[tabNames.length - 1]; //gets name of the image
        String pathToFirstSlice=tabNames[0];
        Object[] vecs = TransformUsingMetadata.getTrMatrices(inputDir+pathToFirstSlice, inputDir+pathToLastSlice);
        double[]tabPosFirst=(double[])(vecs[0]);//Getmetadatafrom file in pathToFirstSlice
        double[]tabPosLast=(double[])(vecs[3]);//Getmetadatafrom file in pathToLastSlice
        double[]vectDelta=TransformUtils.vectorialSubstraction(tabPosLast, tabPosFirst);
        double[]vectDeltaPerSlice=TransformUsingMetadata.divideVector(vectDelta, nbSlices);  
        double voxelSizeZ=TransformUtils.norm(vectDeltaPerSlice);
        return voxelSizeZ;
    }
}



        
