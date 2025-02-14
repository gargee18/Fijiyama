package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import inra.ijpb.morphology.Morphology;
import io.github.rocsg.fijiyama.common.ItkImagePlusInterface;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.fijiyamaplugin.RegistrationAction;
import io.github.rocsg.fijiyama.registration.BlockMatchingRegistration;
import io.github.rocsg.fijiyama.registration.ItkTransform;
import io.github.rocsg.fijiyama.registration.Transform3DType;


public class TestRegistrationWithMask {
    
    public static void main(String[] args) {
        ImageJ ij = new ImageJ();

//        test1();
        test1();
        // Images
    }
        public static void test2(){
            ImagePlus imgRef=IJ.openImage("/home/gargee/Bureau/test/B_201_J077_normalized.tif");
            ImagePlus imgMov=IJ.openImage("/home/gargee/Bureau/test/B_201_J141_normalized.tif");
            imgRef.show();
            imgMov.show();
            ItkTransform trGood = ItkTransform.readTransformFromFile("/home/gargee/Bureau/test/J077_J141_TR_MAT_RigidReg.txt");
            ItkTransform trBad = ItkTransform.readTransformFromFile("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/Cuttings/Processing/02_RigidRegistration/B_201_J077_J141_TR_MAT_Rigid_Reg.txt");
    
    
            // Inoc Alignment/^C
            ItkTransform trInocRef= ItkTransform.readTransformFromFile("/home/gargee/Bureau/test/B_201_J077_TR_MAT_Inoculation_Alignment.txt");
            ItkTransform trInocMov= ItkTransform.readTransformFromFile("/home/gargee/Bureau/test/B_201_J141_TR_MAT_Inoculation_Alignment.txt");
    
            ImagePlus imgRefToInoc=trInocRef.transformImage(imgRef,imgRef);
            imgRefToInoc.setTitle("Ref");
            imgRefToInoc.show();
            ImagePlus imgMovToInoc=trInocMov.transformImage(imgMov,imgMov);

            ImagePlus imgMovToInocBad=trBad.transformImage(imgRefToInoc, imgMovToInoc);
            ImagePlus imgMovToInocGood=trGood.transformImage(imgRefToInoc, imgMovToInoc);

            imgMovToInocBad.setTitle("Bad");
            imgMovToInocBad.show();
            imgMovToInocGood.setTitle("Good");
            imgMovToInocGood.show();
        }
    


        public static void test1(){
        ImagePlus imgRef=IJ.openImage("/home/gargee/Bureau/test/B_201_J029_sub.tif");
        ImagePlus imgMov=IJ.openImage("/home/gargee/Bureau/test/B_201_J077_sub.tif");
        ImagePlus imgMask=IJ.openImage("/home/gargee/Bureau/test/B_201_mask.tif");
        imgMask=Morphology.dilation(imgMask, inra.ijpb.morphology.strel.DiskStrel.fromRadius(10));
        imgRef.show();
        imgMov.show();
        imgMask.show();


        // Inoc Alignment
        ItkTransform trInocRef= ItkTransform.readTransformFromFile("/home/gargee/Bureau/test/B_201_J077_TR_MAT_Inoculation_Alignment.txt");
        ItkTransform trInocMov= ItkTransform.readTransformFromFile("/home/gargee/Bureau/test/B_201_J141_TR_MAT_Inoculation_Alignment.txt");

        ImagePlus imgRefToInoc=trInocRef.transformImage(imgRef,imgRef);
        ImagePlus imgMovToInoc=trInocMov.transformImage(imgMov,imgMov);

        imgRefToInoc.show();
        imgMovToInoc.show();
        imgRefToInoc.setTitle("Ref");
        imgMovToInoc.setTitle("Mov");

        // Rigid reg
        RegistrationAction regAct = new RegistrationAction();
        regAct.defineSettingsFromTwoImages(imgRefToInoc, imgMovToInoc, null, false);
        regAct.typeTrans=Transform3DType.RIGID;
        regAct.typeAutoDisplay=2;   
        regAct.higherAcc=0; 
        regAct.levelMaxLinear=2;
        regAct.levelMinLinear= 1;
        regAct.bhsX=3;
        regAct.bhsY=3;
        regAct.bhsZ=3;
        regAct.strideX=3;
        regAct.strideY=3;
        regAct.strideZ=3;
        regAct.iterationsBMLin= 8;
        regAct.neighX=2;
        regAct.neighX=2;
        regAct.neighZ=2;

        BlockMatchingRegistration bmRegistration = BlockMatchingRegistration.setupBlockMatchingRegistration(imgRefToInoc, imgMovToInoc, regAct);
        bmRegistration.mask=imgMask;
        VitimageUtils.printImageResume(imgRefToInoc);
        VitimageUtils.printImageResume(imgMovToInoc);
        VitimageUtils.printImageResume(bmRegistration.mask);
        BlockMatchingRegistration.setRandomSelection();
        ItkImagePlusInterface.setSpeedupOnForFiji();
        ItkTransform trFinal=bmRegistration.runBlockMatching(null, false);
        bmRegistration.closeLastImages();
        bmRegistration.freeMemory();
        ImagePlus imgMovTransformed=trFinal.transformImage(imgRefToInoc, imgMovToInoc);
        ImagePlus composed=VitimageUtils.compositeNoAdjustOf(imgRefToInoc, imgMovTransformed,"composite");
        trFinal.writeMatrixTransformToFile("/home/gargee/Bureau/test/J077_J141_TR_MAT_RigidReg.txt");
        composed.show();


        //Test on high res
        ImagePlus imgRefHiRes=IJ.openImage("/home/gargee/Bureau/test/B_201_J077_normalized.tif");
        ImagePlus imgMovHiRes=IJ.openImage("/home/gargee/Bureau/test/B_201_J141_normalized.tif");

        ImagePlus imgRefToInocHiRes=trInocRef.transformImage(imgRefHiRes,imgRefHiRes);
        ImagePlus imgMovToInocHiRes=trInocMov.transformImage(imgMovHiRes,imgMovHiRes);


        ItkTransform tr = ItkTransform.readTransformFromFile("/home/gargee/Bureau/test/J077_J141_TR_MAT_RigidReg.txt");
        ImagePlus imgFinal =tr.transformImage(imgRefToInocHiRes,imgMovToInocHiRes);
        ImagePlus composedFinal=VitimageUtils.compositeNoAdjustOf(imgRefToInoc, imgFinal,"composite");
        composedFinal.show();
    }   
}
