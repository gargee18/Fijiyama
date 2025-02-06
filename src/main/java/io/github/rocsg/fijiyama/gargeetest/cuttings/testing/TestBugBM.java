package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.fijiyamaplugin.RegistrationAction;
import io.github.rocsg.fijiyama.registration.BlockMatchingRegistration;
import io.github.rocsg.fijiyama.registration.ItkTransform;
import io.github.rocsg.fijiyama.registration.Transform3DType;


public class TestBugBM {
    public static void main(String[] args) {
        ImageJ ij = new ImageJ();
        String pathToMask="/home/gargee/Bureau/TestBugBM/Mask.tif";
        ImagePlus imgRef=IJ.openImage("/home/gargee/Bureau/TestBugBM/Ref.tif");
        ImagePlus imgMov=IJ.openImage("/home/gargee/Bureau/TestBugBM/Mov.tif");
        ImagePlus imgMask=IJ.openImage(pathToMask);
        ImagePlus imgRefCopy=VitimageUtils.imageCopy(imgRef);
        boolean testing=true;
        boolean tweaking=false;

        if(tweaking){
            VitimageUtils.adjustImageCalibration(imgRef, new double[]{1,1,1}, "mm");
            VitimageUtils.adjustImageCalibration(imgMov, new double[]{1,1,1}, "mm");
            VitimageUtils.adjustImageCalibration(imgMask, new double[]{1,1,1}, "mm");
        }
        RegistrationAction regAct = new RegistrationAction();
        regAct.defineSettingsFromTwoImages(imgRef, imgMov, null, false);
        regAct.typeTrans=Transform3DType.RIGID;
        regAct.typeAutoDisplay=(testing ? 2 : 0);   
        regAct.higherAcc=0; 
        regAct.levelMaxLinear=3;
        regAct.levelMinLinear=(testing ? 2 : 2);
        regAct.bhsX=3;
        regAct.bhsY=3;
        regAct.bhsZ=3;
        regAct.strideX=7;
        regAct.strideY=7;
        regAct.strideZ=7;
        regAct.iterationsBMLin=(testing ? 2 : 8);
        regAct.neighX=2;
        regAct.neighX=2;
        regAct.neighZ=2;

        BlockMatchingRegistration bmRegistration = BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef, imgMov, regAct);
        bmRegistration.mask=imgMask;
        VitimageUtils.printImageResume(imgRef);
        VitimageUtils.printImageResume(imgMov);
        VitimageUtils.printImageResume(bmRegistration.mask);
        ItkTransform trFinal=bmRegistration.runBlockMatching(null, false);
        bmRegistration.closeLastImages();
        bmRegistration.freeMemory();
        ImagePlus imgMovTransformed=trFinal.transformImage(imgRef, imgMov);
        if(tweaking){
            VitimageUtils.copyImageCalibrationAndRange(imgMovTransformed,imgRefCopy);
            VitimageUtils.copyImageCalibrationAndRange(imgMov,imgRefCopy);
            VitimageUtils.copyImageCalibrationAndRange(imgRef,imgRefCopy);
        }
        ImagePlus composed=VitimageUtils.compositeNoAdjustOf(imgRef, imgMovTransformed,"composite");
        composed.show();
    }   
}
