
/* Actual comment :
 * 
 *  The purpose of this step is to perform BM registration from time i to time i+1 (1 reg on 0 , 2 on 1... )
 *  and save the corresponding matrices
 * 
 * 
 */


package io.github.rocsg.fijiyama.gargeetest.cuttings.processing;


//specific libraries
import ij.IJ;
// import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.fijiyamaplugin.RegistrationAction;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.registration.BlockMatchingRegistration;
import io.github.rocsg.fijiyama.registration.ItkTransform;
import io.github.rocsg.fijiyama.registration.Transform3DType;



public class Step_3_RegistrationRigid implements PipelineStep{
    public static void main(String[] args) throws Exception{
        Specimen spec= new Specimen("B_201");
        new Step_3_RegistrationRigid().execute(spec,true); 
        // seeResultsOfRigidRegistration(spec, 0);
    }
 

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }
    
    
    public void execute(Specimen specimen,boolean testing) throws Exception {
        String[] days=Config.timestamps;

        int N=days.length;

        ItkTransform trRigid0PlusInoc = null;
        ItkTransform trRigid1PlusRigid0PlusInoc = null;
        ItkTransform trRigid2PlusRigid1PlusRigid0PlusInoc = null;
       
        for(int i=0;i<N-1;i++){
            int indexRef=i;
            int indexMov=i+1;
            
            //Open images
            ImagePlus imgRef=IJ.openImage(Config.getPathToSubsampledImage(specimen,indexRef));
            ImagePlus imgMov=IJ.openImage(Config.getPathToSubsampledImage(specimen,indexMov));
            ImagePlus imgRefToInoc=imgRef;
            // Get path to mask
            String mask = Config.getPathToMask(specimen, indexRef);

            if (i == 0) {
                //Get Inoc Matrix
                imgRef.show();
                ItkTransform trInoc= ItkTransform.readTransformFromFile(Config.getPathToInoculationAlignmentTransformation(specimen,indexRef));
                //Apply the inoc transform to the reference one.
                imgRefToInoc=trInoc.transformImage(imgRef,imgRefToInoc);
                // Register and obtain transformation matrix of J029 on J001
                ItkTransform trRigid0 = autoLinearRegistrationWithPith(imgRefToInoc, imgMov, trInoc, specimen, mask);
                trRigid0PlusInoc = trInoc.addTransform(trRigid0);
                trRigid0PlusInoc.writeMatrixTransformToFile(Config.getPathToRigidRegistrationMatrix(specimen,indexRef, indexMov));
            }else if (i == 1) {
                //Get Inoc Matrix
                imgRef.show();
                ItkTransform trInoc= ItkTransform.readTransformFromFile(Config.getPathToInoculationAlignmentTransformation(specimen,indexRef));
                imgRefToInoc = trRigid0PlusInoc.transformImage(imgRef,imgRefToInoc);
                // Register and obtain transformation matrix of J077 on J029
                ItkTransform trRigid1 = autoLinearRegistrationWithPith(imgRefToInoc, imgMov, trInoc, specimen, mask);
                // Compose transformation matrices (J077 + J029 + J001)
                trRigid1PlusRigid0PlusInoc = trRigid0PlusInoc.addTransform(trRigid1);
                trRigid1PlusRigid0PlusInoc.writeMatrixTransformToFile(Config.getPathToRigidRegistrationMatrix(specimen,indexRef, indexMov));
            }else if (i == 2) {
                //Get Inoc Matrix
                imgRef.show();
                ItkTransform trInoc= ItkTransform.readTransformFromFile(Config.getPathToInoculationAlignmentTransformation(specimen,indexRef));
                imgRefToInoc = trRigid1PlusRigid0PlusInoc.transformImage(imgRef,imgRefToInoc);
                imgRefToInoc.show();
                // Register and obtain transformation matrix of J141 on J077
                ItkTransform trRigid2 = autoLinearRegistrationWithPith(imgRefToInoc, imgMov, trInoc, specimen, mask);
                // Compose transformation matrices (J141 + J077 + J029 + J001)
                trRigid2PlusRigid1PlusRigid0PlusInoc = trRigid1PlusRigid0PlusInoc.addTransform(trRigid2);
                ImagePlus fi = trRigid2PlusRigid1PlusRigid0PlusInoc.transformImage(imgRefToInoc,imgMov);
                fi.show();
                trRigid2PlusRigid1PlusRigid0PlusInoc.writeMatrixTransformToFile(Config.getPathToRigidRegistrationMatrix(specimen,indexRef, indexMov));
            } 
            
        }

    }
    

    public static ItkTransform autoLinearRegistrationWithPith(ImagePlus imgRef, ImagePlus imgMov, ItkTransform trInit, Specimen specimen, String pathToMask){
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
        regAct.iterationsBMLin=1;
        regAct.neighX=2;
        regAct.neighX=2;
        regAct.neighZ=2;
        BlockMatchingRegistration bmRegistration = BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef, imgMov, regAct);
        bmRegistration.mask=IJ.openImage(pathToMask);
        bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=false;
        ItkTransform trFinal=bmRegistration.runBlockMatching(trInit, false);
        bmRegistration.closeLastImages();
        bmRegistration.freeMemory();
        return trFinal;
    }


  //What is done : we have ItkTransform for each single time of each single specimen. 
        //By applying this, everybody have its inoculation point at 256,380, and that's so cool !
        //But :
        //TODO 1 : The z alignment is not perfect. We can leverage the inner structures (in the pith for running an automatic registration)
        //TODO 2 : The "inoculation plane" is not horizontal. One could make a click stuff in order to realign this. 
        //TODO 3 : For having a proper alignment of inner tissues, some dense registration is required (with some big sigma, and by ckecking
        //  carefuillly that it does not corrupt the evolving part, including disappearing area and growing cambium)
        //Caution : the transformations have to be computed in an order that makes it possible to compose them



    public static ItkTransform autoLinearRegistrationOlder(ImagePlus imgRef, ImagePlus imgMov, String pathToTrInit, String specimen, String pathToMask){
        RegistrationAction regAct = new RegistrationAction();
        regAct.defineSettingsFromTwoImages(imgRef, imgMov, null, false);
        regAct.typeTrans=Transform3DType.RIGID;
        regAct.typeAutoDisplay=0;   
        regAct.higherAcc=0; 
        regAct.levelMaxLinear=3;
        regAct.levelMinLinear=1;
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
        BlockMatchingRegistration bmRegistration = BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef, imgMov, regAct);
        bmRegistration.mask=IJ.openImage(pathToMask);
        ItkTransform trInit=ItkTransform.readTransformFromFile(pathToTrInit);
        ItkTransform trFinal=bmRegistration.runBlockMatching(trInit, false);
        bmRegistration.closeLastImages();
        bmRegistration.freeMemory();
        return trFinal;
    }


    public static void seeResultsOfRigidRegistration(Specimen specimen, int step){
        int indexRef = step;
        int indexMov = step+1;
        ImagePlus imgRef = IJ.openImage(Config.getPathToRawImage(specimen,indexRef));
        ImagePlus imgRefToInoc=imgRef;
        ItkTransform trInoc = ItkTransform.readTransformFromFile(Config.getPathToInoculationAlignmentTransformation(specimen,indexRef));
        imgRefToInoc=trInoc.transformImage(imgRef,imgRefToInoc);
        ImagePlus imgMov = IJ.openImage(Config.getPathToRawImage(specimen,indexMov));
        ItkTransform tr=ItkTransform.readTransformFromFile(Config.getPathToRigidRegistrationMatrix(specimen, indexRef, indexMov));
        ImagePlus movRegistered=tr.transformImage(imgRefToInoc, imgMov);
        VitimageUtils.compositeNoAdjustOf(imgRefToInoc, movRegistered, "composite").show();
    }
}


         