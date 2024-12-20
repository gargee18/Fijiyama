
/* Actual comment :
 * 
 *  The purpose of this step is to perform BM registration from time i to time i+1 (1 reg on 0 , 2 on 1... )
 *  and save the corresponding matrices
 * 
 * 
 */


//CRAPPY SHIT
/* The intend of this class is to automate actions of registration (rigid) over a set of specimens
    * corresponding to 20 grapevine trunk, observed in 2022 in both MRI(7) and Xray(20)
    * The registration will be computed with Xray as fixed image an MRI as moving image or Xray of year 1 as 
    * fixed and Xray of year 2 as moving or MRI T1 as fixed and MRI T2 as moving for the same year 
    * autoLinearRegistration performs Rigid registration and autoNonLinearRegistration performs Dense field registration. Both return the transformation matrix. 
    * The result (the ItkTransform) is saved in a designated directory.
    * runRigidRegistration and runDenseRegistration calls its respective registration methods and saves the transformed moving image
    * and composite image to a designated directory. 
    */
//END OF CRAP


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



public class Step_2_RegistrationRigid implements PipelineStep{
    public static void main(String[] args) {
        test();
 
    }
    
    public static void test(){

    }


    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }


    public void execute(Specimen specimen,boolean testing) throws Exception {
        String[] days=Config.timestamps;

        int N=days.length;
        int nbReg=N-1;

        //What is done : we have ItkTransform for each single time of each single specimen. 
        //By applying this, everybody have its inoculation point at 256,380, and that's so cool !
        //But :
        //TODO 1 : The z alignment is not perfect. We can leverage the inner structures (in the pith for running an automatic registration)
        //TODO 2 : The "inoculation plane" is not horizontal. One could make a click stuff in order to realign this. 
        //TODO 3 : For having a proper alignment of inner tissues, some dense registration is required (with some big sigma, and by ckecking
        //  carefuillly that it does not corrupt the evolving part, including disappearing area and growing cambium)
        //Caution : the transformations have to be computed in an order that makes it possible to compose them


        for(int i=0;i< N-2;i++){
            int indexRef=i;
            int indexMov=i+1;

            //Open images
            ImagePlus imgRef=IJ.openImage(Config.getPathToSpecimenImageAtTimeStamp(specimen,i));
            imgRef=VitimageUtils.Sub222(imgRef);
            ImagePlus imgMov=IJ.openImage(Config.getPathToSpecimenImageAtTimeStamp(specimen,i));
            imgMov=VitimageUtils.Sub222(imgMov);


            //Apply the inoc transform to the reference one.
            ItkTransform tr=null;//TODO

            //Run a rigid body blockMatching limiting to the pith area from Imgmov to imgreftoinoc , by initializing with transform mov to inoc
            //And by setting the output transform to be of no composed style

            //Save the "blue" transform.
            
            //execute linear Reg and keep the resulting transform
            autoLinearRegistrationWithPith(imgRef, imgMov, String pathToTrInit, String specimen, String pathToMask){
        
        

            //save it in a proper place, with a proper name.
        }

    }
    

    int X=image.getWidth();
    int Y=image.getHeight();
    int Z=image.getNSlices();
    //Resize to 256x256x512 and convert to 8bit
    image = image.resize(X/subsampleRatioStandard, Y/subsampleRatioStandard, Z/subsampleRatioStandard, "bilinear");
    ImageConverter.setDoScaling(true);
    IJ.run(image, "8-bit", "");







    public static ItkTransform autoLinearRegistrationWithPith(ImagePlus imgRef, ImagePlus imgMov, String pathToTrInit, String specimen, String pathToMask){
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
       bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=false;
        ItkTransform trFinal=bmRegistration.runBlockMatching(trInit, false);
        bmRegistration.closeLastImages();
        bmRegistration.freeMemory();
        return trFinal;
    }





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


    public static void seeResultsOfRigidRegistration(String specimen){
        ImagePlus imgRef = IJ.openImage(Config.getPathToImageHighRes(specimen,0));
        ImagePlus imgMov = IJ.openImage(Config.getPathToImageHighRes(specimen,3));
        ItkTransform tr=ItkTransform.readTransformFromFile(Config.mainDir+specimen+"/transforms_corrected/Automatic_Transform_3to0.txt");
        ImagePlus movRegistered=tr.transformImage(imgRef, imgMov);
        VitimageUtils.compositeNoAdjustOf(imgRef, movRegistered, "composite").show();
        

    }
}


             /* 
        for (int i = 201; i <= 201; i++) {
            if (i == 204) {
                continue;
            }
            String specimen = "B_" + i;
            ImagePlus imgRef = IJ.openImage(Config.getPathToImageLowRes(specimen,0));
            ImagePlus imgMov = IJ.openImage(Config.getPathToImageLowRes(specimen,3));
            String pathToTrInit = Config.getPathToTrInit(specimen);
            String pathToMask = Config.getPathToMask(specimen);
            // runRigidRegistration(specimen, imgRef, imgMov, pathToTrInit, pathToMask);
            autoLinearRegistration(imgRef,imgMov, pathToTrInit,specimen, pathToMask);
            // save tr
            // IJ.saveAsTiff(finalResult, mainDir+specimen+"/outputs/"+specimen+"_001_029_TransformedImage.tif");
            // System.out.println(mainDir+specimen+"/outputs/"+specimen+"_001_029_TransformedImage.tif");

            // IJ.saveAsTiff(finalResult,mainDir+specimen+"_XR_MRI/result_automate/"+specimen+"_TransformedImageXRMRI.tif");
            // IJ.saveAsTiff(imgComposite,mainDir+specimen+"/outputs/"+specimen+"_001_029_Composite.tif");
            // System.out.println(mainDir+specimen+"/outputs/"+specimen+"_001_029_Composite.tif");
            // IJ.saveAsTiff(imgComposite,mainDir+specimen+"_XR_MRI/result_automate/"+specimen+"_CompositeImageXRMRI.tif");        
            // imgComposite.show();
        }
    */
    // public static ItkTransform runRigidRegistration(String specimen, ImagePlus imgRef, ImagePlus imgMov, String pathToTrInit, String pathToMask){
    //     System.out.println("Now running Rigid Registration for specimen: " + specimen);
    //     System.out.println("-----------------------------------------------");
    //     // Open images
    //     imgRef.show();
    //     imgMov.show();
    //     // Get the global matrix from manual registration and run automatic registration step (rigid body)
    //     ItkTransform trInit=ItkTransform.readTransformFromFile(pathToTrInit);
    //     ItkTransform trFinal = autoLinearRegistration(imgRef,imgMov, trInit,specimen, pathToMask);
    //     ImagePlus finalResult=trFinal.transformImage(imgRef,imgMov);
    //     ImagePlus imgComposite=VitimageUtils.compositeNoAdjustOf(imgRef, finalResult);
    //     imgComposite.close();
    //     imgRef.close();
    //     imgMov.close();
    //     return trFinal;
    // }
    

