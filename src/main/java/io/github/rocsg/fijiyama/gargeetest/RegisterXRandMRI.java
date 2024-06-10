package io.github.rocsg.fijiyama.gargeetest;

//specific libraries
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.registration.ItkTransform;


public class RegisterXRandMRI {
    
    // public static String pathToTrInit = "/home/phukon/Desktop/MRI/02_ceps/2022_CEPS_suivi_clinique/transMatrices_2022_XR_MRI/";
    public static String mainDir = "/home/phukon/Desktop/registration/CEP_2022_XR_2023_XR_2024_XR/";
    public static String trMatPath = "/home/phukon/Desktop/CEP_Registration/CEP_2022_XR_to_2023_XR_HIGHRES/";
    public static String year1 = GeneralUtils.years[0];
    public static String year2 = GeneralUtils.years[1];


    public static void main(String[] args) {
 
        ImageJ ij = new ImageJ();
        
        String[] specimenList = GeneralUtils.specimenListXR; 
        // get the tr matrix
        for(String specimen: specimenList){
            System.out.println(mainDir);
            String imgRefpath = mainDir+"cep_"+specimen+"/raw/"+specimen+"_"+year1+"_crop_sub_z.tif"; 
            System.out.println(imgRefpath);
            String imgMovpath = mainDir+"cep_"+specimen+"/raw/"+specimen+"_"+year2+"_crop_sub_z.tif"; 

            ImagePlus imgRef =IJ.openImage(imgRefpath);
            ImagePlus imgMov = IJ.openImage(imgMovpath);


            ItkTransform trInit=ItkTransform.readTransformFromFile(trMatPath+specimen+"tRigHigh.txt");

            //ItkTransform trFinal = RegistrationAutomate.autoLinearRegistration(imgRef,imgMov,trInit,specimen);
            // ImagePlus finalResult=trFinal.transformImage(imgRef,imgMov);
            // IJ.saveAsTiff(finalResult,mainDir+"_XR_MRI/result_automate/"+specimen+"_TransformedImage.tif");
            // ImagePlus imgComposite=VitimageUtils.compositeNoAdjustOf(imgRef, finalResult);
            // IJ.saveAsTiff(imgComposite,mainDir+"_XR_MRI/result_automate/"+specimen+"_TransformedImage.tif");
            // imgComposite.show();
            // imgComposite.close();
            // imgRef.close();
            // imgMov.close();



            ItkTransform trMat=ItkTransform.readTransformFromFile(trMatPath+specimen+"_tRigHigh.txt");
            
            ImagePlus finalResult=trMat.transformImage(imgRef,imgMov);
            ImagePlus imgComposite=VitimageUtils.compositeNoAdjustOf(imgRef, imgMov);
            IJ.saveAsTiff(finalResult,trMatPath+ specimen+"_TransformedImage.tif");
            IJ.saveAsTiff(imgComposite,trMatPath+ specimen+"_CompositeImage.tif");

        }
    }

}
