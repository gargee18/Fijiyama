package io.github.rocsg.fijiyama.gargeetest.ceps;

//specific libraries
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.registration.ItkTransform;


public class RegisterXRandMRI {
    
    // public static String pathToTrInit = "/home/phukon/Desktop/MRI/02_ceps/2022_CEPS_suivi_clinique/transMatrices_2022_XR_MRI/";
    public static String mainDir = "/home/phukon/Desktop/";
    public static String trMatPath = "/home/phukon/Desktop/master_repo/xray_2024_mrit1_2024/";
    public static String outPath = "/home/phukon/Desktop/";
    public static String year1 = GeneralUtils.years[2];
    //public static String year2 = GeneralUtils.years[2];

  
    public static void main(String[] args) {
 
        ImageJ ij = new ImageJ();
        
        String[] specimen = GeneralUtils.specimenListMRI; 
        // get the tr matrix
        for (int i = 2; i < 3/*specimen.length*/; i ++){
            String imgRefpath = mainDir+"master_repo/xray_2024_mrit1_2024/cep_323/raw/"+specimen[i]+"_imgRef.tif"; 
            System.out.println(imgRefpath);
            String imgMovpath =  mainDir+"master_repo/xray_2024_mrit1_2024/cep_323/raw/"+specimen[i]+"_imgMov.tif"; 
            System.out.println(imgMovpath);
            ImagePlus imgRef =IJ.openImage(imgRefpath);
            imgRef.show();
            ImagePlus imgMov = IJ.openImage(imgMovpath);
            imgMov.show();
            String tr = trMatPath+"cep_"+specimen[i]+"/res/Exported_data/"+specimen[i]+"_trMatrix_2023_2024.txt"; 
            System.out.println(tr);

            ItkTransform trMat=ItkTransform.readTransformFromFile(tr);
            ImagePlus finalResult=trMat.transformImage(imgRef,imgMov);
            IJ.saveAsTiff(finalResult,outPath+ specimen[i]+"_TransformedImage_XR_T1_"+year1+".tif");
            trMat.writeMatrixTransformToFile("/home/phukon/Desktop/"+specimen[i]+"_img_moving_after_registration_highres_T1.txt");

     
        }
    }

}











            //ItkTransform trFinal = RegistrationAutomate.autoLinearRegistration(imgRef,imgMov,trInit,specimen);
            // ImagePlus finalResult=trFinal.transformImage(imgRef,imgMov);
            // IJ.saveAsTiff(finalResult,mainDir+"_XR_MRI/result_automate/"+specimen+"_TransformedImage.tif");
            // ImagePlus imgComposite=VitimageUtils.compositeNoAdjustOf(imgRef, finalResult);
            // IJ.saveAsTiff(imgComposite,mainDir+"_XR_MRI/result_automate/"+specimen+"_TransformedImage.tif");
            // imgComposite.show();
            // imgComposite.close();
            // imgRef.close();
            // imgMov.close();