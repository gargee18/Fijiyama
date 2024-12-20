package io.github.rocsg.fijiyama.gargeetest.cuttings.processing;

import ij.IJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.registration.ItkTransform;

public class Step_3_1_ObtainTransformationMatrixCombinedWithAligningMatrix {
    
    public static ItkTransform getFinalTrMatrix(String specimen, int timeStamp, String mainDir, String pathToRegistrationMatrix, String pathToAlignmentMatrix , String pathToSaveFinalTransform){
        ItkTransform tr_0=ItkTransform.readTransformFromFile(pathToAlignmentMatrix);
        ItkTransform tr_t_to_0=ItkTransform.readTransformFromFile(pathToRegistrationMatrix);
        ItkTransform tr_t_to_0_aligned=(new ItkTransform(tr_t_to_0)).addTransform(tr_0);
        tr_t_to_0_aligned.writeMatrixTransformToFile(pathToSaveFinalTransform);
        return tr_t_to_0_aligned;
    }

    public static ImagePlus imgVerfiyAlignment(String specimen, int timeStamp, String pathToRef, String pathToMov, ItkTransform trMat){
        ImagePlus img0 = IJ.openImage(pathToRef);
        ImagePlus imgt = IJ.openImage(pathToMov);//mainDir+specimen+"/outputs/"+specimen+"_001_141_TransformedImage_Reg_Aligned.tif");
        ImagePlus imgtTo0= trMat.transformImage(img0, imgt);
        ImagePlus imgtTo0_composite = VitimageUtils.compositeNoAdjustOf(img0, imgtTo0);
        imgtTo0_composite.setTitle(specimen + " J001 "+ timeStamp);
        imgtTo0_composite.show();imgtTo0_composite.setSlice(512);
        VitimageUtils.waitFor(3000);
        imgtTo0.close();
        return imgtTo0;
    }


       // public static void main(String[] args) {
    //     ImageJ ij = new ImageJ();
    //     for (int i = 201; i <= 201; i++) {
    //         if (i == 204) {
    //             continue;
    //         }
    //         String specimen = "B_" + i;
    //         String pathToRegistrationMatrix =  Config.getPathToAutomaticRegistrationTransformation(specimen,3);
    //         String pathToAlignmentMatrix = Config.getPathtoSaveAlignedTransformationMatrix(specimen);
    //         String pathToSaveFinalTransform = Config.getPathToSaveFinalTransform(specimen, 3);
    //         ItkTransform tr_3_0 = getFinalTrMatrix(specimen,3, Config.mainDir, pathToRegistrationMatrix, pathToAlignmentMatrix, pathToSaveFinalTransform);
                      
    //         String pathToRef = Config.getPathToImageRef(specimen, 0);
    //         String pathToMov = Config.getPathToImageMov(specimen, 3);
    //         ImagePlus img_3_0 = imgVerfiyAlignment(specimen, 3, pathToRef,pathToMov, tr_3_0);
    //         img_3_0.show();

        
    //     }
    //     System.out.println("done");
    // }

}
