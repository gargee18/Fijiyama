/*
 * This module performs subsampling of high-resolution MRI images of specimens. It reads raw images, resizes them to a standard size (256x256x512), converts them to 8-bit format for reduced file size and compatibility, and saves the processed images to a designated directory.
 * Key Features:
 *  Input: Reads raw .tif images from the specified directory.
 *  Processing:
 *      Resizes images to a fixed resolution (256x256x512) using bilinear interpolation.
 *      Converts images to 8-bit format.
 *  Output: Saves the processed (subsampled) images in a separate directory.
 * 
 */
package io.github.rocsg.fijiyama.gargeetest.cuttings.processing;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import ij.process.ImageConverter;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;

public class Step_1_Subsample implements PipelineStep {
    private int subsampleRatioStandard;
    
    public Step_1_Subsample(int subRatio){
        subsampleRatioStandard=subRatio;
    }
    public Step_1_Subsample(){
        subsampleRatioStandard=2;
    }
    public static void main(String[] args) throws Exception{
        Specimen spec= new Specimen("B_230");
        new Step_1_Subsample().execute(spec,true); 
    }

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,true);
    }

    public void execute(Specimen specimen,boolean testing) throws Exception {
        String[] days=Config.timestamps;

        int N=days.length;
        for(int i=0;i<N;i++){

            //Open input image
            ImagePlus image=IJ.openImage(Config.getPathToNormalizedImage(specimen, i)); 
            int X=image.getWidth();
            int Y=image.getHeight();
            int Z=image.getNSlices();

            //Resize to 256x256x512 and convert to 8bit
            image = image.resize(X/subsampleRatioStandard, Y/subsampleRatioStandard, Z/subsampleRatioStandard, "bilinear");
            ImageConverter.setDoScaling(true);
            IJ.run(image, "8-bit", "");

            // Adjust voxel sizes
            Calibration cal = image.getCalibration();
            System.out.println(cal.pixelWidth);
            cal.setUnit("mm");
            image.setCalibration(cal);

            if(testing){
                image.setTitle(specimen + ".tif");
                image.show();
            }

            // Save the subsampled image
            IJ.saveAsTiff(image,Config.getPathToSubsampledImage(specimen, i));
            System.out.println("Saved subsampled image to: " + Config.getPathToSubsampledImage(specimen, i));
            if(testing){
                VitimageUtils.waitFor(3000);                
                image.close();
            }
        }
    }

   
}
