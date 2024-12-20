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
import ij.process.ImageConverter;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;

public class Step_2_Subsample implements PipelineStep {
    private int subsampleRatioStandard;
    
    public Step_2_Subsample(int subRatio){
        subsampleRatioStandard=subRatio;
    }
    public Step_2_Subsample(){
        subsampleRatioStandard=2;
    }

    public void execute(Specimen specimen,boolean testing) throws Exception {
        String[]days=Config.timestamps;
        for(String day : days){
            // Open the images
            String pathToInputImage=specimen.getSpecimenDirectory()+"raw/"+specimen.getName()+"_"+day+".tif";
            String pathToOutputImage=specimen.getSpecimenDirectory()+"raw_subsampled/"+specimen.getName()+"_"+day+".tif";

            //Open input image
            ImagePlus image=IJ.openImage(pathToInputImage);
            int X=image.getWidth();
            int Y=image.getHeight();
            int Z=image.getNSlices();
            //Resize to 256x256x512 and convert to 8bit
            image = image.resize(X/subsampleRatioStandard, Y/subsampleRatioStandard, Z/subsampleRatioStandard, "bilinear");
            ImageConverter.setDoScaling(true);
            IJ.run(image, "8-bit", "");

            if(testing){
                image.setTitle(specimen + ".tif");
                image.show();
            }

            // Save the subsampled image
            IJ.saveAsTiff(image,pathToOutputImage);
            System.out.println("Saved subsampled image to: " + pathToOutputImage);
            if(testing){
                VitimageUtils.waitFor(3000);                
                image.close();
            }
        }
    }

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }
}
