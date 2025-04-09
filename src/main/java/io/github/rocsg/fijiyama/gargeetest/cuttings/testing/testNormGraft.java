package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;
import java.util.Arrays;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;

public class testNormGraft {
    public static void main(String[] args) {
       
        ImageJ ij = new ImageJ();           
        // Get hyperstack
        ImagePlus imgToNormalize = IJ.openImage("/home/phukon/Desktop/ge3d_MMC_1_20190410_01/test.tif");
        int numChannels = imgToNormalize.getNChannels(); // Number of channels
        int width = imgToNormalize.getWidth();
        int height = imgToNormalize.getHeight();
        int depth = imgToNormalize.getNSlices(); // Number of Z slices
        int frames = imgToNormalize.getNFrames(); // Number of time points

        double capillarySize = 14; // Capillary size parameter
        int[] capCentre = {146,197,512};//29={265,132,20};//77={224,124,20};// 141={194,139,20}; // Manually set the capillary central coordinates (x, y, z)
      
       // Create an array to store normalized ImagePlus for each channel
        ImagePlus[] normalizedChannels = new ImagePlus[numChannels];  // This is the array you want

        for (int c = 1; c <= numChannels; c++) {
            ImageStack channelStack = new ImageStack(width, height); // Stack for the current channel

            // Loop through slices for the current channel
            for (int z = 1; z <= depth; z++) {
                imgToNormalize.setPosition(c, z, 1); // Set channel, slice, frame
                ImageProcessor ip = imgToNormalize.getProcessor().duplicate();
                channelStack.addSlice(ip); // Add the slice for this channel
            }

            // Create ImagePlus from the stack of slices for the current channel
            ImagePlus imgChannel = new ImagePlus("Channel " + c, channelStack);

            // Normalize the channel
            double[] capValues = VitimageUtils.capillaryValuesAlongZStatic(imgChannel, capCentre, capillarySize);
            double medianCap = VitimageUtils.MADeStatsDoubleSided(capValues, capValues)[0];
            double[] offset = VitimageUtils.caracterizeBackgroundOfImage(imgChannel);

            // Apply normalization
            ImagePlus imgNorm = processNormalizationWithRespectToCapillary(imgChannel, offset[0], medianCap, 1.0);

            // Store the normalized ImagePlus in the array at the correct index (c-1)
            normalizedChannels[c - 1] = imgNorm;  // Store in the array
        }

        // After the loop, you have the `normalizedChannels[]` array filled with ImagePlus objects for each channel

        // Now pass this array to VitimageUtils.hyperStackingChannels()
        ImagePlus imgNormFinal = VitimageUtils.hyperStackingChannels(normalizedChannels);

        // Now `imgNormFinal` is the final hyperstack that you can display or save
        imgNormFinal.setDisplayRange(0, Config.max_display_val);
        VitimageUtils.setLutToFire(imgNormFinal);

        imgNormFinal.show();
        // IJ.saveAsTiff(imgNormFinal,"/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/hypermap_registered_normalized/206_J001_hyperMap_reg_normalized.tif");
        
    }
    

    public static ImagePlus processNormalizationWithRespectToCapillary(ImagePlus image, double backgroundVal,double capillaryVal, double targetCapVal) {
        ImageStack stack = image.getStack();
        int stackSize = stack.getSize();
        int width = stack.getWidth();
        int height = stack.getHeight();

        ImageStack normalizedStack = new ImageStack(width, height);

        for (int z = 1; z <= stackSize; z++) {
            FloatProcessor sliceProcessor = (FloatProcessor) stack.getProcessor(z);
            FloatProcessor normalizedSlice = (FloatProcessor) sliceProcessor.duplicate();

            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    double intensity = sliceProcessor.getPixelValue(x, y);
                    double normalizedIntensity = ((intensity - 0*backgroundVal) / (capillaryVal - 0*backgroundVal))* targetCapVal;
                    normalizedSlice.putPixelValue(x, y, (float) normalizedIntensity);
                }
            }
            normalizedStack.addSlice(stack.getSliceLabel(z), normalizedSlice);
        }
        ImagePlus normalizedImage = new ImagePlus(image.getTitle() + "_normalized", normalizedStack);

        return normalizedImage;


    } 
    
}
   
   

    



