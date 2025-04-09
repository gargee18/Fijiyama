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

public class NormalizeMRIT1T2 {
    public static void main(String[] args) {
       
        ImageJ ij = new ImageJ();           
        // Get hyperstack
        ImagePlus imgToNormalize = IJ.openImage("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/T1_T2_Specimen_Analysis/processing/stacksRegistered");
        ImagePlus hypermap = normalizeMRImages(imgToNormalize);

        // ImagePlus hypermap = IJ.openImage("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/raw_registered_normalized/206_hypermap.tif");
        ImagePlus correctedHypermap = sigmaCorrectionPostNormalization(hypermap);
        correctedHypermap.show();
        
    }
// }
    
    public static ImagePlus normalizeMRImages(ImagePlus imgToNormalize) {

        int numChannels = imgToNormalize.getNChannels(); // Number of channels
        int width = imgToNormalize.getWidth();
        int height = imgToNormalize.getHeight();
        int depth = imgToNormalize.getNSlices(); // Number of Z slices

        // double capillarySize = 14; // Capillary size parameter
        // int[] capCentre = {270,129,20};//29=  {265,132,20};//77={224,124,20};// 141={194,139,20}; // Manually set the capillary central coordinates (x, y, z)

        // Create an array to store normalized ImagePlus for each channel
        ImagePlus[] normalizedChannels = new ImagePlus[numChannels];  // This is the array you want

        for (int c = 1; c <= numChannels; c++) {
            // System.out.println("SLice Label: " + imgToNormalize.getStack().getSliceLabel(c));
            ImageStack channelStack = new ImageStack(width, height); // Stack for the current channel
            // Loop through slices for the current channel
            for (int z = 1; z <= depth; z++) {
                imgToNormalize.setPosition(c, z, 1); // Set channel, slice, frame
                ImageProcessor ip = imgToNormalize.getProcessor().duplicate(); 
                String sliceLabel = imgToNormalize.getStack().getSliceLabel(imgToNormalize.getCurrentSlice());
                channelStack.addSlice(sliceLabel,ip); // Add the slice for this channel
            }

            imgToNormalize.getStack().setSliceLabel("", c);
            // Create ImagePlus from the stack of slices for the current channel
            ImagePlus imgChannel = new ImagePlus("Channel " + c, channelStack);

            // Normalize the channel
            // double[] capValues = VitimageUtils.capillaryValuesAlongZStatic(imgChannel, capCentre, capillarySize);
            // double medianCap = VitimageUtils.MADeStatsDoubleSided(capValues, capValues)[0];
            // double[] offset = VitimageUtils.caracterizeBackgroundOfImage(imgChannel);
            // VitimageUtils.writeDoubleArray1DInFile(capValues, "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/hypermap_registered_normalized/capillaryValues.txt");
            // Apply normalization
            ImagePlus imgNorm = processNormalizationWithRespectToCapillary(imgChannel, 0, 57040, 1.0);
            // imgChannel.setTitle("T="+c);
            // Store the normalized ImagePlus in the array at the correct index (c-1)
            normalizedChannels[c - 1] = imgNorm;  // Store in the array
        }
            
        // After the loop, you have the `normalizedChannels[]` array filled with ImagePlus objects for each channel

        // Now pass this array to VitimageUtils.hyperStackingChannels()
        ImagePlus imgNormFinal = VitimageUtils.hyperStackingChannels(normalizedChannels);
        // Now `imgNormFinal` is the final hyperstack that you can display or save
        imgNormFinal.setDisplayRange(0, Config.max_display_val);
        VitimageUtils.setLutToFire(imgNormFinal);

        // Adjust voxel sizes
        Calibration cal = imgToNormalize.getCalibration();
        // System.out.println(cal.pixelWidth);
        cal.setUnit("mm");
        imgNormFinal.setCalibration(cal);
        imgNormFinal.show();
        IJ.saveAsTiff(imgNormFinal,"/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/raw_registered_normalized/206_J141_hyperMap_reg_normalized.tif");
        return imgNormFinal;
    }


    public static ImagePlus processNormalizationWithRespectToCapillary(ImagePlus image, double backgroundVal,double capillaryVal, double targetCapVal) {
        ImageStack stack = image.getStack();
        int stackSize = stack.getSize();
        int width = stack.getWidth();
        int height = stack.getHeight();
        System.out.println(capillaryVal);
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

    public static ImagePlus sigmaCorrectionPostNormalization(ImagePlus image) {
        int C = image.getNChannels();
        int Z = image.getNSlices();
        int T = image.getNFrames();

        double[]normFactors=new double[]{65956,61754,63562,57040};

        for(int c=0; c<C;c++){
            for(int z=0; z<Z;z++){
                for(int t=0; t<T;t++){
                    int sli=VitimageUtils.getCorrespondingSliceInHyperImage(image,c, z, t);
                    String chain=image.getStack().getSliceLabel(sli);
                    System.out.println(chain);
                    String[]strs=chain.split("_");
                    String totSig=strs[strs.length-1];
                    String valueString=totSig.split("=")[1];
                    double d=Double.parseDouble(valueString);
                    System.out.println(d);
                    d/=normFactors[t];
                    String totFutureLabel="";
                    for(int s=0;s<strs.length-1;s++)totFutureLabel=totFutureLabel+=strs[s]+"_";
                    totFutureLabel+="_SIGMARICE="+VitimageUtils.dou(d);
                    image.getStack().setSliceLabel(totFutureLabel, sli);
                }
            }
        }

        IJ.saveAsTiff(image, "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/raw_registered_normalized/206_hypermap_sigma_corrected.tif");
        return image;
    }
    
}
   
   

    

