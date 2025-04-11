package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.T1T2;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;

public class Step_6_NormalizeMRIT1T2 implements PipelineStep{


    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }
    public static void main(String[] args) throws Exception {
        ImageJ ij = new ImageJ();           
        Specimen spec= new Specimen("B_223");
        new Step_6_NormalizeMRIT1T2 ();
        Step_6_NormalizeMRIT1T2.execute(spec,true); 
    }
   
   
    /**
     * Executes the normalization process for a given specimen, applying channel-wise normalization
     * and combining the results into a hyperstack.
     *
     * @param specimen The specimen to process.
     * @param testing  Flag to indicate if testing mode is enabled.
     * @return The final normalized and processed ImagePlus object.
     * @throws Exception If an error occurs during processing.
     */
    public static ImagePlus execute(Specimen specimen, boolean testing) throws Exception {
        // Open the image to be normalized
        ImagePlus imgToNormalize = IJ.openImage("/home/phukon/Desktop/B_223_J029_registered.tif");

        // Get image dimensions and channel count
        int numChannels = imgToNormalize.getNChannels();
        int width = imgToNormalize.getWidth();
        int height = imgToNormalize.getHeight();
        int depth = imgToNormalize.getNSlices();

        // Array to hold normalized ImagePlus objects for each channel
        ImagePlus[] normalizedChannels = new ImagePlus[numChannels];

        // Process each channel
        for (int c = 1; c <= numChannels; c++) {
            // Create a new ImageStack for the current channel
            ImageStack channelStack = new ImageStack(width, height);

            // Loop through slices for the current channel
            for (int z = 1; z <= depth; z++) {
                imgToNormalize.setPosition(c, z, 1); // Set channel, slice, frame
                ImageProcessor ip = imgToNormalize.getProcessor().duplicate();
                String sliceLabel = imgToNormalize.getStack().getSliceLabel(imgToNormalize.getCurrentSlice());
                channelStack.addSlice(sliceLabel, ip); // Add the slice for this channel
            }

            // Clear the slice label (not used further)
            imgToNormalize.getStack().setSliceLabel("", c);

            // Create ImagePlus from the stack of slices for the current channel
            ImagePlus imgChannel = new ImagePlus("Channel " + c, channelStack);

            // Apply normalization to the channel image
            ImagePlus imgNorm = processNormalizationWithRespectToCapillary(imgChannel, 0, 57040, 1.0);

            // Store the normalized ImagePlus in the array at the correct index (c-1)
            normalizedChannels[c - 1] = imgNorm;
        }

        // Combine all normalized channels into a single hyperstack
        ImagePlus imgNormFinal = VitimageUtils.hyperStackingChannels(normalizedChannels);

        // Set display range and apply LUT
        imgNormFinal.setDisplayRange(0, Config.max_display_val);
        VitimageUtils.setLutToFire(imgNormFinal);

        // Adjust voxel sizes and calibration
        Calibration cal = imgToNormalize.getCalibration();
        cal.setUnit("mm");
        imgNormFinal.setCalibration(cal);

        // Perform sigma correction on the normalized image
        imgNormFinal = sigmaCorrectionPostNormalization(imgNormFinal);

        // Display the final hyperstack
        imgNormFinal.show();

        // Return the final processed ImagePlus object
        return imgNormFinal;
    }


    /**
     * This method normalizes an image with respect to a given capillary value
     * 
     * @param image The image to be normalized
     * @param backgroundVal The minimum intensity value in the image
     * @param capillaryVal The intensity value corresponding to the capillary
     * @param targetCapVal The desired intensity value for the capillary in the normalized image
     * @return The normalized image
     */
    public static ImagePlus processNormalizationWithRespectToCapillary(ImagePlus image, double backgroundVal,double capillaryVal, double targetCapVal) {
        // Get the stack of the image
        ImageStack stack = image.getStack();
        int stackSize = stack.getSize();
        int width = stack.getWidth();
        int height = stack.getHeight();
        System.out.println("Capillary value: "+capillaryVal);

        // Create an empty stack for the normalized image
        ImageStack normalizedStack = new ImageStack(width, height);

        // Loop through each slice in the stack
        for (int z = 1; z <= stackSize; z++) {
            // Get the current slice
            FloatProcessor sliceProcessor = (FloatProcessor) stack.getProcessor(z);
            // Create a duplicate of the slice to store the normalized values
            FloatProcessor normalizedSlice = (FloatProcessor) sliceProcessor.duplicate();

            // Loop through each pixel in the slice
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    // Get the intensity value at the current pixel position
                    double intensity = sliceProcessor.getPixelValue(x, y);
                    // Calculate the normalized value
                    double normalizedIntensity = ((intensity - 0*backgroundVal) / (capillaryVal - 0*backgroundVal))* targetCapVal;
                    // Store the normalized value in the normalized slice
                    normalizedSlice.putPixelValue(x, y, (float) normalizedIntensity);
                }
            }
            // Add the normalized slice to the normalized stack
            normalizedStack.addSlice(stack.getSliceLabel(z), normalizedSlice);
        }
        // Create a new image from the normalized stack
        ImagePlus normalizedImage = new ImagePlus(image.getTitle() + "_normalized", normalizedStack);

        // Return the normalized image
        return normalizedImage;

    }

    /**
     * This function performs sigma correction on a normalized image.
     * It adjusts the sigma values of the image slices based on predefined normalization factors.
     *
     * @param image The ImagePlus object representing the image to be corrected.
     * @return The sigma-corrected ImagePlus object.
     */
    public static ImagePlus sigmaCorrectionPostNormalization(ImagePlus image) {
        // Get the number of channels, slices, and frames
        int C = image.getNChannels();
        int Z = image.getNSlices();
        int T = image.getNFrames();

        // Predefined normalization factors for each frame
        double[] normFactors = new double[]{65956, 61754, 63562, 57040};

        // Iterate through each channel, slice, and frame
        for (int c = 0; c < C; c++) {
            for (int z = 0; z < Z; z++) {
                for (int t = 0; t < T; t++) {
                    // Get the corresponding slice index in the hyper image
                    int sli = VitimageUtils.getCorrespondingSliceInHyperImage(image, c, z, t);
                    // Retrieve the slice label
                    String chain = image.getStack().getSliceLabel(sli);
                    System.out.println(chain);

                    // Parse the sigma value from the slice label
                    String[] strs = chain.split("_");
                    String totSig = strs[strs.length - 1];
                    String valueString = totSig.split("=")[1];
                    double d = Double.parseDouble(valueString);
                    System.out.println(d);

                    // Normalize the sigma value using the corresponding normalization factor
                    d /= normFactors[t];

                    // Construct the new slice label with the corrected sigma value
                    String totFutureLabel = "";
                    for (int s = 0; s < strs.length - 1; s++) totFutureLabel += strs[s] + "_";
                    totFutureLabel += "_SIGMARICE=" + VitimageUtils.dou(d);

                    // Update the slice label in the image stack
                    image.getStack().setSliceLabel(totFutureLabel, sli);
                }
            }
        }

        // Save the corrected image as a TIFF file
        IJ.saveAsTiff(image, "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/raw_registered_normalized/206_hypermap_sigma_corrected.tif");
        
        // Return the corrected image
        return image;
    }

    
    
}
   
   

    

