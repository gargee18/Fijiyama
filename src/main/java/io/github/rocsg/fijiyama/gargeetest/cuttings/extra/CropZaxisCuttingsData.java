package io.github.rocsg.fijiyama.gargeetest.cuttings.extra;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;


public class CropZaxisCuttingsData {

    public static void main(String[] args) {
            // Define the base directory
            String baseDir = "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/Cuttings_MRI_registration/";
    
            // Loop through the desired range of images
            for (int i = 201; i <= 240; i++) {
                // Generate folder and file names dynamically
                String folderName = "B_" + i + "_J141";
                String rawPath = baseDir + "B_" + i + "/raw/" + folderName + ".tif";
                String savePath = baseDir + "B_" + i + "/raw/" + folderName + "_cropped_z.tif";
    
                System.out.println("Raw: " + rawPath);
                System.out.println("Crop: " + savePath);
    
                // Open the sequence
                ImagePlus imp = IJ.openImage(rawPath);
                if (imp == null) {
                    System.err.println("Failed to open image at: " + rawPath);
                    continue;
                }
    
                // Remove the first 256 slices and the last 256 slices
                ImageStack stack = imp.getStack();
                int totalSlices = stack.getSize();
    
                if (totalSlices <= 512) {
                    System.err.println("Not enough slices to crop: " + totalSlices);
                    imp.close();
                    continue;
                }
    
                // Keep slices in the range [257, totalSlices - 256]
                ImageStack croppedStack = new ImageStack(stack.getWidth(), stack.getHeight());
                for (int slice = 257; slice <= totalSlices - 256; slice++) {
                    croppedStack.addSlice(stack.getSliceLabel(slice), stack.getProcessor(slice));
                }
    
                // Create a new ImagePlus with the cropped stack
                ImagePlus croppedImp = new ImagePlus(imp.getTitle() + "_cropped", croppedStack);
    
                // Save the processed image
                IJ.saveAsTiff(croppedImp, savePath);
    
                // Close images
                imp.close();
                croppedImp.close();
            }
    
            System.out.println("Processing completed!");
        }
    }

