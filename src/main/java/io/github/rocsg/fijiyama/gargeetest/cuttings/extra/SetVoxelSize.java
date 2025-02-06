package io.github.rocsg.fijiyama.gargeetest.cuttings.extra;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;

public class SetVoxelSize {
    public static void main(String[] args) {
        int start = 201;
        int end = 240; // Only processing one image

        String baseDir = "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/Cuttings/Data/02_Normalized/";

        for (int i = start; i <= end; i++) {
            // Construct file name
            if( i == 204) continue;
            String fileName = "B_" + i + "_J141_normalized.tif";
            String filePath = baseDir + fileName;

            // Open image
            ImagePlus imp = IJ.openImage(filePath);
            if (imp == null) {
                System.err.println("Error: Could not open file " + filePath);
                continue;
            }

            // Set calibration
            Calibration cal = imp.getCalibration();
            cal.setUnit("mm");
            // cal.pixelWidth = 0.0703126;
            // cal.pixelHeight = 0.0703126;
            // cal.pixelDepth = 0.0703126;
            cal.pixelWidth = 0.0351563;
            cal.pixelHeight = 0.0351563;
            cal.pixelDepth = 0.0351563;
            
            // Save the image
            IJ.save(imp, filePath);

            // Close the image
            imp.close();
        }
    }
    
}
