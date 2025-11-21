package io.github.rocsg.fijiyama.gargeetest.cuttings.helpers;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;

public class shit {
    public static void main(String[] args) throws Exception {
        new ImageJ();

        int nT = Config.timestamps.length;  // should be 4: J001, J029, J077, J141

        // --- Open the first negative stack to get dimensions ---
        String firstNegPath = Config.mainDir
                + "Results/04_EllipsoidFitting/01_MirroredEllipsoids/Negative_Z_B_202_"
                + Config.timestamps[0] + "_mirrored_ellipsoid.tif";

        ImagePlus firstNeg = IJ.openImage(firstNegPath);
        if (firstNeg == null) {
            throw new RuntimeException("Could not open " + firstNegPath);
        }

        int w        = firstNeg.getWidth();
        int h        = firstNeg.getHeight();
        int nZ       = firstNeg.getNSlices();
        int bitDepth = firstNeg.getBitDepth();

        System.out.println("Dims: w=" + w + " h=" + h + " z=" + nZ + " bit=" + bitDepth);

        // --- Create output: C=2 (neg,pos), Z=nZ, T=nT ---
        ImagePlus negPosHyper = IJ.createHyperStack(
                "B_202_neg_pos_2C",
                w, h,
                2,      // channels: 1=negative, 2=positive
                nZ,
                nT,
                bitDepth
        );

        // --- Fill the hyperstack ---
        for (int t = 0; t < nT; t++) {
            String ts = Config.timestamps[t];

            String negPath = Config.mainDir
                    + "Results/04_EllipsoidFitting/01_MirroredEllipsoids/Negative_Z_B_202_"
                    + ts + "_mirrored_ellipsoid.tif";

            String posPath = Config.mainDir
                    + "Results/04_EllipsoidFitting/01_MirroredEllipsoids/Positive_Z_B_202_"
                    + ts + "_mirrored_ellipsoid.tif";

            ImagePlus neg = IJ.openImage(negPath);
            ImagePlus pos = IJ.openImage(posPath);

            if (neg == null) throw new RuntimeException("Could not open " + negPath);
            if (pos == null) throw new RuntimeException("Could not open " + posPath);

            // sanity check: dimensions must match the first one
            if (neg.getWidth() != w || neg.getHeight() != h || neg.getNSlices() != nZ)
                throw new RuntimeException("Negative stack at " + ts + " has wrong dims");
            if (pos.getWidth() != w || pos.getHeight() != h || pos.getNSlices() != nZ)
                throw new RuntimeException("Positive stack at " + ts + " has wrong dims");

            int tOut = t + 1; // T is 1-based in ImageJ

            for (int z = 1; z <= nZ; z++) {

                // ---- C=1: negative ----
                negPosHyper.setPosition(1, z, tOut);
                negPosHyper.getProcessor().setPixels(
                        neg.getStack().getProcessor(z).getPixels()
                );

                // ---- C=2: positive ----
                negPosHyper.setPosition(2, z, tOut);
                negPosHyper.getProcessor().setPixels(
                        pos.getStack().getProcessor(z).getPixels()
                );
            }
        }

        negPosHyper.setCalibration(firstNeg.getCalibration());
        negPosHyper.show();
    }
        
}
