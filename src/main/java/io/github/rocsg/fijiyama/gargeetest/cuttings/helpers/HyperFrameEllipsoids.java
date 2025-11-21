package io.github.rocsg.fijiyama.gargeetest.cuttings.helpers;
import ij.CompositeImage;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import ij.plugin.HyperStackConverter;
import ij.plugin.ImageCalculator;
import ij.plugin.RGBStackMerge;
import ij.process.ImageProcessor;
import ij.process.LUT;
import java.awt.Color;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;


public class HyperFrameEllipsoids {

    public static void main(String[] args) throws Exception {
        new ImageJ();

        String specimenId = "B_206";  // or "B_202", etc.

        ImagePlus hyper2C = buildNegPosCompositeHyperstack(specimenId);
        // hyper2C.show();
        IJ.saveAsTiff(hyper2C, Config.mainDir + "Results/04_EllipsoidFitting/03_EllipsoidHyperStack/02_MirroredEllipsoidsandContour/"+specimenId+"_Ellipsoid_Hyperstack.tif");
    }
        
    public static ImagePlus buildNegPosCompositeHyperstack(String specimenId) {

        String[] timestamps = Config.timestamps;
        int nT = timestamps.length;

        // --- Open first negative mirrored stack to get dimensions ---
        String ts0 = timestamps[0];

        String negMirPath0 = Config.mainDir +
                "Results/04_EllipsoidFitting/01_MirroredEllipsoids/Negative_Z_" +
                specimenId + "_" + ts0 + "_mirrored_ellipsoid.tif";

        String negConPath0 = Config.mainDir +
                "Results/04_EllipsoidFitting/03_EllipsoidContours/Negative_Z_" +
                specimenId + "_" + ts0 + "_ellipsoid_contour.tif";

        ImagePlus negMir0 = IJ.openImage(negMirPath0);
        ImagePlus negCon0 = IJ.openImage(negConPath0);

        if (negMir0 == null) throw new RuntimeException("Cannot open " + negMirPath0);
        if (negCon0 == null) throw new RuntimeException("Cannot open " + negConPath0);

        int w   = negMir0.getWidth();
        int h   = negMir0.getHeight();
        int nZ  = negMir0.getNSlices();
        int bit = negMir0.getBitDepth();

        // Sanity check: first contour must match first mirrored
        if (negCon0.getWidth()  != w ||
            negCon0.getHeight() != h ||
            negCon0.getNSlices() != nZ) {
            throw new RuntimeException("Negative contour " + ts0 + " does not match mirrored in size");
        }

        // --- Create final hyperstack: C=2 (neg,pos), Z=nZ, T=nT ---
        String title = specimenId + "_NegPos_composite_2C";
        ImagePlus out = IJ.createHyperStack(title, w, h, 2, nZ, nT, bit);

        // --- Fill per timepoint and slice ---
        for (int ti = 0; ti < nT; ti++) {
            String ts = timestamps[ti];

            // Paths for this timepoint
            String negMirPath = Config.mainDir +
                    "Results/04_EllipsoidFitting/01_MirroredEllipsoids/Negative_Z_" +
                    specimenId + "_" + ts + "_mirrored_ellipsoid.tif";

            String negConPath = Config.mainDir +
                    "Results/04_EllipsoidFitting/03_EllipsoidContours/Negative_Z_" +
                    specimenId + "_" + ts + "_ellipsoid_contour.tif";

            String posMirPath = Config.mainDir +
                    "Results/04_EllipsoidFitting/01_MirroredEllipsoids/Positive_Z_" +
                    specimenId + "_" + ts + "_mirrored_ellipsoid.tif";

            String posConPath = Config.mainDir +
                    "Results/04_EllipsoidFitting/03_EllipsoidContours/Positive_Z_" +
                    specimenId + "_" + ts + "_ellipsoid_contour.tif";

            ImagePlus negMir = IJ.openImage(negMirPath);
            ImagePlus negCon = IJ.openImage(negConPath);
            ImagePlus posMir = IJ.openImage(posMirPath);
            ImagePlus posCon = IJ.openImage(posConPath);

            if (negMir == null) throw new RuntimeException("Cannot open " + negMirPath);
            if (negCon == null) throw new RuntimeException("Cannot open " + negConPath);
            if (posMir == null) throw new RuntimeException("Cannot open " + posMirPath);
            if (posCon == null) throw new RuntimeException("Cannot open " + posConPath);

            // Check dimensions match reference
            if (negMir.getWidth() != w || negMir.getHeight() != h || negMir.getNSlices() != nZ)
                throw new RuntimeException("Negative mirrored " + ts + " has wrong size");
            if (negCon.getWidth() != w || negCon.getHeight() != h || negCon.getNSlices() != nZ)
                throw new RuntimeException("Negative contour " + ts + " has wrong size");
            if (posMir.getWidth() != w || posMir.getHeight() != h || posMir.getNSlices() != nZ)
                throw new RuntimeException("Positive mirrored " + ts + " has wrong size");
            if (posCon.getWidth() != w || posCon.getHeight() != h || posCon.getNSlices() != nZ)
                throw new RuntimeException("Positive contour " + ts + " has wrong size");

            int t = ti + 1; // ImageJ T index (1-based)

            for (int z = 1; z <= nZ; z++) {

                // ---- NEGATIVE composite: max(mirrored, contour) → channel 1 ----
                ImageProcessor ipNegMir = negMir.getStack().getProcessor(z);
                ImageProcessor ipNegCon = negCon.getStack().getProcessor(z);

                out.setPosition(1, z, t);
                ImageProcessor ipOutNeg = out.getProcessor();

                for (int y = 0; y < h; y++) {
                    for (int x = 0; x < w; x++) {
                        int vM = ipNegMir.get(x, y);
                        int vC = ipNegCon.get(x, y);
                        int v  = (vM > vC) ? vM : vC;
                        ipOutNeg.set(x, y, v);
                    }
                }

                out.getStack().setSliceLabel(
                        specimenId + " Negative_Z " + ts + " Z=" + z,
                        out.getStackIndex(1, z, t)
                );

                // ---- POSITIVE composite: max(mirrored, contour) → channel 2 ----
                ImageProcessor ipPosMir = posMir.getStack().getProcessor(z);
                ImageProcessor ipPosCon = posCon.getStack().getProcessor(z);

                out.setPosition(2, z, t);
                ImageProcessor ipOutPos = out.getProcessor();

                for (int y = 0; y < h; y++) {
                    for (int x = 0; x < w; x++) {
                        int vM = ipPosMir.get(x, y);
                        int vC = ipPosCon.get(x, y);
                        int v  = (vM > vC) ? vM : vC;
                        ipOutPos.set(x, y, v);
                    }
                }

                out.getStack().setSliceLabel(
                        specimenId + " Positive_Z " + ts + " Z=" + z,
                        out.getStackIndex(2, z, t)
                );
            }
        }

        out.setCalibration(negMir0.getCalibration());

        // Wrap into CompositeImage just to manage LUTs
        CompositeImage ci = new CompositeImage(out, CompositeImage.COMPOSITE);
          
        ci.setDisplayMode(IJ.GRAYSCALE); // no color overlay by default

        ci.setTitle(title);

        return ci;
    }
            
    
    public static ImagePlus buildTwoChannelHyperstack(ImagePlus neg, ImagePlus pos, String title) {

        int w = neg.getWidth();
        int h = neg.getHeight();
        int nZ = neg.getNSlices();
        int nT = neg.getNFrames();
        int bit = neg.getBitDepth();

        ImagePlus out = IJ.createHyperStack(title, w, h, 2, nZ, nT, bit);

        for (int t = 1; t <= nT; t++) {
            for (int z = 1; z <= nZ; z++) {

                
            int idxNZ = neg.getStackIndex(1, z, t);
            out.setPosition(1, z, t);
            out.getProcessor().setPixels(neg.getStack().getProcessor(idxNZ).getPixels());
            out.getStack().setSliceLabel("Negative_Z  T=" + t + "  Z=" + z,
                                         out.getStackIndex(1, z, t));

            // ---- Channel 2: POSITIVE ----
            int idxPZ = pos.getStackIndex(1, z, t);
            out.setPosition(2, z, t);
            out.getProcessor().setPixels(pos.getStack().getProcessor(idxPZ).getPixels());
            out.getStack().setSliceLabel("Positive_Z  T=" + t + "  Z=" + z,
                                         out.getStackIndex(2, z, t));
            }
        }

        // VERY IMPORTANT: force grayscale mode (no overlay)
        out.setDisplayMode(IJ.GRAYSCALE);

        // optional: turn off composite LUTs
        CompositeImage ci = new CompositeImage(out, CompositeImage.GRAYSCALE);
        ci.setTitle(title);

        return ci;
    }

    public static ImagePlus compositeMirroredAndContour(ImagePlus mirrored, ImagePlus contour, String title) {
    if (mirrored == null || contour == null) {
        throw new IllegalArgumentException("One of the inputs is null");
    }

    if (mirrored.getWidth()  != contour.getWidth() ||
        mirrored.getHeight() != contour.getHeight() ||
        mirrored.getNSlices() != contour.getNSlices() ||
        mirrored.getNFrames() != contour.getNFrames()) {
        throw new IllegalArgumentException("Mirrored and contour do not match in size");
    }

    int w   = mirrored.getWidth();
    int h   = mirrored.getHeight();
    int nZ  = mirrored.getNSlices();
    int nT  = mirrored.getNFrames();
    int bit = mirrored.getBitDepth();

    ImagePlus out = IJ.createHyperStack(title, w, h, 1, nZ, nT, bit);

    for (int t = 1; t <= nT; t++) {
        for (int z = 1; z <= nZ; z++) {

            int idxM = mirrored.getStackIndex(1, z, t);
            int idxC = contour.getStackIndex(1, z, t);

            ImageProcessor ipM = mirrored.getStack().getProcessor(idxM);
            ImageProcessor ipC = contour.getStack().getProcessor(idxC);

            out.setPosition(1, z, t);
            ImageProcessor ipOut = out.getProcessor();

            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    int v = Math.max(ipM.get(x,y), ipC.get(x,y));
                    ipOut.set(x, y, v);
                }
            }
        }
    }

    out.setCalibration(mirrored.getCalibration());
    return out;
}
    public static ImagePlus buildTwoChannelHyperstack(
            ImagePlus neg,
            ImagePlus pos,
            String title,
            String specimenId,
            String[] timestamps
    ) {
        int w = neg.getWidth();
        int h = neg.getHeight();
        int nZ = neg.getNSlices();
        int nT = neg.getNFrames();
        int bit = neg.getBitDepth();

        ImagePlus out = IJ.createHyperStack(title, w, h, 2, nZ, nT, bit);

        for (int t = 1; t <= nT; t++) {
            String ts = timestamps[t - 1];

            for (int z = 1; z <= nZ; z++) {

                // NEGATIVE channel
                int idxNZ = neg.getStackIndex(1, z, t);
                out.setPosition(1, z, t);
                out.getProcessor().setPixels(
                    neg.getStack().getProcessor(idxNZ).getPixels()
                );
                out.getStack().setSliceLabel(
                    specimenId + " Negative_Z " + ts + " Z=" + z,
                    out.getStackIndex(1, z, t)
                );

                // POSITIVE channel
                int idxPZ = pos.getStackIndex(1, z, t);
                out.setPosition(2, z, t);
                out.getProcessor().setPixels(
                    pos.getStack().getProcessor(idxPZ).getPixels()
                );
                out.getStack().setSliceLabel(
                    specimenId + " Positive_Z " + ts + " Z=" + z,
                    out.getStackIndex(2, z, t)
                );
            }
        }

        // Make CompositeImage so LUTs can be applied
        CompositeImage ci = new CompositeImage(out, CompositeImage.COMPOSITE);
        ci.setTitle(title);

        return ci;
    }

   
    
   

    
}
