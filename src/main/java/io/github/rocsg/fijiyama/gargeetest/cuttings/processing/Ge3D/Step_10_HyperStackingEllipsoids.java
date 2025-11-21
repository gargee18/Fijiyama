package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D;

import ij.CompositeImage;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.common.VitimageUtils;


public class Step_10_HyperStackingEllipsoids implements PipelineStep {
     @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen, true);
    }

    public static void main(String[] args) throws Exception {
        Specimen spec = new Specimen("B_206");
        ImageJ ij = new ImageJ();
        new  Step_10_HyperStackingEllipsoids().execute(spec, true);
    }

    public void execute(Specimen specimen, boolean testing) throws Exception {
           // ImagePlus[] negativeZTab = new ImagePlus[4];
        // ImagePlus[] positiveZTab = new ImagePlus[4];
        // for (int i = 0; i<4; i++){
        //     ImagePlus maskNegativeZ = IJ.openImage(Config.mainDir + "Results/04_EllipsoidFitting/01_MirroredEllipsoids/Negative_Z_"+specimen.getName()+"_"+Config.timestamps[i]+"_mirrored_ellipsoid.tif");
        //     ImagePlus maskPositiveZ = IJ.openImage(Config.mainDir + "Results/04_EllipsoidFitting/01_MirroredEllipsoids/Positive_Z_"+specimen.getName()+"_"+Config.timestamps[i]+"_mirrored_ellipsoid.tif");
        //     negativeZTab[i] = maskNegativeZ;
        //     positiveZTab[i] = maskPositiveZ;
        //     System.out.println(negativeZTab[i]);
        //     System.out.println(positiveZTab[i]);
        //     }
        // ImagePlus hyperStackNegativeZ = VitimageUtils.hyperStackingChannels(negativeZTab);
        // ImagePlus hyperStackPositiveZ = VitimageUtils.hyperStackingChannels(positiveZTab);
        // hyperStackNegativeZ.setTitle("HyperStackNegativeZ");
        // hyperStackPositiveZ.setTitle("HyperStackPositiveZ");


        
        // ImagePlus hyperFrameNZ = VitimageUtils.hyperStackChannelToHyperStackFrame(hyperStackNegativeZ);
        // ImagePlus hyperFramePZ = VitimageUtils.hyperStackChannelToHyperStackFrame(hyperStackPositiveZ);
        // hyperFrameNZ.setTitle("HyperFrameNegativeZ");
        // hyperFramePZ.setTitle("HyperFramePositiveZ");

        // ImagePlus mergedHyperFrameNZPZ = buildTwoChannelHyperstack(specimen,Config.timestamps,hyperFrameNZ, hyperFramePZ, "MergedHyperFrameNZPZ");
        // // mergedHyperFrameNZPZ.show();
        // IJ.saveAsTiff(mergedHyperFrameNZPZ, Config.mainDir + "Results/04_EllipsoidFitting/03_EllipsoidHyperStack/"+specimen.getName()+"Mirrored_Ellipsoid_Hyperstack.tif");

        ImagePlus hyper2C = buildNegPosCompositeHyperstack(specimen);
        IJ.saveAsTiff(hyper2C, Config.mainDir + "Results/04_EllipsoidFitting/03_EllipsoidHyperStack/02_MirroredEllipsoidsandContour/"+specimen.getName()+"_Ellipsoid_Hyperstack.tif");
        
    }


    public static ImagePlus buildTwoChannelHyperstack(Specimen specimen, String[] timestamps, ImagePlus neg, ImagePlus pos, String title) {

        int w = neg.getWidth();
        int h = neg.getHeight();
        int nZ = neg.getNSlices();
        int nT = neg.getNFrames();
        int bit = neg.getBitDepth();

        ImagePlus out = IJ.createHyperStack(title, w, h, 2, nZ, nT, bit);

        for (int t = 1; t <= nT; t++) {
            for (int z = 1; z <= nZ; z++) {

            // ---- Channel 1: NEGATIVE ----
            int idxNZ = neg.getStackIndex(1, z, t);
            out.setPosition(1, z, t);
            out.getProcessor().setPixels(neg.getStack().getProcessor(idxNZ).getPixels());
            out.getStack().setSliceLabel(specimen.getName() + "_Negative_Z_T=" + timestamps[t-1] + "_Z=" + z,out.getStackIndex(1, z, t));

            // ---- Channel 2: POSITIVE ----
            int idxPZ = pos.getStackIndex(1, z, t);
            out.setPosition(2, z, t);
            out.getProcessor().setPixels(pos.getStack().getProcessor(idxPZ).getPixels());
            out.getStack().setSliceLabel(specimen.getName() + "_Positive_Z_T=" + timestamps[t-1] + "_Z=" + z,out.getStackIndex(2, z, t));
            }
        }

        // VERY IMPORTANT: force grayscale mode (no overlay)
        out.setDisplayMode(IJ.GRAYSCALE);

        // optional: turn off composite LUTs
        CompositeImage ci = new CompositeImage(out, CompositeImage.GRAYSCALE);
        ci.setTitle(title);

        return ci;
    }


     public static ImagePlus buildNegPosCompositeHyperstack(Specimen specimen) {

        String[] timestamps = Config.timestamps;
        int nT = timestamps.length;

        // --- Open first negative mirrored stack to get dimensions ---
        String ts0 = timestamps[0];

        String negMirPath0 = Config.mainDir + "Results/04_EllipsoidFitting/01_MirroredEllipsoids/Negative_Z_" + specimen.getName() + "_" + ts0 + "_mirrored_ellipsoid.tif";

        String negConPath0 = Config.mainDir + "Results/04_EllipsoidFitting/03_EllipsoidContours/Negative_Z_" + specimen.getName() + "_" + ts0 + "_ellipsoid_contour.tif";

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
        String title = specimen.getName() + "_NegPos_composite_2C";
        ImagePlus out = IJ.createHyperStack(title, w, h, 2, nZ, nT, bit);

        // --- Fill per timepoint and slice ---
        for (int ti = 0; ti < nT; ti++) {
            String ts = timestamps[ti];

            // Paths for this timepoint
            String negMirPath = Config.mainDir + "Results/04_EllipsoidFitting/01_MirroredEllipsoids/Negative_Z_" + specimen.getName() + "_" + ts + "_mirrored_ellipsoid.tif";

            String negConPath = Config.mainDir + "Results/04_EllipsoidFitting/03_EllipsoidContours/Negative_Z_" + specimen.getName() + "_" + ts + "_ellipsoid_contour.tif";

            String posMirPath = Config.mainDir + "Results/04_EllipsoidFitting/01_MirroredEllipsoids/Positive_Z_" + specimen.getName() + "_" + ts + "_mirrored_ellipsoid.tif";

            String posConPath = Config.mainDir + "Results/04_EllipsoidFitting/03_EllipsoidContours/Positive_Z_" + specimen.getName() + "_" + ts + "_ellipsoid_contour.tif";

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

                out.getStack().setSliceLabel(specimen.getName() + "_Negative_Z_T=" + ts + "_Z=" + z, out.getStackIndex(1, z, t)
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

                out.getStack().setSliceLabel(specimen.getName() + "_Positive_Z_T= " + ts + "_Z=" + z, out.getStackIndex(2, z, t)
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
            


    
}




