package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D;
import java.io.File;
import java.util.ArrayList;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.Duplicator;
import ij.plugin.LutLoader;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.gargeetest.cuttings.helpers.ImgProUtils;
import trainableSegmentation.WekaSegmentation;

public class Step_7_ProbabilisticAtlas implements PipelineStep {
    static final int cond_PCH = 0;
    static final int cond_CONTROL = 1;

    static final int var_CHARD = 0;
    static final int var_MERLOT = 1;
    static final int var_TEMPRA = 2;
    static final int var_UGNI = 3;

    static final String[] VAR_NAMES = {
        "var_CHARD", "var_MERLOT", "var_TEMPRA", "var_UGNI"
    };

    public static String getVarietyName(int value) {
        if (value >= 0 && value < VAR_NAMES.length) {
            return VAR_NAMES[value];
        }
        return "Unknown";
    }

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen, true);
    }

    public static void main (String[] args) throws Exception {
        Specimen spec = new Specimen("B_201");
        
        ImageJ ij = new ImageJ();
        new Step_7_ProbabilisticAtlas().execute(spec,true);
    }

    public void execute(Specimen specimen, boolean testing) throws Exception {
           
        // for(int t = 1; t <5; t++){
        //      ImagePlus m = runFillTopToFullRowOrBest(specimen, t); 
        // }   
       
    String[] timestamps = Config.timestamps;
        
    for(int t = 1; t <5; t++){
        ImagePlus atlas = buildProbabilisticAtlasAllVar(cond_CONTROL, t);
        atlas.show();
        System.out.println(Config.mainDir+"Results/05_ProbabilisticAtlas/ProbabilisticAtlas_allVar_CT_"+timestamps[t-1]+".tif");

        IJ.saveAsTiff(atlas, Config.mainDir+"Results/05_ProbabilisticAtlas/ProbabilisticAtlas_allVar_CT_"+timestamps[t-1]+".tif");
    }
        
    }


   
    public static ImagePlus buildProbabilisticAtlas(int condition, int variety, int t) {
        String[] timestamps = Config.timestamps;
        String[] specimens = Config.getSpecimensName(condition, variety);
        ImagePlus mask1 = IJ.openImage(Config.mainDir + "Results/04_MaskPithFilled/" + specimens[0] +"_"+ timestamps[t-1]+"_mask.tif");
        mask1.show();
        int h = mask1.getHeight();
        int w = mask1.getWidth();
        int d = mask1.getNSlices();

        // accumulate masks
        ImageStack sumStack = new ImageStack(w, h);
        for (int i = 0; i < d; i++) {
            sumStack.addSlice(new FloatProcessor(w, h));
        }

         // add all masks
        for (String specimen : specimens) {
            String path =Config.mainDir + "Results/04_MaskPithFilled/" + specimen +"_"+ timestamps[t-1]+"_mask.tif";
            ImagePlus mask = IJ.openImage(path);
            for (int z = 1; z <= d; z++) {
                float[] src = (float[]) mask.getStack().getProcessor(z).convertToFloat().getPixels();
                float[] dst = (float[]) sumStack.getProcessor(z).getPixels();
                for (int i = 0; i < src.length; i++) {
                    if (src[i] > 0) dst[i] += 1f; // treat any nonzero as foreground
                }
            }
        }

        // normalize by number of specimens → probability [0,1]
        int n = specimens.length;
        for (int z = 1; z <= d; z++) {
            float[] dst = (float[]) sumStack.getProcessor(z).getPixels();
            for (int i = 0; i < dst.length; i++) {
                dst[i] /= n;
            }
        }

    return new ImagePlus(getVarietyName(variety) + "_probAtlas_t" + timestamps[t-1], sumStack);

 
    }

   public static ImagePlus buildProbabilisticAtlasAllVar(int cond,int t) {
        String ts = Config.timestamps[t - 1];

        // collect all specimen IDs in PCH across all varieties
        java.util.List<String> allSpecs = new java.util.ArrayList<>();
        for (int v = 0; v < VAR_NAMES.length; v++) {
            for (String s : Config.getSpecimensName(cond, v)) allSpecs.add(s);
        }
         // >>> PRINT WHICH SPECIMENS YOU'RE USING <<<
        System.out.println("[ProbAtlas] Condition=PCH  Timestamp=" + ts);
        System.out.println("[ProbAtlas] N specimens = " + allSpecs.size());
        System.out.println("[ProbAtlas] Specimens   = " + String.join(", ", allSpecs));
        // open first to get size
        ImagePlus first = IJ.openImage(Config.mainDir + "Results/04_MaskPithFilled/" + allSpecs.get(0) + "_" + ts + "_mask.tif");
        int w = first.getWidth(), h = first.getHeight(), d = first.getNSlices();

        // accumulator
        ImageStack sum = new ImageStack(w, h);
        for (int z = 0; z < d; z++) sum.addSlice(new FloatProcessor(w, h));

        // add all masks (nonzero → +1)
        int n = 0;
        for (String s : allSpecs) {
            ImagePlus m = IJ.openImage(Config.mainDir + "Results/04_MaskPithFilled/" + s + "_" + ts + "_mask.tif");
            for (int z = 1; z <= d; z++) {
                float[] src = (float[]) m.getStack().getProcessor(z).convertToFloat().getPixels();
                float[] dst = (float[]) sum.getProcessor(z).getPixels();
                for (int i = 0; i < src.length; i++) if (src[i] > 0f) dst[i] += 1f;
            }
            n++;
        }

        // normalize to [0,1]
        for (int z = 1; z <= d; z++) {
            float[] px = (float[]) sum.getProcessor(z).getPixels();
            for (int i = 0; i < px.length; i++) px[i] /= (float) n;
        }

        return new ImagePlus("ProbAtlas_allVar_PCH_" + ts, sum);
    }
    
    public static void binarizeAndFillTopAuto(ImagePlus imp, double coverage) {
        int w = imp.getWidth(), h = imp.getHeight(), d = imp.getNSlices();
        coverage = Math.max(0.0, Math.min(1.0, coverage));
        int needed = Math.max(1, (int) Math.ceil(w * coverage)); // columns needed to count a row as foreground

        for (int z = 1; z <= d; z++) {
            ImageProcessor ip = imp.getStack().getProcessor(z).convertToByte(false); // work in 8-bit
            byte[] px = (byte[]) ip.getPixels();

            // 1) binarize to 0/1
            for (int i = 0; i < px.length; i++) px[i] = (byte) (((px[i] & 0xFF) > 0) ? 1 : 0);

            // 2) scan from top to find first "real" foreground row (robust to noise)
            int stopY = -1;
            for (int y = 0; y < h; y++) {
                int count = 0;
                int off = y * w;
                for (int x = 0; x < w; x++) {
                    if (px[off + x] != 0) {
                        if (++count >= needed) { stopY = y; break; }
                    }
                }
                if (stopY >= 0) break;
            }
            System.out.println("Slice " + z + " needed=" + needed + " stopY=" + stopY + " (h=" + h + ")");

            // 3) set all rows above (and including) stopY to 1
            if (stopY >= 0) {
                for (int y = 0; y <= stopY; y++) {
                    int off = y * w;
                    for (int x = 0; x < w; x++) px[off + x] = 1;
                }
            }

            // write back
            ip.setPixels(px);
            imp.getStack().setProcessor(ip, z);
        }
        imp.updateAndDraw();

    }

    public static void fillTopToFullRowOrBest(ImagePlus imp, double minCoverageGate) {
    int w = imp.getWidth(), h = imp.getHeight(), d = imp.getNSlices();
    int minNeeded = Math.max(1, (int)Math.ceil(w * Math.max(0.0, Math.min(1.0, minCoverageGate)))); // gate

    for (int z = 1; z <= d; z++) {
        ImageProcessor ip = imp.getStack().getProcessor(z).convertToByte(false);
        byte[] px = (byte[]) ip.getPixels();

        // binarize to 0/1
        for (int i = 0; i < px.length; i++) px[i] = (byte)(((px[i] & 0xFF) > 0) ? 1 : 0);

        int stopY = -1;

        // pass 1: look for first 100% row
        for (int y = 0; y < h; y++) {
            int off = y * w;
            int count = 0;
            for (int x = 0; x < w; x++) if (px[off + x] != 0) count++;
            if (count == w) { stopY = y; break; }
        }

        // pass 2 (fallback): pick best-coverage row if no full row exists
        if (stopY < 0) {
            int bestY = -1, bestCount = -1;
            for (int y = 0; y < h; y++) {
                int off = y * w;
                int count = 0;
                for (int x = 0; x < w; x++) if (px[off + x] != 0) count++;
                if (count > bestCount) { bestCount = count; bestY = y; }
            }
            // only use best row if it beats the minimum coverage gate
            if (bestCount >= minNeeded) stopY = bestY;
        }

        // fill top if we have a valid stopY
        if (stopY >= 0) {
            for (int y = 0; y <= stopY; y++) {
                int off = y * w;
                for (int x = 0; x < w; x++) px[off + x] = 1; // use 255 if you prefer 8-bit white
            }
        }

        ip.setPixels(px);
        imp.getStack().setProcessor(ip, z);
    }
    imp.updateAndDraw();
}
    
    public static ImagePlus runFillTopToFullRowOrBest(Specimen specimen, int frame) {
        String[] timestamps = Config.timestamps;
        String specimenName = specimen.getName();
        System.out.println(Config.mainDir+ "Results/03_MaskforPA/" + specimenName + "_" + timestamps[frame-1] + "_mask.tif");
        
        ImagePlus img =  IJ.openImage(Config.mainDir+ "Results/03_MaskforPA/" + specimenName + "_" + timestamps[frame-1] + "_mask.tif");
        fillTopToFullRowOrBest(img, 0.40);
        IJ.saveAsTiff(img, Config.mainDir+ "Results/04_MaskPithFilled/" + specimenName + "_" + timestamps[frame-1] + "_mask.tif");     
        // System.out.println("done!");
        return img;
        
    }
   

    

}

 //     ImagePlus atlas_chard = buildProbabilisticAtlas(cond_PCH, var_CHARD, t);
    //     atlas_chard.show();
    //     IJ.saveAsTiff(atlas_chard, Config.mainDir + "Results/05_ProbabilisticAtlas/" + specimen.getName() + "_ProbabilisticAtlas_CHARD_PCH_"+timestamps[t-1]+".tif");
    //     ImagePlus atlas_mer = buildProbabilisticAtlas(cond_PCH, var_MERLOT, t);
    //     atlas_mer.show();
    //     IJ.saveAsTiff(atlas_mer, Config.mainDir + "Results/05_ProbabilisticAtlas/" + specimen.getName() + "_ProbabilisticAtlas_MER_PCH_"+timestamps[t-1]+".tif");
    //     ImagePlus atlas_temp = buildProbabilisticAtlas(cond_PCH, var_TEMPRA, t);
    //     atlas_temp.show();
    //     IJ.saveAsTiff(atlas_temp, Config.mainDir + "Results/05_ProbabilisticAtlas/" + specimen.getName() + "_ProbabilisticAtlas_TEMP_PCH_"+timestamps[t-1]+".tif");
    //     ImagePlus atlas_ug = buildProbabilisticAtlas(cond_PCH, var_UGNI, t);
    //     atlas_ug.show();
    //     IJ.saveAsTiff(atlas_ug, Config.mainDir + "Results/05_ProbabilisticAtlas/" + specimen.getName() + "_ProbabilisticAtlas_UGNI_PCH_"+timestamps[t-1]+".tif");

      
        
    // }
