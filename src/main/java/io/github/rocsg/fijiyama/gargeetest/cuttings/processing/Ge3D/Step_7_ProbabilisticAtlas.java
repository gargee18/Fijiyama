package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.process.BinaryProcessor;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.gargeetest.cuttings.helpers.ImgProUtils;

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
        Specimen spec = new Specimen("B_213");
        
        ImageJ ij = new ImageJ();
        new Step_7_ProbabilisticAtlas().execute(spec,true);
    }

    public void execute(Specimen specimen, boolean testing) throws Exception {
        String[] timestamps = Config.timestamps;
        
        
       
     
        // ALL VARIETIES     
        // for(int t = 1; t <5; t++){
        //     ImagePlus atlas = buildProbabilisticAtlasAllVar(cond_CONTROL, t);
        //     VitimageUtils.setLutToFire(atlas);
        //     atlas.setDisplayRange(0,1);
        //     // atlas.show();
        //     System.out.println(Config.mainDir+"Results/03_ProbabilisticAtlas/04_ProbAtlasInvertedZeroBelow100/ProbabilisticAtlas_allVar_CT_"+timestamps[t-1]+".tif");
        //     IJ.saveAsTiff(atlas, Config.mainDir+"Results/03_ProbabilisticAtlas/04_ProbAtlasInvertedZeroBelow100/ProbabilisticAtlas_allVar_CT_"+timestamps[t-1]+".tif");
        // }

        // Specific Varieties
        // for(int t = 1; t <5; t++){
        //     ImagePlus atlas_chard = buildProbabilisticAtlas(cond_CONTROL, var_CHARD, t);
        //     VitimageUtils.setLutToFire(atlas_chard);
        //     atlas_chard.setDisplayRange(0,1);
        //     IJ.saveAsTiff(atlas_chard, Config.mainDir + "Results/03_ProbabilisticAtlas/04_ProbAtlasInvertedZeroBelow100/ProbabilisticAtlas_CHARD_CT_"+timestamps[t-1]+".tif");
        //     ImagePlus atlas_mer = buildProbabilisticAtlas(cond_CONTROL, var_MERLOT, t);
        //     VitimageUtils.setLutToFire(atlas_mer);
        //     atlas_mer.setDisplayRange(0,1);
        //     IJ.saveAsTiff(atlas_mer, Config.mainDir + "Results/03_ProbabilisticAtlas/04_ProbAtlasInvertedZeroBelow100/ProbabilisticAtlas_MER_CT_"+timestamps[t-1]+".tif");
        //     ImagePlus atlas_temp = buildProbabilisticAtlas(cond_CONTROL, var_TEMPRA, t);
        //     VitimageUtils.setLutToFire(atlas_temp);
        //     atlas_temp.setDisplayRange(0,1);
        //     IJ.saveAsTiff(atlas_temp, Config.mainDir + "Results/03_ProbabilisticAtlas/04_ProbAtlasInvertedZeroBelow100/ProbabilisticAtlas_TEMP_CT_"+timestamps[t-1]+".tif");
        //     ImagePlus atlas_ug = buildProbabilisticAtlas(cond_CONTROL, var_UGNI, t);
        //     VitimageUtils.setLutToFire(atlas_ug);
        //     atlas_ug.setDisplayRange(0,1);
        //     IJ.saveAsTiff(atlas_ug, Config.mainDir + "Results/03_ProbabilisticAtlas/04_ProbAtlasInvertedZeroBelow100/ProbabilisticAtlas_UGNI_CT_"+timestamps[t-1]+".tif");
        // }
        // System.out.println("done!");
        // for(int t = 1; t <5; t++){
        //     ImagePlus img = new ImagePlus(Config.mainDir + "Processing/04_Masks/03_MaskPithFilled/"+specimen.getName()+"_"+timestamps[t-1]+"_mask.tif");
        //     ImagePlus invert = invertBinaryMask(img);
        //     invert.setDisplayRange(0,1);
        //     IJ.saveAsTiff(invert, Config.mainDir + "Processing/04_Masks/04_MaskPithFilledInverted/"+specimen.getName()+"_"+timestamps[t-1]+"_mask.tif");
        // }
         
    }


   
    public static ImagePlus buildProbabilisticAtlas(int condition, int variety, int t) {
        String[] timestamps = Config.timestamps;
        String[] specimens = Config.getSpecimensName(condition, variety);
        ImagePlus mask1 = IJ.openImage(Config.mainDir + "Processing/04_Masks/06_MaskPithFilledInvertedZeroBelow100/" + specimens[0] +"_"+ timestamps[t-1]+"_mask.tif");
        
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
            String path =Config.mainDir + "Processing/04_Masks/06_MaskPithFilledInvertedZeroBelow100/" + specimen +"_"+ timestamps[t-1]+"_mask.tif";
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
        ImagePlus first = IJ.openImage(Config.mainDir + "Processing/04_Masks/06_MaskPithFilledInvertedZeroBelow100/" + allSpecs.get(0) + "_" + ts + "_mask.tif");
        int w = first.getWidth(), h = first.getHeight(), d = first.getNSlices();

        // accumulator
        ImageStack sum = new ImageStack(w, h);
        for (int z = 0; z < d; z++) sum.addSlice(new FloatProcessor(w, h));

        // add all masks (nonzero → +1)
        int n = 0;
        for (String s : allSpecs) {
            ImagePlus m = IJ.openImage(Config.mainDir + "Processing/04_Masks/06_MaskPithFilledInvertedZeroBelow100/" + s + "_" + ts + "_mask.tif");
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
    


}

