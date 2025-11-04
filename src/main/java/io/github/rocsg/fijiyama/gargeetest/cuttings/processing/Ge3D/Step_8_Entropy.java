package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.Duplicator;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.common.VitimageUtils;


import java.util.Arrays;

public class Step_8_Entropy implements PipelineStep {
     @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen, true);
    }

    public static void main (String[] args) throws Exception {
        Specimen spec = new Specimen("B_201");
        
        ImageJ ij = new ImageJ();
        new  Step_8_Entropy().execute(spec,true);
    }

    public void execute (Specimen specimen, boolean testing) throws Exception {
        String[] spec = Config.getSpecimensName(specimen.getCondition(), specimen.getVariety());
        ImagePlus img = IJ.openImage(Config.mainDir + "Processing/03_PolarTransform/" + spec[0] + "_GeneralizedPolarTransform.tif");
        img.show();
        ImagePlus imgP = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, 2, 2);
        int nBins = 500;
        double[] mmP = findMinMax(imgP);
        double minValP = mmP[0], maxValP = mmP[1];
        // System.out.println("min=" + minValP + " max=" + maxValP);
        // ImagePlus entropyPolar = entropyStack(imgP, nBins, minValP, maxValP, 2);
        // entropyPolar.show();
        ImagePlus imp = IJ.openImage(Config.mainDir + "Results/01_Hyperstack/" + spec[0] + "_Hyperstack.tif");
        imp.show();
        ImagePlus imgC = new Duplicator().run(imp, 1, 1, 256, img.getNSlices() - 256, 2, 2);
        double[] mmC = findMinMax(imgC);
        double minValC = mmC[0], maxValC = mmC[1];
        // System.out.println("min=" + minValC + " max=" + maxValC);
        // ImagePlus entropyCart = entropyStack3D(imgC, nBins, minValC, maxValC, 2);
        // entropyCart.show();
        double invp = 1.0 / (maxValP - minValP);
        double invc = 1.0 / (maxValC - minValC);
        double Hc = entropyGlobal(imgC,  nBins, minValC, maxValC, invc);
        double Hp = entropyGlobal(imgP, nBins, minValP, maxValP, invp);
        double d  = Hc - Hp;
        System.out.println("Hc=" + Hc + " Hp=" + Hp + " d=" + d);
    }


    public static ImagePlus entropyStack(ImagePlus img, int nBins, double minVal, double maxVal, int radius) {
        int W = img.getWidth();
        int H = img.getHeight();
        int Z = img.getStackSize();

        ImageStack outStack = new ImageStack(W, H);

        for (int z = 1; z <= Z; z++) {
            ImageProcessor ip = img.getStack().getProcessor(z).convertToFloatProcessor();
            FloatProcessor fpOut = new FloatProcessor(W, H);

            long[] hist = new long[nBins];

            for (int y = 0; y < H; y++) {
                int y0 = Math.max(0, y - radius), y1 = Math.min(H - 1, y + radius);
                for (int x = 0; x < W; x++) {
                    int x0 = Math.max(0, x - radius), x1 = Math.min(W - 1, x + radius);

                    Arrays.fill(hist, 0L);
                    long count = 0;

                    for (int yy = y0; yy <= y1; yy++) {
                        for (int xx = x0; xx <= x1; xx++) {
                            double v = ip.getf(xx, yy);
                            int bin = (int) ((v - minVal) / (maxVal - minVal) * (nBins - 1));
                            if (bin >= 0 && bin < nBins) {
                                hist[bin]++;
                                count++;
                            }
                        }
                    }
                    fpOut.setf(x, y, (float) entropy(hist, count));
                }
            }
            outStack.addSlice(fpOut);
        }
        System.out.println("Returning entropy stack...");
        return new ImagePlus("Entropy", outStack);
    }
    private static double entropy(long[] hist, long total) {
        if (total <= 0) return 0.0;
        double H = 0.0;
        double log2 = Math.log(2);
        for (long h : hist) {
            if (h > 0) {
                double p = (double) h / total;
                H -= p * (Math.log(p) / log2);
            }
        }
        return H;

      
    }
      private static double[] findMinMax(ImagePlus img) {
        double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
        for (int z = 1; z <= img.getStackSize(); z++) {
            FloatProcessor fp = img.getStack().getProcessor(z).convertToFloatProcessor().toFloat(0, null);
            float[] a = (float[]) fp.getPixels();
            for (float v : a) if (Float.isFinite(v)) {
                if (v < min) min = v;
                if (v > max) max = v;
            }
        }
        if (!(min < max)) { max = min + 1e-6; } 
        return new double[]{min, max};
    }

    public static ImagePlus entropyStack3D(ImagePlus img, int nBins, double minVal, double maxVal, int radius) {

        final int W = img.getWidth();
        final int H = img.getHeight();
        final int Z = img.getStackSize();
        final double invRange = 1.0 / (maxVal - minVal);

        // Pre-convert all slices to float arrays for speed
        final float[][] pix = new float[Z][];
        for (int z = 0; z < Z; z++) {
            ImageProcessor ip = img.getStack().getProcessor(z + 1).convertToFloatProcessor();
            pix[z] = (float[]) ((FloatProcessor) ip).getPixels();
        }

        ImageStack outStack = new ImageStack(W, H);

        // Reusable histogram
        final long[] hist = new long[nBins];

        for (int z = 0; z < Z; z++) {
            int z0 = Math.max(0, z - radius), z1 = Math.min(Z - 1, z + radius);

            FloatProcessor fpOut = new FloatProcessor(W, H);
            float[] out = (float[]) fpOut.getPixels();

            for (int y = 0; y < H; y++) {
                int y0 = Math.max(0, y - radius), y1 = Math.min(H - 1, y + radius);
                final int rowBase = y * W;

                for (int x = 0; x < W; x++) {
                    int x0 = Math.max(0, x - radius), x1 = Math.min(W - 1, x + radius);

                    Arrays.fill(hist, 0L);
                    long count = 0L;
                    for (int zz = z0; zz <= z1; zz++) {
                        float[] pz = pix[zz];
                        for (int yy = y0; yy <= y1; yy++) {
                            int base = yy * W;
                            for (int xx = x0; xx <= x1; xx++) {
                                float v = pz[base + xx];
                                if (!Float.isFinite(v)) continue;

                                int bin = (int) Math.floor((v - minVal) * invRange * (nBins - 1));
                                if (bin < 0) bin = 0;
                                else if (bin >= nBins) bin = nBins - 1;

                                hist[bin]++; count++;
                            }
                        }
                    }

                    out[rowBase + x] = (float) entropy(hist, count);
                }
            }
            outStack.addSlice(fpOut);
        }

        ImagePlus outImg = new ImagePlus("Entropy (3D)", outStack);
        outImg.setCalibration(img.getCalibration());
        return outImg;
    }
    static double entropyGlobal(ImagePlus img, int bins, double min, double max, double invRange){
        long[] hist=new long[bins]; long total=0; int Z=img.getStackSize(), W=img.getWidth(), H=img.getHeight();
        for(int z=1; z<=Z; z++){
            ImageProcessor ip=img.getStack().getProcessor(z).convertToFloatProcessor();
            float[] p=(float[])((FloatProcessor)ip).getPixels();
            for(float v: p){ if(!Float.isFinite(v)) continue;
                int b=(int)Math.floor((v-min)*invRange*(bins-1));
                if(b<0) b=0; else if(b>=bins) b=bins-1; hist[b]++; total++;
            }
        }
        if(total==0) return 0.0;
        double Hh=0, log2=Math.log(2.0);
        for(long h: hist) if(h>0){ double p=(double)h/total; Hh -= p*(Math.log(p)/log2); }
        return Hh;
    }
   
   
}
