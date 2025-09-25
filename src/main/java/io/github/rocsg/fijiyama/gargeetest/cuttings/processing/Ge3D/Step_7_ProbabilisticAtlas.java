package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D;
import java.util.ArrayList;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.Duplicator;
import ij.process.AutoThresholder;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
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
        String[] timestamps = Config.timestamps;
        String specimenName = specimen.getName();
        
        // String path = Config.mainDir + "Processing/03_PolarTransform/" + specimenName + "_GeneralizedPolarTransform.tif";  
         
        // for (int t = 1; t < 5; t++) {
        //     ImagePlus mask = null;
        //     try {
        //         // String path = Config.mainDir + "Results/02_Atlas/PolarAtlas/test_fullpop/mean_all_var_PCH_"+timestamps[t-1]+".tif";  
        //         System.out.println(path);
        //         mask = getMaskWithConnectedComponents(t, path); 
        //         String outPath = Config.getPathToMask()+ specimenName + "_" + timestamps[t-1] + "_mask.tif";
        //         IJ.saveAsTiff(mask, outPath);
        //     } finally {
        //         if (mask != null) mask.close();
        //     }
        // }
        ImagePlus img = getMaskUsingWekaSegmentation(cond_PCH, var_CHARD, 4);
        img.show();
                
        System.out.println("done!");
        
    }
    

    public static ImagePlus getMaskWithConnectedComponents(int frame, String path) {
        System.out.println(path);
        ImagePlus img = IJ.openImage(path);
        // Get image
        if (img.isHyperStack()) {
            img = new Duplicator().run(img, 1, 1, 256, 768, frame, frame);
        } else {
            img = new Duplicator().run(img);
        }

       
        // img.show();
        // Enhance contrast
        IJ.run(img, "Enhance Contrast", "saturated=0.50");
        // Add blur
        IJ.run(img, "Gaussian Blur...", "sigma=1 stack");
        // Compute connected components to get mask
        ImagePlus connectedComponents = VitimageUtils.connexeNoFuckWithVolume(img, 0.2, 2, 400, 2000, 6, 1, true);
        // Convert to 8 bit
        VitimageUtils.convertToGray8(connectedComponents);
        IJ.run(connectedComponents, "Invert", "stack");
        IJ.run(connectedComponents, "Fill Holes", "stack");
        IJ.run(connectedComponents, "Invert", "stack");
        IJ.run(connectedComponents, "Divide...", "value=255 stack");
        connectedComponents.setDisplayRange(0, 1);
        // connectedComponents.show();
        return connectedComponents;
    }

    public static ImagePlus getMaskOfPolarTransforms(int condition, int variety){
        String[] spec = Config.getSpecimensName(condition, variety);
        if (spec == null || spec.length == 0) {
            throw new IllegalStateException("No specimens for condition=" + condition + ", variety=" + variety);
        }

        String path = Config.mainDir + "Processing/03_PolarTransform/" + spec[0] + "_GeneralizedPolarTransform.tif";
        System.out.println(path);
        ImagePlus img = IJ.openImage(path);

        ImagePlus imgToMask = new Duplicator().run(img, 1, 1, 256, 768, 2,2 );
         imgToMask.show();
         imgToMask.setTitle(spec[0] + "_Image");
        // adjust 0.04_0.96
        imgToMask.setDisplayRange(0.39, 1.037);
        imgToMask.show();
        imgToMask.setTitle(spec[0] + "_Adjust_Brightness");
        // gaussian blur 20 rad
        IJ.run(imgToMask, "Gaussian Blur...", "sigma=2 stack");
        imgToMask.show();
        imgToMask.setTitle(spec[0] + "_Gaussian_Blur");
        // 8 bit
        VitimageUtils.convertToGray8(imgToMask);
        imgToMask.show();
        imgToMask.setTitle(spec[0] + "_8_bit");


        // auto local threshold 
        // IJ.run(imgToMask, "Auto Local Threshold", "method=otsu radius=15 white stack" );
        int nSlices = imgToMask.getStackSize();
       
       
        // AutoThresholder thresholder = new AutoThresholder();
        // AutoThresholder.Method method = AutoThresholder.Method.Otsu;
        // for (int i = 1; i<= nSlices; i++){
        //     imgToMask.setSlice(i);
        //     ImageProcessor ip = imgToMask.getProcessor();
        //     int[] histogram = ip.getHistogram();
        //     int threshold = thresholder.getThreshold(method, histogram);
        //     ip.threshold(threshold);
            
        // }
        for (int i = 1; i<= nSlices; i++){
            imgToMask.setSlice(i);
            IJ.runPlugIn(imgToMask, "Auto Local Threshold", "method=otsu radius=15 parameter_1=0 parameter_2=0 white");
        }

        imgToMask.show();
        imgToMask.setTitle(spec[0] + "_Otsu_Threshold");
        
        
        // fill holes
        IJ.run(imgToMask, "Fill Holes", "stack");
        imgToMask.show();
        imgToMask.setTitle(spec[0] + "_Fill_Holes");
        // close
        IJ.run(imgToMask, "Close...", "stack");
        imgToMask.show();
        return imgToMask;
    }
    


    public static ImagePlus createOtsuMask(int condition, int variety){
        String[] spec = Config.getSpecimensName(condition, variety);
       
        if (spec == null || spec.length == 0) {
            throw new IllegalStateException("No specimens for condition=" + condition + ", variety=" + variety);
        }
        String path = Config.mainDir + "Processing/03_PolarTransform/" + spec[0] + "_GeneralizedPolarTransform.tif";   
        ImagePlus hyperFrame = IJ.openImage(path);
        // Get image
        ImagePlus img = new Duplicator().run(hyperFrame, 1, 1, 256, 768,2,2);
        // ImagePlus img  = getMaskWithConnectedComponents(cond_PCH, var_CHARD);
        if (spec == null || spec.length == 0) {
            throw new IllegalStateException("No specimens for condition=" + condition + ", variety=" + variety);
        }

     
       
        // ImagePlus hyperFrame = IJ.openImage(path);
        // ImagePlus img = new Duplicator().run(hyperFrame, 1, 1, 256, 768, 1, 1);
        // if (img == null) {
        //     throw new IllegalStateException("Not Found: " + path);
        // }

        ImageStack stack = img.getStack();
        double min = img.getProcessor().getMin();
        double max = img.getProcessor().getMax();
        long[] hist = new long[256];


        for(int z= 1; z<stack.getSize(); z++){
            ImageProcessor ip = stack.getProcessor(z);
            // ip.blurGaussian(2.0);
            for(int y = 0; y<ip.getHeight(); y++){
                for(int x=0; x<ip.getWidth(); x++){
                    double v = ip.getf(x,y);
                    int b = (int)((v-min)*255.0/(max-min));
                    if(b<0)b=0; if (b>255) b=255;
                    hist[b]++;
                }
            }
        }
        //Otsu cutoff
        int t = otsuThreshold(hist);
        double cutoff = min +(max+min)*(t/255.0);

        ImageStack mask = new ImageStack(stack.getWidth(), stack.getHeight());
        for (int z=1; z<=stack.getSize();z++){
            ImageProcessor ip = stack.getProcessor(z);
            ByteProcessor bp = new ByteProcessor(ip.getWidth(),ip.getHeight());
            for(int y = 0; y <ip.getHeight(); y++){
                for(int x=0; x< ip.getWidth(); x++){
                    bp.set(x,y,ip.getf(x,y)>= cutoff ? 255 : 0);
                }
            }
            mask.addSlice(bp);
        }
        ImagePlus maskOtsu = new ImagePlus("OtsuMask",mask);
        maskOtsu.show();
        return maskOtsu;
    }

    // Compute Otsu threshold index from a histogram
    public static int otsuThreshold(long[] hist) {
        long total = 0, sum = 0;
        for (int i = 0; i < hist.length; i++) {
            total += hist[i];
            sum += (long) i * hist[i];
        }
        long sumB = 0, wB = 0;
        double maxBetween = -1.0;
        int threshold = 0;
        for (int t = 0; t < hist.length; t++) {
            wB += hist[t];
            if (wB == 0) continue;
            long wF = total - wB;
            if (wF == 0) break;
            sumB += (long) t * hist[t];
            double mB = (double) sumB / wB;
            double mF = (double) (sum - sumB) / wF;
            double between = (double) wB * (double) wF * (mB - mF) * (mB - mF);
            if (between > maxBetween) {
                maxBetween = between;
                threshold = t;
            }
        }
        return threshold;
    }

    
    public static ImagePlus computeOtsuMask(int condition, int variety) {
        int t = 2;
        String[] spec = Config.getSpecimensName(condition, variety);
        String path = Config.mainDir + "Processing/03_PolarTransform/" + spec[0] + "_GeneralizedPolarTransform.tif";
        ImagePlus imgT1 = getMaskWithConnectedComponents(t,path);
       
        if (spec == null || spec.length == 0) {
            throw new IllegalStateException("No specimens for condition=" + condition + ", variety=" + variety);
        }

        // // for (int cond = condition; cond <= condition; cond++) {
        // //     for (int i = 0; i < spec.length; i++) {
        //         String path = Config.mainDir + "Processing/03_PolarTransform/" + spec[0] + "_GeneralizedPolarTransform.tif";
        //         System.out.println(path);

        //         ImagePlus img = IJ.openImage(path);
        //         if (img == null) {
        //             throw new IllegalStateException("Not Found: " + path);
        //         }
        //         int n = img.getNSlices();
        //         int start = 256;
        //         int end   = Math.max(start, n - 256);

        //         ImagePlus imgT1 = new Duplicator().run(img, t, t, start, end, t, t);

                // IJ.run(imgT1, "Gaussian Blur...", "sigma=2 stack");
                // IJ.run(imgT1, "Median...", "radius=1 stack");
                // imgT1.setDisplayRange(0.00, 1.80);
                // imgT1.updateAndDraw();

                IJ.setAutoThreshold(imgT1, "Default dark no-reset");
                double threshold = imgT1.getProcessor().getMinThreshold();

                for (int z = 1; z <= imgT1.getNSlices(); z++) {
                    ImageProcessor ip = imgT1.getStack().getProcessor(z).convertToFloatProcessor();
                    float[] px = (float[]) ((FloatProcessor) ip).getPixels();
                    for (int j = 0; j < px.length; j++) {
                        px[j] = (px[j] > threshold) ? 1.0f : 0.0f;
                    }
                    imgT1.getStack().setProcessor(ip, z);
                }
            
            imgT1.setTitle(spec[0] + "_Otsu_Mask");
            return imgT1; 
            }

    public static ImagePlus getMaskUsingWekaSegmentation(int condition, int variety, int frame) {
        int t = 2;
        String[] spec = Config.getSpecimensName(condition, variety);
        String path = Config.mainDir + "Processing/03_PolarTransform/" + spec[0] + "_GeneralizedPolarTransform.tif";
        System.out.println(path);
        ImagePlus hyperframe = IJ.openImage(path);
        ImagePlus img = new Duplicator().run(hyperframe, 1, 1, 256, 768, frame, frame);
        ImagePlus imgPti = new Duplicator().run(hyperframe, 1, 1, 256, 257, frame, frame);

        System.out.println("Toto 1");
        WekaSegmentation weka = new WekaSegmentation(imgPti);
        System.out.println("Toto 2");
        System.out.println(Config.mainDir +"classifier.model");
        System.out.println("Toto 3");
        weka.loadClassifier(Config.mainDir +"classifier.model" );
        System.out.println("Toto 4");
        ImagePlus proba = weka.applyClassifier(img, 0, true);
        System.out.println("Toto 5");
        ImagePlus imgMask=new Duplicator().run(proba,1,1,1,proba.getNSlices(),1,1);
        System.out.println("Toto 6");
        imgMask.setDisplayRange(0.5, 0.5);
        return imgMask;
    }

}
