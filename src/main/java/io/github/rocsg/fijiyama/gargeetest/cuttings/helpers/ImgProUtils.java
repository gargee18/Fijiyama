package io.github.rocsg.fijiyama.gargeetest.cuttings.helpers;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.Duplicator;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import trainableSegmentation.WekaSegmentation;

public class ImgProUtils {

    public static void runWeka(Specimen specimen) {
        String[] timestamps = Config.timestamps;
        String specimenName = specimen.getName();
        System.out.println( Config.mainDir + "Processing/03_PolarTransform/" + specimenName + "_GeneralizedPolarTransform.tif");
        String path =  Config.mainDir + "Processing/03_PolarTransform/" + specimenName + "_GeneralizedPolarTransform.tif";
        ImagePlus img = trainWekaToGetMask(path, 3);
        IJ.saveAsTiff(img, Config.mainDir+ "Results/04_MaskWeka/" + specimenName + "_" + timestamps[3-1] + "_mask.tif");     
        System.out.println("done!");
        
    }
    
  
    public static ImagePlus trainWekaToGetMask(String path, int frame){
        // String[] spec = Config.getSpecimensName(condition, variety);
        // ImagePlus hyperframe = IJ.openImage(Config.mainDir + "Processing/03_PolarTransform/" + spec[0] + "_GeneralizedPolarTransform.tif");
        ImagePlus hyperframe = IJ.openImage(path);
        ImagePlus img = new Duplicator().run(hyperframe, 1, 1, 256, 768, frame, frame);
        WekaSegmentation weka = new WekaSegmentation(img);
        weka.loadClassifier(Config.mainDir +"segment_pch_all_var_t3.model" );
        ImagePlus proba = weka.applyClassifier(img, 0, true);
        ImagePlus imgMask=new Duplicator().run(proba,1,1,1,proba.getNSlices(),1,1);
        imgMask.setDisplayRange(0.5, 0.5);
        VitimageUtils.convertToGray8(imgMask);
        IJ.run(imgMask, "Invert", "stack");
        IJ.run(imgMask, "Fill Holes", "stack");
        IJ.run(imgMask, "Invert", "stack");
        IJ.run(imgMask, "Divide...", "value=255 stack");      
        imgMask.setDisplayRange(0, 1);
        IJ.run(imgMask,"Median...", "radius=1 stack");
        return imgMask;
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

    public static ImagePlus buildWekaTrainingStack(int condition, int variety, int frame) {
        String[] specimens = Config.getSpecimensName(condition, variety);

        ImageStack out = null;

        for (int i = 0; i < specimens.length; i++) {

            String specimen = specimens[i];
            String path = Config.mainDir + "Processing/03_PolarTransform/" + specimen + "_GeneralizedPolarTransform.tif";

            IJ.log("Opening: " + path);
            ImagePlus hyper = IJ.openImage(path);
            if (hyper == null) {
                IJ.log("WARNING: couldn't open " + path + " (skipping).");
                continue;
            }

            try {
                // Sanity checks
                if (frame < 1 || frame > hyper.getNFrames()) {
                    IJ.log("WARNING: " + specimen + " has only " + hyper.getNFrames() + " frames; requested T=" + frame + " (skipping).");
                    continue;
                }
                if (hyper.getNSlices() < 627) {
                    IJ.log("WARNING: " + specimen + " has only " + hyper.getNSlices() + " Z-slices; need 627 (skipping).");
                    continue;
                }

                // Extract Z=512, T=frame (channel 1); ImageJ indexing is 1-based
                ImagePlus z512 = new Duplicator().run(hyper, 1, 1, 512, 512, frame, frame);
                ImagePlus z627 = new Duplicator().run(hyper, 1, 1, 627, 627, frame, frame);

                // Initialize output stack with correct dimensions & type
                if (out == null) {
                    out = new ImageStack(z512.getWidth(), z512.getHeight());
                }

                // Add slices with informative labels
                ImageProcessor ip512 = z512.getProcessor();
                ImageProcessor ip627 = z627.getProcessor();

                out.addSlice(specimen + "_Z512_T" + frame, ip512.duplicate());
                out.addSlice(specimen + "_Z627_T" + frame, ip627.duplicate());

            } finally {
                hyper.close();
            }
        }

        if (out == null || out.getSize() == 0) {
            IJ.log("No slices were added. Check file paths and indices.");
            return null;
        }

        ImagePlus result = new ImagePlus("Weka_TrainStack_T" + frame + "_Z512_627", out);
        return result;
    }
    
    public static ImagePlus combineStacks(String name, ImagePlus... stacks) {
        if (stacks == null || stacks.length == 0) return null;

        // Use dimensions of the first stack
        int w = stacks[0].getWidth();
        int h = stacks[0].getHeight();
        ImageStack out = new ImageStack(w, h);

        for (ImagePlus imp : stacks) {
            if (imp == null) continue;
            ImageStack s = imp.getStack();

            for (int i = 1; i <= s.getSize(); i++) {
                ImageProcessor ip = s.getProcessor(i);
                String label = imp.getTitle() + "_" + s.getSliceLabel(i);
                out.addSlice(label, ip.duplicate());
            }
        }

        return new ImagePlus(name, out);
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

}
