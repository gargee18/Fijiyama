package io.github.rocsg.fijiyama.gargeetest.cuttings.helpers;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.plugin.Duplicator;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import trainableSegmentation.WekaSegmentation;

public class ImgProUtils {

    public static ImagePlus runWeka(Specimen specimen, int frame) {
        String[] timestamps = Config.timestamps;
        String specimenName = specimen.getName();
        System.out.println( Config.mainDir + "Processing/03_PolarTransform/" + specimenName + "_GeneralizedPolarTransform.tif");
        String path =  Config.mainDir + "Processing/03_PolarTransform/" + specimenName + "_GeneralizedPolarTransform.tif";
        ImagePlus img = trainWekaToGetMask(path, frame);   
        System.out.println("done!");
        return img;
        
    }
    
  
    public static ImagePlus trainWekaToGetMask(String path, int frame){
        String[] timestamps = Config.timestamps;
        // String[] spec = Config.getSpecimensName(condition, variety);
        // ImagePlus hyperframe = IJ.openImage(Config.mainDir + "Processing/03_PolarTransform/" + spec[0] + "_GeneralizedPolarTransform.tif");
        ImagePlus hyperframe = IJ.openImage(path);
        ImagePlus img = new Duplicator().run(hyperframe, 1, 1, 256, 767, frame, frame);
        WekaSegmentation weka = new WekaSegmentation(img);
        weka.loadClassifier(Config.mainDir +"Processing/Trained_Models/All_Var_PCH_z512_545_610_t"+timestamps[frame-1]+".model" );
        ImagePlus proba = weka.applyClassifier(img, 0, true);
        ImagePlus imgMask=new Duplicator().run(proba,1,1,1,proba.getNSlices(),1,1);
        imgMask.setDisplayRange(0.5, 0.5);
        VitimageUtils.convertToGray8(imgMask);      
        IJ.run(imgMask,"Median...", "radius=1 stack");
        IJ.run(imgMask, "Divide...", "value=255 stack"); 
        imgMask.setDisplayRange(0, 1);
        return imgMask;
    }

    

    public static ImagePlus getMaskUsingWekaSegmentation(int condition, int variety, int frame) {
        int t = 2;
        String[] spec = Config.getSpecimensName(condition, variety);
        String path = Config.mainDir + "Processing/03_PolarTransform/" + spec[0] + "_GeneralizedPolarTransform.tif";
        System.out.println(path);
        ImagePlus hyperframe = IJ.openImage(path);
        ImagePlus img = new Duplicator().run(hyperframe, 1, 1, 256, 767, frame, frame);
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
        for (int i = 0; i < specimens.length; i++) { String specimen = specimens[i]; String path = Config.mainDir + "Processing/03_PolarTransform/" + specimen + "_GeneralizedPolarTransform.tif"; 
        IJ.log("Opening: " + path); ImagePlus hyper = IJ.openImage(path); if (hyper == null) { IJ.log("WARNING: couldn't open " + path + " (skipping)."); continue; }
        try { // Sanity checks 
            if (frame < 1 || frame > hyper.getNFrames()) { IJ.log("WARNING: " + specimen + " has only " + hyper.getNFrames() + " frames; requested T=" + frame + " (skipping).");
                continue; 
                } 
                if (hyper.getNSlices() < 610) { IJ.log("WARNING: " + specimen + " has only " + hyper.getNSlices() + " Z-slices; need 610 (skipping).");
                continue; 
                } // Extract Z=512, T=frame (channel 1); ImageJ indexing is 1-based 
                ImagePlus z512 = new Duplicator().run(hyper, 1, 1, 512, 512, frame, frame); 
                ImagePlus z545 = new Duplicator().run(hyper, 1, 1, 545, 545, frame, frame); 
                ImagePlus z610 = new Duplicator().run(hyper, 1, 1, 610, 610, frame, frame); 
            
                // Initialize output stack with correct dimensions & type 
                if (out == null) { out = new ImageStack(z512.getWidth(), z512.getHeight()); } 
                // Add slices with informative labels 
                ImageProcessor ip512 = z512.getProcessor(); 
                ImageProcessor ip545 = z545.getProcessor(); 
                ImageProcessor ip610 = z610.getProcessor(); 
                out.addSlice(specimen + "_Z512_T" + frame, ip512.duplicate()); 
                out.addSlice(specimen + "_Z545_T" + frame, ip545.duplicate()); 
                out.addSlice(specimen + "_Z610_T" + frame, ip610.duplicate()); 
            } 
            finally { hyper.close(); } 
        } 
        if (out == null || out.getSize() == 0) { IJ.log("No slices were added. Check file paths and indices."); 
        return null; 
        } 
        ImagePlus result = new ImagePlus("Weka_TrainStack_T" + frame + "_Z512_545_610", out); return result; 
    }

    public static ImagePlus buildWekaTrainingFullStack(int condition, int variety, int frame) {

        /* DEBUG
         * for(int t = 1; t <5; t++){
            ImagePlus stack_CHARD = ImgProUtils.buildWekaTrainingStack(cond_PCH, var_CHARD, t);
            ImagePlus stack_MER = ImgProUtils.buildWekaTrainingStack(cond_PCH, var_MERLOT, t);
            ImagePlus stack_TEMP= ImgProUtils.buildWekaTrainingStack(cond_PCH, var_TEMPRA, t);
            ImagePlus stack_UG = ImgProUtils.buildWekaTrainingStack(cond_PCH, var_UGNI, t);
            ImagePlus result = ImgProUtils.combineStacks("All_Var_PCH_z512_627_t"+timestamps[t-1]+".tif",  stack_CHARD   ,stack_MER, stack_TEMP, stack_UG); 
            result.show(); 
            IJ.saveAsTiff(result, Config.mainDir + "Processing/04_Masks/01_StacksForWeka/" + "All_Var_PCH_z350_650_t"+timestamps[t-1]+".tif");
        }
         */
        String[] specimens = Config.getSpecimensName(condition, variety);

        ImageStack out = null;
        final int zStart = 350; // inclusive
        final int zEnd   = 650; // inclusive
        final int cStart = 1, cEnd = 1; // channel range (modify if needed)

        for (String specimen : specimens) {

            String path = Config.mainDir + "Processing/03_PolarTransform/" +
                        specimen + "_GeneralizedPolarTransform.tif";
            IJ.log("Opening: " + path);

            ImagePlus hyper = IJ.openImage(path);
            if (hyper == null) {
                IJ.log("WARNING: couldn't open " + path + " (skipping).");
                continue;
            }

            // Sanity checks
            int nT = hyper.getNFrames();
            int nZ = hyper.getNSlices();

            if (frame < 1 || frame > nT) {
                IJ.log("WARNING: " + specimen + " has only " + nT +
                    " frames; requested T=" + frame + " (skipping).");
                continue;
            }
            if (nZ < zEnd) {
                IJ.log("WARNING: " + specimen + " has only " + nZ +
                    " Z-slices; need >= " + zEnd + " (skipping).");
                continue;
            }

            // Duplicate full Z-band [zStart..zEnd] at requested T
            ImagePlus zRange = new Duplicator().run(
                    hyper, cStart, cEnd, zStart, zEnd, frame, frame
            );

            // Init output stack once
            if (out == null) {
                out = new ImageStack(zRange.getWidth(), zRange.getHeight());
            }

            // Append slices into training stack
            ImageStack zs = zRange.getStack();
            int size = zs.getSize();

            for (int i = 1; i <= size; i++) {
                int zIdx = zStart + i - 1;
                out.addSlice(specimen + "_Z" + zIdx + "_T" + frame,
                            zs.getProcessor(i).duplicate());
            }

            // free original stack from RAM
            hyper.flush();
            zRange.flush();
        }

        if (out == null || out.getSize() == 0) {
            IJ.log("No slices were added. Check file paths and indices.");
            return null;
        }

        return new ImagePlus("Weka_TrainStack_T" + frame +"_Z" + zStart + "-" + zEnd,out);
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

    public static long countForegroundVoxels(ImagePlus mask, int foregroundValue) {
        if (mask == null) return 0;

        long count = 0;
        int W = mask.getWidth();
        int H = mask.getHeight();
        int Z = mask.getNSlices();

        ImageStack stack = mask.getStack();

        for (int z = 0; z < Z; z++) {
            byte[] pixels = (byte[]) stack.getProcessor(z + 1).getPixels();
            for (int i = 0; i < W * H; i++) {
                if ((pixels[i] & 0xFF) == foregroundValue) {
                    count++;
                }
            }
        }
        return count;
    }


    public static ImagePlus getMaskWithConnectedComponents(int frame, String path) {
        System.out.println(path);
        ImagePlus img = IJ.openImage(path);
        // Get image
        if (img.isHyperStack()) {
            img = new Duplicator().run(img, 1, 1, 256, 767, frame, frame);
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

        ImagePlus imgToMask = new Duplicator().run(img, 1, 1, 256, 767, 2,2 );
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

        int nSlices = imgToMask.getStackSize();
    
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

        /**
         * Fills the top rows of the given image to full coverage rows.
         * If no full coverage row exists, it will pick the best coverage row
         * if it beats the minimum coverage gate.
         * 
         * @param imp the image to fill
         * @param minCoverageGate the minimum coverage gate for the best coverage row
         * @return the filled image
         */
    public static ImagePlus fillTopToFullRowOrBest(ImagePlus imp, double minCoverageGate) {
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
        return imp;
    }

    /**
     * Crops non-tissue at the given y-coordinate from the given image.
     * 
     * @param imp the image to crop
     * @param y the y-coordinate of the top row
     * @return the cropped image
     */
    public static ImagePlus cropNonTissueAtGivenCoords(ImagePlus imp, int y){
        // imp.show();
        int w = imp.getWidth(), h = Math.max(0, Math.min(y, imp.getHeight()));

        ImageStack croppedStack = new ImageStack(w,h);

        for (int z = 0; z < imp.getNSlices(); z++) {
            ImageProcessor ip = imp.getStack().getProcessor(z+1);
            ip.setRoi(new Roi(0,0,w,h));
            ImageProcessor cropped = ip.crop();
            croppedStack.addSlice(cropped);

        }
        
        ImagePlus croppedMask = new ImagePlus("croppedMask", croppedStack);
        croppedMask.copyScale(imp);
        // croppedMask.show();

        // ImagePlus keepBigBlob = largest3DBlob(croppedMask);
        // keepBigBlob.setTitle("largest3D_binary");
        // keepBigBlob.show();

        return croppedMask;
    }


    /**
     * Set all pixels in the region [y .. H-1] to zero for all slices of the given image.
     * The region is defined by the given y-coordinate (exclusive) and spans the full width of the image.
     * 
     * @param imp the image to modify
     * @param yCutExclusive the y-coordinate of the top row to zero (exclusive)
     * @return the modified image
     */
    public static ImagePlus zeroBelowY(ImagePlus imp, int yCutExclusive) {
        if (imp == null || imp.getStackSize() < 1)
            throw new IllegalArgumentException("zeroBelowY: input has no slices");

        int W = imp.getWidth(), H = imp.getHeight(), N = imp.getStackSize();
        int y = Math.max(0, Math.min(yCutExclusive, H)); // clamp

        ImageStack out = new ImageStack(W, H);

        for (int s = 1; s <= N; s++) {
            ImageProcessor ip = imp.getStack().getProcessor(s).duplicate();

            // Zero the region [y .. H-1] across full width
            if (y < H) {
                ip.setRoi(new Roi(0, y, W, H - y));
                ip.setValue(0);
                ip.fill();
                ip.resetRoi();
            }

            out.addSlice(ip);
        }

        ImagePlus res = new ImagePlus(imp.getShortTitle() + "_zeroBelow_" + y, out);
        res.copyScale(imp);
        return res;


    }
    /**
     * Given a binary 3D image (mask), find the largest 3D blob and return it as a new binary 3D image.
     * The output image will have the same dimensions and calibration as the input image.
     * The largest blob is defined as the connected component with the maximum number of pixels.
     * 
     * @param mask the input binary 3D image
     * @return the largest 3D blob as a binary 3D image
     */
    public static ImagePlus largest3DBlob(ImagePlus mask){
        int W=mask.getWidth(), H=mask.getHeight(), Z=mask.getNSlices();
        int N=W*H*Z;
        byte[][] sl=new byte[Z][]; for(int z=0;z<Z;z++) sl[z]=(byte[])mask.getStack().getProcessor(z+1).getPixels();
        int[] lab=new int[N], st=new int[N]; int cur=0,bLab=0,bCnt=0;

        int[] dz={0,0,0,0,1,-1}, dy={0,1,0,-1,0,0}, dx={1,0,-1,0,0,0}; // 6-neigh

        for(int z=0;z<Z;z++) for(int y=0;y<H;y++) for(int x=0;x<W;x++){
            int i=x+W*(y+H*z); if(lab[i]!=0 || (sl[z][x+W*y]&0xff)==0) continue;
            int cnt=0, sp=0; lab[i]=++cur; st[sp++]=i;
            while(sp>0){
                int p=st[--sp], zz=p/(W*H), rem=p%(W*H), yy=rem/W, xx=rem%W;
                for(int k=0;k<6;k++){
                    int xn=xx+dx[k], yn=yy+dy[k], zn=zz+dz[k];
                    if(xn<0||yn<0||zn<0||xn>=W||yn>=H||zn>=Z) continue;
                    int q=xn+W*(yn+H*zn);
                    if(lab[q]==0 && (sl[zn][xn+W*yn]&0xff)>0){ lab[q]=cur; st[sp++]=q; }
                }
                cnt++;
            }
            if(cnt>bCnt){ bCnt=cnt; bLab=cur; }
        }

        ImageStack out=new ImageStack(W,H);
        for(int z=0;z<Z;z++){
            byte[] dst=new byte[W*H];
            for(int y=0;y<H;y++) for(int x=0;x<W;x++){
                int i=x+W*(y+H*z); if(lab[i]==bLab) dst[x+W*y]=(byte)255;
            }
            out.addSlice(new ByteProcessor(W,H,dst));
        }
        ImagePlus res=new ImagePlus("largest3D",out); res.copyScale(mask); return res;
    }
      

/**
 * Returns a new binary 3D image that is the inverse of the input.
 *
 * @param imp The input image.
 * @return The inverse of the input image.
 */
    public static ImagePlus invertBinaryMask(ImagePlus imp) {
        int W = imp.getWidth();
        int H = imp.getHeight();
        int Z = imp.getNSlices();

        ImageStack out = new ImageStack(W, H);

        for (int z = 0; z < Z; z++) {
            ImageProcessor ip = imp.getStack().getProcessor(z + 1).duplicate();

            for (int y = 0; y < H; y++) {
                for (int x = 0; x < W; x++) {
                    float v = ip.getf(x, y);
                    ip.setf(x, y, (v == 0f) ? 1f : 0f);
                }
            }
            out.addSlice(ip);
        }

        ImagePlus inverted = new ImagePlus("inverted_mask", out);
        inverted.copyScale(imp);

        ImagePlus keepBigBlob = largest3DBlob(inverted);
        keepBigBlob.setTitle("largest3D_binary");
    
        return keepBigBlob;
    }


    public static ImagePlus surfaceVoxels3D(ImagePlus mask) {
        int W = mask.getWidth(), H = mask.getHeight(), Z = mask.getStackSize();
        byte[][] sl = new byte[Z][];
        for (int z = 0; z < Z; z++)
            sl[z] = (byte[]) mask.getStack().getProcessor(z + 1).convertToByte(false).getPixels();

        ImageStack out = new ImageStack(W, H);

        for (int z = 0; z < Z; z++) {
            byte[] dst = new byte[W * H];

            for (int y = 0; y < H; y++) {
                for (int x = 0; x < W; x++) {
                    if ((sl[z][x + W * y] & 0xFF) == 0) continue;

                    boolean edge = false;

                    // Check only IN-BOUNDS neighbors; ignore out-of-bounds so borders aren't marked just for touching the frame.
                    if (x > 0     && (sl[z][x - 1 + W * y] & 0xFF) == 0) edge = true;
                    else if (x < W-1 && (sl[z][x + 1 + W * y] & 0xFF) == 0) edge = true;
                    else if (y > 0     && (sl[z][x + W * (y - 1)] & 0xFF) == 0) edge = true;
                    else if (y < H-1   && (sl[z][x + W * (y + 1)] & 0xFF) == 0) edge = true;
                    // else if (z > 0     && (sl[z - 1][x + W * y] & 0xFF) == 0) edge = true;
                    // else if (z < Z-1   && (sl[z + 1][x + W * y] & 0xFF) == 0) edge = true;

                    if (edge) dst[x + W * y] = (byte) 255;
                }
            }

            out.addSlice(new ByteProcessor(W, H, dst));
        }

        ImagePlus surf = new ImagePlus(mask.getShortTitle() + "_surface3D_noBorderJoin", out);
        surf.copyScale(mask);
        return surf;
    }

    /*
     * DEBUG: mask extraction
     *  for(int t = 1; t <5; t++){  
            ImagePlus imgMask = IJ.openImage(Config.mainDir + "Processing/04_Masks/08_MaskWeka/"+specimen.getName()+"_"+timestamps[t-1]+"_mask.tif");
            ImagePlus imgFilled = ImgProUtils.fillTopToFullRowOrBest(imgMask, 0.40); 
            IJ.saveAsTiff(imgFilled, Config.mainDir + "Processing/04_Masks/03_MaskPithFilled/"+specimen.getName()+"_"+timestamps[t-1]+"_mask.tif");
            ImagePlus imgInv = ImgProUtils.invertBinaryMask(imgFilled);
            IJ.saveAsTiff(imgInv, Config.mainDir + "Processing/04_Masks/04_MaskPithFilledInverted/"+specimen.getName()+"_"+timestamps[t-1]+"_mask.tif");
            ImagePlus zeroBelowY = ImgProUtils.zeroBelowY(imgInv, 100);
            IJ.saveAsTiff(zeroBelowY, Config.mainDir + "Processing/04_Masks/06_MaskPithFilledInvertedZeroBelow100/"+specimen.getName()+"_"+timestamps[t-1]+"_mask.tif");  
            ImagePlus imgCrop = ImgProUtils.cropNonTissueAtGivenCoords(imgInv, 100);
            IJ.saveAsTiff(imgCrop, Config.mainDir + "Processing/04_Masks/05_MaskPithFilledCropped/"+specimen.getName()+"_"+timestamps[t-1]+"_mask.tif");
            ImagePlus imgContour = ImgProUtils.surfaceVoxels3D(imgCrop);
            IJ.saveAsTiff(imgContour, Config.mainDir + "Processing/04_Masks/07_MaskSurface3D/"+specimen.getName()+"_"+timestamps[t-1]+"_mask_contour.tif");
        }
     */




    

}
