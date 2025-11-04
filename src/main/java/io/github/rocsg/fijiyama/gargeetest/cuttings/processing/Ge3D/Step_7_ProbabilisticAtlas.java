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
        Specimen spec = new Specimen("B_206");
        
        ImageJ ij = new ImageJ();
        new Step_7_ProbabilisticAtlas().execute(spec,true);
    }

    public void execute(Specimen specimen, boolean testing) throws Exception {
        String[] timestamps = Config.timestamps;
        
        for(int t = 1; t <5; t++){
                ImagePlus mask = IJ.openImage(Config.mainDir + "Processing/04_Masks/05_MaskPithFilledCropped/"+specimen.getName()+"_"+timestamps[t-1]+"_mask.tif");
                ImagePlus contour = surfaceVoxels3D(mask);
                // contour.show();
                IJ.saveAsTiff(contour, Config.mainDir + "Processing/04_Masks/07_MaskSurface3D/"+specimen.getName()+"_"+timestamps[t-1]+"_mask_contour.tif");
            }

        // for(int t = 1; t <5; t++){
        // //  
        //     ImagePlus img = IJ.openImage(Config.mainDir + "Processing/04_Masks/04_MaskPithFilledInverted/"+specimen.getName()+"_"+timestamps[t-1]+"_mask.tif");
        //     ImagePlus m =zeroBelowY(img,100);
        //     // m.show();
        //     IJ.saveAsTiff(m, Config.mainDir + "Processing/04_Masks/06_MaskPithFilledInvertedZeroBelow100/"+specimen.getName()+"_"+timestamps[t-1]+"_mask.tif");
            
        // } 
      
       
      
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
    
    // public static void binarizeAndFillTopAuto(ImagePlus imp, double coverage) {
    //     int w = imp.getWidth(), h = imp.getHeight(), d = imp.getNSlices();
    //     coverage = Math.max(0.0, Math.min(1.0, coverage));
    //     int needed = Math.max(1, (int) Math.ceil(w * coverage)); // columns needed to count a row as foreground

    //     for (int z = 1; z <= d; z++) {
    //         ImageProcessor ip = imp.getStack().getProcessor(z).convertToByte(false); // work in 8-bit
    //         byte[] px = (byte[]) ip.getPixels();

    //         // 1) binarize to 0/1
    //         for (int i = 0; i < px.length; i++) px[i] = (byte) (((px[i] & 0xFF) > 0) ? 1 : 0);

    //         // 2) scan from top to find first "real" foreground row (robust to noise)
    //         int stopY = -1;
    //         for (int y = 0; y < h; y++) {
    //             int count = 0;
    //             int off = y * w;
    //             for (int x = 0; x < w; x++) {
    //                 if (px[off + x] != 0) {
    //                     if (++count >= needed) { stopY = y; break; }
    //                 }
    //             }
    //             if (stopY >= 0) break;
    //         }
    //         System.out.println("Slice " + z + " needed=" + needed + " stopY=" + stopY + " (h=" + h + ")");

    //         // 3) set all rows above (and including) stopY to 1
    //         if (stopY >= 0) {
    //             for (int y = 0; y <= stopY; y++) {
    //                 int off = y * w;
    //                 for (int x = 0; x < w; x++) px[off + x] = 1;
    //             }
    //         }

    //         // write back
    //         ip.setPixels(px);
    //         imp.getStack().setProcessor(ip, z);
    //     }
    //     imp.updateAndDraw();

    // }

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

    public static ImagePlus cropNonTissueAtGivenCoords(ImagePlus imp, int y){
        // imp.show();
        int w = imp.getWidth(), h = Math.max(0, Math.min(y, imp.getHeight()));

        ImageStack croppedStack = new ImageStack(w,h);

        for (int z = 0; z < imp.getNSlices(); z++) {
            ImageProcessor ip = imp.getStack().getProcessor(z+1);
            ip.setRoi(new Roi(0,0,w,h));
            ImageProcessor cropped = ip.crop();
            croppedStack.addSlice(cropped);

               // Invert binary: 0 ↔ 1 
        //     if (cropped.getBitDepth() == 8) {
        //         for (int yy = 0; yy < h; yy++) {
        //             for (int xx = 0; xx < w; xx++) {
        //                 float v = cropped.getf(xx, yy);
        //                 cropped.setf(xx, yy, (v == 0f) ? 1f : 0f); // for binary {0,1}
        //                 // If your mask uses {0,255}, use:
        //                 // cropped.setf(xx, yy, (v == 0f) ? 255f : 0f);
        //             }
        //         }
        //     } else {
        //         // Fallback: generic invert for float processors
        //         for (int yy = 0; yy < h; yy++) {
        //             for (int xx = 0; xx < w; xx++) {
        //                 float v = cropped.getf(xx, yy);
        //                 cropped.setf(xx, yy, (v == 0f) ? 1f : 0f);
        //             }
        //         }
        //     }

            }
        
        ImagePlus croppedMask = new ImagePlus("croppedMask", croppedStack);
        croppedMask.copyScale(imp);
        // croppedMask.show();

        // ImagePlus keepBigBlob = largest3DBlob(croppedMask);
        // keepBigBlob.setTitle("largest3D_binary");
        // keepBigBlob.show();

        return croppedMask;
    }


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

        // ImagePlus keepBigBlob = largest3DBlob(res);
        // keepBigBlob.setTitle("largest3D_binary");
        // // keepBigBlob.show();

        return res;


    }
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
    
    public static ImagePlus runFillTopToFullRowOrBest(Specimen specimen, int frame) {
        String[] timestamps = Config.timestamps;
        String specimenName = specimen.getName();
        System.out.println(Config.mainDir+ "Processing/04_Masks/02_MasksforPA/" + specimenName + "_" + timestamps[frame-1] + "_mask.tif");
        
        ImagePlus img =  IJ.openImage(Config.mainDir+ "Processing/04_Masks/02_MaskforPA/" + specimenName + "_" + timestamps[frame-1] + "_mask.tif");
        fillTopToFullRowOrBest(img, 0.40);
        // IJ.saveAsTiff(img, Config.mainDir+ "Results/04_Masks/03_MaskPithFilled/" + specimenName + "_" + timestamps[frame-1] + "_mask.tif");     
        // System.out.println("done!");
        return img;
        
    }
   

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
        return inverted;
    }

        // --- 2D inner contour per slice (8-connectivity), binary {0,255} ---
    public static ImagePlus innerContour2D(ImagePlus mask) {
        int W=mask.getWidth(), H=mask.getHeight(), Z=mask.getStackSize();
        ImageStack out = new ImageStack(W,H);

        for (int s=1; s<=Z; s++) {
            // Work on 8-bit binary slice
            ByteProcessor src = (ByteProcessor) mask.getStack().getProcessor(s).convertToByte(false);
            ByteProcessor er  = (ByteProcessor) src.duplicate();
            new BinaryProcessor(er).erode();               // 3x3, 8-connected

            byte[] a = (byte[]) src.getPixels();
            byte[] b = (byte[]) er.getPixels();
            byte[] c = new byte[W*H];

            for (int i=0; i<c.length; i++) c[i] = (byte)( (a[i]&0xFF)>0 && (b[i]&0xFF)==0 ? 255 : 0 );
            out.addSlice(new ByteProcessor(W,H,c));
        }
        ImagePlus contour = new ImagePlus(mask.getShortTitle()+"_contour2D", out);
        contour.copyScale(mask);
        return contour;
    }

    // --- 3D surface voxels (6-neighborhood), binary {0,255} ---
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

}
 //ImagePlus img = IJ.openImage(Config.mainDir + "Processing/04_Masks/04_MaskPithFilledInverted/"+specimen.getName()+"_"+timestamps[t-1]+"_mask.tif");
