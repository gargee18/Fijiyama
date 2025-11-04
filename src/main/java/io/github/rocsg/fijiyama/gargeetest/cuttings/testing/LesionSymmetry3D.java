package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
public class LesionSymmetry3D {
    // ==== PUBLIC API =========================================================
    public static class SymmetryResult {
        public final ImagePlus fullFromMinus; // built from Z-
        public final ImagePlus fullFromPlus;  // built from Z+
        public SymmetryResult(ImagePlus a, ImagePlus b){ fullFromMinus=a; fullFromPlus=b; }
    }

    /** Build two double-symmetric volumes from a 3D mask given center (0-based). */
    public static SymmetryResult buildSymmetricEllipsoids(ImagePlus mask, int cx, int cy, int cz) {
        // 1) split
        ImagePlus minus = zMinus(mask, cz);
        ImagePlus plus  = zPlus(mask,  cz);

        // 2) boxes
        int[] boxMinus = bbox3D(minus); // {xmin,xmax,ymin,ymax,zmin,zmax}
        int[] boxPlus  = bbox3D(plus);

        // 3) crops
        ImagePlus cropMinus = crop3D(minus, boxMinus);
        ImagePlus cropPlus  = crop3D(plus,  boxPlus);

        // 4) local centers (global -> local) + CLAMP to crop bounds
        int cyMinusLocal = clamp(cy - boxMinus[2], 0, cropMinus.getHeight()-1);
        int czMinusLocal = clamp(cz - boxMinus[4], 0, cropMinus.getStackSize()-1);
        int cyPlusLocal  = clamp(cy - boxPlus[2],  0, cropPlus.getHeight()-1);
        int czPlusLocal  = clamp(cz - boxPlus[4],  0, cropPlus.getStackSize()-1);

        // 5) mirror: Z (XY plane) then Y (XZ plane)
        ImagePlus fullFromMinus = mirrorY(mirrorZ(cropMinus, czMinusLocal), cyMinusLocal);
        ImagePlus fullFromPlus  = mirrorY(mirrorZ(cropPlus,  czPlusLocal ), cyPlusLocal );

        fullFromMinus.setTitle("Ellipsoid_Zminus_full"); fullFromMinus.setCalibration(mask.getCalibration());
        fullFromPlus.setTitle("Ellipsoid_Zplus_full");   fullFromPlus.setCalibration(mask.getCalibration());
        return new SymmetryResult(fullFromMinus, fullFromPlus);
    }

    // ==== SPLIT / BOX / CROP =================================================
    static ImagePlus zMinus(ImagePlus mask, int cz0){
        int W=mask.getWidth(), H=mask.getHeight(), Z=mask.getStackSize();
        ImageStack out=new ImageStack(W,H);
        for(int z=0; z<Z; z++){
            FloatProcessor fp=new FloatProcessor(W,H);
            if(z<cz0) fp.setPixels(mask.getStack().getProcessor(z+1).toFloat(0,null).getPixels());
            out.addSlice(fp);
        }
        return new ImagePlus(mask.getTitle()+"_Zminus", out);
    }
    static ImagePlus zPlus(ImagePlus mask, int cz0){
        int W=mask.getWidth(), H=mask.getHeight(), Z=mask.getStackSize();
        ImageStack out=new ImageStack(W,H);
        for(int z=0; z<Z; z++){
            FloatProcessor fp=new FloatProcessor(W,H);
            if(z>=cz0) fp.setPixels(mask.getStack().getProcessor(z+1).toFloat(0,null).getPixels());
            out.addSlice(fp);
        }
        return new ImagePlus(mask.getTitle()+"_Zplus", out);
    }

    /** Return {xmin,xmax,ymin,ymax,zmin,zmax}; empty -> {-1,-1,-1,-1,-1,-1}. */
    static int[] bbox3D(ImagePlus img){
        int W=img.getWidth(), H=img.getHeight(), Z=img.getStackSize();
        int xmin=W, xmax=-1, ymin=H, ymax=-1, zmin=Z, zmax=-1;
        for(int z=0; z<Z; z++){
            ImageProcessor ip=img.getStack().getProcessor(z+1); boolean any=false;
            for(int y=0; y<H; y++) for(int x=0; x<W; x++){
                if(ip.getf(x,y)>0f){
                    if(x<xmin) xmin=x; if(x>xmax) xmax=x;
                    if(y<ymin) ymin=y; if(y>ymax) ymax=y;
                    any=true;
                }
            }
            if(any){ if(z<zmin) zmin=z; if(z>zmax) zmax=z; }
        }
        if(xmax<0) return new int[]{-1,-1,-1,-1,-1,-1};
        return new int[]{xmin,xmax,ymin,ymax,zmin,zmax};
    }

    static ImagePlus crop3D(ImagePlus img, int[] b){
        if(b[0]<0) return new ImagePlus(img.getTitle()+"_empty", new ImageStack(1,1));
        int x0=b[0], x1=b[1], y0=b[2], y1=b[3], z0=b[4], z1=b[5];
        int w=x1-x0+1, h=y1-y0+1;
        ImageStack out=new ImageStack(w,h);
        for(int z=z0; z<=z1; z++){
            ImageProcessor ip=img.getStack().getProcessor(z+1);
            ip.setRoi(new java.awt.Rectangle(x0,y0,w,h));
            out.addSlice(ip.crop());
        }
        return new ImagePlus(img.getTitle()+"_crop", out);
    }

    // ==== MIRRORING ==========================================================
    static int clamp(int v, int lo, int hi){ return v<lo?lo : (v>hi?hi : v); }
    static int reflectIdx(int i, int c, int N){
        int r = 2*c - i;
        if(r<0) r=0; else if(r>=N) r=N-1;
        return r;
    }

    /** Mirror across XY plane at z=cz (OR-combine). */
    static ImagePlus mirrorZ(ImagePlus img, int cz){
        int W=img.getWidth(), H=img.getHeight(), Z=img.getStackSize();
        ImageStack out=new ImageStack(W,H);
        for(int z=0; z<Z; z++){
            int zm=reflectIdx(z, cz, Z);
            ImageProcessor a=img.getStack().getProcessor(z+1);
            ImageProcessor m=img.getStack().getProcessor(zm+1);
            FloatProcessor fp=new FloatProcessor(W,H);
            for(int y=0; y<H; y++) for(int x=0; x<W; x++){
                boolean v = a.getf(x,y)>0f || m.getf(x,y)>0f;
                fp.setf(x,y, v?1f:0f);
            }
            out.addSlice(fp);
        }
        ImagePlus r=new ImagePlus(img.getTitle()+"_mirZ", out);
        r.setCalibration(img.getCalibration());
        return r;
    }

    /** Mirror across XZ plane at y=cy (OR-combine). */
    static ImagePlus mirrorY(ImagePlus img, int cy){
        int W=img.getWidth(), H=img.getHeight(), Z=img.getStackSize();
        ImageStack out=new ImageStack(W,H);
        for(int z=0; z<Z; z++){
            ImageProcessor a=img.getStack().getProcessor(z+1);
            FloatProcessor fp=new FloatProcessor(W,H);
            for(int y=0; y<H; y++){
                int ym = reflectIdx(y, cy, H);
                for(int x=0; x<W; x++){
                    boolean v = a.getf(x,y)>0f || a.getf(x,ym)>0f;
                    fp.setf(x,y, v?1f:0f);
                }
            }
            out.addSlice(fp);
        }
        ImagePlus r=new ImagePlus(img.getTitle()+"_mirY", out);
        r.setCalibration(img.getCalibration());
        return r;
    }


    public static void demo(){
        ImagePlus mask = IJ.openImage(Config.mainDir + "/Data/test/B_206_J141_mask_lesion.tif"); // 3D, 0/1 or 0/255
        int cx=160, cy=100, cz=256; // 0-based indices
        SymmetryResult res = buildSymmetricEllipsoids(mask, cx, cy, cz);
        res.fullFromMinus.show();
        res.fullFromPlus.show();
        System.out.println("done!");
    }
    
}
