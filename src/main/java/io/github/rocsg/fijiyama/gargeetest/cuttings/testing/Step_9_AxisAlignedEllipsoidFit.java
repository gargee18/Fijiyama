package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;


// AxisAlignedEllipsoidFit.java
// Fit a 3D, axis-aligned ellipsoid to SURFACE points (x,y,z) with equation:
//   ((x-x0)/rx)^2 + ((y-y0)/ry)^2 + ((z-z0)/rz)^2 = 1
// Unknowns: x0,y0,z0, rx,ry,rz. Axis-aligned (no rotations).
// Method: linear least-squares on the quadric with NO cross terms
//         A x^2 + B y^2 + C z^2 + D x + E y + F z + G = 0,
// then complete the square to recover center and radii.
//
// Also includes helpers to read an ImageJ ImagePlus (0=background, 255=object),
// extract SURFACE points in pixel units, and run the fit.
//

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import org.apache.commons.math3.linear.*;

import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.common.VitimageUtils;

public class Step_9_AxisAlignedEllipsoidFit implements PipelineStep {

    // === Configuration ===
    private static final String DEFAULT_PATH = Config.mainDir + "/Processing/04_Masks/07_MaskSurface3D/B_231_J077_mask_contour.tif";
    private static final int FOREGROUND = 255;   
    private static final double SUBSAMPLE = 1.0; 

    // === Entry Points ===
    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen, true);
    }

    public static void main(String[] args) throws Exception {
        ImageJ ij = new ImageJ();
        String path = (args.length > 0) ? args[0] : DEFAULT_PATH;
        ImagePlus imp = IJ.openImage(path);
        if (imp == null) { System.err.println("Cannot open: " + path); return; }
        Specimen spec = new Specimen("B_231");
        new Step_9_AxisAlignedEllipsoidFit().execute(spec, true);
    }

    public void execute(Specimen specimen, boolean testing) throws Exception {
        ImagePlus mask = IJ.openImage(DEFAULT_PATH);
        int cx = 160, cy = 100, cz = 255;  
        String[] timestamps = Config.timestamps;
        // runEllipsoidBatch(specimen, timestamps, cx, cy, cz);
        buildSymmetricEllipsoids(mask, cx, cy, cz);
    }

    // === Verification / Visualization ===
    public static void verifyResults(Specimen specimen, ImagePlus mask, String[] timestamps, int cx, int cy, int cz) {
        SymmetryResult res = buildSymmetricEllipsoids(mask, cx, cy, cz);

        ImagePlus res_minus = res.fullFromMinus;
        res_minus.setTitle("negative_z"); res_minus.show();
        ImagePlus res_plus = res.fullFromPlus;
        res_plus.setTitle("positive_z"); res_plus.show();
        
        Result rm = fitFromImagePlus(res_minus, FOREGROUND, SUBSAMPLE);
        System.out.printf(java.util.Locale.US, "ok=%s note=%s\n", rm.ok, rm.note);
        System.out.printf(java.util.Locale.US, "center_negative_z_ellipsoid = [%.6f, %.6f, %.6f]\n", rm.center[0], rm.center[1], rm.center[2]);
        System.out.printf(java.util.Locale.US, "radii_negative_z_ellipsoid = [%.6f, %.6f, %.6f]\n", rm.radii[0], rm.radii[1], rm.radii[2]);
        System.out.printf(java.util.Locale.US, "mean residual_negative_z_ellipsoid = %.6g  std = %.6g  (N=%d)\n",
                          rm.meanResidual, rm.stdResidual, rm.nPoints);

        Result rp = fitFromImagePlus(res_plus, FOREGROUND, SUBSAMPLE);
        System.out.printf(java.util.Locale.US, "ok=%s note=%s\n", rp.ok, rp.note);
        System.out.printf(java.util.Locale.US, "center_positive_z_ellipsoid = [%.6f, %.6f, %.6f]\n", rp.center[0], rp.center[1], rp.center[2]);
        System.out.printf(java.util.Locale.US, "radii_positive_z_ellipsoid  = [%.6f, %.6f, %.6f]\n", rp.radii[0], rp.radii[1], rp.radii[2]);
        System.out.printf(java.util.Locale.US, "mean residual_positive_z_ellipsoid = %.6g  std = %.6g  (N=%d)\n",
                          rp.meanResidual, rp.stdResidual, rp.nPoints);

        if (rm.ok) {
                double pw = res_minus.getCalibration().pixelWidth;  if (pw<=0) pw=1;
                double ph = res_minus.getCalibration().pixelHeight; if (ph<=0) ph=1;
                double pd = res_minus.getCalibration().pixelDepth;  if (pd<=0) pd=1;

                double[] center_mm_m = new double[]{ rm.center[0]*pw, rm.center[1]*ph, rm.center[2]*pd };
                double[] radii_mm_m  = new double[]{ rm.radii[0]*pw,  rm.radii[1]*ph,  rm.radii[2]*pd  };

                double shellThickness = Math.max(pw, Math.max(ph, pd));
                ImagePlus contour = generateEllipsoidContour(center_mm_m, radii_mm_m,
                        res_minus.getWidth(), res_minus.getHeight(), res_minus.getNSlices(),
                        pw, ph, pd, shellThickness);
                ImagePlus composite = VitimageUtils.compositeNoAdjustOf(res_minus, contour, "Composite");
                composite.setTitle("Negative Z Ellipsoid fit");
                composite.show();

        }

        if (rp.ok) {
                double pw = res_plus.getCalibration().pixelWidth;  if (pw<=0) pw=1;
                double ph = res_plus.getCalibration().pixelHeight; if (ph<=0) ph=1;
                double pd = res_plus.getCalibration().pixelDepth;  if (pd<=0) pd=1;

                double[] center_mm_p = new double[]{ rp.center[0]*pw, rp.center[1]*ph, rp.center[2]*pd };
                double[] radii_mm_p  = new double[]{ rp.radii[0]*pw,  rp.radii[1]*ph,  rp.radii[2]*pd  };

                double shellThickness = Math.max(pw, Math.max(ph, pd));
                ImagePlus contour = generateEllipsoidContour(center_mm_p, radii_mm_p,
                        res_plus.getWidth(), res_plus.getHeight(), res_plus.getNSlices(),
                        pw, ph, pd, shellThickness);
                ImagePlus composite = VitimageUtils.compositeNoAdjustOf(res_plus, contour, "Composite");
                composite.setTitle("Positive Z Ellipsoid fit");
                composite.show();
        }
    }

    // === Symmetry Construction ===
    public static SymmetryResult buildSymmetricEllipsoids(ImagePlus mask, int cx, int cy, int cz) {
        //Split the mask into two halves
        ImagePlus minus = zMinus(mask, cz); minus.show();  
        ImagePlus plus  = zPlus(mask,  cz); plus.show();

        ImagePlus cropMinus = cropZMinus(minus, cz);  cropMinus.show();
        ImagePlus cropPlus  = cropZPlus(plus, cz);  cropPlus.show();


        ImagePlus fullFromMinus = expandMirrorY(expandMirrorZ(cropMinus, false)); fullFromMinus.show();
        ImagePlus fullFromPlus  = expandMirrorY(expandMirrorZ(cropPlus, true)); fullFromPlus.show();

        VitimageUtils.convertToGray8(fullFromMinus);
        VitimageUtils.convertToGray8(fullFromPlus);

        fullFromMinus.setTitle("Ellipsoid_Zminus_full");
        fullFromPlus.setTitle("Ellipsoid_Zplus_full");
        fullFromMinus.setCalibration(mask.getCalibration());
        fullFromPlus.setCalibration(mask.getCalibration());

        SymmetryResult r = new SymmetryResult();
        r.fullFromMinus = fullFromMinus;
        r.fullFromPlus  = fullFromPlus;
        return r;
    }

    public static class SymmetryResult {
        public ImagePlus fullFromMinus;
        public ImagePlus fullFromPlus;
    }

    public static int reflectIdx(int i, int center, int size) {
        int r = 2 * center - i;
        if (r < 0) r = 0;
        if (r >= size) r = size - 1;
        return r;
    }

    // === Batch Processing and CSV Output ===
    public static void runEllipsoidBatch(Specimen specimen, String[] timestamps, int cx, int cy, int cz) {
        String csv = Config.mainDir + "/Results/04_EllipsoidFitting/test_Ellipsoid_fit_results.csv";
        ensureCsvHeader(csv);

        String[] regions = {"negative_z", "positive_z"};

        for (String region : regions) {
        for (int t = 1; t < 5; t++) {
            String maskPath = Config.mainDir + "Processing/04_Masks/07_MaskSurface3D/"+ specimen.getName() + "_" + timestamps[t-1] + "_mask_contour.tif";
            ImagePlus mask = IJ.openImage(maskPath);
            if (mask == null) continue;

            double pw = mask.getCalibration().pixelWidth;  if (pw <= 0) pw = 1;
            double ph = mask.getCalibration().pixelHeight; if (ph <= 0) ph = 1;
            double pd = mask.getCalibration().pixelDepth;  if (pd <= 0) pd = 1;

            SymmetryResult res = buildSymmetricEllipsoids(mask, cx, cy, cz);
            ImagePlus img = region.equals("negative_z") ? res.fullFromMinus : res.fullFromPlus;

            Result r = fitFromImagePlus(img, 255, SUBSAMPLE);

            double[][] pts = extractSurfacePoints(img, FOREGROUND, SUBSAMPLE,img.getCalibration().pixelWidth, img.getCalibration().pixelHeight, img.getCalibration().pixelDepth);

            // use your existing CSV appender (or the distances version if you added it)
            appendEllipsoidResults(csv, specimen.getName(), timestamps[t-1], region, r, pw, ph, pd,pts );

            System.out.println("Saved " + region + " for " + specimen.getName() + " at " + timestamps[t-1]);
            }
        }
    }

    private static void ensureCsvHeader(String csvPath) {
        java.io.File f = new java.io.File(csvPath);
        if (!f.exists()) {
            f.getParentFile().mkdirs();
            IJ.saveString("specimen,timestamp,region," +
                "center_x_px,center_y_px,center_z_px," +
                "radii_x_px,radii_y_px,radii_z_px," +
                "mean_residual_mh,std_residual_mh," +
                "mean_residual_euc,std_residual_euc," +
                "mean_residual_rmse_mm,mean_residual_hausdorff_mm\n", csvPath);
        }
    }

    public static void saveEllipsoidResultsToCSV(Specimen specimen, Result rm, Result rp, String outputPath) {
        try (FileWriter fw = new FileWriter(outputPath)) {
            fw.write("specimen,region,ok,note,center_x,center_y,center_z," +
                     "radii_x,radii_y,radii_z,mean_residual,std_residual,n_points\n");

            fw.write(String.format(java.util.Locale.US,
                "%s,negative_z,%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6g,%.6g,%d\n",
                specimen.getName(), rm.ok, rm.note, rm.center[0], rm.center[1], rm.center[2],
                rm.radii[0], rm.radii[1], rm.radii[2], rm.meanResidual, rm.stdResidual, rm.nPoints));

            fw.write(String.format(java.util.Locale.US,
                "%s,positive_z,%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6g,%.6g,%d\n",
                specimen.getName(), rp.ok, rp.note, rp.center[0], rp.center[1], rp.center[2],
                rp.radii[0], rp.radii[1], rp.radii[2], rp.meanResidual, rp.stdResidual, rp.nPoints));

            fw.flush();
            System.out.println("Ellipsoid results saved to: " + outputPath);
        } catch (IOException e) { e.printStackTrace(); }
    }

    // === Ellipsoid Fitting ===
    public static class Result {
        public double[] center;
        public double[] radii;
        public double meanResidual;
        public double stdResidual;
        public int nPoints;
        public boolean ok;
        public String note;
    }

    public static Result fitFromPoints(double[][] pts) {
        // System.out.println("Number of points = " + (pts == null ? "null" : pts.length));
        if (pts==null || pts.length < 6) throw new IllegalArgumentException("Need >= 6 points");

        // 1) center + scale (use stddev or range)
        int n = pts.length;
        double mx=0,my=0,mz=0;
        for (double[] p : pts){ mx+=p[0]; my+=p[1]; mz+=p[2]; }
        mx/=n; my/=n; mz/=n;

        double sx=0, sy=0, sz=0;
        for (double[] p : pts){
            sx += (p[0]-mx)*(p[0]-mx);
            sy += (p[1]-my)*(p[1]-my);
            sz += (p[2]-mz)*(p[2]-mz);
        }
        sx = Math.sqrt(sx/n); sy = Math.sqrt(sy/n); sz = Math.sqrt(sz/n);
        if (sx==0) sx=1; if (sy==0) sy=1; if (sz==0) sz=1;  // avoid divide-by-zero

        double[][] u = new double[n][3];
        for (int i=0;i<n;i++){
            u[i][0] = (pts[i][0]-mx)/sx;
            u[i][1] = (pts[i][1]-my)/sy;
            u[i][2] = (pts[i][2]-mz)/sz;
        }

        // 2) build design matrix on normalized coords
        double[][] D = new double[n][7];
        for (int i=0;i<n;i++){
            double x=u[i][0], y=u[i][1], z=u[i][2];
            D[i][0] = x*x;  // A
            D[i][1] = y*y;  // B
            D[i][2] = z*z;  // C
            D[i][3] = x;    // D
            D[i][4] = y;    // E
            D[i][5] = z;    // F
            D[i][6] = 1.0;  // G
        }
        RealVector p = new SingularValueDecomposition(MatrixUtils.createRealMatrix(D))
                        .getV().getColumnVector(6);

        double A = p.getEntry(0), B = p.getEntry(1), C = p.getEntry(2);
        double Dp = p.getEntry(3), Ep = p.getEntry(4), Fp = p.getEntry(5), G = p.getEntry(6);

        // 3) normalize sign so A,B,C > 0 (helps a lot)
        if (A<0 || B<0 || C<0){
            A=-A; B=-B; C=-C; Dp=-Dp; Ep=-Ep; Fp=-Fp; G=-G;
        }

        // 4) complete the square in normalized space
        Result R = new Result();
        if (A==0 || B==0 || C==0){ R.ok=false; R.note="degenerate A/B/C=0"; return R; }

        double x0n = -Dp/(2*A), y0n = -Ep/(2*B), z0n = -Fp/(2*C);
        double Gc  = G - (Dp*Dp)/(4*A) - (Ep*Ep)/(4*B) - (Fp*Fp)/(4*C);
        double S   = -Gc;
        if (S <= 0){ R.ok=false; R.note="S<=0 (not an ellipsoid)"; return R; }

        double rxn2 = S/A, ryn2 = S/B, rzn2 = S/C;
        if (rxn2<=0 || ryn2<=0 || rzn2<=0){ R.ok=false; R.note="non-positive radii^2"; return R; }

        // 5) unscale back to original units
        double x0 = mx + sx*x0n;
        double y0 = my + sy*y0n;
        double z0 = mz + sz*z0n;

        double rx = sx*Math.sqrt(rxn2);
        double ry = sy*Math.sqrt(ryn2);
        double rz = sz*Math.sqrt(rzn2);

        R.ok=true;
        R.center=new double[]{x0,y0,z0};
        R.radii =new double[]{rx,ry,rz};

        // residuals (use original pts)
        double sum=0,sum2=0; int m=0;
        for (double[] q : pts){
            double dx=(q[0]-x0)/rx, dy=(q[1]-y0)/ry, dz=(q[2]-z0)/rz;
            double val = Math.abs(dx*dx + dy*dy + dz*dz - 1.0);
            sum += val; sum2 += val*val; m++;
        }
        double mean=sum/m, var=Math.max(0, sum2/m - mean*mean);
        R.meanResidual=mean; R.stdResidual=Math.sqrt(var); R.nPoints=m;
        R.note="normalized fit";
        return R;
    
    }

    public static Result fitFromImagePlus(ImagePlus imp, int foreground, double subsample) {
        double pw = imp.getCalibration().pixelWidth;  if (pw <= 0) pw = 1;
        double ph = imp.getCalibration().pixelHeight; if (ph <= 0) ph = 1;
        double pd = imp.getCalibration().pixelDepth;  if (pd <= 0) pd = 1;
        double[][] pts = extractSurfacePoints(imp, foreground, subsample, pw, ph, pd);
        if (pts == null || pts.length < 6) return noContour(pts == null ? 0 : pts.length);
        return fitFromPoints(pts);
    }

    private static double[][] extractSurfacePoints(ImagePlus imp, int fg, double keepRate, double pw, double ph, double pd) {
        List<double[]> list = new ArrayList<>();
        ImageStack stack = imp.getStack();
        int w = imp.getWidth(), h = imp.getHeight(), d = imp.getNSlices();
        for (int z = 0; z < d; z++) {
            ImageProcessor ip = stack.getProcessor(z + 1);
            for (int y = 1; y < h - 1; y++)
                for (int x = 1; x < w - 1; x++) {
                    int v = ip.get(x, y) & 0xFF;
                    if (v != fg) continue;
                    list.add(new double[]{x, y, z});
                }
        }
        double[][] pts = new double[list.size()][];
        for (int i = 0; i < list.size(); i++) pts[i] = list.get(i);
        return pts;
    }

    private static Result noContour(int n) {
        Result r = new Result();
        r.ok = false;
        r.note = "no_contour (N=" + n + ")";
        r.center = new double[]{Double.NaN, Double.NaN, Double.NaN};
        r.radii  = new double[]{Double.NaN, Double.NaN, Double.NaN};
        r.meanResidual = Double.NaN;
        r.stdResidual  = Double.NaN;
        r.nPoints = n;
        return r;
    }

    // === Z/Y Split, Crop, Mirror Utilities ===
    public static ImagePlus zMinus(ImagePlus mask, int cz0) {
        final int W = mask.getWidth(), H = mask.getHeight(), Z = mask.getStackSize();
        ImageStack out = new ImageStack(W, H);
        for (int z = 0; z < Z; z++) {
            ImageProcessor ip = mask.getStack().getProcessor(z+1); 
            FloatProcessor fp = (z <= cz0) ? (FloatProcessor) ip.convertToFloat() : new FloatProcessor(W, H); 
            out.addSlice(mask.getStack().getSliceLabel(z + 1), fp);
        }
        ImagePlus res = new ImagePlus(mask.getTitle() + "_Zminus", out); 
        res.copyScale(mask);
        res.setCalibration(mask.getCalibration().copy());
        return res;
    }
    public static ImagePlus zPlus(ImagePlus mask, int cz0) { 
        final int W = mask.getWidth(), H = mask.getHeight(), Z = mask.getStackSize();
        cz0 = Math.max(0, Math.min(cz0, Z - 1)); 
        ImageStack out = new ImageStack(W, H);
        for (int z = 0; z < Z; z++) {
            ImageProcessor ip = mask.getStack().getProcessor(z + 1);
            FloatProcessor fp = (z >= cz0) ? (FloatProcessor) ip.convertToFloat() : new FloatProcessor(W, H); 
            out.addSlice(mask.getStack().getSliceLabel(z + 1), fp);
        }
        ImagePlus res = new ImagePlus(mask.getTitle() + "_Zplus", out);
        res.copyScale(mask);
        res.setCalibration(mask.getCalibration().copy());
        return res;
    }
    public static ImagePlus cropZMinus(ImagePlus img, int cz) {
        int W = img.getWidth(), H = img.getHeight(), Z = img.getStackSize();
        cz = Math.max(0, Math.min(cz, Z - 1));
        ImageStack out = new ImageStack(W, H);

        for (int z = 0; z <= cz; z++) { // include cz
            out.addSlice(img.getStack().getProcessor(z + 1).duplicate());
        }

        ImagePlus res = new ImagePlus(img.getTitle() + "_Zminus", out);
        res.setCalibration(img.getCalibration());
        return res;
    }
    public static ImagePlus cropZPlus(ImagePlus img, int cz) {
        int W = img.getWidth(), H = img.getHeight(), Z = img.getStackSize();
        cz = Math.max(0, Math.min(cz, Z - 1)); // clamp cz to valid range
        ImageStack out = new ImageStack(W, H);

        for (int z = cz; z < Z-1; z++) { // start at cz, include it
            out.addSlice(img.getStack().getProcessor(z + 1).duplicate());
        }

        ImagePlus res = new ImagePlus(img.getTitle() + "_Zplus", out);
        res.setCalibration(img.getCalibration());
        return res;
    }
    public static ImagePlus expandMirrorZ(ImagePlus halfZ, boolean midAtStart) {
          int W = halfZ.getWidth(), H = halfZ.getHeight(), Z = halfZ.getStackSize();
        ImageStack out = new ImageStack(W, H);

        if (midAtStart) {
            // mirror about the *first* slice
            for (int z = Z - 1; z >= 1; z--) {
                float[] src = (float[]) halfZ.getStack().getProcessor(z + 1).convertToFloat().getPixels();
                out.addSlice(new FloatProcessor(W, H, src.clone()));
            }
            for (int z = 0; z < Z; z++) {
                float[] src = (float[]) halfZ.getStack().getProcessor(z + 1).convertToFloat().getPixels();
                out.addSlice(new FloatProcessor(W, H, src.clone()));
            }
        } else {
            // mirror about the *last* slice (your current behavior)
            for (int z = 0; z < Z; z++) {
                float[] src = (float[]) halfZ.getStack().getProcessor(z + 1).convertToFloat().getPixels();
                out.addSlice(new FloatProcessor(W, H, src.clone()));
            }
            for (int z = Z - 2; z >= 0; z--) {
                float[] src = (float[]) halfZ.getStack().getProcessor(z + 1).convertToFloat().getPixels();
                out.addSlice(new FloatProcessor(W, H, src.clone()));
            }
        }

        ImagePlus res = new ImagePlus(halfZ.getTitle() + "_expandZ", out);
        res.setCalibration(halfZ.getCalibration());
        return res;
    }
    public static ImagePlus expandMirrorY(ImagePlus halfY) {
           int W = halfY.getWidth(), H = halfY.getHeight(), Z = halfY.getStackSize();
        int Hnew = 2*H - 1;
        ImageStack out = new ImageStack(W, Hnew);

        for (int z=0; z<Z; z++) {
            float[] src = (float[]) halfY.getStack().getProcessor(z+1).convertToFloat().getPixels();
            float[] dst = new float[W * Hnew];

            // copy original rows 0..H-1
            System.arraycopy(src, 0, dst, 0, src.length);
            // append mirrored rows H..2H-2 from rows H-2..0
            for (int y=0; y<H-1; y++) {
                int ySrc = (H-2) - y;
                int yDst = H + y;
                System.arraycopy(src, ySrc*W, dst, yDst*W, W);
            }
            out.addSlice(new FloatProcessor(W, Hnew, dst));
        }
        ImagePlus res = new ImagePlus(halfY.getTitle()+"_expandY", out);
        res.setCalibration(halfY.getCalibration());
        return res;
    }
    public static ImagePlus mirrorZ(ImagePlus img, int cz) {
        int W = img.getWidth();
        int H = img.getHeight();
        int Z = img.getStack().getSize();

        ImageStack out = new ImageStack(W, H);

        for (int z = 0; z < Z; z++) {
            int zm = reflectIdx(z, cz, Z);
            ImageProcessor ip  = img.getStack().getProcessor(z+1);
            ImageProcessor ipm = img.getStack().getProcessor(zm+1);

            FloatProcessor fp = new FloatProcessor(W, H);
            for (int y = 0; y < H; y++) {
                for (int x = 0; x < W; x++) {
                    boolean v = ip.getf(x,y)  > 0f;
                    boolean m = ipm.getf(x,y) > 0f;
                    fp.setf(x,y, (v || m) ? 1f : 0f);
                }
            }
            out.addSlice(fp);
        }

    ImagePlus result = new ImagePlus(img.getTitle()+"_mirZ", out);
        result.setCalibration(img.getCalibration());
        return result;
    }
    
    public static ImagePlus mirrorY(ImagePlus img, int cy) {
        int W = img.getWidth(),H = img.getHeight(), Z = img.getStack().getSize();
        ImageStack out = new ImageStack(W, H);
        for (int z = 0; z < Z; z++) {
            ImageProcessor ip = img.getStack().getProcessor(z+1);
            FloatProcessor fp = new FloatProcessor(W, H);

            for (int y = 0; y < H; y++) {
                int ym = reflectIdx(y, cy, H);
                for (int x = 0; x < W; x++) {
                    boolean v = ip.getf(x,y)  > 0f;
                    boolean m = ip.getf(x,ym) > 0f;
                    fp.setf(x,y, (v || m) ? 1f : 0f);
                }
            }
            out.addSlice(fp);
        }
        ImagePlus result = new ImagePlus(img.getTitle()+"_mirY", out);
        result.setCalibration(img.getCalibration());
        return result;
    }

    // === Bounding Box & Cropping ===
    public static int[] bbox3D(ImagePlus img) {
        int W=img.getWidth(), H=img.getHeight(), Z=img.getStackSize();
        int xmin=W, xmax=-1, ymin=H, ymax=-1, zmin=Z, zmax=-1;
        for (int z=0; z<Z; z++) {
            ImageProcessor ip = img.getStack().getProcessor(z+1);
            boolean any=false;
            for (int y=0; y<H; y++) for (int x=0; x<W; x++) {
                if (ip.getf(x,y) > 0f) {
                    if (x<xmin) xmin=x; if (x>xmax) xmax=x;
                    if (y<ymin) ymin=y; if (y>ymax) ymax=y;
                    any = true;
                }
            }
            if (any) { if (z<zmin) zmin=z; if (z>zmax) zmax=z; }
        }
        return new int[]{xmin,xmax,ymin,ymax,zmin,zmax};
    }
    public static ImagePlus crop3D(ImagePlus img, int[] b) {
        if (b[1]<b[0] || b[3]<b[2] || b[5]<b[4]) return new ImagePlus(img.getTitle()+"_empty", new ImageStack(1,1));
        int x0=b[0], y0=b[2], w=b[1]-b[0]+1, h=b[3]-b[2]+1;
        ImageStack out = new ImageStack(w,h);
        for (int z=b[4]; z<=b[5]; z++) {
            ImageProcessor ip = img.getStack().getProcessor(z+1);
            ip.setRoi(new java.awt.Rectangle(x0,y0,w,h));
            out.addSlice(ip.crop());
        }
        ImagePlus res = new ImagePlus(img.getTitle()+"_crop", out);
        res.setCalibration(img.getCalibration());    // <-- keep mm units
        return res;
    }
    public static void doBoxes(ImagePlus mask, int cx0, int cy0, int cz0) {
        ImagePlus minus = zMinus(mask, cz0);
        ImagePlus plus  = zPlus(mask,  cz0);
        // 2) boxes
        int[] boxMinus = bbox3D(minus);
        int[] boxPlus  = bbox3D(plus);

        System.out.println(String.format("Z- box: x[%d..%d] y[%d..%d] z[%d..%d]", boxMinus[0],boxMinus[1],boxMinus[2],boxMinus[3],boxMinus[4],boxMinus[5]));
        System.out.println(String.format("Z+ box: x[%d..%d] y[%d..%d] z[%d..%d]", boxPlus[0], boxPlus[1], boxPlus[2], boxPlus[3], boxPlus[4], boxPlus[5]));

        // 3) crop
        ImagePlus cropMinus = cropZMinus(mask, cz0);
        VitimageUtils.convertToGray8(cropMinus);
        cropMinus.show();
        ImagePlus cropPlus  = cropZPlus(mask, cz0);
        VitimageUtils.convertToGray8(cropPlus);
        cropPlus.show();
        cropMinus.setTitle("Lesion_Zminus_crop"); cropMinus.show();
        cropPlus.setTitle("Lesion_Zplus_crop");   cropPlus.show();

        // (for later mirroring, convert global centers to local-in-crop)
        int cyMinusLocal = cy0 - boxMinus[2];
        int czMinusLocal = cz0 - boxMinus[4];
        int cyPlusLocal  = cy0 - boxPlus[2];
        int czPlusLocal  = cz0 - boxPlus[4];
        System.out.println("locals (for mirroring later)  Z-: cy="+cyMinusLocal+" cz="+czMinusLocal+"   Z+: cy="+cyPlusLocal+" cz="+czPlusLocal);
    }

    // === Ellipsoid Visualization ===
    public static ImagePlus generateFilledEllipsoid(double[] center, double[] radii,
                                                     int width, int height, int depth,
                                                     double pixelWidth, double pixelHeight, double pixelDepth) {
                                                     ImageStack stack = ImageStack.create(width, height, depth, 8);
        double x0 = center[0], y0 = center[1], z0 = center[2];
        double rx = radii[0], ry = radii[1], rz = radii[2];
        double rx2 = rx*rx, ry2 = ry*ry, rz2 = rz*rz;
        
        for (int z = 0; z < depth; z++) {
            ImageProcessor ip = stack.getProcessor(z + 1);
            double Z = z  * pixelDepth;
            // double Z = (z + 0.5) * pixelDepth;
            double dz = Z - z0;
            
            for (int y = 0; y < height; y++) {
                double Y = y * pixelHeight;
                // double Y = (y + 0.5) * pixelHeight;
                double dy = Y - y0;
                
                for (int x = 0; x < width; x++) {
                    double X = x * pixelWidth;
                    // double X = (x + 0.5) * pixelWidth;
                    double dx = X - x0;
                    
                    // Check if point is inside ellipsoid
                    double f = (dx*dx)/rx2 + (dy*dy)/ry2 + (dz*dz)/rz2;
                    if (f <= 1.0) {
                        ip.set(x, y, 255);
                    }
                }
            }
        }
        
        ImagePlus imp = new ImagePlus("Ellipsoid_Filled", stack);
        imp.getCalibration().pixelWidth = pixelWidth;
        imp.getCalibration().pixelHeight = pixelHeight;
        imp.getCalibration().pixelDepth = pixelDepth;
        return imp;   
    }
    public static ImagePlus generateEllipsoidContour(double[] center, double[] radii,
                                                      int width, int height, int depth,
                                                      double pixelWidth, double pixelHeight, double pixelDepth,
                                                      double shellThickness) {
                                                    ImageStack stack = ImageStack.create(width, height, depth, 8);
        double x0 = center[0], y0 = center[1], z0 = center[2];
        double rx = radii[0], ry = radii[1], rz = radii[2];
        double rx2 = rx*rx, ry2 = ry*ry, rz2 = rz*rz;
        
        for (int z = 0; z < depth; z++) {
            ImageProcessor ip = stack.getProcessor(z + 1);
            double Z = z  * pixelDepth;
            double dz = Z - z0;
            
            for (int y = 0; y < height; y++) {
                double Y = y  * pixelHeight;
                double dy = Y - y0;
                
                for (int x = 0; x < width; x++) {
                    double X = x  * pixelWidth;
                    double dx = X - x0;
                    
                    // Compute normalized ellipsoid function value
                    double f = (dx*dx)/rx2 + (dy*dy)/ry2 + (dz*dz)/rz2;
                    
                    // Compute the gradient magnitude at this point to get proper distance
                    // gradient of f = [2*dx/rx^2, 2*dy/ry^2, 2*dz/rz^2]
                    double gradX = 2*dx/rx2;
                    double gradY = 2*dy/ry2;
                    double gradZ = 2*dz/rz2;
                    double gradNorm = Math.sqrt(gradX*gradX + gradY*gradY + gradZ*gradZ);
                    
                    // Approximate geometric distance to surface: |f - 1| / ||grad(f)||
                    double dist = Math.abs(f - 1.0) / (gradNorm + 1e-12);
                    
                    // If within shell thickness of the surface, mark as foreground
                    if (dist <= shellThickness) {
                        ip.set(x, y, 255);
                    }
                }
            }
        }
        
        ImagePlus imp = new ImagePlus("Ellipsoid_Contour", stack);
        imp.getCalibration().pixelWidth = pixelWidth;
        imp.getCalibration().pixelHeight = pixelHeight;
        imp.getCalibration().pixelDepth = pixelDepth;
        return imp;
    }

    // === Internal Helpers ===
    private static void appendEllipsoidResults(String csvPath, String specimen, String ts, String region, Result r, double pw, double ph, double pd, double[][]pts) {
        double cx_px = (r.center == null) ? Double.NaN : r.center[0];
        double cy_px = (r.center == null) ? Double.NaN : r.center[1];
        double cz_px = (r.center == null) ? Double.NaN : r.center[2];

        double rx_px = (r.radii  == null) ? Double.NaN : r.radii[0];
        double ry_px = (r.radii  == null) ? Double.NaN : r.radii[1];
        double rz_px = (r.radii  == null) ? Double.NaN : r.radii[2];

        // mm = px * (mm per px)
        double cx_mm = cx_px * pw, cy_mm = cy_px * ph, cz_mm = cz_px * pd;
        double rx_mm = rx_px * pw, ry_mm = ry_px * ph, rz_mm = rz_px * pd;

        // residuals are in px
        double mean_px = r.meanResidual;
        double std_px  = r.stdResidual;

        // geometric distances in mm
        double[] mm = residualMetricsForEllipsoidFit(pts, r, pw, ph, pd);
        double mean_mm = mm[0], std_mm = mm[1], rmse_mm = mm[2], haus_mm = mm[3];

        String line = String.format(java.util.Locale.US,
            "%s,%s,%s," +                 // specimen, timestamp, region
            "%.6f,%.6f,%.6f," +           // center px
            "%.6f,%.6f,%.6f," +           // radii px
            "%.6g,%.6g," +                // mean/std residual px (unitless)
            "%.6g,%.6g," +                // mean/std residual mm (geom)
            "%.6g,%.6g",                  // rmse_mm, hausdorff_mm
            specimen, ts, region,
            nz(cx_px), nz(cy_px), nz(cz_px),
            nz(rx_px), nz(ry_px), nz(rz_px),
            nz(mean_px), nz(std_px),
            nz(mean_mm), nz(std_mm),
            nz(rmse_mm), nz(haus_mm)
        );
    IJ.append(line, csvPath);
    }
    private static double nz(double v) { return (Double.isNaN(v) || Double.isInfinite(v)) ? 0.0 : v; }

    private static double[] residualMetricsForEllipsoidFit(double[][] pts, Result r, double pw, double ph, double pd) {

        if (pts == null || pts.length == 0 || r == null || !r.ok) {
            return new double[]{Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        }

        // center/radii to mm
        final double x0 = r.center[0] * pw, y0 = r.center[1] * ph, z0 = r.center[2] * pd;
        final double rx = r.radii[0]  * pw, ry = r.radii[1]  * ph, rz = r.radii[2]  * pd;
        final double rx2 = rx*rx, ry2 = ry*ry, rz2 = rz*rz;

        double sum = 0, sum2 = 0, sumsq = 0, haus = 0;
        int m = 0;

        for (double[] p : pts) {
            // voxel -> mm
            double X = p[0] * pw, Y = p[1] * ph, Z = p[2] * pd;
            double dx = X - x0, dy = Y - y0, dz = Z - z0;

            double f = (dx*dx)/rx2 + (dy*dy)/ry2 + (dz*dz)/rz2;  // implicit ellipsoid
            double gx = 2*dx/rx2, gy = 2*dy/ry2, gz = 2*dz/rz2;  // gradient components
            double g  = Math.sqrt(gx*gx + gy*gy + gz*gz) + 1e-12;
            double dist = Math.abs(f - 1.0) / g;  // ≈ Euclidean distance to surface (mm)
            sum   += dist;
            sum2  += dist * dist;
            sumsq += dist * dist;                 // for RMSE
            if (dist > haus) haus = dist;          // Hausdorff
            m++;
        }
        double mean = sum / m;
        double var  = Math.max(0, sum2 / m - mean * mean);
        double std  = Math.sqrt(var);
        double rmse = Math.sqrt(sumsq / m);

        return new double[]{ mean, std, rmse, haus};
    }

}
