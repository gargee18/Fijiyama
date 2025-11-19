package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.gargeetest.cuttings.helpers.ImgProUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;


public class Step_9_BoundingBox implements PipelineStep {
     @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen, true);
    }

    public static void main (String[] args) throws Exception {
        Specimen spec = new Specimen("B_201");
        
        ImageJ ij = new ImageJ();
        new  Step_9_BoundingBox().execute(spec,true);
    }

    public void execute (Specimen specimen, boolean testing) throws Exception {
        String[] timestamps = Config.timestamps;
        int cx = 150, cy = 100, cz = 255;  // Center coordinates for symmetry
        String csv = Config.mainDir + "/Results/05_BoundingBox/BoundingBox_fit_results.csv";
        ensureCsvHeader(csv);
        
        runBoundingBoxBatch(specimen, timestamps, cx, cy, cz, csv);
    }

     // === Main Batch Processing ===
    
   public static void runBoundingBoxBatch(Specimen specimen, String[] timestamps, int cx, int cy, int cz, String csvPath) {
        String[] regions = {"negative_z", "positive_z"};
        
        for (String region : regions) {
            // Process ALL timestamps, not just 1 to 4
            for (int t = 0; t < timestamps.length; t++) {
                String maskPath = Config.mainDir + "Processing/04_Masks/07_MaskSurface3D/" + 
                    specimen.getName() + "_" + timestamps[t] + "_mask_contour.tif";
                
                ImagePlus mask = IJ.openImage(maskPath);
                if (mask == null) {
                    System.out.println(String.format("%s %s %s: Missing mask - writing zeros",
                        specimen.getName(), timestamps[t], region));
                    appendEmptyBBoxResults(csvPath, specimen.getName(), timestamps[t], region);
                    continue;
                }
                
                Step_8_AxisAlignedEllipsoidFit.SymmetryResult res = Step_8_AxisAlignedEllipsoidFit.buildSymmetricEllipsoids(mask, cx, cy, cz);
                ImagePlus img = region.equals("negative_z") ? res.fullFromMinus : res.fullFromPlus;
                
                BoundingBoxResult bbox = computeBoundingBox(img);
                
                if (!bbox.valid) {
                    System.out.println(String.format("%s %s %s: Empty image - writing zeros",
                        specimen.getName(), timestamps[t], region));
                    appendEmptyBBoxResults(csvPath, specimen.getName(), timestamps[t], region);
                    continue;
                }
                
                printBoundingBoxResults(specimen.getName(), timestamps[t], region, bbox);
                appendBBoxResults(csvPath, specimen.getName(), timestamps[t], region, bbox);
            }
        }
    }


    // === Bounding Box Computation ===
    
    public static class BoundingBoxResult {
        public int[] bounds;      // [xmin, xmax, ymin, ymax, zmin, zmax]
        public int width;
        public int height;
        public int depth;
        public int volume;        // in voxels
        public boolean valid;
        
        public BoundingBoxResult(int[] bounds) {
            this.bounds = bounds;
            this.valid = (bounds[1] >= bounds[0] && bounds[3] >= bounds[2] && bounds[5] >= bounds[4]);
            
            if (valid) {
                this.width = bounds[1] - bounds[0] + 1;
                this.height = bounds[3] - bounds[2] + 1;
                this.depth = bounds[5] - bounds[4] + 1;
                this.volume = width * height * depth;
            } else {
                this.width = 0;
                this.height = 0;
                this.depth = 0;
                this.volume = 0;
            }
        }
        
        public double getHalfWidth() { return width / 2.0; }
        public double getHalfHeight() { return height / 2.0; }
        public double getHalfDepth() { return depth / 2.0; }
    }
    private static void appendEmptyBBoxResults(String csv, String specimen, String timestamp, String region) {
        try (FileWriter fw = new FileWriter(csv, true);
            PrintWriter pw = new PrintWriter(fw)) {
            // All zeros, including ratio columns
            pw.printf(java.util.Locale.US,
                    "%s,%s,%s," +
                    "0,0,0,0,0,0," +      // bounds
                    "0,0,0," +            // width/height/depth
                    "0.0,0.0,0.0," +      // half-width/half-height/half-depth
                    "0," +                // volume
                    "0.0,0.0,0.0%n",      // ratios
                    specimen, timestamp, region);
            } catch (Exception e) {
            e.printStackTrace();
        }
    }
     // === Output Formatting ===
    
    private static void printBoundingBoxResults(String specimen, String timestamp, String region, BoundingBoxResult bbox) {
            System.out.println(String.format("%s %s %s:", specimen, timestamp, region));
            System.out.println(String.format("  BBox bounds: x[%d-%d] y[%d-%d] z[%d-%d]",
                bbox.bounds[0], bbox.bounds[1], bbox.bounds[2], 
                bbox.bounds[3], bbox.bounds[4], bbox.bounds[5]));
            System.out.println(String.format("  BBox full lengths: %d × %d × %d (volume=%d voxels)",
                bbox.width, bbox.height, bbox.depth, bbox.volume));
            System.out.println(String.format("  BBox half-lengths: %.1f, %.1f, %.1f",
                bbox.getHalfWidth(), bbox.getHalfHeight(), bbox.getHalfDepth()));
            System.out.println();
        }

     public static BoundingBoxResult computeBoundingBox(ImagePlus img) {
        int[] bounds = ImgProUtils.bbox3D(img);
        return new BoundingBoxResult(bounds);
    }

      // === CSV Management ===
    
    private static void ensureCsvHeader(String csv) {
        File f = new File(csv);
        if (f.exists()) return;
        
        f.getParentFile().mkdirs();
        
        try (FileWriter fw = new FileWriter(csv);
             PrintWriter pw = new PrintWriter(fw)) {
            pw.println("specimen,timestamp,region," +
                      "xmin,xmax,ymin,ymax,zmin,zmax," +
                      "width_px,height_px,depth_px," +
                      "half_width_px,half_height_px,half_depth_px," +
                      "volume_voxels,"+
                      "ratio_width_height,ratio_width_depth,ratio_height_depth");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private static void appendBBoxResults(String csv, String specimen, String timestamp,String region, BoundingBoxResult bbox) {
        try (FileWriter fw = new FileWriter(csv, true);
            PrintWriter pw = new PrintWriter(fw)) {

            // Ratios (protect against division by zero)
            double ratioWH = (bbox.height == 0) ? 0.0 : (double) bbox.width  / (double) bbox.height;
            double ratioWD = (bbox.depth  == 0) ? 0.0 : (double) bbox.width  / (double) bbox.depth;
            double ratioHD = (bbox.depth  == 0) ? 0.0 : (double) bbox.height / (double) bbox.depth;

            pw.printf(java.util.Locale.US,
                    "%s,%s,%s," +
                    "%d,%d,%d,%d,%d,%d," +     // bounds
                    "%d,%d,%d," +             // width, height, depth
                    "%.3f,%.3f,%.3f," +       // half-width, half-height, half-depth
                    "%d," +                   // volume
                    "%.6f,%.6f,%.6f%n",       // ratios
                    specimen, timestamp, region,
                    bbox.bounds[0], bbox.bounds[1], bbox.bounds[2],
                    bbox.bounds[3], bbox.bounds[4], bbox.bounds[5],
                    bbox.width, bbox.height, bbox.depth,
                    bbox.getHalfWidth(), bbox.getHalfHeight(), bbox.getHalfDepth(),
                    bbox.volume,
                    ratioWH, ratioWD, ratioHD);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}