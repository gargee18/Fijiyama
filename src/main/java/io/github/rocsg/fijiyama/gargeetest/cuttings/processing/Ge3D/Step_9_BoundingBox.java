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
import java.util.Locale;

public class Step_9_BoundingBox implements PipelineStep {

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen, true);
    }

    public static void main(String[] args) throws Exception {
        Specimen spec = new Specimen("B_201");
        ImageJ ij = new ImageJ();
        new Step_9_BoundingBox().execute(spec, true);
    }

    public void execute(Specimen specimen, boolean testing) throws Exception {
        String[] timestamps = Config.timestamps;
        int cx = 150, cy = 100, cz = 255;  // Center coordinates for symmetry
        
        // // Measure surface masks (whole lesion)
        // String csvSurface = Config.mainDir + "/Results/05_BoundingBox/BoundingBox_surface_masks.csv";
        // ensureCsvHeader(csvSurface, false);  // no region column
        // runBoundingBoxBatch(specimen, timestamps, csvSurface, false, cx, cy, cz);
        
        // Measure symmetric regions (z- and z+)
        String csvRegions = Config.mainDir + "/Results/05_BoundingBox/BoundingBox_symmetric_regions.csv";
        // ensureCsvHeader(csvRegions, true);   // with region column
        runBoundingBoxBatch(specimen, timestamps, csvRegions, true, cx, cy, cz);
    }

    // === Main Batch Processing ===

    public static void runBoundingBoxBatch(Specimen specimen, String[] timestamps, 
                                          String csvPath, boolean withRegions,
                                          int cx, int cy, int cz) {
        if (withRegions) {
            // Process z- and z+ regions separately
            String[] regions = {"negative_z", "positive_z"};
            
            for (String region : regions) {
                for (int t = 0; t < timestamps.length; t++) {
                    String ts = timestamps[t];
                    String maskPath = Config.mainDir + "Processing/04_Masks/07_MaskSurface3D/" +  specimen.getName() + "_" + ts + "_mask_contour.tif";

                    ImagePlus mask = IJ.openImage(maskPath);
                    if (mask == null) {
                        System.out.printf("%s %s %s: Missing mask - writing zeros%n", 
                            specimen.getName(), ts, region);
                        appendEmptyBBoxResults(csvPath, specimen.getName(), ts, region);
                        continue;
                    }

                    // Build symmetric mirrored region
                    
                    ImagePlus zMinus = ImgProUtils.zMinus(mask, cz);
                    zMinus.show();
                    ImagePlus zPlus  = ImgProUtils.zPlus(mask,  cz);
                    zPlus.show();
                    ImagePlus img = region.equals("negative_z") ? zMinus : zPlus;

                    BoundingBoxResult bbox = computeBoundingBox(img);

                    if (!bbox.valid) {
                        System.out.printf("%s %s %s: Empty image - writing zeros%n", 
                            specimen.getName(), ts, region);
                        appendEmptyBBoxResults(csvPath, specimen.getName(), ts, region);
                        continue;
                    }

                    printBoundingBoxResults(specimen.getName(), ts, bbox, region);
                    // appendBBoxResults(csvPath, specimen.getName(), ts, bbox, region);
                }
            }
        } else {
            // Process whole surface masks (no regions)
            for (int t = 0; t < timestamps.length; t++) {
                String ts = timestamps[t];
                String maskPath = Config.mainDir + "Processing/04_Masks/07_MaskSurface3D/" +   specimen.getName() + "_" + ts + "_mask_contour.tif";

                ImagePlus mask = IJ.openImage(maskPath);
                if (mask == null) {
                    System.out.printf("%s %s: Missing mask - writing zeros%n", 
                        specimen.getName(), ts);
                    appendEmptyBBoxResults(csvPath, specimen.getName(), ts, null);
                    continue;
                }

                BoundingBoxResult bbox = computeBoundingBox(mask);

                if (!bbox.valid) {
                    System.out.printf("%s %s: Empty image - writing zeros%n", 
                        specimen.getName(), ts);
                    appendEmptyBBoxResults(csvPath, specimen.getName(), ts, null);
                    continue;
                }

                printBoundingBoxResults(specimen.getName(), ts, bbox, null);
                appendBBoxResults(csvPath, specimen.getName(), ts, bbox, null);
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
            this.valid = (bounds[1] >= bounds[0] &&
                          bounds[3] >= bounds[2] &&
                          bounds[5] >= bounds[4]);

            if (valid) {
                this.width  = bounds[1] - bounds[0] + 1;
                this.height = bounds[3] - bounds[2] + 1;
                this.depth  = bounds[5] - bounds[4] + 1;
                this.volume = width * height * depth;
            } else {
                this.width = 0;
                this.height = 0;
                this.depth = 0;
                this.volume = 0;
            }
        }
    }

    public static BoundingBoxResult computeBoundingBox(ImagePlus img) {
        int[] bounds = ImgProUtils.bbox3D(img);
        return new BoundingBoxResult(bounds);
    }

    // === Output Formatting ===

    private static void printBoundingBoxResults(String specimen, String timestamp, 
                                               BoundingBoxResult bbox, String region) {
        if (region != null) {
            System.out.println(String.format("%s %s %s:", specimen, timestamp, region));
        } else {
            System.out.println(String.format("%s %s:", specimen, timestamp));
        }
        System.out.println(String.format("  BBox bounds: x[%d-%d] y[%d-%d] z[%d-%d]",
                bbox.bounds[0], bbox.bounds[1],
                bbox.bounds[2], bbox.bounds[3],
                bbox.bounds[4], bbox.bounds[5]));
        System.out.println(String.format("  BBox lengths: %d × %d × %d (volume=%d)",
                bbox.width, bbox.height, bbox.depth, bbox.volume));
        System.out.println();
    }

    // === CSV Management ===

    private static void ensureCsvHeader(String csv, boolean includeRegion) {
        File f = new File(csv);
        if (f.exists()) return;

        f.getParentFile().mkdirs();

        try (FileWriter fw = new FileWriter(csv);
             PrintWriter pw = new PrintWriter(fw)) {
            if (includeRegion) {
                pw.println("specimen,timestamp,region," +
                          "xmin,xmax,ymin,ymax,zmin,zmax," +
                          "width_px,height_px,depth_px," +
                          "volume_voxels");
            } else {
                pw.println("specimen,timestamp," +
                          "xmin,xmax,ymin,ymax,zmin,zmax," +
                          "width_px,height_px,depth_px," +
                          "volume_voxels");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void appendEmptyBBoxResults(String csv, String specimen, 
                                              String timestamp, String region) {
        try (FileWriter fw = new FileWriter(csv, true);
             PrintWriter pw = new PrintWriter(fw)) {
            if (region != null) {
                pw.printf(Locale.US,
                        "%s,%s,%s," +           // specimen,timestamp,region
                        "0,0,0,0,0,0," +        // bounds
                        "0,0,0," +              // dimensions
                        "0%n",                  // volume
                        specimen, timestamp, region);
            } else {
                pw.printf(Locale.US,
                        "%s,%s," +              // specimen,timestamp
                        "0,0,0,0,0,0," +        // bounds
                        "0,0,0," +              // dimensions
                        "0%n",                  // volume
                        specimen, timestamp);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void appendBBoxResults(String csv, String specimen, String timestamp, 
                                         BoundingBoxResult bbox, String region) {
        try (FileWriter fw = new FileWriter(csv, true);
             PrintWriter pw = new PrintWriter(fw)) {
            if (region != null) {
                pw.printf(Locale.US,
                        "%s,%s,%s," +           // specimen,timestamp,region
                        "%d,%d,%d,%d,%d,%d," +  // xmin,xmax,ymin,ymax,zmin,zmax
                        "%d,%d,%d," +           // width,height,depth
                        "%d%n",                 // volume
                        specimen, timestamp, region,
                        bbox.bounds[0], bbox.bounds[1],
                        bbox.bounds[2], bbox.bounds[3],
                        bbox.bounds[4], bbox.bounds[5],
                        bbox.width, bbox.height, bbox.depth,
                        bbox.volume);
            } else {
                pw.printf(Locale.US,
                        "%s,%s," +              // specimen,timestamp
                        "%d,%d,%d,%d,%d,%d," +  // xmin,xmax,ymin,ymax,zmin,zmax
                        "%d,%d,%d," +           // width,height,depth
                        "%d%n",                 // volume
                        specimen, timestamp,
                        bbox.bounds[0], bbox.bounds[1],
                        bbox.bounds[2], bbox.bounds[3],
                        bbox.bounds[4], bbox.bounds[5],
                        bbox.width, bbox.height, bbox.depth,
                        bbox.volume);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}








// package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D;

// import ij.IJ;
// import ij.ImageJ;
// import ij.ImagePlus;

// import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
// import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
// import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
// import io.github.rocsg.fijiyama.gargeetest.cuttings.helpers.ImgProUtils;

// import java.io.File;
// import java.io.FileWriter;
// import java.io.PrintWriter;
// import java.util.Locale;

// public class Step_9_BoundingBox implements PipelineStep {

//     @Override
//     public void execute(Specimen specimen) throws Exception {
//         execute(specimen, true);
//     }

//     public static void main(String[] args) throws Exception {
//         Specimen spec = new Specimen("B_202");
//         ImageJ ij = new ImageJ();
//         new Step_9_BoundingBox().execute(spec, true);
//     }

//     public void execute(Specimen specimen, boolean testing) throws Exception {
//         String[] timestamps = Config.timestamps;
//         String csv = Config.mainDir + "/Results/05_BoundingBox/BoundingBox_fit_results_of_mask.csv";
//         ensureCsvHeader(csv);
//         runBoundingBoxBatch(specimen, timestamps, csv);
//     }

//     // === Main Batch Processing ===

//     public static void runBoundingBoxBatch(Specimen specimen, String[] timestamps, String csvPath) {

//         for (int t = 0; t < timestamps.length; t++) {
//             String ts = timestamps[t];

//             String maskPath = Config.mainDir + "Processing/04_Masks/07_MaskSurface3D/" + specimen.getName() + "_" + ts + "_mask_contour.tif";

//             ImagePlus mask = IJ.openImage(maskPath);
//             if (mask == null) {
//                 System.out.printf("%s %s: Missing mask - writing zeros%n", specimen.getName(), ts);
//                 appendEmptyBBoxResults(csvPath, specimen.getName(), ts);
//                 continue;
//             }

//             BoundingBoxResult bbox = computeBoundingBox(mask);

//             if (!bbox.valid) {
//                 System.out.printf("%s %s: Empty image - writing zeros%n", specimen.getName(), ts);
//                 appendEmptyBBoxResults(csvPath, specimen.getName(), ts);
//                 continue;
//             }

//             printBoundingBoxResults(specimen.getName(), ts, bbox);
//             appendBBoxResults(csvPath, specimen.getName(), ts, bbox);
//         }
//     }

//     // === Bounding Box Computation ===

//     public static class BoundingBoxResult {
//         public int[] bounds;      // [xmin, xmax, ymin, ymax, zmin, zmax]
//         public int width;
//         public int height;
//         public int depth;
//         public int volume;        // in voxels
//         public boolean valid;

//         public BoundingBoxResult(int[] bounds) {
//             this.bounds = bounds;
//             this.valid = (bounds[1] >= bounds[0] &&
//                           bounds[3] >= bounds[2] &&
//                           bounds[5] >= bounds[4]);

//             if (valid) {
//                 this.width  = bounds[1] - bounds[0] + 1;
//                 this.height = bounds[3] - bounds[2] + 1;
//                 this.depth  = bounds[5] - bounds[4] + 1;
//                 this.volume = width * height * depth;
//             } else {
//                 this.width = 0;
//                 this.height = 0;
//                 this.depth = 0;
//                 this.volume = 0;
//             }
//         }
//     }

//     public static BoundingBoxResult computeBoundingBox(ImagePlus img) {
//         int[] bounds = ImgProUtils.bbox3D(img);
//         return new BoundingBoxResult(bounds);
//     }

//     // === Output Formatting ===

//     private static void printBoundingBoxResults(String specimen, String timestamp, BoundingBoxResult bbox) {
//         System.out.println(String.format("%s %s:", specimen, timestamp));
//         System.out.println(String.format("  BBox bounds: x[%d-%d] y[%d-%d] z[%d-%d]",
//                 bbox.bounds[0], bbox.bounds[1],
//                 bbox.bounds[2], bbox.bounds[3],
//                 bbox.bounds[4], bbox.bounds[5]));
//         System.out.println(String.format("  BBox lengths: %d, %d, %d",
//                 bbox.width, bbox.height, bbox.depth));
//         System.out.println();
//     }

//     // === CSV Management ===

//     private static void ensureCsvHeader(String csv) {
//         File f = new File(csv);
//         if (f.exists()) return;

//         f.getParentFile().mkdirs();

//         try (FileWriter fw = new FileWriter(csv);
//              PrintWriter pw = new PrintWriter(fw)) {
//             pw.println(
//                     "specimen,timestamp," +
//                     "xmin,xmax,ymin,ymax,zmin,zmax," +
//                     "width_px,height_px,depth_px," +
//                     "volume_voxels," 
//                 );
//         } catch (Exception e) {
//             e.printStackTrace();
//         }
//     }

//     private static void appendEmptyBBoxResults(String csv, String specimen, String timestamp) {
//         try (FileWriter fw = new FileWriter(csv, true);
//              PrintWriter pw = new PrintWriter(fw)) {

//             pw.printf(Locale.US,
//                     "%s,%s," +              // specimen,timestamp
//                     "0,0,0,0,0,0," +        // bounds
//                     "0,0,0," +        // half-width/half-height/half-depth
//                     "0%n",                 // volume
//                     specimen, timestamp);
//         } catch (Exception e) {
//             e.printStackTrace();
//         }
//     }

//     private static void appendBBoxResults(String csv, String specimen, String timestamp, BoundingBoxResult bbox) {
//         try (FileWriter fw = new FileWriter(csv, true);
//             PrintWriter pw = new PrintWriter(fw)) {

//             pw.printf(java.util.Locale.US,
//                     "%s,%s," +              // specimen,timestamp
//                     "%d,%d,%d,%d,%d,%d," +  // xmin,xmax,ymin,ymax,zmin,zmax
//                     "%d,%d,%d," +           // width_px,height_px,depth_px
//                     "%d%n",                 // volume_voxels
//                     specimen, timestamp,
//                     bbox.bounds[0], bbox.bounds[1],
//                     bbox.bounds[2], bbox.bounds[3],
//                     bbox.bounds[4], bbox.bounds[5],
//                     bbox.width, bbox.height, bbox.depth,
//                     bbox.volume);
//         } catch (Exception e) {
//             e.printStackTrace();
//         }
//     }
// }
