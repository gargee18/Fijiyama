// package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;
// import java.io.File;

// import ij.IJ;
// import ij.ImagePlus;
// import ij.gui.Roi;
// import ij.io.RoiDecoder;
// import ij.plugin.Duplicator;
// import ij.plugin.frame.RoiManager;
// import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
// import java.io.File;

// public class TestCartesianVsPolar {
//    public static void main(String[] args) {
//         String base1 = Config.mainDir + "/Results/01_Hyperstack/";
//         String base2 = Config.mainDir +"/Processing/03_PolarTransform/";
//         String save  = Config.mainDir + "/Test_cartesian_vs_polar/";

//         RoiManager rm = RoiManager.getInstance2();
//         rm.reset();
// 		   Roi loaded = RoiDecoder.open(new File(roiPath));   // ✅ loads a single .roi
//         if (loaded == null) {
//             IJ.error("Could not load ROI: " + roiPath);
//             return;
//         }
//         rm.addRoi(loaded);

        

//         for (int i = 201; i <= 240; i++) {
//             if (i == 204) continue;
//             String id = String.format("B_%03d", i);
//             System.out.println("Processing " + id);

//             ImagePlus imp1 = IJ.openImage(base1 + id + "_Hyperstack.tif");
//             imp1.show();
//             if (imp1 != null) {
//                 if (rm.getCount() > 0) rm.select(imp1, 0);
//                 new Duplicator().run(imp1, 1, 1, 1, 1024, 3, 3);
//                 imp1.show();
//             }

//             ImagePlus imp2 = IJ.openImage(base2 + id + "_GeneralizedPolarTransform.tif");
//             if (imp2 != null) {
//                 if (rm.getCount() > 0) rm.select(imp2, 0);
//                 ImagePlus dup = new Duplicator().run(imp2, 1, 1, 1, 1024, 3, 3);
//                 System.out.println( save + id + "_GeneralizedPolarTransform_duplicated.tif");
//                 // IJ.saveAsTiff(dup, save + id + "_GeneralizedPolarTransform_duplicated.tif");
//                 dup.show();
//                 imp2.close();
//             }
//         }
//         System.out.println("Done.");
//     }
    
// }
