package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.registration.ItkTransform;

public class HyperMapforFijirelax {

    public static void main(String[] args) {
        ImageJ ij = new ImageJ();
        ImagePlus img001 = IJ.openImage("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/hypermap_registered_normalized/206_J001_hyperMap_reg_normalized.tif");
        ImagePlus img029 = IJ.openImage("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/hypermap_registered_normalized/206_J029_hyperMap_reg_normalized.tif");
        ImagePlus img077 = IJ.openImage("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/hypermap_registered_normalized/206_J077_hyperMap_reg_normalized.tif");
        ImagePlus img141 = IJ.openImage("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/hypermap_registered_normalized/206_J141_hyperMap_reg_normalized.tif");

        // ItkTransform tr001 = ItkTransform.readTransformFromFile("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/test/output/Exported_data/transform_global_206_J001_TR2400_TE12.txt");
        // ItkTransform tr029 = ItkTransform.readTransformFromFile("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/test/output/Exported_data/transform_global_206_J029_TR2400_TE12.txt");
        // ItkTransform tr077 = ItkTransform.readTransformFromFile("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/test/output/Exported_data/transform_global_206_J077_TR2400_TE12.txt");
        // ItkTransform tr141= ItkTransform.readTransformFromFile("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/test/output/Exported_data/transform_global_206_J141_TR2400_TE12.txt");
        


        ImagePlus imgTransformed001 = img001;
        ImagePlus imgTransformed029 = img029;
        ImagePlus imgTransformed077 = img077;
        ImagePlus imgTransformed141 = img141;


        // imgTransformed001 = tr001.transformImage(img001,imgTransformed001);
        // imgTransformed029 = tr029.transformImage(img029,imgTransformed029);
        // imgTransformed077 = tr077.transformImage(img077,imgTransformed077);
        // imgTransformed141 = tr141.transformImage(img141,imgTransformed141);



        ImagePlus[] tabImgover4setsofdays = new ImagePlus[4];
      
        tabImgover4setsofdays[0] = imgTransformed001;
        tabImgover4setsofdays[1] = imgTransformed029;
        tabImgover4setsofdays[2] = imgTransformed077;
        tabImgover4setsofdays[3] = imgTransformed141;
        for (int i = 0; i < tabImgover4setsofdays.length; i++) {
            ImagePlus imgTransformed = tabImgover4setsofdays[i];
            // save to tab
            tabImgover4setsofdays[i] = imgTransformed;
            
            System.out.println(tabImgover4setsofdays[i]);
        }
       
        ImagePlus hyperStackXR = VitimageUtils.hyperStackingChannels(tabImgover4setsofdays);
        ImagePlus hyperFrame = VitimageUtils.hyperStackChannelToHyperStackFrame(hyperStackXR);
        hyperFrame.show();
        hyperFrame.setTitle("B_206_hypermap");
        // IJ.saveAsTiff(hyperFrame,"/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/test/output_highres/B_206_hypermap.tif");


    }
    
}
