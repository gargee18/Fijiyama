package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;
//specific libraries
import ij.IJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.registration.ItkTransform;
public class ApplyTrMatrixToAllSequence {

    public static void main(String[] args) {
        test1();
    }

    public static ImagePlus test1() {
        ImagePlus openImgSeq = IJ.openImage("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/206_J141_hyperMap.tif");
        openImgSeq.show();
        ItkTransform tr001 = ItkTransform.readTransformFromFile("/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/test/output/Exported_data/transform_global_206_J141_TR2400_TE12.txt");
        ImagePlus imgTransformed001 = openImgSeq;
        openImgSeq.show();
        imgTransformed001 = tr001.transformImage(openImgSeq,imgTransformed001);
        imgTransformed001.show();
        IJ.saveAsTiff(imgTransformed001, "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/raw_registered/206_J141_hyperMap_reg.tif");
        return null;
}
    
}
