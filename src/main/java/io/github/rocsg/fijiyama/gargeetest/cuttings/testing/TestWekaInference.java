package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import trainableSegmentation.WekaSegmentation;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
public class TestWekaInference {
    
    public static void main(String[] args) {
        String timestamps = "J141";
        String specimenName = "B_206";
        int frame = 4;
        String path = Config.mainDir+"/Processing/03_PolarTransform/B_206_GeneralizedPolarTransform.tif";
        ImagePlus img = trainWekaToGetMask(path, frame); 
        System.out.println("done!");
        img.show();
    }
    public static ImagePlus trainWekaToGetMask(String path, int frame){
        ImagePlus hyperframe = IJ.openImage(path);
        System.out.println("image opened"+hyperframe);
        ImagePlus img = new Duplicator().run(hyperframe, 1, 1, 256, 767, frame, frame);
        WekaSegmentation weka = new WekaSegmentation(img);
        weka.loadClassifier( Config.mainDir+"/Processing/Trained_Models/All_Var_PCH_z512_545_610_tJ141.model");
        ImagePlus proba = weka.applyClassifier(img, 0, true);
        IJ.saveAsTiff(proba, Config.mainDir+"Processing/test_trained_model/res_206_141_proba.tif");
        ImagePlus imgMask=new Duplicator().run(proba,1,1,1,proba.getNSlices(),1,1);
        imgMask.setDisplayRange(0.5, 0.5);
        VitimageUtils.convertToGray8(imgMask);
        IJ.run(imgMask,"Median...", "radius=1 stack");
        IJ.run(imgMask, "Divide...", "value=255 stack");   
        imgMask.setDisplayRange(0, 1); 
        IJ.saveAsTiff(imgMask, Config.mainDir+"Processing/test_trained_model/res_206_141_mask.tif");
        return imgMask;
    }
}