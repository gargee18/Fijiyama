package io.github.rocsg.fijiyama.gargeetest.ceps;

//specific libraries
import ij.IJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
public class MakeCompositeQuick {

    //public static String[] spec = GeneralUtils.specimenListXR;
    public static String[] spec = {"1195"};
    public static String commondir = "/home/phukon/Desktop/xray_2023_2024_hires/";

    public static void main(String[] args){
        
        for(String s: spec){
            String dirRefImage = "/home/phukon/Desktop/xray_2023_2024_hires/cep_1195/raw/1195_2023.tif";
            String dirRegImage = "/home/phukon/Desktop/xray_2023_2024_hires/cep_1195/res/1195_TransformedImage.tif";
            String dirsaveimg = "home/phukon/Desktop/xray_2023_2024_hires/cep_1195/res/";
            ImagePlus refimg = IJ.openImage(dirRefImage);
            ImagePlus regimg = IJ.openImage(dirRegImage);
            ImagePlus compositeimg = VitimageUtils.compositeNoAdjustOf(refimg, regimg);
            compositeimg.show();
            IJ.saveAsTiff(compositeimg, "/home/phukon/Desktop/xray_2023_2024_hires/cep_1195/res/1195_Composite.tif");
        }


    }

}
