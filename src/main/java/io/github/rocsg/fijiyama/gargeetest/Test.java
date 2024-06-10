package io.github.rocsg.fijiyama.gargeetest;

import org.apache.commons.math3.transform.TransformUtils;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;

public class Test {
    public static void main(String[] args) {
        ImageJ ij=new ImageJ();
        String dir="/home/phukon/Desktop/";
        ImagePlus img1=IJ.openImage(dir+"Bla.tif");
        ImagePlus img2=IJ.openImage(dir+"Bla2.tif");
        img1.setTitle("img1");
        img1.show();
        img2.setTitle("img2");
        img2.show();
        ImagePlus img1_switched=VitimageUtils.switchAxis(img1,1);
        ImagePlus img2_switched=VitimageUtils.switchAxis(img2,1);
        img1_switched.setTitle("img1_switched");
        img1_switched.show();
        img2_switched.setTitle("img2_switched");
        img2_switched.show();
    }
        
}
