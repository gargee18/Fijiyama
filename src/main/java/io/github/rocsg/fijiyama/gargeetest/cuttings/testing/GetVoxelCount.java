package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.helpers.ImgProUtils;

public class GetVoxelCount {

    public static void main(String[] args) {
        ImageJ ij = new ImageJ();

        ImagePlus img = IJ.openImage(Config.mainDir+"Processing/04_Masks/05_MaskPithFilledCropped/B_201_J141_mask.tif");
        long num_surface_voxels = ImgProUtils.countForegroundVoxels(img,255);
        System.out.println(num_surface_voxels);
    }
    
}
