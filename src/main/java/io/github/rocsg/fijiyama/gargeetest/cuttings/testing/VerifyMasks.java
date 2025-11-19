package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;

import java.util.List;

import ij.ImageJ;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.gargeetest.cuttings.helpers.ImgProUtils;

public class VerifyMasks {
    public static void main(String[] args) {
        ImageJ ij = new ImageJ();
        List<Specimen> specimens = Specimen.loadAllSpecimens();
        for (Specimen specimen : specimens) {
            ImgProUtils.verfifyProcessing(specimen, Config.timestamps);
        }
    }
    
}
