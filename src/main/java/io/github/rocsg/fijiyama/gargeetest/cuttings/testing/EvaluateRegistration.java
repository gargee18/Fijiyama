package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;

public class EvaluateRegistration implements PipelineStep{

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }
    public static void main(String[] args) throws Exception {

        Specimen spec= new Specimen("B_202");
        ImageJ ij=new ImageJ();//Needed for testing
        new  EvaluateRegistration().execute(spec,true);
        System.out.println("Saved!");
    }
    
    public void execute(Specimen specimen,boolean testing) throws Exception {
    
        // Load images (modify paths to your actual image files)
        String[] timestamps = Config.timestamps;
        ImagePlus img = IJ.openImage("/home/phukon/Desktop/Cuttings/Results/01_Hyperstack/"+specimen.getName()+"_Hyperstack.tif");
        System.out.println(specimen);
        for (int i = 1; i<4; i++){
            ImagePlus refImage = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, i    , i);
            refImage.show();
            ImagePlus movingImage = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, i+1, i+1);
            movingImage.show();
            double nccScore = computeNCC(refImage, movingImage);
            System.out.println("3D NCC Score of image "+timestamps[i-1]+" and "+timestamps[i]+" : " + nccScore);
            VitimageUtils.waitFor(10000);
            refImage.close();
            movingImage.close();
        }
        
       
        
    }
    public static double computeNCC(ImagePlus refImage, ImagePlus movingImage) {

        int depth = refImage.getStackSize(); // Number of slices
        int width = refImage.getWidth();
        int height = refImage.getHeight();

        // Compute means
        double meanA = computeMean(refImage);
        double meanB = computeMean(movingImage);

        double numerator = 0, denominatorA = 0, denominatorB = 0;

        // Iterate through all voxels
        for (int z = 1; z <= depth; z++) {  // ImageJ slice indices start at 1
            ImageProcessor refSlice = refImage.getStack().getProcessor(z);
            ImageProcessor movSlice = movingImage.getStack().getProcessor(z);

            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    double voxelA = refSlice.getPixelValue(x, y);
                    double voxelB = movSlice.getPixelValue(x, y);

                    double diffA = voxelA - meanA;
                    double diffB = voxelB - meanB;

                    numerator += diffA * diffB;
                    denominatorA += diffA * diffA;
                    denominatorB += diffB * diffB;
                }
            }
        }

        double denominator = Math.sqrt(denominatorA) * Math.sqrt(denominatorB);
        return (denominator == 0) ? 0 : numerator / denominator;
    }

    private static double computeMean(ImagePlus image) {
        double sum = 0;
        int count = 0;

        int depth = image.getStackSize();
        int width = image.getWidth();
        int height = image.getHeight();

        for (int z = 1; z <= depth; z++) {
            ImageProcessor slice = image.getStack().getProcessor(z);
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    sum += slice.getPixelValue(x, y);
                    count++;
                }
            }
        }
        return sum / count;
    }

    

    
}
