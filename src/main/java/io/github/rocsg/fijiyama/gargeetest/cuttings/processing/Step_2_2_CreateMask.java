package io.github.rocsg.fijiyama.gargeetest.cuttings.processing;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import inra.ijpb.morphology.Morphology;


public class Step_2_2_CreateMask implements PipelineStep{

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }

    public static void main(String[] args) throws Exception {
        ImageJ ij = new ImageJ();
        Specimen spec = Specimen.getSpecimen("B_238");
        System.out.println(spec);
        new Step_2_2_CreateMask().execute(spec);

    }

    public void execute(Specimen specimen,boolean testing) throws Exception {
        String[]timestamps=Config.timestamps;
        // int N=timestamps.length;
        for(int n=0;n<1;n++){
        
            // Load reference image
            // System.out.println(Config.getPathToInocAlignedImageSub(specimen, n));
            String imagePath = Config.getPathToInocAlignedImageSub(specimen, n);

            ImagePlus imgRef = IJ.openImage(imagePath);

            // Get connected components
            ImagePlus connectedComponents = VitimageUtils.connexeNoFuckWithVolume(imgRef, 0, 8, 2000, 36864, 6, 1, true);

            // Convert to 8-bit grayscale
            VitimageUtils.convertToGray8(connectedComponents);

            // Apply median filter
            IJ.run(connectedComponents, "Median...", "radius=10 stack");

            // Duplicate image and fill holes
            ImagePlus filledImage = connectedComponents.duplicate();
            IJ.run(filledImage, "Fill Holes", "stack");

            // Subtract images to extract only pith mask
            ImagePlus pithMask = VitimageUtils.makeOperationBetweenTwoImages(connectedComponents, filledImage, 4, false);

            // Define the outer and inner radii for morphological operations  
            int radiusOutside=20;
            int radiusInside=15;

            // Perform dilation on the pithMask image using a disk-shaped structuring element of radius 'radiusOutside'
            ImagePlus dilationTest = Morphology.dilation(pithMask, inra.ijpb.morphology.strel.DiskStrel.fromRadius(radiusOutside)); 

            // Perform erosion on the pithMask image using a disk-shaped structuring element of radius 'radiusInside'
            ImagePlus erosionTest = Morphology.erosion(pithMask, inra.ijpb.morphology.strel.DiskStrel.fromRadius(radiusInside));

            // Subtract the eroded image from the dilated image to get a specific region of interest
            ImagePlus imgSub = VitimageUtils.makeOperationBetweenTwoImages(dilationTest, erosionTest, 4, false);

            // Get the dimensions (width, height, and stack size) of the resulting image
            int X=imgSub.getWidth();
            int Y=imgSub.getHeight();
            int Z=imgSub.getStackSize();

            // Draw a parallelepiped shape on the image 
            imgSub=VitimageUtils.drawParallepipedInImage(imgSub,0,0,0,X-1,Y-1,110, 0);
            imgSub=VitimageUtils.drawParallepipedInImage(imgSub,0,0,Z-110,X-1,Y-1,Z-1, 0);

            IJ.saveAsTiff(imgSub, Config.getPathToMask(specimen,n));

            // Create a composite image by overlaying imgSub onto imgRef with the label "Mask on ref"
            ImagePlus composite= VitimageUtils.compositeNoAdjustOf(imgRef,imgSub,"Mask on ref");
            composite.show(); 
            VitimageUtils.waitFor(20000);
            composite.close(); 
        }
    }

    
    
}
