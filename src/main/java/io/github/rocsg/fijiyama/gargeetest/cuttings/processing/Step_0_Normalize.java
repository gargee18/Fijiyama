package io.github.rocsg.fijiyama.gargeetest.cuttings.processing;

import java.util.Arrays;

import ij.IJ;
import ij.ImageJ;
// import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.FloatProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;

public class Step_0_Normalize implements PipelineStep{

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }

    public static void main(String[] args) throws Exception {
        ImageJ ij = new ImageJ();
        Specimen spec = Specimen.getSpecimen("B_239");
        System.out.println(spec);
        new Step_0_Normalize().execute(spec);

    }

    public void execute(Specimen specimen,boolean testing) throws Exception {
        String[]timestamps=Config.timestamps;
        int N=timestamps.length;
        for(int n=0;n<(testing ? 1 : N);n++){
            ImagePlus imgToGetCapillary = IJ.openImage( Config.getPathToCroppedImage(specimen, n));
            System.out.println( Config.getPathToCroppedImage(specimen, n));
            ImagePlus imgToNormalize = IJ.openImage( Config.getPathToRawImage(specimen, n));
            System.out.println(Config.getPathToRawImage(specimen, n));
            double capillarySize = 14;// Care : when use with pixel size, divide by 28 (35µm), making 0.5 (diameter)

            int[] capCentre = VitimageUtils.findCapillaryCenterInSlice(imgToGetCapillary,capillarySize);
            // int[] capCentre = {437,284,785}; //use when you need to manualy enter the capillary central coordinates
            System.out.println("Capillary Centre: " + Arrays.toString(capCentre));
            double[] capValues = VitimageUtils.capillaryValuesAlongZStatic(imgToGetCapillary, capCentre, capillarySize);

            double medianCap = VitimageUtils.MADeStatsDoubleSided(capValues, capValues)[0];

            double[] offset = VitimageUtils.caracterizeBackgroundOfImage(imgToGetCapillary);
            System.out.println("Mean Background: " + offset[0] + " Std: " + offset[1]);

            ImagePlus imgNorm = processNormalizationWithRespectToCapillary(imgToNormalize, offset[0], medianCap, 1.0);

            Calibration cal = imgToNormalize.getCalibration();
            System.out.println(cal.pixelWidth);
            cal.setUnit("mm");
            imgNorm.setCalibration(cal);

            if(testing)imgNorm.show();
            imgNorm.setDisplayRange(0, Config.max_display_val);
            VitimageUtils.setLutToFire(imgNorm);
            IJ.saveAsTiff(imgNorm,Config.getPathToNormalizedImage(specimen,n));
        }
    }

    public static ImagePlus processNormalizationWithRespectToCapillary(ImagePlus image, double backgroundVal,double capillaryVal, double targetCapVal) {
        ImageStack stack = image.getStack();
        int stackSize = stack.getSize();
        int width = stack.getWidth();
        int height = stack.getHeight();

        ImageStack normalizedStack = new ImageStack(width, height);

        for (int z = 1; z <= stackSize; z++) {
            FloatProcessor sliceProcessor = (FloatProcessor) stack.getProcessor(z);
            FloatProcessor normalizedSlice = (FloatProcessor) sliceProcessor.duplicate();

            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    double intensity = sliceProcessor.getPixelValue(x, y);
                    double normalizedIntensity = ((intensity - backgroundVal) / (capillaryVal - backgroundVal))
                            * targetCapVal;
                    normalizedSlice.putPixelValue(x, y, (float) normalizedIntensity);
                }
            }
            normalizedStack.addSlice(stack.getSliceLabel(z), normalizedSlice);
        }
        ImagePlus normalizedImage = new ImagePlus(image.getTitle() + "_normalized", normalizedStack);

        return normalizedImage;

    }
   
   
}
