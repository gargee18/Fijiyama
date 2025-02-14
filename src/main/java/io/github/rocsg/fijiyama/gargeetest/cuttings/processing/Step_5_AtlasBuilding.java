/**
 * This class, `step_5_AtlasBuilding`, is part of the image processing pipeline
 * for analyzing plant specimen hyperstacks. It provides methods to compute 
 * average and standard deviation images based on timepoint differences (e.g., T1-T0, T2-T0, T2-T1) and
 * for various conditions and plant varieties. 
 */


package io.github.rocsg.fijiyama.gargeetest.cuttings.processing;

import java.util.ArrayList;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;

public class Step_5_AtlasBuilding {

    static final int cond_PCH = 0;
    static final int cond_CONTROL = 1;

    static final int var_CHARD = 0;
    static final int var_MERLOT = 1;
    static final int var_TEMPRA = 2;
    static final int var_UGNI = 3;


    static final String[] VAR_NAMES = {
        "var_CHARD", "var_MERLOT", "var_TEMPRA", "var_UGNI"
    };

    public static String getVarietyName(int value) {
        if (value >= 0 && value < VAR_NAMES.length) {
            return VAR_NAMES[value];
        }
        return "Unknown";
    }

    public static void main(String[] args) {
        
        ImageJ ij=new ImageJ();
   //    computeAllHyperStacks();

    for ( int var = 0; var < 4; var ++){
            ImagePlus[] atlasCtAllVar = computeAverageAndStdForSpecificVariety( cond_CONTROL, var);
            atlasCtAllVar[1].setTitle("stdCT");
            atlasCtAllVar[1].show();
            atlasCtAllVar[1].setDisplayRange(0, 70);
            atlasCtAllVar[0].setTitle("meanCT");
            atlasCtAllVar[0].show();
            atlasCtAllVar[0].setDisplayRange(0, 70);
               IJ.saveAsTiff( atlasCtAllVar[0]    ,Config.mainDir+"PolarAtlas/mean_"+getVarietyName(var)+"_CT_day_141.tif");
               IJ.saveAsTiff( atlasCtAllVar[1]    ,Config.mainDir+"PolarAtlas/std_"+getVarietyName(var)+"_CT_day_141.tif");
            ImagePlus[] atlasPchAllVar = computeAverageAndStdForSpecificVariety( cond_PCH, var);
            atlasPchAllVar[1].setTitle("stdPCH");
            atlasPchAllVar[1].show();
            atlasPchAllVar[1].setDisplayRange(0, 70);
            atlasPchAllVar[0].setTitle("meanPCH");
            atlasPchAllVar[0].show();
            atlasPchAllVar[0].setDisplayRange(0, 70);
            IJ.saveAsTiff( atlasPchAllVar[0]    ,Config.mainDir+"PolarAtlas/mean_"+getVarietyName(var)+"_PCH_day_141.tif");
            IJ.saveAsTiff( atlasPchAllVar[1]    ,Config.mainDir+"PolarAtlas/std_"+getVarietyName(var)+"_PCH_day_141.tif");
        }
       
    }





    /**
     * getSpecimensName:
     * - Retrieves the names of specimen identifiers based on the provided condition and variety.
     * - Input: 
     *   - `condition`: An integer representing the condition (e.g., PCH, CONTROL).
     *   - `variety`: An integer representing the variety (e.g., CHARD, MERLOT).
     * - Output: A string array of specimen identifiers for the specified condition and variety.
     */
    public static String[] getSpecimensName(int condition, int variety) {
        String[][][] specs = new String[][][] {
                {
                        { "B_201", "B_202", "B_203", "B_205", "B_206" }/* CHARD */,
                        { "B_211", "B_212", "B_213", "B_214", "B_215", "B_216", "B_219" }/* MERLOT */,
                        { "B_221", "B_222", "B_223", "B_224", "B_225", "B_226" }/* TEMPRA */,
                        { "B_231", "B_232", "B_233", "B_234", "B_235", "B_236" }/* UGNI */
                }/* PCH */,
                {
                        { "B_207", "B_208", "B_209", "B_210" }/* CHARD */,
                        { "B_217", "B_218", "B_220" }/* MERLOT */,
                        { "B_227", "B_228", "B_229", "B_230" }/* TEMPRA */,
                        { "B_237", "B_238", "B_239", "B_240" }/* UGNI */
                }/* CT */
        };

        return specs[condition][variety];
    }

    public static String getDirOfSpecimen(String specimenName) {
        return Config.mainDir + specimenName;
    }
    /**
     * computeAverageAndStdIndividualAtT1MinusT0:
     * - Processes hyperstack images to compute the average and standard deviation of T1-T0 images.
     * - Duplicates the hyperstack to isolate timepoints 0 and 1, calculates the difference, 
     *   and aggregates the results for multiple specimens.
     * - Input: 
     *   - `condition`: An integer representing the condition.
     * - Output: An array of two `ImagePlus` objects (mean and standard deviation images).
     */


    public static ImagePlus[] computeAverageAndStdIndividualAtT1MinusT0(int condition) {
        ArrayList<ImagePlus> stacks = new ArrayList<ImagePlus>();
        for (int cond = condition; cond <= condition; cond++) {
            for (int var = 0; var < 1; var++) {
                String[] spec = getSpecimensName(cond, var);
                for (int i = 0; i < spec.length; i++) {
                    System.out.println("Specimen: " + spec[i] + " Condition: " + cond + " Variety: " + var);
                    ImagePlus img = IJ
                            .openImage(getDirOfSpecimen(spec[i]) + "/hyperimage/" + spec[i] + "_Hyperstack.tif");
                    
                    // I would like to display if img is null, and then print the argument of
                    // openImage
                    System.out.println("Image: " + img);
                    System.out.println(getDirOfSpecimen(spec[i]));

                    ImagePlus imgT0 = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, 1, 1);
                    ImagePlus imgT1 = new Duplicator().run(img, 2, 2, 256, img.getNSlices() - 256, 1, 1);
                    img = VitimageUtils.makeOperationBetweenTwoImages(imgT1, imgT0, 4, true);
                    imgT0 = null;
                    imgT1 = null;
                    stacks.add(img);
                }
            }
        }

        ImagePlus[] imgTab = new ImagePlus[stacks.size()];
        for (int i = 0; i < stacks.size(); i++) {
            System.out.println("Image: " + stacks.get(i));
            imgTab[i] = stacks.get(i);
        }
        ImagePlus[] res = meanAndStdOfImageArrayFloatFreshNew(imgTab);
        return res;
    }
    /**
     * computeAverageAndStdIndividualAtT2MinusT0:
     * - Similar to `computeAverageAndStdIndividualAtT1MinusT0`, but processes T2-T0 images.
     * - Input: 
     *   - `condition`: An integer representing the condition.
     * - Output: An array of two `ImagePlus` objects (mean and standard deviation images).
     */
    public static ImagePlus[] computeAverageAndStdIndividualAtT2MinusT0(int condition) {
        ArrayList<ImagePlus> stacks = new ArrayList<ImagePlus>();
        for (int cond = condition; cond <= condition; cond++) {
            for (int var = 0; var < 4; var++) {
                String[] spec = getSpecimensName(cond, var);
                for (int i = 0; i < spec.length; i++) {
                    System.out.println("Specimen: " + spec[i] + " Condition: " + cond + " Variety: " + var);
                    ImagePlus img = IJ
                            .openImage(getDirOfSpecimen(spec[i]) + "/hyperimage/" + spec[i] + "_Hyperstack.tif");
                    // I would like to display if img is null, and then print the argument of
                    // openImage
                    System.out.println("Image: " + img);
                    System.out.println(getDirOfSpecimen(spec[i]));

                    ImagePlus imgT0 = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, 1, 1);
                    ImagePlus imgT2 = new Duplicator().run(img, 3, 3, 256, img.getNSlices() - 256, 1, 1);
                    img = VitimageUtils.makeOperationBetweenTwoImages(imgT2, imgT0, 4, true);
                    imgT0 = null;
                    imgT2 = null;
                    stacks.add(img);
                }
            }
        }

        ImagePlus[] imgTab = new ImagePlus[stacks.size()];
        for (int i = 0; i < stacks.size(); i++) {
            System.out.println("Image: " + stacks.get(i));
            imgTab[i] = stacks.get(i);
        }
        // ImagePlus[] res=meanAndStdOfImageArrayByteFreshNew(imgTab);
        ImagePlus[] res = meanAndStdOfImageArrayFloatFreshNew(imgTab);
        return res;
    }
    /**
     * computeAverageAndStdIndividualAtT2MinusT1:
     * - Similar to `computeAverageAndStdIndividualAtT1MinusT0`, but processes T2-T1 images.
     * - Input: 
     *   - `condition`: An integer representing the condition.
     * - Output: An array of two `ImagePlus` objects (mean and standard deviation images).
     */
    public static ImagePlus[] computeAverageAndStdIndividualAtT2MinusT1(int condition) {
        ArrayList<ImagePlus> stacks = new ArrayList<ImagePlus>();
        for (int cond = condition; cond <= condition; cond++) {
            for (int var = 0; var < 4; var++) {
                String[] spec = getSpecimensName(cond, var);
                for (int i = 0; i < spec.length; i++) {
                    System.out.println("Specimen: " + spec[i] + " Condition: " + cond + " Variety: " + var);
                    ImagePlus img = IJ
                            .openImage(getDirOfSpecimen(spec[i]) + "/hyperimage/" + spec[i] + "_Hyperstack.tif");
                    // I would like to display if img is null, and then print the argument of
                    // openImage
                    System.out.println("Image: " + img);
                    System.out.println(getDirOfSpecimen(spec[i]));

                    ImagePlus imgT1 = new Duplicator().run(img, 2, 2, 256, img.getNSlices() - 256, 1, 1);
                    ImagePlus imgT2 = new Duplicator().run(img, 3, 3, 256, img.getNSlices() - 256, 1, 1);
                    img = VitimageUtils.makeOperationBetweenTwoImages(imgT2, imgT1, 4, true);
                    imgT1 = null;
                    imgT2 = null;
                    stacks.add(img);
                }
            }
        }

        ImagePlus[] imgTab = new ImagePlus[stacks.size()];
        for (int i = 0; i < stacks.size(); i++) {
            System.out.println("Image: " + stacks.get(i));
            imgTab[i] = stacks.get(i);
        }
        // ImagePlus[] res=meanAndStdOfImageArrayByteFreshNew(imgTab);
        ImagePlus[] res = meanAndStdOfImageArrayFloatFreshNew(imgTab);
        return res;
    }
    /**
     * computeAverageAndStdForSpecificVariety:
     * - Processes hyperstack images for a specific variety (e.g., UGNI) and computes 
     *   the average and standard deviation for the selected timepoints.
     * - Input: 
     *   - `condition`: An integer representing the condition.
     * - Output: An array of two `ImagePlus` objects (mean and standard deviation images).
     */
    public static ImagePlus[] computeAverageAndStdForSpecificVariety(int condition, int var) {
        ArrayList<ImagePlus> stacks = new ArrayList<ImagePlus>();
        for (int cond = condition; cond <= condition; cond++) {
            // for (int var = 0; var < variety; var++) {
                String[] spec = getSpecimensName(cond, var);
                for (int i = 0; i < spec.length; i++) {
                    System.out.println("Specimen: " + spec[i] + " Condition: " + cond + " Variety: " + var);
                    // ImagePlus img = IJ.openImage(getDirOfSpecimen(spec[i]) + "/hyperimage/" + spec[i] + "_Hyperstack.tif");
                    ImagePlus img = IJ.openImage(getDirOfSpecimen(spec[i]) +"/Polar_Transform_Results/GeneralizedPolarTransform.tif");
                    // I would like to display if img is null, and then print the argument of
                    // openImage
                    System.out.println("Image: " + img);
                    System.out.println(getDirOfSpecimen(spec[i]));
                    ImagePlus imgVar = new Duplicator().run(img, 3, 3, 256, img.getNSlices() - 256, 3, 3);
                    stacks.add(imgVar);

                // }
            }
        }

        ImagePlus[] imgTab = new ImagePlus[stacks.size()];
        for (int i = 0; i < stacks.size(); i++) {
            System.out.println("Image: " + stacks.get(i));
            imgTab[i] = stacks.get(i);
        }
        // ImagePlus[] res = meanAndStdOfImageArrayByteFreshNew(imgTab);
        ImagePlus[] res=meanAndStdOfImageArrayFloatFreshNew(imgTab);
        return res;
    }
    /**
     * computeAverageAndStdAtDifferentT:
     * - Computes the average and standard deviation across hyperstacks at a specified timepoint (e.g., T3).
     * - Processes multiple specimens for different conditions and varieties.
     * - Input: 
     *   - `condition`: An integer representing the condition.
     * - Output: An array of two `ImagePlus` objects (mean and standard deviation images).
     */
    public static ImagePlus[] computeAverageAndStdAtDifferentT(int condition) {
        ArrayList<ImagePlus> stacks = new ArrayList<ImagePlus>();
        for (int cond = condition; cond <= condition; cond++) {
            for (int var = 0; var < 4; var++) {
                String[] spec = getSpecimensName(cond, var);
                for (int i = 0; i < spec.length; i++) {
                    System.out.println("Specimen: " + spec[i] + " Condition: " + cond + " Variety: " + var);
                    // ImagePlus img = IJ.openImage(getDirOfSpecimen(spec[i]) + "/hyperimage/" + spec[i] + "_Hyperstack.tif");
                    ImagePlus img = IJ.openImage("/home/phukon/Desktop/Cuttings_MRI_registration/"+ spec[i] +"/Polar_Transform_Results/GeneralizedPolarTransform.tif");
                    // I would like to display if img is null, and then print the argument of
                    // openImage
                    System.out.println("Image: " + img);
                    System.out.println(getDirOfSpecimen(spec[i]));
                    ImagePlus imgT = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, 0, 0);
                    stacks.add(imgT);

                }
            }
        }

        ImagePlus[] imgTab = new ImagePlus[stacks.size()];
        for (int i = 0; i < stacks.size(); i++) {
            System.out.println("Image: " + stacks.get(i));
            imgTab[i] = stacks.get(i);
        }
        // ImagePlus[] res = meanAndStdOfImageArrayByteFreshNew(imgTab);
        ImagePlus[] res=meanAndStdOfImageArrayFloatFreshNew(imgTab);
        return res;
    }
    /**
     * meanAndStdOfImageArrayByteFreshNew:
     * - Calculates the mean and standard deviation images for a set of `ImagePlus` objects
     *   containing 8-bit image data.
     * - Input: 
     *   - `imgs`: An array of `ImagePlus` objects to process.
     * - Output: An array of two `ImagePlus` objects (mean and standard deviation images).
     */
    public static ImagePlus[] meanAndStdOfImageArrayByteFreshNew(ImagePlus[] imgs) {
        int nImg = imgs.length;
        ImagePlus imgMean = imgs[0];
        imgMean = VitimageUtils.convertByteToFloatWithoutDynamicChanges(imgMean);
        imgMean = VitimageUtils.nullImage(imgMean);
        ImagePlus imgStd = VitimageUtils.nullImage(imgMean);
        // int xM = imgMean.getWidth();
        // int yM = imgMean.getHeight();
        int zM = imgMean.getStackSize();

        for (int z = 0; z < zM; z++) {
            float[] valsImgMean = (float[]) imgMean.getStack().getProcessor(z + 1).getPixels();
            float[] valsImgStd = (float[]) imgStd.getStack().getProcessor(z + 1).getPixels();
            byte[][] valsImgSpec = new byte[nImg][];
            for (int i = 0; i < imgs.length; i++) {
                valsImgSpec[i] = (byte[]) imgs[i].getStack().getProcessor(z + 1).getPixels();
            }
            int nP = valsImgMean.length;
            double[] tabValuesForThisPixel = new double[nImg];
            for (int p = 0; p < nP; p++) {
                for (int i = 0; i < imgs.length; i++) {
                    tabValuesForThisPixel[i] = (double) (valsImgSpec[i][p] & 0xff);
                }
                double[] tabValues = VitimageUtils.statistics1D(tabValuesForThisPixel);
                valsImgMean[p] += tabValues[0];
                valsImgStd[p] += tabValues[1];
            }
        }
        return new ImagePlus[] { imgMean, imgStd };
    }
    /**
     * meanAndStdOfImageArrayFloatFreshNew:
     * - Calculates the mean and standard deviation images for a set of `ImagePlus` objects
     *   containing 32-bit float image data.
     * - Input: 
     *   - `imgs`: An array of `ImagePlus` objects to process.
     * - Output: An array of two `ImagePlus` objects (mean and standard deviation images).
     */
    public static ImagePlus[] meanAndStdOfImageArrayFloatFreshNew(ImagePlus[] imgs) {
        int nImg = imgs.length;
        ImagePlus imgMean = imgs[0];
        imgMean = VitimageUtils.nullImage(imgMean);
        ImagePlus imgStd = VitimageUtils.nullImage(imgMean);
        // int xM = imgMean.getWidth();
        // int yM = imgMean.getHeight();
        int zM = imgMean.getStackSize();

        for (int z = 0; z < zM; z++) {
            float[] valsImgMean = (float[]) imgMean.getStack().getProcessor(z + 1).getPixels();
            float[] valsImgStd = (float[]) imgStd.getStack().getProcessor(z + 1).getPixels();
            float[][] valsImgSpec = new float[nImg][];
            for (int i = 0; i < imgs.length; i++) {
                valsImgSpec[i] = (float[]) imgs[i].getStack().getProcessor(z + 1).getPixels();
            }
            int nP = valsImgMean.length;
            double[] tabValuesForThisPixel = new double[nImg];
            for (int p = 0; p < nP; p++) {
                for (int i = 0; i < imgs.length; i++) {
                    tabValuesForThisPixel[i] = (double) (valsImgSpec[i][p]);
                }
                double[] tabValues = VitimageUtils.statistics1D(tabValuesForThisPixel);
                valsImgMean[p] = (float) tabValues[0];
                valsImgStd[p] = (float) tabValues[1];
            }
        }
        return new ImagePlus[] { imgMean, imgStd };
    }


    // public static void main(String[] args) {

    //     ImageJ ij = new ImageJ();
    //     ImagePlus[] atlas = computeAverageAndStdForSpecificVariety(1);
    //     atlas[0].show();
    // }

}
