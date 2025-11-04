/**
 * This class, `step_5_AtlasBuilding`, is part of the image processing pipeline
 * for analyzing plant specimen hyperstacks. It provides methods to compute 
 * average and standard deviation images based on timepoint differences (e.g., T1-T0, T2-T0, T2-T1) and
 * for various conditions and plant varieties. 
 */


package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D;

import java.util.ArrayList;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;


public class Step_6_MeanIndividual implements PipelineStep {

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

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }

    public static void main(String[] args) throws Exception {

        Specimen spec= new Specimen("B_202");
        ImageJ ij=new ImageJ();//Needed for testing
        new  Step_6_MeanIndividual().execute(spec,true);
        System.out.println("Saved!");
    }
    
    public void execute(Specimen specimen,boolean testing) throws Exception {
        String[] timestamps = Config.timestamps;
        
        int step = 4;
        ImagePlus[] allVarCT = computeAverageAndStdAllVar(cond_CONTROL, step);
        allVarCT[0].setDisplayRange(0,1);
        VitimageUtils.setLutToFire(allVarCT[0]);
        allVarCT[1].setDisplayRange(0.1,1);
        VitimageUtils.setLutToFire(allVarCT[1]);
        IJ.saveAsTiff( allVarCT[0]    ,Config.getPathToCartesianAtlas()+"/01_FullPop/mean_all_var_CT_"+timestamps[step -1]+".tif");
        IJ.saveAsTiff( allVarCT[1]    ,Config.getPathToCartesianAtlas()+"/01_FullPop/std_all_var_CT_"+timestamps[step -1]+".tif");
        ImagePlus[] allVarPCH = computeAverageAndStdAllVar(cond_PCH, step);
        allVarPCH[0].setDisplayRange(0,1);
        VitimageUtils.setLutToFire(allVarPCH[0]);
        allVarPCH[1].setDisplayRange(0,1);
        VitimageUtils.setLutToFire(allVarPCH[1]);
        IJ.saveAsTiff( allVarPCH[0]    ,Config.getPathToCartesianAtlas()+"/01_FullPop/mean_all_var_PCH_"+timestamps[step -1]+".tif");
        IJ.saveAsTiff( allVarPCH[1]    ,Config.getPathToCartesianAtlas()+"/01_FullPop/std_all_var_PCH_"+timestamps[step -1]+".tif");
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

    public static void computeAndSaveDiffMaps( int initialTime, int finalTime) {
        String[] timestamps = Config.timestamps;
        // int initialTime = 0; // 0 = J_001
        // int finalTime = 1; // 1 = J_029, 2 = J_077, 3 = J_141 
        int initTimeFrame = initialTime+1;
        int finTimeFrame = finalTime+1;
        ImagePlus[] atlasCtAllVar = computeAverageAndStdIndividualAtT1MinusT0( cond_CONTROL, initTimeFrame, finTimeFrame);
        atlasCtAllVar[0].setDisplayRange(-0.2,1);
        VitimageUtils.setLutToFire(atlasCtAllVar[0]);
        atlasCtAllVar[1].setDisplayRange(-0.2,1);
        VitimageUtils.setLutToFire(atlasCtAllVar[1]);
        IJ.saveAsTiff( atlasCtAllVar[0]    ,Config.getPathToPolarAtlas()+"/03_DiffMap/mean_diff_UGNI_CT_"+timestamps[finalTime]+"_"+timestamps[initialTime]+".tif");
        IJ.saveAsTiff( atlasCtAllVar[1]    ,Config.getPathToPolarAtlas()+"/03_DiffMap/std_diff_UGNI_CT_"+timestamps[finalTime]+"_"+timestamps[initialTime]+".tif");
        ImagePlus[] atlasPchAllVar = computeAverageAndStdIndividualAtT1MinusT0( cond_PCH, initTimeFrame, finTimeFrame);
        atlasPchAllVar[0].setDisplayRange(-0.2,1);
        VitimageUtils.setLutToFire(atlasPchAllVar[0]);
        atlasPchAllVar[1].setDisplayRange(-0.2,1);
        VitimageUtils.setLutToFire(atlasPchAllVar[1]);
        IJ.saveAsTiff( atlasPchAllVar[0]    ,Config.getPathToPolarAtlas()+"/03_DiffMap/mean_diff_UGNI_PCH_"+timestamps[finalTime]+"_"+timestamps[initialTime]+".tif");
        IJ.saveAsTiff( atlasPchAllVar[1]    ,Config.getPathToPolarAtlas()+"/03_DiffMap/std_diff_UGNI_PCH_"+timestamps[finalTime]+"_"+timestamps[initialTime]+".tif");
        

    }


    public static ImagePlus[] computeAverageAndStdIndividualAtT1MinusT0(int condition, int initialTime, int finalTime) {
        ArrayList<ImagePlus> stacks = new ArrayList<ImagePlus>();
        for (int cond = condition; cond <= condition; cond++) {
            for (int var = 3; var < 4; var++) {
                String[] spec = getSpecimensName(cond, var);
                for (int i = 0; i < spec.length; i++) {
                    System.out.println("Specimen: " + spec[i] + " Condition: " + cond + " Variety: " + var);
                    ImagePlus img = IJ.openImage(Config.mainDir +"/Processing/03_PolarTransform/"+spec[i]+"_GeneralizedPolarTransform.tif");
                    // openImage
                    System.out.println("Image: " + img);
                    System.out.println(getDirOfSpecimen(spec[i]));

                    ImagePlus imgT1 = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, initialTime, initialTime);
                    ImagePlus imgT2 = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, finalTime, finalTime);
        
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
    public static ImagePlus[] computeAverageAndStdIndividualAtT3MinusT0(int condition) {
        ArrayList<ImagePlus> stacks = new ArrayList<ImagePlus>();
        for (int cond = condition; cond <= condition; cond++) {
            for (int var = 0; var < 4; var++) {
                String[] spec = getSpecimensName(cond, var);
                for (int i = 0; i < spec.length; i++) {
                    System.out.println("Specimen: " + spec[i] + " Condition: " + cond + " Variety: " + var);
                    ImagePlus img = IJ.openImage(Config.mainDir +"/Processing/03_PolarTransform/"+spec[i]+"_GeneralizedPolarTransform.tif");
                    System.out.println("Image: " + img);
                    System.out.println(getDirOfSpecimen(spec[i]));

                    ImagePlus imgT0 = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, 0, 0);
                    ImagePlus imgT3 = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, 3, 3);
                    img = VitimageUtils.makeOperationBetweenTwoImages(imgT3, imgT0, 4, true);
                    imgT0 = null;
                    imgT3 = null;
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
    public static ImagePlus[] computeAverageAndStdIndividualAtT2MinusT0(int condition) {
        ArrayList<ImagePlus> stacks = new ArrayList<ImagePlus>();
        for (int cond = condition; cond <= condition; cond++) {
            for (int var = 0; var < 4; var++) {
                String[] spec = getSpecimensName(cond, var);
                for (int i = 0; i < spec.length; i++) {
                    System.out.println("Specimen: " + spec[i] + " Condition: " + cond + " Variety: " + var);
                    ImagePlus img = IJ.openImage(Config.mainDir +"/Processing/03_PolarTransform/"+spec[i]+"_GeneralizedPolarTransform.tif");
                    System.out.println("Image: " + img);
                    System.out.println(getDirOfSpecimen(spec[i]));

                    ImagePlus imgT0 = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, 0, 0);
                    ImagePlus imgT2 = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, 2    , 2);
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
     * computeAverageAndStdForSpecificVariety:
     * - Processes hyperstack images for a specific variety (e.g., UGNI) and computes 
     *   the average and standard deviation for the selected timepoints.
     * - Input: 
     *   - `condition`: An integer representing the condition.
     * - Output: An array of two `ImagePlus` objects (mean and standard deviation images).
     */

     public static void computeAndSaveSpecificVariety(){
        String[] timestamps = Config.timestamps;

        for (int t = 1; t < 5; t++) {
            // ImagePlus[] atlasCtAllSpecVar = computeAverageAndStdForSpecificVariety(cond_CONTROL, var_UGNI, t);
            // IJ.saveAsTiff( atlasCtAllSpecVar[0]    ,Config.getPathToPolarAtlas()+"/02_SpecificVariety/mean_var_UGNI_CT_"+timestamps[t-1]+".tif");
            // IJ.saveAsTiff( atlasCtAllSpecVar[1]    ,Config.getPathToPolarAtlas()+"/02_SpecificVariety/std_var_UGNI_CT_"+timestamps[t-1]+".tif");
            // ImagePlus[] atlasPchAllSpecVar = computeAverageAndStdForSpecificVariety(cond_PCH, var_UGNI, t);
            // IJ.saveAsTiff( atlasPchAllSpecVar[0]    ,Config.getPathToPolarAtlas()+"/02_SpecificVariety/mean_var_UGNI_PCH_"+timestamps[t-1]+".tif");
            // IJ.saveAsTiff( atlasPchAllSpecVar[1]    ,Config.getPathToPolarAtlas()+"/02_SpecificVariety/std_var_UGNI_PCH_"+timestamps[t-1]+".tif");
            ImagePlus[] atlasCtAllSpecVar = computeAverageAndStdForSpecificVariety(cond_CONTROL, var_UGNI, t);
             atlasCtAllSpecVar[0].setDisplayRange(-0.2,2);
            VitimageUtils.setLutToFire(atlasCtAllSpecVar[0]);
            IJ.saveAsTiff( atlasCtAllSpecVar[0]    ,Config.getPathToPolarAtlas()+"/test_spec_var/mean_var_UGNI_CT_"+timestamps[t-1]+".tif");
            IJ.saveAsTiff( atlasCtAllSpecVar[1]    ,Config.getPathToPolarAtlas()+"/test_spec_var/std_var_UGNI_CT_"+timestamps[t-1]+".tif");
            ImagePlus[] atlasPchAllSpecVar = computeAverageAndStdForSpecificVariety(cond_PCH, var_UGNI, t);
            atlasPchAllSpecVar[1].setDisplayRange(-0.2,2);
            VitimageUtils.setLutToFire(atlasPchAllSpecVar[1]);
            IJ.saveAsTiff( atlasPchAllSpecVar[0]    ,Config.getPathToPolarAtlas()+"/test_spec_var/mean_var_UGNI_PCH_"+timestamps[t-1]+".tif");
            IJ.saveAsTiff( atlasPchAllSpecVar[1]    ,Config.getPathToPolarAtlas()+"/test_spec_var/std_var_UGNI_PCH_"+timestamps[t-1]+".tif");
        }
    }
    public static ImagePlus[] computeAverageAndStdForSpecificVariety(int condition, int var,  int step) {
        ArrayList<ImagePlus> stacks = new ArrayList<ImagePlus>();
        for (int cond = condition; cond <= condition; cond++) {
            // for (int var = 0; var < variety; var++) {
                String[] spec = getSpecimensName(cond, var);
                for (int i = 0; i < spec.length; i++) {
                    System.out.println("Specimen: " + spec[i] + " Condition: " + cond + " Variety: " + var);
                    ImagePlus img = IJ.openImage(Config.mainDir +"/Processing/03_PolarTransform/"+spec[i]+"_GeneralizedPolarTransform.tif");
                    // I would like to display if img is null, and then print the argument of
                    // openImage
                    System.out.println("Image: " + img);
                    System.out.println(Config.mainDir +"/Processing/03_PolarTransform/"+spec[i]+"_GeneralizedPolarTransform.tif");
                    ImagePlus imgVar = new Duplicator().run(img, 2, 2, 256, img.getNSlices() - 256, step, step);
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

    
    public static ImagePlus[] computeAverageAndStdAtDifferentT(int condition, int var,  int step) {
        ArrayList<ImagePlus> stacks = new ArrayList<ImagePlus>();
        for (int cond = condition; cond <= condition; cond++) {
            // for (int var = 0; var < 4; var++) {
                String[] spec = getSpecimensName(cond, var);
                for (int i = 0; i < spec.length; i++) {
                    System.out.println("Specimen: " + spec[i] + " Condition: " + cond + " Variety: " + var);
                    // ImagePlus img = IJ.openImage(getDirOfSpecimen(spec[i]) + "/hyperimage/" + spec[i] + "_Hyperstack.tif");
                    ImagePlus img = IJ.openImage(Config.mainDir +"/Processing/03_PolarTransform/"+spec[i]+"_GeneralizedPolarTransform.tif");
                    System.out.println("Image: " + img);
                    System.out.println(getDirOfSpecimen(spec[i]));
                    ImagePlus imgT = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, step, step);
                    stacks.add(imgT);

                }
            // }
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


    public static ImagePlus[] computeAverageAndStdAllVar(int condition, int step) {
        ArrayList<ImagePlus> stacks = new ArrayList<ImagePlus>();
        for (int cond = condition; cond <= condition; cond++) {
            for (int var = 0; var < 4; var++) {
                String[] spec = getSpecimensName(cond, var);
                for (int i = 0; i < spec.length; i++) {
                    System.out.println("Specimen: " + spec[i] + " Condition: " + cond + " Variety: " + var);
                    ImagePlus img = IJ.openImage(Config.mainDir + "Results/01_Hyperstack/" + spec[i] + "_Hyperstack.tif");
                    // ImagePlus img = IJ.openImage(Config.mainDir +"/Processing/03_PolarTransform/"+spec[i]+"_GeneralizedPolarTransform.tif");
                    System.out.println("Image: " + img);
                    System.out.println(getDirOfSpecimen(spec[i]));
                    ImagePlus imgT = new Duplicator().run(img, 1, 1, 256, img.getNSlices() - 256, step, step);
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

    
}
// int initialTime = 0; // 0 = J_001
        // int finalTime = 1; // 1 = J_029, 2 = J_077, 3 = J_141 
        // int initTimeFrame = initialTime+1;
        // int finTimeFrame = finalTime+1;
        // ImagePlus[] atlasCtAllVar = computeAverageAndStdIndividualAtT1MinusT0( cond_CONTROL, initTimeFrame, finTimeFrame);
        // atlasCtAllVar[0].setDisplayRange(-0.2,1);
        // VitimageUtils.setLutToFire(atlasCtAllVar[0]);
        // atlasCtAllVar[1].setDisplayRange(-0.2,1);
        // VitimageUtils.setLutToFire(atlasCtAllVar[1]);
        // IJ.saveAsTiff( atlasCtAllVar[0]    ,Config.getPathToPolarAtlas()+"/03_DiffMap/mean_CHARD_CT_"+timestamps[finalTime]+"_"+timestamps[initialTime]+".tif");
        // IJ.saveAsTiff( atlasCtAllVar[1]    ,Config.getPathToPolarAtlas()+"/03_DiffMap/std_CHARD_CT_"+timestamps[finalTime]+"_"+timestamps[initialTime]+".tif");
        // ImagePlus[] atlasPchAllVar = computeAverageAndStdIndividualAtT1MinusT0( cond_PCH, initTimeFrame, finTimeFrame);
        // atlasPchAllVar[0].setDisplayRange(-0.2,1);
        // VitimageUtils.setLutToFire(atlasPchAllVar[0]);
        // atlasPchAllVar[1].setDisplayRange(-0.2,1);
        // VitimageUtils.setLutToFire(atlasPchAllVar[1]);
        // IJ.saveAsTiff( atlasPchAllVar[0]    ,Config.getPathToPolarAtlas()+"/03_DiffMap/mean_CHARD_PCH_"+timestamps[finalTime]+"_"+timestamps[initialTime]+".tif");
        // IJ.saveAsTiff( atlasPchAllVar[1]    ,Config.getPathToPolarAtlas()+"/03_DiffMap/std_CHARD_PCH_"+timestamps[finalTime]+"_"+timestamps[initialTime]+".tif");