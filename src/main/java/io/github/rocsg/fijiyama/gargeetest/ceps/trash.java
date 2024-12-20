package io.github.rocsg.fijiyama.gargeetest.ceps;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.Plot;
import ij.process.ImageProcessor;
import ij.plugin.frame.RoiManager;
import io.github.rocsg.fijiyama.common.VitimageUtils;

public class trash {
    public static String specimen = GeneralUtils.getXRSpecimenNameWithIndex(0);
    public static String[] specimenList = GeneralUtils.getSpecimenListXR();
    public static String year1 = "2022";
    public static String year2 = "2023";

    
    public static void main(String args[]){
        ImageJ ij = new ImageJ();

        for (int i = 0; i < 3/*specimenList.length*/; i++) {
            String imgPathR = "/home/phukon/Desktop/Results_2024/Cep_Monitoring/XR/4_Registered_Images/8bit_Registered_2022_2023_2024/CEP_"
                    + specimenList[i] + "_" + year1 + "_XR.tif";// Red channel
            String imgPathG = "/home/phukon/Desktop/Results_2024/Cep_Monitoring/XR/4_Registered_Images/8bit_Registered_2022_2023_2024/CEP_"
                    + specimenList[i] + "_" + year2 + "_XR.tif";// Green channel

            ImagePlus imgR = IJ.openImage(imgPathR);
            ImagePlus imgG = IJ.openImage(imgPathG);
            // ImagePlus imgComposite = VitimageUtils.compositeNoAdjustOf(imgR, imgG) ;
            ImagePlus[] jointHistograms = computeAndDisplayJointHistogramsWithNoBinSize(imgR, imgG);

            jointHistograms[0].show();
            jointHistograms[1].show();
            jointHistograms[2].show();

            // IJ.save(jointHistograms[0],
            //         "/home/phukon/Desktop/Results_2024/Cep_Monitoring/XR/5_Histogram/jointHistogram/" + specimenList[i]
            //                 + "_" + year1 + "_" + year2 + "_Probability.tif");
            // IJ.save(jointHistograms[1],"/home/phukon/Desktop/jointHistoandROI/"+specimenList[i]+"_"+year1+"_"+year2+"_LogProbability.tif");
            // System.out.println("Saved
            // to"+"/home/phukon/Desktop/jointHistoandROI/"+specimenList[i]+"_"+year1+"_"+year2+"_LogProbability.tif");
            // IJ.save(jointHistograms[2],"/home/phukon/Desktop/jointHistoandROI/"+specimenList[i]+"_"+year1+"_"+year2+"_MutualProbability.tif");
            // System.out.println("Saved
            // to"+"/home/phukon/Desktop/jointHistoandROI/"+specimenList[i]+"_"+year1+"_"+year2+"_MutualProbability.tif");
        }
        System.out.println("fini");

    }

    public static ImagePlus makeSegmentationBasedOnRoiInputFromTheUser(ImagePlus imgR, ImagePlus imgG, ImagePlus histo,
            ImagePlus histoLog, String specimen) {

        // Ask for three rois
        RoiManager roiManager = RoiManager.getRoiManager();
        int nbAreas = 4;
        boolean finished = false;

        do {
            VitimageUtils.waitFor(1000);
            if (roiManager.getCount() == nbAreas)
                finished = true;
            System.out.println("Waiting " + nbAreas + ". Current number=" + roiManager.getCount());
        } while (!finished);

        Roi roi1 = roiManager.getRoi(0); // index for first polygon
        Roi roi2 = roiManager.getRoi(1); // index for second polygon
        Roi roi3 = roiManager.getRoi(2); // index for third polygon
        Roi roi4 = roiManager.getRoi(3); // index for third polygon

        // save roi
        roiManager.runCommand("Save", "/home/phukon/Desktop/jointHistoandROI/roi_" + specimen + ".zip");

        // Get access to original images (imgR and imgG) as dataR[] and dataG[]
        int X = imgR.getWidth();
        int Y = imgR.getHeight();
        int Z = imgR.getStackSize();

        // Initialize arrays to store pixel data
        byte[] dataR;
        byte[] dataG;

        // Create a destination image (float), of the size of the original image
        ImagePlus segmentation = VitimageUtils.nullImage(imgG);
        segmentation = VitimageUtils.convertByteToFloatWithoutDynamicChanges(segmentation);

        // Take access to the pixels of segmentation as dataSegmentation[]
        for (int z = 0; z < Z; z++) {
            float[] dataSegmentation = (float[]) segmentation.getStack().getProcessor(z + 1).getPixels();
            // For each pixl (x,y,z)
            // Collect the value of imgR and imgG. It yields equivalent coordinates within
            // the histogram (x,y)
            dataR = (byte[]) imgR.getStack().getProcessor(z + 1).getPixels();
            dataG = (byte[]) imgG.getStack().getProcessor(z + 1).getPixels();
            for (int y = 0; y < Y; y++) {
                for (int x = 0; x < X; x++) {
                    // Collect the pixel values of imgR and imgG
                    int valR = ((int) (dataR[x + X * y] & 0xff));
                    int valG = ((int) (dataG[x + X * y] & 0xff));

                    int xHisto = valR;
                    int yHisto = 255 - valG;

                    // For each roi, check if the point (x,y) is contained by the Roi
                    // Check if the pixel is contained by ROI 1
                    if (roi1.contains(xHisto, yHisto)) {
                        dataSegmentation[x + X * y] = 1;
                    }
                    // Check if the pixel is contained by ROI 2
                    else if (roi2.contains(xHisto, yHisto)) {
                        dataSegmentation[x + X * y] = 2;
                    }
                    // Check if the pixel is contained by ROI 3
                    else if (roi3.contains(xHisto, yHisto)) {
                        dataSegmentation[x + X * y] = 3;
                    }
                    // Check if the pixel is contained by ROI 4
                    else if (roi4.contains(xHisto, yHisto)) {
                        dataSegmentation[x + X * y] = 4;
                    }
                    // If the pixel is not contained by any ROI
                    else {
                        dataSegmentation[x + X * y] = 0;
                    }
                }
            }

        }
        segmentation.setDisplayRange(0, 5);

        return segmentation;
    }

    public static ImagePlus[] computeAndDisplayJointHistogramsWithNoBinSize(ImagePlus imgR, ImagePlus imgG) {
        ImagePlus imgProba = ij.gui.NewImage.createImage("Proba", 256, 256, 1, 32, ij.gui.NewImage.FILL_BLACK);
        ImagePlus imgLogProba = ij.gui.NewImage.createImage("Log-Proba", 256, 256, 1, 32, ij.gui.NewImage.FILL_BLACK);
        ImagePlus imgMutualProba = ij.gui.NewImage.createImage("Mutual Probability", 256, 256, 1, 32,
                ij.gui.NewImage.FILL_BLACK);
        float[] dataProba = (float[]) (imgProba.getStack().getProcessor(1).getPixels());
        float[] dataLogProba = (float[]) (imgLogProba.getStack().getProcessor(1).getPixels());
        float[] dataMutualProba = (float[]) (imgMutualProba.getStack().getProcessor(1).getPixels());

        // Normally we should verify that they have the same size
        int X = imgR.getWidth();
        int Y = imgR.getHeight();
        int Z = imgR.getStackSize();
        int Npix = X * Y * Z;
        double fuzzyValue = Math.log(1 / Npix) - 5;

        // We know that values are between 0 and 255
        int[][] jointHistogramSumOfPixels = new int[256][256];
        double[][] jointHistogramProbability = new double[256][256];
        double[][] jointHistogramLogProbability = new double[256][256];
        double[] marginalDistributionX = new double[256];
        double[] marginalDistributionY = new double[256];
        double[][] mutualProba = new double[256][256];

        int[] tabFirstNoFuck = new int[256];
        // Gain access to pixels of image
        byte[] dataR;
        byte[] dataG;
        for (int z = 0; z < Z; z++) {
            dataR = (byte[]) (imgR.getStack().getProcessor(z + 1).getPixels());
            dataG = (byte[]) (imgG.getStack().getProcessor(z + 1).getPixels());

            for (int x = 0; x < X; x++) {
                for (int y = 0; y < Y; y++) {
                    int val = (int) (dataR[x + X * y] & 0xff);
                    tabFirstNoFuck[val]++;
                    int valR = ((int) (dataR[x + X * y] & 0xff));
                    int valG = ((int) (dataG[x + X * y] & 0xff));
                    jointHistogramSumOfPixels[valR][valG] += 1;
                    
                }
            }
        }
        for (int i = 1; i < 255; i++)
            //System.out.println((tabFirstNoFuck[i]) * 2.0 / (tabFirstNoFuck[i - 1] + tabFirstNoFuck[i + 1]));
        // Normalize for getting a probability
        for (int r = 0; r < (256); r++)
            for (int g = 0; g < (256); g++) {
                jointHistogramProbability[r][g] = jointHistogramSumOfPixels[r][g] * 1.0 / Npix;

                // write things in the proba image
                dataProba[r + 256 * ((256 - 1) - g)] = (float) jointHistogramProbability[r][g];

                if (jointHistogramSumOfPixels[r][g] != 0) {
                    jointHistogramLogProbability[r][g] = Math.log(jointHistogramSumOfPixels[r][g] * 1.0 / Npix);
                    // Marginal distribution calculation for mutual information
                    marginalDistributionX[r] += jointHistogramProbability[r][g];
                    marginalDistributionY[g] += jointHistogramProbability[r][g];

                    // Calculate mutual information
                    if (jointHistogramSumOfPixels[r][g] != 0) {
                        double pxy = jointHistogramProbability[r][g];
                        double px = marginalDistributionX[r];
                        double py = marginalDistributionY[g];

                        mutualProba[r][g] = pxy * Math.log(pxy / (px * py));
                    }
                } else {
                    jointHistogramLogProbability[r][g] = fuzzyValue;
                }

                // write things in the LOGPROBA image
                dataLogProba[r + 256 * ((256 - 1) - g)] = (float) jointHistogramLogProbability[r][g];

                // write things in MUTUALPROBA image
                // dataMutualProba[r+256*(255-g)]=(float)(jointHistogramProbability[r][g]
                // *(jointHistogramLogProbability[r][g])/(marginalDistributionX[r] *
                // marginalDistributionY[g]));
                dataMutualProba[r + 256 * ((256 - 1) - g)] = (float) (mutualProba[r][g]);

            }
        
        
        double[] regressionLine = calcRegressionLine(jointHistogramLogProbability[0], jointHistogramLogProbability[1]);
        

        // // Draw lines on imgLogProba
        drawRegressionLine(imgLogProba, regressionLine[0], regressionLine[1]);

        // // Draw lines on imgMutualProba
        drawRegressionLine(imgMutualProba, regressionLine[0], regressionLine[1]);

        // imgProba.show();
        // imgLogProba.show();
        // imgMutualProba.show();
        
        
        
        IJ.run(imgProba, "Fire", "");
        IJ.run(imgLogProba, "Fire", "");
        IJ.run(imgMutualProba, "Fire", "");

        imgProba.setDisplayRange(10e-14, 10e-5);
        imgLogProba.setDisplayRange(-18, 4);
        imgMutualProba.setDisplayRange(10e-14, 10e-5);
        return new ImagePlus[] { imgProba, imgLogProba, imgMutualProba };
    }


    public static double[] calcRegressionLine(double[] x, double[] y){

        double a = 0;
        double b = 0;
        
        double xMean = calcMean(x);
        double yMean = calcMean(y);

        double numerator = 0;
        double denominator = 0;

        for(int i = 0;i<255;i++){
            numerator += (x[i] - xMean) * (yMean - y[255 - i]);
            denominator += Math.pow((x[i] - xMean),2);
        }
        a = numerator/denominator;
        b = yMean - a * xMean;
        System.out.println("Slope: "+a);
        System.out.println("Intercept: "+b);
        return new double[]{a,b};
    }

    private static double calcMean(double[] data) {
        double sum = 0;
        for (double value : data) {
            sum += value;
        }
        return sum / data.length;
    }

     // Method to draw a straight line on an ImagePlus object
     private static void drawRegressionLine(ImagePlus image, double a, double b) {
        ImageProcessor ip = image.getProcessor();
        int width = ip.getWidth();
        int height = ip.getHeight();

        int x1 = 0; // Starting x-coordinate (usually 0 or another appropriate value)
        int x2 = width - 1; // Ending x-coordinate (width-1 or another appropriate value)
        int y1 = width - (int) Math.round(a * x1 + b); // y-coordinate corresponding to x1
        int y2 = width - (int) Math.round(a * x2 + b); // y-coordinate corresponding to x2

        // Ensure the line is within the bounds of the image
        y1 = Math.max(0, Math.min(height - 1, y1));
        y2 = Math.max(0, Math.min(height - 1, y2));

        // Draw the line on the image
        ip.drawLine(x1, y1, x2, y2);
    }

}
