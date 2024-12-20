/*
 * This Java class, ImageStandardization, provides methods for standardizing and normalizing image data.
 * It uses the ImageJ library to handle image processing tasks such as normalizing pixel values to a specified range,
 * converting between different image bit depths, and calculating statistical properties of the image.
 * 
 * Methods:
 * - meanPixelIntensityCalc: Calculates the mean pixel intensity of a given image.
 * - normalizeImageMinMax8bit: Normalizes an image to 8-bit depth using min-max scaling.
 * - normalizeImageMinMax16bit: Normalizes a 16-bit image using min-max scaling.
 * - convertTo8bit: Converts a 16-bit image to 8-bit depth.
 * - normalizeImageMinMax32bit: Normalizes an image to 32-bit depth using min-max scaling.
 * - normalizeImageZScore: Normalizes an image using Z-score standardization.
 * - calculateMean: Calculates the mean pixel intensity for each slice in the image stack.
 * - calculateStdDev: Calculates the standard deviation of pixel intensities for each slice in the image stack.
 * 
 * This class utilizes the ImageJ library for image manipulation and processing tasks.
 * 
 * Input Parameters:
 * - img: An instance of ImagePlus representing the image to be processed.
 * - img16bit: An instance of ImagePlus representing a 16-bit image.
 * 
 * Output Parameters:
 * - normImg2022, normImg2023, normImg2024: Normalized images for different years.
 * - mean22, mean23, mean24: Mean pixel intensities for images from different years.
 * - imgNorm: A normalized version of the input image.
 * - img8bit: A converted 8-bit image from a 16-bit image.
 * - mean: Array of mean pixel intensities for each slice.
 * - std: Array of standard deviations for pixel intensities for each slice.
 * 
 * (This is an experimental class; in the end did not use this method for normalizing data)
 */

package io.github.rocsg.fijiyama.gargeetest.ceps;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.NewImage;
import io.github.rocsg.fijiyama.common.VitimageUtils;

public class ImageStandardization {
    public static void main(String[] args) {
        ImageJ ij = new ImageJ();
        // Load images
        ImagePlus img2022 = IJ.openImage(
                "/home/phukon/Desktop/registration/CEP_2022_XR_2023_XR_2024_XR/cep_1266A/raw/1266A_2022_crop_sub_z.tif");
        ImagePlus img2023 = IJ.openImage(
                "/home/phukon/Desktop/CEP_Registration/CEP_2022_XR_to_2023_XR_LOWRES/1266A_img_moving_after_registration.tif");
        ImagePlus img2024 = IJ
                .openImage("/home/phukon/Desktop/xray_2023_2024/cep_1266A/res/1266A_TransformedImage.tif");
        // Normalize images to 8-bit depth
        ImagePlus normImg2022 = normalizeImageMinMax8bit(img2022);
        normImg2022.show();
        ImagePlus normImg2023 = normalizeImageMinMax8bit(img2023);
        normImg2023.show();
        ImagePlus normImg2024 = normalizeImageMinMax8bit(img2024);
        normImg2024.show();
        // Calculate and print mean pixel intensities
        double mean22 = meanPixelIntensityCalc(normImg2022);
        double mean23 = meanPixelIntensityCalc(normImg2023);
        double mean24 = meanPixelIntensityCalc(normImg2024);
        System.out.println("Mean for year 2022: " + mean22 + "\n" + "Mean for year 2023: " + mean23 + "\n"
                + "Mean for year 2024: " + mean24);
    }
    // Calculates the mean pixel intensity for a given image
    public static double meanPixelIntensityCalc(ImagePlus img) {
        double sum = 0.0;
        int count = 0;
        for (int z = 0; z < img.getStackSize(); z++) {
            byte[] data = (byte[]) (img.getStack().getProcessor(z + 1).getPixels());
            for (int x = 0; x < img.getWidth(); x++) {
                for (int y = 0; y < img.getHeight(); y++) {
                    int val = (data[x + img.getWidth() * y] & 0xff);
                    if (val != 0) {
                        sum += val;
                        count++;
                    }
                }
            }
        }
        if (count == 0) {
            return 0; // Avoid division by zero
        } else {
            return sum / count;
        }
    }
    // Normalizes an image to 8-bit depth using min-max scaling
    public static ImagePlus normalizeImageMinMax8bit(ImagePlus img) {
        // Normally we should verify that they have the same size
        int X = img.getWidth();
        int Y = img.getHeight();
        int Z = img.getStackSize();
        ImagePlus imgNorm = ij.gui.NewImage.createImage("Normalized Image", X, Y, Z, 8, ij.gui.NewImage.FILL_BLACK);
        VitimageUtils.adjustImageCalibration(imgNorm, img);
        int min = 255;
        int max = 0;
        int val = 0;
        // Find the min and max pixel values
        for (int z = 0; z < Z; z++) {
            byte[] data = (byte[]) (img.getStack().getProcessor(z + 1).getPixels());
            for (int x = 0; x < X; x++) {
                for (int y = 0; y < Y; y++) {
                    val = (data[x + X * y] & 0xff);
                    if (val < min) {
                        min = val;
                    }
                    if (val > max) {
                        max = val;
                    }
                }
            }
        }
        System.out.println(min);
        System.out.println(max);
        int range = max - min;
        double thresholdPercentage = 0.15;
        int threshold = min + (int) (range * thresholdPercentage);
        System.out.println(threshold);
        // Normalize the pixel values and set background to 0
        for (int z = 0; z < Z; z++) {
            byte[] imgData = (byte[]) (imgNorm.getStack().getProcessor(z + 1).getPixels());
            byte[] data = (byte[]) (img.getStack().getProcessor(z + 1).getPixels());
            for (int x = 0; x < X; x++) {
                for (int y = 0; y < Y; y++) {
                    val = (data[x + X * y] & 0xff);
                    if (val < threshold) {
                        imgData[x + X * y] = 0;
                    } else {
                        imgData[x + X * y] = (byte) ((float) (val - min) / (max - min) * 255);
                    }
                }
            }
        }
        return imgNorm;
    }
    // Normalizes a 16-bit image using min-max scaling
    public static ImagePlus normalizeImageMinMax16bit(ImagePlus img) {
        int X = img.getWidth();
        int Y = img.getHeight();
        int Z = img.getStackSize();

        // Create an output image with 16-bit depth
        ImagePlus imgNorm = NewImage.createImage("Normalized Image", X, Y, Z, 16, NewImage.FILL_BLACK);
        VitimageUtils.adjustImageCalibration(imgNorm, img);
        int min = 65535; // Maximum value for 16-bit
        int max = 0;
        int val;

        // Find the min and max pixel values
        for (int z = 0; z < Z; z++) {
            short[] data = (short[]) (img.getStack().getProcessor(z + 1).getPixels());
            for (int i = 0; i < data.length; i++) {
                val = data[i] & 0xffff; // Convert signed short to unsigned int
                if (val < min) {
                    min = val;
                }
                if (val > max) {
                    max = val;
                }
            }
        }
        System.out.println(min);
        System.out.println(max);
        double thresholdPercentage = 0.1;
        int range = max - min;
        int threshold = min + (int) (range * thresholdPercentage);
        System.out.println(threshold);
        // Normalize the pixel values and set background to 0
        for (int z = 0; z < Z; z++) {
            short[] imgData = (short[]) (imgNorm.getStack().getProcessor(z + 1).getPixels());
            short[] data = (short[]) (img.getStack().getProcessor(z + 1).getPixels());
            for (int i = 0; i < data.length; i++) {
                val = data[i] & 0xffff; // Convert signed short to unsigned int
                if (val < threshold) {
                    imgData[i] = 0;
                } else {
                    imgData[i] = (short) ((float) val * 1 + (-1024));
                }
            }
        }
        return imgNorm;
    }
    // Converts a 16-bit image to 8-bit depth
    public static ImagePlus convertTo8bit(ImagePlus img16bit) {
        int X = img16bit.getWidth();
        int Y = img16bit.getHeight();
        int Z = img16bit.getStackSize();
        // Create a new 8-bit image for the converted result
        ImagePlus img8bit = NewImage.createByteImage("Converted to 8-bit", X, Y, Z, NewImage.FILL_BLACK);
        VitimageUtils.adjustImageCalibration(img8bit, img16bit);

        // Convert each slice from 16-bit to 8-bit
        for (int z = 0; z < Z; z++) {
            short[] imgData16bit = (short[]) (img16bit.getStack().getProcessor(z + 1).getPixels());
            byte[] imgData8bit = (byte[]) (img8bit.getStack().getProcessor(z + 1).getPixels());
            for (int i = 0; i < imgData16bit.length; i++) {
                // Scale 16-bit value to 8-bit range
                int value8bit = (int) (imgData16bit[i] / 65535.0 * 255);
                imgData8bit[i] = (byte) value8bit;
            }
        }
        return img8bit;
    }
    // Normalizes an image to 32-bit depth using min-max scaling
    
    public static ImagePlus normalizeImageMinMax32bit(ImagePlus img) {
        // Normally we should verify that they have the same size
        int X = img.getWidth();
        int Y = img.getHeight();
        int Z = img.getStackSize();
        ImagePlus imgNorm = ij.gui.NewImage.createImage("Normalized Image", X, Y, Z, 32, ij.gui.NewImage.FILL_BLACK);
        int min = 255;
        int max = 0;
        int val = 0;
        // Find the min and max pixel values
        for (int z = 0; z < Z; z++) {
            byte[] data = (byte[]) (img.getStack().getProcessor(z + 1).getPixels());
            for (int x = 0; x < X; x++) {
                for (int y = 0; y < Y; y++) {
                    val = (data[x + X * y] & 0xff);
                    if (val < min) {
                        min = val;
                    }
                    if (val > max) {
                        max = val;
                    }
                }
            }
        }
        System.out.println(min);
        System.out.println(max);
        int range = max - min;
        double thresholdPercentage = 0.15;
        int threshold = min + (int) (range * thresholdPercentage);
        System.out.println(threshold);
        // Normalize the pixel values and set background to 0
        for (int z = 0; z < Z; z++) {
            float[] imgData = (float[]) (imgNorm.getStack().getProcessor(z + 1).getPixels());
            byte[] data = (byte[]) (img.getStack().getProcessor(z + 1).getPixels());
            for (int x = 0; x < X; x++) {
                for (int y = 0; y < Y; y++) {
                    val = (data[x + X * y] & 0xff);
                    if (val < threshold) {
                        imgData[x + X * y] = 0;
                    } else {
                        imgData[x + X * y] = (float) (val - min) / (max - min);
                    }
                }
            }
        }
        return imgNorm;
    }
    // Normalizes an image using Z-score standardization
    public static ImagePlus normalizeImageZScore(ImagePlus img) {
        // Normally we should verify that they have the same size
        int X = img.getWidth();
        int Y = img.getHeight();
        int Z = img.getStackSize();
        ImagePlus imgNorm = ij.gui.NewImage.createImage("Normalized Image", X, Y, Z, 32, ij.gui.NewImage.FILL_BLACK);
        imgNorm.show();
        double[] mean = calculateMean(img);
        double[] std = calculateStdDev(img, mean);
        // Gain access to pixels of image
        for (int z = 0; z < Z; z++) {
            float[] imgData = (float[]) (imgNorm.getStack().getProcessor(z + 1).getPixels());
            byte[] data = (byte[]) (img.getStack().getProcessor(z + 1).getPixels());
            for (int x = 0; x < X; x++) {
                for (int y = 0; y < Y; y++) {
                    int val = ((int) (data[x + X * y] & 0xff));
                    double normalizedVal = (val - mean[z]) / std[z];
                    imgData[x + X * y] = (float) normalizedVal;
                }
            }
        }
        return imgNorm;
    }
    // Calculates the mean pixel intensity for each slice in the image stack
    private static double[] calculateMean(ImagePlus img) {
        int depth = img.getStackSize();
        double[] mean = new double[depth];
        int npix = img.getWidth() * img.getHeight();
        for (int z = 0; z < depth; z++) {
            double sum = 0;
            byte[] data = (byte[]) img.getStack().getProcessor(z + 1).getPixels();
            for (int y = 0; y < img.getHeight(); y++) {
                for (int x = 0; x < img.getWidth(); x++) {
                    sum += (double) (data[x + img.getWidth() * y] & 0xff); // Convert byte to int (0-255)
                }
            }
            mean[z] = sum / npix;
            System.out.println("Mean at slice " + z + " =" + mean[z]);
        }
        return mean;
    }
    // Calculates the standard deviation of pixel intensities for each slice in the image stack
    private static double[] calculateStdDev(ImagePlus img, double[] mean) {
        int depth = img.getStackSize();
        double[] std = new double[depth];
        int npix = img.getWidth() * img.getHeight();
        for (int z = 0; z < depth; z++) {
            double sum = 0;
            byte[] data = (byte[]) img.getStack().getProcessor(z + 1).getPixels();
            for (int y = 0; y < img.getHeight(); y++) {
                for (int x = 0; x < img.getWidth(); x++) {
                    double diff = (double) (data[x + img.getWidth() * y] & 0xff) - mean[z]; // Convert byte to int
                                                                                            // (0-255)
                    sum += diff * diff;
                }
            }
            std[z] = Math.sqrt(sum / npix);
            System.out.println("Std Dev at slice " + z + " =" + std[z]);
        }
        return std;
    }
}
