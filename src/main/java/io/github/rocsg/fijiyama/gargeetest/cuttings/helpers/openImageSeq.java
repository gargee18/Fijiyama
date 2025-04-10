package io.github.rocsg.fijiyama.gargeetest.cuttings.helpers;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import ij.plugin.FolderOpener;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;


public class openImageSeq {

    public static void main(String[] args) {
        ImageJ ij = new ImageJ();
        String[] timestamp = Config.timestamps;
        int N = timestamp.length;
        String pathToImage = "/home/phukon/Desktop/201_T1_T2/B201_HyperStack.tif";   
        for (int t = 1; t<N+1;t++){
            ImagePlus img = getSquenceFromAHypermap(pathToImage, 3,t);
            img.show();
        } 
    }


    /**
     * Opens an image from the specified path and extracts a sequence from a hypermap.
     * @param pathToImage The path to the image file.
     * @return The extracted sequence as an ImagePlus object.
     */
    public static ImagePlus getSquenceFromAHypermap(String pathToImage,  int channel, int time) {
        // Open the image from the given path
        ImagePlus img = IJ.openImage(pathToImage);
        // Initialize an empty ImagePlus object for the sequence
        ImagePlus seqT1 = new ImagePlus();

        // Check if the image is a hypermap with multiple channels or frames
        if (img.getNChannels() > 1 || img.getNFrames() > 1) {
            System.out.println("It's a hypermap");
            // Duplicate the third channel ( TR=2400 TE=12 )of the current frame to create a new sequence
            seqT1 = new Duplicator().run(img, channel, channel, 1, img.getNSlices(), time, time);
            seqT1.show();
        } else {
            // Print a message if the image is not a hypermap
            System.out.println("Not an image");
        }

        // Return the last extracted sequence
        return seqT1;
    }

    /**
     * Given a directory path, this method opens the sequence of images in that directory as a stack.
     * @param pathToDir The path to the directory containing the image sequence.
     * @return The opened image sequence as an ImagePlus stack.
     */
    public static ImagePlus getSequenceFromDir(String pathToDir){
        // Open the sequence of images in the given directory as a stack
        ImagePlus imp = FolderOpener.open(pathToDir, "");
        imp.show();
        return imp;
    }
    
}
