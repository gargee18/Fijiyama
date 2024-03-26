/*
 * 
 * 
 * 
 */

 
package io.github.rocsg.fijiyama.gargeetest;


import org.apache.commons.math3.transform.TransformUtils;

import java.util.ArrayList;
import java.util.Random;

import org.itk.simple.ComposeImageFilter;
import org.itk.simple.DisplacementFieldJacobianDeterminantFilter;
import org.itk.simple.DisplacementFieldTransform;
import org.itk.simple.Image;
import org.itk.simple.ImageFileReader;
import org.itk.simple.ImageFileWriter;
import org.itk.simple.MultiplyImageFilter;
import org.itk.simple.ResampleImageFilter;
import org.itk.simple.SimpleITK;
import org.itk.simple.Transform;
import org.itk.simple.VectorDouble;
import org.itk.simple.VectorIndexSelectionCastImageFilter;
import org.scijava.java3d.Transform3D;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.ChannelSplitter;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.plugin.HyperStackConverter;
import ij.plugin.RGBStackMerge;
import ij.process.ImageProcessor;
import ij.process.StackConverter;
import io.github.rocsg.fijiyama.common.ItkImagePlusInterface;
import io.github.rocsg.fijiyama.common.VitiDialogs;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.registration.ItkTransform;
import math3d.Point3d;
import vib.FastMatrix;

public class TransformUsingMetadata {
    
    static String mainDir="/home/phukon/Desktop/MRI/02_ceps/2022-03_CEPS_suivi_clinique/1_STACK_16-bits/";
    static String specimen = "318";
       
    public static void main(String[] args) {
        ImageJ ij=new ImageJ();
        System.out.println("Starting test...");
        twoMaricesAlignmentWithMetadata();
    }

    // public static String pathT1Data(){
    //     return mainDir+"/CEP_"+specimen+"_2022_MRI_T1.tif";
    // }

    // public static String pathT2Data(){
    //     return mainDir+"/CEP_"+specimen+"_2022_MRI_T2.tif";
    // }
    
    public static void twoMaricesAlignmentWithMetadata(){
        double[] tRefX=  {-0.649988,0.75979,0.0153097};
        double[] tRefY= {-0.0620298,-0.0329653,-0.99753}; 
        double []vzz= TransformUtils.vectorialProduct(tRefX,tRefY);
        double[] tMovX=  {-0.649987,0.759791,0.0153221};
        double[] tMovY= {-0.0620379,-0.0329559,-0.99753};

        ItkTransform tT1= ItkTransform.itkTransformFromDICOMVectors(null, null, null, null);
    }
    
    
}
