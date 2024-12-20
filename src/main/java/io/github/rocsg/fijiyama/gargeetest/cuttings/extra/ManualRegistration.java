package io.github.rocsg.fijiyama.gargeetest.cuttings.extra;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.fijiyamaplugin.Fijiyama_GUI;
import io.github.rocsg.fijiyama.fijiyamaplugin.RegistrationManager;
import io.github.rocsg.fijiyama.registration.ItkTransform;



// This class is intended to run batch manual registration for the cutting experiences.
//The number of manual registration to do is approximately 120, meaning there is a need to automate the gathering of files, opening of 3d universe, saving and 
//organizing of the transformations, etc.
public class ManualRegistration {
    static String testSpecimen = "B_240";   // An-Information-with-elements-separated-with-minus_while-different-fields-will-be-separated-by-underscores
    static String[] timestamps = new String[]{"J001", "J029", "J077", "J141"};
    static String mainDir = "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/Cuttings_MRI_registration/";

    public static String getPathToReferenceImage(String specimen,int step){
        System.out.println(mainDir+specimen+"/raw_subsampled/"+specimen+"_"+timestamps[step]+"_sub222.tif");
        return mainDir+specimen+"/raw_subsampled/"+specimen+"_"+timestamps[step]+"_sub222.tif";
        
    }

    public static String getPathToImage(String specimen,int step){
  //      System.out.println("/home/phukon/Desktop/Cuttings_MRI_registration/"+specimen+"/raw_subsampled/"+specimen+"_"+timestamps[step+1]+"_sub222.tif");
        return mainDir+specimen+"/raw_subsampled/"+specimen+"_"+timestamps[step]+"_sub222.tif";
        
    }

    public static String destinationPath(String specimen){
        return mainDir+specimen+"/transforms_corrected/Transform_Step_2to0.txt";
        
    }

    public static String getPathToManualRegistrationTransformation(String specimen,int step){
        return (mainDir+specimen+"/transforms_corrected/Transform_Step_"+step+".txt");
    }

    public static ItkTransform  manualReg (ImagePlus refImage, ImagePlus movImage){
        //Fijiyama_GUI fg=new Fijiyama_GUI();
        RegistrationManager regManager = new RegistrationManager(new Fijiyama_GUI());
        regManager.start3dManualRegistration(refImage,movImage);
        ItkTransform tr=regManager.finish3dManualRegistration(); 
        return tr;		
    }



    public static void main(String[] args) {
        ImageJ ij=new ImageJ();
/*         System.setProperty("jogl.forceGL3", "true");
        System.setProperty("jogl.x11.display", ":0");
        GLProfile.initSingleton();
        GLProfile glp = GLProfile.get(GLProfile.GL3);*/
        //test1_testStep1With();
//        decode_transform_nofuck();
        // String[]specs=new String[]{"B_212"};//"B_227",
        // for(String sp : specs){
        //     test2_testAllsteps(sp);
        //     VitimageUtils.waitFor(300000000);
        // }
        for (int i = 201; i <= 240; i++) {
            if (i == 204) {
                continue;
            }
            String specimen = "B_" + i;
            //decode_transform_nofuck(specimen);
            test2_testAllsteps(specimen);
            //VitimageUtils.waitFor(30000);
        }

        System.out.println("I have finished");
    }


    public static void decode_transform_nofuck(String specimen){
        // for(String specimen : new String[]{"B_202","B_203"}){
            for(int step:new int[]{0,1}){

                String originalPath=getPathToManualRegistrationTransformation(specimen,step);
                // String destinationPath=originalPath+"butgood.txt";
                String destinationPath="/home/phukon/Desktop/Cuttings_MRI_registration/"+specimen+"/transforms_corrected/Transform_Step_"+step+".txt";
                ItkTransform trRepaired=new ItkTransform(ItkTransform.quickFixBrokenItkMatrixTo4x4Matrix(originalPath));
                trRepaired.writeMatrixTransformToFile(destinationPath);
            }
        // }
    }

    public static void test2_testAllsteps(String specimen){
        ImagePlus img0 = IJ.openImage(getPathToImage(specimen,0));
        ImagePlus img1 = IJ.openImage(getPathToImage(specimen,1));
        ImagePlus img2 = IJ.openImage(getPathToImage(specimen,2));
        // img0.setTitle("img0");img0.show();img0.setSlice(255);
        // img1.setTitle("img1");img1.show();img1.setSlice(255);
        // img2.setTitle("img2");img2.show();img2.setSlice(255);

        //Gather the corresponding transformation
        //And compute the composed one
        ItkTransform tr_1to0=ItkTransform.readTransformFromFile(getPathToManualRegistrationTransformation(specimen,0));
        ItkTransform tr_2to1=ItkTransform.readTransformFromFile(getPathToManualRegistrationTransformation(specimen,1));
        ItkTransform tr_2to0=(new ItkTransform(tr_2to1)).addTransform(tr_1to0);
        tr_2to0.writeMatrixTransformToFile(destinationPath(specimen));
        // System.err.println("New matrix saved to: "+destinationPath(specimen));
        // System.out.println("Tr_1to0="+tr_1to0);
        // System.out.println("Tr_2to1="+tr_2to1);
        // System.out.println("Tr_2to0="+tr_2to0);

        //Apply it to the moving image
        ImagePlus img1To0=tr_1to0.transformImage(img0, img1);
        ImagePlus img2To0=tr_2to0.transformImage(img0, img2);
        ImagePlus img2To1=tr_2to1.transformImage(img1, img2);
        
        //Make superimposition (composite)
        ImagePlus img_001_029 = VitimageUtils.compositeNoAdjustOf(img0, img1To0);
        img_001_029.setTitle(specimen + " J001 J029reg");
        img_001_029.show();img_001_029.setSlice(255);
        ImagePlus img_001_077 = VitimageUtils.compositeNoAdjustOf(img0, img2To0);
        img_001_077.setTitle(specimen + " J001 J077");
        img_001_077.show();img_001_077.setSlice(255);

        ImagePlus img_029_077 = VitimageUtils.compositeNoAdjustOf(img1, img2To1);
        img_029_077.setTitle(specimen + " J029 J077");
        img_029_077.show();img_029_077.setSlice(255);


        img_001_029.close();
        img_001_077.close();
        img_029_077.close();
        //System.exit(0);
    }


    public static void test1_testStep1With(){
        String specimen=testSpecimen;
        int step=1;
        ImagePlus imgRef = IJ.openImage(getPathToReferenceImage(specimen,step));
        ImagePlus imgMov = IJ.openImage(getPathToImage(specimen,step));
        //Gather the corresponding transformation
        ItkTransform tr=ItkTransform.readTransformFromFile(getPathToManualRegistrationTransformation(specimen,step));

        //Apply it to the moving image
        ImagePlus imgManuallyRegistered=tr.transformImage(imgRef, imgMov);

        imgRef.show();
        imgMov.show();

        //Make superimposition (composite)
        VitimageUtils.compositeNoAdjustOf(imgRef, imgManuallyRegistered).show();
        VitimageUtils.waitFor(10000);
        System.exit(0);
    }
        

    
}
