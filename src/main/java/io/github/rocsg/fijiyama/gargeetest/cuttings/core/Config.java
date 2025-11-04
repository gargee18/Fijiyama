package io.github.rocsg.fijiyama.gargeetest.cuttings.core;

public class Config {

    public static final double max_display_val = 1.8;

    // Aqcuition time for cuttings data
    public final static String[] timestamps = new String[]{"J001", "J029", "J077", "J141"};
   
    // Base directory for this project
    // public final static String mainDir = "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/Cuttings_MRI_registration/";
    public final static String mainDir = "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/Cuttings/"; //If phenodrone
    // public final static String mainDir = "/home/phukon/Desktop/Cuttings/"; //if local
    // public final static String mainDir = "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/T1_T2_Specimen_Analysis/"; 

    // Path to get high res raw data
    public static String getPathToRawImageDir() {
        return mainDir+"Data/01_Raw/";
    }

    // Get Raw data 
    public static String getPathToRawImage(Specimen specimen, int step) {
        return mainDir+"Data/01_Raw/"+specimen.getName()+"_"+timestamps[step]+".tif";
    }

    // Path to get high res normalized data
    public static String getPathToNormalizedImagesDir(){
        return mainDir+"Data/02_Normalized/";
    }

    // Get high res normalized data
    public static String getPathToNormalizedImage(Specimen specimen, int step){
        return mainDir+"Data/02_Normalized/"+specimen.getName()+"_"+timestamps[step]+"_normalized.tif";
    }

    // Get path to aligned images
    public static String getPathToInocAlignedImage(Specimen specimen, int step) {
        return mainDir+"Data/03_InocAligned/"+specimen.getName()+"_"+timestamps[step]+"_aligned.tif";
    }

    // Path to get low res data 
    public static String getPathToSubsampledImageDir(){
        return mainDir+"Data/04_Subsampled/";
    }

     // Get low res data 
     public static String getPathToSubsampledImage(Specimen specimen, int step){
        return mainDir+"Data/04_Subsampled/"+specimen.getName()+"_"+timestamps[step]+"_sub.tif";
    }

    public static String getPathToMask(Specimen specimen, int step){
        return mainDir+"Data/05_Mask/"+specimen.getName()+"_mask.tif";
        // return mainDir+"Data/05_Mask/"+specimen.getName()+"_mask_"+timestamps[step]+".tif";
    }

    // Get Cropped data 
    public static String getPathToCroppedImage(Specimen specimen, int step) {
        return mainDir+"Data/06_Cropped_z/"+specimen.getName()+"_"+timestamps[step]+"_cropped_z.tif";
    }

    // Get path to aligned images sub
    public static String getPathToInocAlignedImageSub(Specimen specimen, int step) {
        return mainDir+"Data/07_InocSub/"+specimen.getName()+"_"+timestamps[step]+".tif";
    }

    //Get contour ROI
    public static String getPathToContourROIForGeneralizedPolarTransform(Specimen specimen){
        return mainDir+"Data/08_ROIContourGPT/"+specimen.getName()+"_contour.roi";
    }

    // Get Transformation Matrix for Inoculation Alignment
    public static String getPathToInoculationAlignmentTransformation(Specimen specimen, int step){
        return mainDir+"Processing/01_InocAlignment/"+specimen.getName()+"_"+timestamps[step]+"_TR_MAT_Inoculation_Alignment.txt";
    }
    // Get Transformation Matrix for Manual Rigid Registration
    public static String getPathToManualRigidRegistrationMatrix(Specimen specimen, int stepRef, int stepMov){
        return mainDir+"Processing/02_RigidRegistration/"+specimen.getName()+"_"+timestamps[stepRef]+"_"+timestamps[stepMov]+"_TR_MAT_Rigid_Reg_man.txt";
    }
    // Get Transformation Matrix for Auto Rigid Registration
    public static String getPathToAutoRigidRegistrationMatrix(Specimen specimen, int stepRef, int stepMov){
        return mainDir+"Processing/02_RigidRegistration/"+specimen.getName()+"_"+timestamps[stepRef]+"_"+timestamps[stepMov]+"_TR_MAT_Rigid_Reg_auto.txt";
    }

    //Get Polar Transformed Images Directory
    public static String getPathtoPolarTransformsDir(Specimen specimen){
        return mainDir+"Processing/03_PolarTransform/";
    }

    //Get hyperstack
    public static String getPathToHyperstack(Specimen specimen){
        return mainDir+"Results/01_Hyperstack/"+specimen.getName()+"_Hyperstack.tif";
    }

    //Get polar atlas
    public static String getPathToPolarAtlas(){
        return mainDir+"Results/02_MeanIndividual/Polar/";
    }

    public static String getPathToCartesianAtlas(){
        return mainDir+"Results/02_MeanIndividual/Cartesian/";
    }


     public static String getPathToMask(){
        return mainDir+"Results/03_MaskforPA/";
    }

    public static String getRawT1T2sequence(Specimen specimen, int time){
        return mainDir + "/processing/input/"+specimen.getName()+"_"+timestamps[time]+"_TR2400_TE12.tif";
    }

    public static String getPathToSerialRegistrationTRMatrix(Specimen specimen, int time){
        return mainDir+"/processing/serialRegistration/Exported_data/transform_global"+specimen.getName()+"_"+timestamps[time]+"_TR2400_TE12.txt";
    }

    public static String getPathToRegisteredT1T2sequence(Specimen specimen, int time){
        return mainDir+"/processing/stacksRegistered/"+specimen.getName()+"_"+timestamps[time]+"_TR2400_TE12_registered.tif";
    }

    public static String getPathToNormalizedT1T2sequence(Specimen specimen, int time){
        return mainDir+"/processing/stacksRegistered/"+specimen.getName()+"_"+timestamps[time]+"_TR2400_TE12_registered_normalized_sigma_corrected.tif";
    }
   
    public static String getPathToHyperMap(Specimen specimen){
        return mainDir+"/results/"+specimen.getName()+"_HyperMap.tif";
    }
   
    
   
    static final int cond_PCH = 0;
    static final int cond_CONTROL = 1;

    static final int var_CHARD = 0;
    static final int var_MERLOT = 1;
    static final int var_TEMPRA = 2;
    static final int var_UGNI = 3;
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



   

   

    public static String getPathToSpecimenImageAtTimeStamp(Specimen specimen,int timestamp){
        // return mainDir+specimen+"/outputs/"+specimen+"_001_141_TransformedImage_Reg_Aligned.tif";
        return mainDir+specimen.getName()+"/raw/"+specimen.getName()+"_"+timestamps[timestamp]+"_normalized.tif";
    }



    //CRAP HERE
    public static String getPathToImageMov(String specimen,int step){
        // return mainDir+specimen+"/outputs/"+specimen+"_001_141_TransformedImage_Reg_Aligned.tif";
        return mainDir+specimen+"/raw/"+specimen+"_"+timestamps[step]+"_normalized.tif";
    }

    public static String getPathToImageRef(String specimen,int step){
        return mainDir+specimen+"/raw/"+specimen+"_"+timestamps[step]+"_aligned.tif";
    }
    //END OF CRAP

    public static String testSpecimen = "B_201";   // An-Information-with-elements-separated-with-minus_while-different-fields-will-be-separated-by-underscores

    
    public static String getPathToImage(Specimen specimen,int step){
        System.out.println(mainDir+specimen+"/raw/"+specimen+"_"+timestamps[step]+"_normalized.tif");
        return mainDir+specimen+"/raw/"+specimen+"_"+timestamps[step]+"_normalized.tif";
    }

    public static String getPathToSaveAlignedImage(String specimen,int step){
        System.out.println(mainDir+specimen+"/raw/"+specimen+"_"+timestamps[step]+"_aligned.tif");
        return mainDir+specimen+"/raw/"+specimen+"_"+timestamps[step]+"_aligned.tif";
    }
    

    public static String getPathToReferenceImage(String specimen){
        return mainDir+specimen+"/raw/"+specimen+"_J001_aligned.tif";
    }

    public static String getPathToMovingImage(String specimen){
        return mainDir+specimen+"/raw/"+specimen+"_J029_normalized.tif";
    }

    public static String getPathToTrInit(String specimen){
        return mainDir+specimen+"/transforms_corrected/Transform_Step_3.txt";
    }

    public static String getPathToTrFinal(String specimen){
        return mainDir+specimen+"/transforms_corrected/Automatic_Transform_3to0.txt";
    }

   

    public static String getPathToImageLowRes(String specimen,int step){
        return mainDir+specimen+"/raw_subsampled/"+specimen+"_"+timestamps[step]+"_sub222.tif";
        
    }

    public static String getPathToImageHighRes(String specimen,int step){
            return mainDir+specimen+"/raw/"+specimen+"_"+timestamps[step]+".tif";
            
    }
    
    public static String getPathtoSaveAlignedTransformationMatrix(String specimen){
        return mainDir+specimen+"/transforms_corrected/Transform_matrix_ref_alignment.txt";
    }  

    public static String getPathToSaveFinalTransform(String specimen, int step){
        return mainDir+specimen+"/transforms_corrected/Final_transform_001_"+step+".txt";
    }

    public static String[] getTimestamps() {
        return timestamps;
    }

    public static String getMaindir() {
        return mainDir;
    }

    public static int getCondPch() {
        return cond_PCH;
    }

    public static int getCondControl() {
        return cond_CONTROL;
    }

    public static int getVarChard() {
        return var_CHARD;
    }

    public static int getVarMerlot() {
        return var_MERLOT;
    }

    public static int getVarTempra() {
        return var_TEMPRA;
    }

    public static int getVarUgni() {
        return var_UGNI;
    }

    public static String getTestSpecimen() {
        return testSpecimen;
    }

    public static void setTestSpecimen(String testSpecimen) {
        Config.testSpecimen = testSpecimen;
    }
        
}
