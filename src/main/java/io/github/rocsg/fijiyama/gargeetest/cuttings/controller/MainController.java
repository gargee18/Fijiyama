package io.github.rocsg.fijiyama.gargeetest.cuttings.controller;
import java.util.ArrayList;
import java.util.List;
import ij.ImageJ;

import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D.Step_0_Normalize;
import io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D.Step_1_Subsample;
import io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D.Step_2_2_CreateMask;
import io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D.Step_2_InocAlignment;
import io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D.Step_3_RegistrationRigid;
import io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D.Step_4_Hyperstack;
import io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D.Step_5_PolarTransform;
import io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D.Step_6_AtlasBuilding;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Pipeline;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;


public class MainController {

    public static void main(String[] args) throws Exception {
       ImageJ ij=new ImageJ();
    //    test();
      run();
   }

   
   public static void test() throws Exception{
       // List<Specimen> specimens = loadPCHSpecimens();
      
       // for(Specimen s:specimens)System.out.println(s);
       Specimen spec= new Specimen("B_201");
       List<PipelineStep> steps = getStepsForDevPipeline();
       Pipeline pipeline = new Pipeline(steps);
       pipeline.run(spec);
   }


   public static void run(){
       try {
           // Load specimen data
           List<Specimen> specimens = loadAllSpecimens();
        //    List<Specimen> specimens = loadTestSpecimens();
           

           // Define the pipeline steps
           List<PipelineStep> steps = getStepsForDevPipeline();

           // Create the pipeline
           Pipeline pipeline = new Pipeline(steps);

           // Run the pipeline for each specimen
           for (Specimen specimen : specimens) {
               System.out.println("Processing specimen: " + specimen.getName());
               pipeline.run(specimen);
           }

           System.out.println("Pipeline completed!");
       } catch (Exception e) {
           e.printStackTrace();
       }
   }

   /**
    * Define the whole pipeline, with all the steps
    */
   public static ArrayList<PipelineStep>getStepsForFullPipeline(){
       ArrayList<PipelineStep> steps = new ArrayList<PipelineStep>();
       steps.add(new Step_0_Normalize());
       steps.add(new Step_1_Subsample());
       steps.add(new Step_2_InocAlignment());
       steps.add(new Step_2_2_CreateMask());
       steps.add(new Step_3_RegistrationRigid());
       steps.add(new Step_4_Hyperstack());
       steps.add(new Step_5_PolarTransform());
       steps.add(new Step_6_AtlasBuilding());
       return steps;
   }

   public static ArrayList<PipelineStep>getStepsForDevPipeline(){
       ArrayList<PipelineStep> steps = new ArrayList<PipelineStep>();
       steps.add(new Step_6_AtlasBuilding());
       return steps;
   }


   public static List<PipelineStep>getSomeSteps(int startingStep,int stoppingStep){
       ArrayList<PipelineStep> steps = getStepsForFullPipeline();
       List<PipelineStep>   someSteps= steps.subList(startingStep, stoppingStep);
       return someSteps;
   }

   /**
    * Define the list of all the specimens
    */
   public static  List<Specimen>loadAllSpecimens(){
       //Aggregate the names of the specimen
       List<Specimen> listSpecimen = new ArrayList<>();
       
        // Specimen names as given in the example
        String[] specimenNames = {
           "B_201", "B_202", "B_203", "B_205", "B_206",   "B_207", "B_208", "B_209", "B_210",
           "B_211", "B_212", "B_213", "B_214", "B_215",   "B_216", "B_217", "B_218", "B_219", "B_220",
           "B_221", "B_222", "B_223", "B_224", "B_225",   "B_226", "B_227", "B_228", "B_229", "B_230",
           "B_231", "B_232", "B_233", "B_234", "B_235",   "B_236", "B_237", "B_238", "B_239", "B_240"
       };

       // Create Specimen objects for each specimen name
       for (String name : specimenNames) {
           Specimen specimen = new Specimen(name); 
           listSpecimen.add(specimen);
       }
       //Transform into Specimen
       return listSpecimen;
   }
   
   public static  List<Specimen>loadTestSpecimens(){ 
    
       List<Specimen> testList=loadAllSpecimens().subList(3, loadAllSpecimens().size());
       return testList;
   }

   public static List<Specimen> loadCTSpecimens() {
       List<Specimen> allSpecimens = loadAllSpecimens();
       List<Specimen> ctSpecimens = new ArrayList<>();
       for (Specimen specimen : allSpecimens) {
           if (specimen.getCondition() == 1) {  // Condition 1 is CT
               ctSpecimens.add(specimen);
           }
       }
       return ctSpecimens;
   }
   public static List<Specimen> loadPCHSpecimens() {
       List<Specimen> allSpecimens = loadAllSpecimens();
       List<Specimen> pchSpecimens = new ArrayList<>();
       for (Specimen specimen : allSpecimens) {
           if (specimen.getCondition() == 0) {  // Condition 0 is PCH
               pchSpecimens.add(specimen);
           }
       }
       return pchSpecimens;
   }
}