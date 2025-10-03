    // @Override
    // public void execute(Specimen specimen) throws Exception {
    //     execute(specimen, true);
    // }

    // public static void main (String[] args) throws Exception {
    //     Specimen spec = new Specimen("B_201");
    //     ImageJ ij = new ImageJ();
    //     new Step_7_ProbabilisticAtlas().execute(spec,true);
    //     System.out.println("Saved!");
    // }

    // public void execute(Specimen specimen, boolean testing) throws Exception {
    //     String[] timestamps = Config.timestamps;
    //     int t = 4;
    //     ImagePlus mask = computeOtsuThreshold(cond_CONTROL, var_CHARD);
    //     mask.show();
    // }

    // public static ImagePlus computeOtsuThreshold(int condition, int variety){
    //     String[] spec = Config.getSpecimensName(condition,variety);
    //     System.out.println(Config.mainDir + "Processing/03_PolarTransform/"+spec[0]+"_GeneralizedPolarTransform.tif");
    //     ImagePlus img = IJ.openImage(Config.mainDir + "Processing/03_PolarTransform/"+spec[0]+"_GeneralizedPolarTransform.tif");
    //      if (img == null) {
    //         throw new IllegalStateException("Not Found!");
    //     }
    //     ImagePlus imgT1 = new Duplicator().run(img,1,1,256,img.getNSlices()-256,1,1);
    //     IJ.run(imgT1, "Gaussian Blur...", "sigma=4 stack");
    //     imgT1.setDisplayRange(0.00,1.80);
    //     imgT1.updateAndDraw();
    //     IJ.setAutoThreshold(imgT1, "Default dark no-reset");
    //     return imgT1;


    // }

     //     ImagePlus stack_CHARD = ImgProUtils.buildWekaTrainingStack(cond_CONTROL, var_CHARD, t);
        //     ImagePlus stack_MER = ImgProUtils.buildWekaTrainingStack(cond_CONTROL, var_MERLOT, t);
        //     ImagePlus stack_TEMP= ImgProUtils.buildWekaTrainingStack(cond_CONTROL, var_TEMPRA, t);
        //     ImagePlus stack_UG = ImgProUtils.buildWekaTrainingStack(cond_CONTROL, var_UGNI, t);
        //     ImagePlus result = ImgProUtils.combineStacks("All_Var_CT_z512_627_t"+t+".tif",  stack_CHARD   ,stack_MER, stack_TEMP, stack_UG); // name of the new combined stack

        //     IJ.saveAsTiff(result, Config.mainDir + "Processing/04_Masks/01_StacksForWeka/" + "All_Var_CT_z512_627_t"+t+".tif");
        // }


        // String[] timestamps = Config.timestamps;

        // ImagePlus mask_chard = ImgProUtils.trainWekaToGetMask(Config.mainDir + "Results/04_MaskWeka/B_201_J141_mask.tif", t);

        // IJ.saveAsTiff(atlas_chard, Config.mainDir + "Results/05_ProbabilisticAtlas/ProbabilisticAtlas_CHARD_CT_"+timestamps[t-1]+".tif");
        // ImagePlus atlas_mer = buildProbabilisticAtlas(cond_CONTROL, var_MERLOT, t);
        // atlas_mer.show();
        // IJ.saveAsTiff(atlas_mer, Config.mainDir + "Results/05_ProbabilisticAtlas/ProbabilisticAtlas_MER_CT_"+timestamps[t-1]+".tif");
        // ImagePlus atlas_temp = buildProbabilisticAtlas(cond_CONTROL, var_TEMPRA, t);
        // atlas_temp.show();
        // IJ.saveAsTiff(atlas_temp, Config.mainDir + "Results/05_ProbabilisticAtlas/ProbabilisticAtlas_TEMP_CT_"+timestamps[t-1]+".tif");
        // ImagePlus atlas_ug = buildProbabilisticAtlas(cond_CONTROL, var_UGNI, t);
        // atlas_ug.show();
        // IJ.saveAsTiff(atlas_ug, Config.mainDir + "Results/05_ProbabilisticAtlas/ProbabilisticAtlas_UGNI_CT_"+timestamps[t-1]+".tif");
        // String path = Config.mainDir + "Results/03_MaskforPA/B_201_J077_mask.tif";
        
        // ImagePlus m = IJ.openImage(path);   // opens the 3D stack
        // fillTopToFullRowOrBest(m, 0.40);  