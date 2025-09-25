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