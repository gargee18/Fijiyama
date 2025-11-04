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

        // public static Result fitFromPoints(double[][] pts){
        
    //     System.out.println("Number of points = " + (pts == null ? "null" : pts.length));
    //     if (pts==null || pts.length < 6) throw new IllegalArgumentException("Need >= 6 points");
    //     // Build linear design matrix for parameters [A,B,C,D,E,F,G]
    //     final int n = pts.length;
    //     double[][] D = new double[n][7];
    //     for (int i=0;i<n;i++){
    //         double x=pts[i][0], y=pts[i][1], z=pts[i][2];
    //         D[i][0] = x*x;  // A
    //         D[i][1] = y*y;  // B
    //         D[i][2] = z*z;  // C
    //         D[i][3] = x;    // D
    //         D[i][4] = y;    // E
    //         D[i][5] = z;    // F
    //         D[i][6] = 1.0;  // G
    //     }
    //     RealMatrix DM = MatrixUtils.createRealMatrix(D);
    //     // Solve homogeneous least squares: find p (||p||=1) minimizing ||D p||
    //     SingularValueDecomposition svd = new SingularValueDecomposition(DM);
    //     RealVector p = svd.getV().getColumnVector(6); // smallest singular vector
    //     double A = p.getEntry(0), B = p.getEntry(1), C = p.getEntry(2);
    //     double Dp = p.getEntry(3), Ep = p.getEntry(4), Fp = p.getEntry(5), G = p.getEntry(6);

    //     // Complete the square to extract center and radii.
    //     // For 1D: A x^2 + D x = A (x - x0)^2 - A x0^2 with x0 = -D/(2A)
    //     // In 3D, after shifting to center, the constant becomes:
    //     //   Gc = G - D^2/(4A) - E^2/(4B) - F^2/(4C)
    //     // Canonical form: A (X)^2 + B (Y)^2 + C (Z)^2 + Gc = 0,
    //     // divide by -Gc (>0) to obtain (X^2/(S/A)) + ... = 1 with S = -Gc
    //     Result R = new Result();
    //     StringBuilder note = new StringBuilder();

    //     // Guard against zeros
    //     if (A==0 || B==0 || C==0) {
    //         R.ok=false; R.note = "A/B/C contains zero -> degenerate"; return R;
    //     }

    //     double x0 = -Dp/(2*A);
    //     double y0 = -Ep/(2*B);
    //     double z0 = -Fp/(2*C);

    //     double Gc = G - (Dp*Dp)/(4*A) - (Ep*Ep)/(4*B) - (Fp*Fp)/(4*C);
    //     // We want S = -Gc > 0; if not, flip the overall sign (p is defined up to scale)
    //     double S = -Gc;
    //     if (S <= 0) {
    //         A = -A; B = -B; C = -C; Dp = -Dp; Ep = -Ep; Fp = -Fp; G = -G;
    //         x0 = -Dp/(2*A); y0 = -Ep/(2*B); z0 = -Fp/(2*C);
    //         Gc = G - (Dp*Dp)/(4*A) - (Ep*Ep)/(4*B) - (Fp*Fp)/(4*C);
    //         S = -Gc; // try again
    //         note.append("(flipped sign)");
    //     }

    //     double rx2 = S / A;
    //     double ry2 = S / B;
    //     double rz2 = S / C;

    //     boolean ok = rx2>0 && ry2>0 && rz2>0;
    //     R.ok = ok;
    //     R.center = new double[]{x0,y0,z0};
    //     if (ok) {
    //         R.radii = new double[]{ Math.sqrt(rx2), Math.sqrt(ry2), Math.sqrt(rz2) };
    //     } else {
    //         R.radii = new double[]{Double.NaN, Double.NaN, Double.NaN};
    //         note.append(" non-positive radii^2");
    //     }

    //     // Residuals: evaluate |(dx^2/rx^2 + dy^2/ry^2 + dz^2/rz^2 - 1)|
    //     double sum=0, sum2=0; int m=0;
    //     if (ok) {
    //         for (double[] q : pts) {
    //             double dx = (q[0]-x0); double dy=(q[1]-y0); double dz=(q[2]-z0);
    //             double val = (dx*dx)/rx2 + (dy*dy)/ry2 + (dz*dz)/rz2 - 1.0;
    //             val = Math.abs(val);
    //             sum += val; sum2 += val*val; m++;
    //         }
    //     }
    //     if (m>0) {
    //         double mean = sum/m; double var = Math.max(0, sum2/m - mean*mean);
    //         R.meanResidual = mean; R.stdResidual = Math.sqrt(var); R.nPoints=m;
    //     } else { R.meanResidual=Double.NaN; R.stdResidual=Double.NaN; R.nPoints=0; }
    //     R.note = note.toString();
    //     return R;
    // }