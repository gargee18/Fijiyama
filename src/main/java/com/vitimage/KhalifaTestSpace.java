package com.vitimage;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import math3d.Point3d;

public class KhalifaTestSpace {
	
	
	public static void main(String[]args) {
		System.out.println("Lancement tests Khalifa");
		testAutomaticRegistration();
		/*ImageJ ij=new ImageJ();
		ImagePlus imgReference=IJ.openImage("/home/fernandr/Bureau/Test/Origine/RX_flip.tif");
		ImagePlus imgMoving=IJ.openImage("/home/fernandr/Bureau/Test/Origine/IRM.tif");
		imgMoving=runManualRegistration(imgReference,imgMoving);
		runAutomaticRegistration(imgReference,imgMoving);
		*/
	}
	
	
	public static void testManualRegistration() {
		
		
	}
	
	public static void testAutomaticRegistration() {
		ImageJ ij=new ImageJ();
		ImagePlus imgReference=IJ.openImage("/home/fernandr/Bureau/Test/Origine/RX_flip_sub.tif");
		ImagePlus imgMoving=IJ.openImage("/home/fernandr/Bureau/Test/Origine/IRM_to_RX_flip_sub.tif");
		runAutomaticRegistration(imgReference,imgMoving);
	}
	
	public static ImagePlus runManualRegistration(ImagePlus imgRef,ImagePlus imgMov) {
		ImagePlus imgReference=new Duplicator().run(imgRef);
		ImagePlus imgMoving=new Duplicator().run(imgMov);
		ImagePlus imgMovingSave=new Duplicator().run(imgMov);
		imgReference.setTitle("Reference");
		imgMoving.setTitle("Moving");
		imgReference.show();
		imgMoving.show();
		boolean computeRealCoordinatesTransformation=true;
		int nbWantedPointsPerImage=5;
		int j=-1;
		ItkTransform transGlobal=new ItkTransform();
		while(j<3) {
			j++;
			Point3d [][]pointTab=VitiDialogs.registrationPointsUI(nbWantedPointsPerImage, imgReference,imgMoving, computeRealCoordinatesTransformation);
			Point3d []pRef=pointTab[0];
			Point3d []pMov=pointTab[1];
			for(int i=0;i<pRef.length;i++)System.out.println("Pref["+i+"]="+pRef[i]);
			for(int i=0;i<pMov.length;i++)System.out.println("PMov["+i+"]="+pMov[i]);
			ItkTransform trans=ItkTransform.estimateBestRigid3D(pMov,pRef);
			transGlobal.addTransform(trans);
			System.out.println("Transformation calculée ="+transGlobal);
			ImagePlus result=transGlobal.transformImage(imgReference, imgMovingSave);	
			result.getProcessor().resetMinAndMax();
			result.setTitle("Result");
			result.show();
			ImagePlus comp=VitimageUtils.compositeOf(imgReference, result);
			VitimageUtils.imageChecking(comp);
			VitimageUtils.waitFor(10000);
			//trans.writeTransform("/home/fernandr/Bureau/Test/Origine/trans_"+j+".itktr");
			imgMoving=result;
		}
		return imgMoving;
	}
	
	
	public static ImagePlus runAutomaticRegistration(ImagePlus imgRef,ImagePlus imgMov) {
		System.out.println("Recalage automatique");
		ItkRegistrationManager manager=new ItkRegistrationManager();
		ImagePlus result=manager.runScenarioInterModal(imgRef, imgMov);
		System.out.println("Recalage ok.");
		result.show();
		ImagePlus result2=VitimageUtils.convertFloatToShortWithoutDynamicChanges(result);
		System.out.println("Conversion  ok");
		result2.setTitle("Resultat recalage auto");
		result2.getProcessor().resetMinAndMax();
		result2.show();
		return result2;
	}
	
	
	
	
	
	
	
	public KhalifaTestSpace() {
		
	}

}
