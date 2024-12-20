package io.github.rocsg.fijiyama.gargeetest.cuttings.testing;

import io.github.rocsg.fijiyama.registration.ItkTransform;
import math3d.Point3d;

public class testTranslationEstimation {
    public static void main(String[] args) {
        
        myFunction();
    }


    public static void myFunction(){

        Point3d[][]correspondancePoints= new Point3d[][]{
            {   new Point3d(0,0,0), new Point3d(1,0,0),new Point3d(0,1,0),new Point3d(0,0,1)}/*This is correspondancePoints[0] */,
            {   new Point3d(10,0,2), new Point3d(11,0,-2),new Point3d(10,1,0),new Point3d(10,0,1)}/*This is correspondancePoints[1] */
        };

        ItkTransform result=ItkTransform.estimateBestTranslation3D(correspondancePoints[1],correspondancePoints[0]);
        System.out.println("Best trans="+result);
        result=ItkTransform.estimateBestRigid3D(correspondancePoints[1],correspondancePoints[0]);
        System.out.println("Best rigid="+result);
    }
}
