package io.github.rocsg.fijiyama.gargeetest.ceps;

import ij.IJ;

public class GeneralUtils {

    public static void main(String[] args) {
//        String a=getSpecimenNameWithMriIndex(1);
        String a= "322";
        int index=getIndexFromMRISpecimenName(a);
        System.out.println(index);
        System.out.println(a);
    }
    // Main Directory



    // Years
    public static String[] years = new String[]{"2022","2023","2024"};

    //private final static String[]specimens=new String[]{"318","322","323","330","335","764B","1181","1193"};

    // Specimens (MRI)
    final static String[]specimenListMRI=new String[]{"318","322","323","330","335","764B","1181","1193"}; 

    // Specimens (XRay)
    final static String[] specimenListXR=new String[]{"313B","318","322","323","330","335","368B","378A","378B","380A","764B","988B","1181","1186A","1189","1191","1193","1195","1266A","2184A"};




    // Modalities for MRI
    static String[]modalities=new String[]{"T1","T2"};


    // Params obtained from the metadata
    static String[]paramsWanted = new String[]{
        "0020,0032", //0: Position
        "0020,0037", //1: Orientation
        "0018,0050", //2: Slice Thickness
        "0028,0030", //3: Pixel Spacing
        "0018,1030", //4: Protocol Name
        "0020,1041", //5: Slice Location
        "0018,1250", //6: Transmitting Coil
        "0018,1251", //7: Receiving Coil
        "0018,1314", //8: Flip Angle
        "0051,100C"  //9: FoV
    };

    /* Getters 
     * 
     * 
     * 
    */
    public static String[]getSpecimenListMRI(){
        String[]ret=new String[specimenListMRI.length];
        for(int i=0;i<ret.length;i++)ret[i]=new String(specimenListMRI[i]);
        return specimenListMRI;
    }

    public static String[]getSpecimenListXR(){
        String[]ret=new String[specimenListXR.length];
        for(int i=0;i<ret.length;i++)ret[i]=new String(specimenListXR[i]);
        return specimenListXR;
    }

    public static String getMRISpecimenNameWithIndex(int index){
        return specimenListMRI[index];
    }

    public static String getXRSpecimenNameWithIndex(int index){
        return specimenListXR[index];
    }

    public static int getIndexFromMRISpecimenName(String str){
        for(int i=0;i<specimenListMRI.length;i++){
            if(specimenListMRI[i].equals(str)){
                return i;
            }
        }
        IJ.showMessage("Big failure : no match for the specimen "+str+" . An error will occur.");
        return -1;
    }

    public static int getIndexFromXRSpecimenName(String str){
        for(int i=0;i<specimenListXR.length;i++){
            if(specimenListXR[i].equals(str)){
                return i;
            }
        }
        IJ.showMessage("Big failure : no match for the specimen "+str+" . An error will occur.");
        return -1;
    }

}
