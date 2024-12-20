package io.github.rocsg.fijiyama.gargeetest.cuttings.core;

public class Specimen {
    private String name;
    private int condition;  // Condition index ( 0 for PCH, 1 for CT)
    private int variety;    // Variety index ( 0 for CHARD, 1 for MERLOT, 2 for TEMPRA, 3 for UGNI)
    
    public Specimen(String name) {
        this.name = name;
        setCondition();
        setVariety();
    }

    public String getName() {
        return name;
    }

    public String getSpecimenDirectory() {
        return Config.mainDir+name+"/";
    }

    public void setCondition(){
        int id = Integer.parseInt(name.split("_")[1]);
        if ((id >= 201 && id <= 206) || (id >= 211 && id <= 216) ||  id == 219 ||  (id >= 221 && id <= 226) || (id >= 231 && id <= 236)) {
            this.condition = 0; // PCH
        } 
        else this.condition = 1; // CT
    }

    public void setVariety(){
        int id = Integer.parseInt(name.split("_")[1]);

        if ((id >= 201 && id <= 210)) {
            this.variety = 0; // CHARD
        } else if ((id >= 211 &&  id <= 220)) {
            this.variety = 1; // MERLOT
        } else if ((id >= 221 &&  id <= 230)) {
            this.variety = 2; // TEMPRA
        } else this.variety = 3; // UGNI
    }

    public int getCondition() {
        return this.condition;
    }

    public int getVariety() {
        return this.variety;
    }

    public static Specimen getSpecimen(String name){
        Specimen spec = new Specimen(name);
        return spec;
    }

    // /* Not in use anymore*/
    // //SHIT HERE
    // public static String[] getSpecimensName(int condition, int variety) {
    //     String[][][] specs = new String[][][] {
    //             {
    //                     { "B_201", "B_202", "B_203", "B_205", "B_206" }/* CHARD */,
    //                     { "B_211", "B_212", "B_213", "B_214", "B_215", "B_216", "B_219" }/* MERLOT */,
    //                     { "B_221", "B_222", "B_223", "B_224", "B_225", "B_226" }/* TEMPRA */,
    //                     { "B_231", "B_232", "B_233", "B_234", "B_235", "B_236" }/* UGNI */
    //             }/* PCH */,
    //             {
    //                     { "B_207", "B_208", "B_209", "B_210" }/* CHARD */,
    //                     { "B_217", "B_218", "B_220" }/* MERLOT */,
    //                     { "B_227", "B_228", "B_229", "B_230" }/* TEMPRA */,
    //                     { "B_237", "B_238", "B_239", "B_240" }/* UGNI */
    //             }/* CT */
    //     };

    //     return specs[condition][variety];
    // }
    // //END OF SHIT

    public String toString(){
        return "Name="+name+" variety="+variety+" condition="+condition;
    }

    public static void main(String[]args){
        test();
    }


    public static void test(){
        Specimen s=new Specimen("B_227");
        System.out.println(s);
        }

}
    


