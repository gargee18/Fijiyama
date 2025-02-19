package io.github.rocsg.fijiyama.gargeetest.cuttings.extra;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

public class FileMove {
    public static void main(String[] args) {
        String sourceDir = "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/Cuttings_MRI_registration/";      
        String destinationDir = "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/Cuttings/Data/08_ROIContourGPT/"; 
        String oldFileName = "contour.roi";  // Original file name
        
        for (int i = 201; i <= 240; i++) {
            String folderName = "B_" + i;
            File sourceFolder = new File(sourceDir + folderName + "/roi_cambium_on_hyperimage/");

            File sourceFile = new File(sourceFolder, oldFileName);

            File destinationFile = new File(destinationDir, folderName + "_" + oldFileName);

            try {
                Files.copy(sourceFile.toPath(), destinationFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                System.out.println("Moved: " + sourceFile.getAbsolutePath() + " -> " + destinationFile.getAbsolutePath());
            } catch (IOException e) {
                System.err.println("Failed to move: " + sourceFile.getAbsolutePath());
                e.printStackTrace();
            }
        }
    }
    
}
