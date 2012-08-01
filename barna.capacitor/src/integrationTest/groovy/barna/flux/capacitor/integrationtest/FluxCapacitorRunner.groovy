package barna.flux.capacitor.integrationtest

import com.martiansoftware.jsap.RequiredParameterMissingException
import barna.io.FileHelper
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings

/**
 *
 * @author  Emilio Palumbo (emiliopalumbo@gmail.com)
 */
class FluxCapacitorRunner {

    static File parFile;
    static File tmpDir;
    static File outFile;

    static void createTestDir(File cwd, Map parameters) throws RequiredParameterMissingException {
        if (!parameters.containsKey("ANNOTATION_FILE"))
            throw new RequiredParameterMissingException("The parameter for annotation file is missing")
        if (!parameters.containsKey("MAPPING_FILE"))
            throw new RequiredParameterMissingException("The parameter for mapping file is missing")
        if (!parameters.containsKey("TMP_DIR")) {
            tmpDir = FileHelper.createTempDir("run_tmp","",cwd);
            parameters.put("TMP_DIR", tmpDir);
        } else {
            tmpDir = parameters["TMP_DIR"];
        }
        File outDir = FileHelper.createTempDir("output","",cwd);
        if (!parameters.containsKey("STDOUT_FILE")) {
            outFile = FileHelper.createTempFile("FluxCapacitor","gtf",outDir);
            parameters.put("STDOUT_FILE", outFile);
            outFile.delete();
        }


        //Writing parameter file
        parFile= File.createTempFile("FluxCapacitorTest", "par", cwd);
        FluxCapacitorSettings settings= new FluxCapacitorSettings();
        for (String s : parameters.keySet()) {
            settings.set(s,parameters[s]);
        }
        BufferedWriter buffy = new BufferedWriter(new FileWriter(parFile));
        buffy.write(settings.toString());
        buffy.close();
    }
}
