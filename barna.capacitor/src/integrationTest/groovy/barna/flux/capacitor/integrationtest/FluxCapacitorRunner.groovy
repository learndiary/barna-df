package barna.flux.capacitor.integrationtest

import com.martiansoftware.jsap.RequiredParameterMissingException
import barna.io.FileHelper
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings

/**
 *
 * @author  Emilio Palumbo (emiliopalumbo@gmail.com)
 */
class FluxCapacitorRunner {

    static void createTestDir(File cwd, Map parameters) throws RequiredParameterMissingException {
        if (!parameters.containsKey("ANNOTATION_FILE"))
            throw new RequiredParameterMissingException("The parameter for annotation file is missing")
        if (!parameters.containsKey("MAPPING_FILE"))
            throw new RequiredParameterMissingException("The parameter for mapping file is missing")
        if (!parameters.containsKey("TMP_DIR")) {
            parameters.put([ "TMP_DIR" : FileHelper.createTempDir("run_tmp","",cwd).path]);
        }
        File outDir = FileHelper.createTempDir("output","",cwd);
        if (!parameters.containsKey("STDOUT_FILE")) {
            parameters.put([ "STDOUT_FILE" : FileHelper.append(outDir.path+File.separator+"FluxCapacitor.gtf")]);
        }


        //Writing parameter file
        File parFile= File.createTempFile("FluxCapacitorTest", "par", cwd);
        FluxCapacitorSettings settings= new FluxCapacitorSettings();
        for (String s : parameters.keySet()) {
            settings.set(s,parameters[s]);
        }
        settings.write(new BufferedOutputStream(parFile));
    }
}
