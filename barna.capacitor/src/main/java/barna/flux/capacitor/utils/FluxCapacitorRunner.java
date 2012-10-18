package barna.flux.capacitor.utils;

import barna.commons.Execute;
import barna.commons.parameters.ParameterException;
import barna.commons.system.OSChecker;
import barna.flux.capacitor.profile.MappingStats;
import barna.flux.capacitor.reconstruction.FluxCapacitor;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.io.rna.UniversalReadDescriptor;
import com.martiansoftware.jsap.RequiredParameterMissingException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.concurrent.Future;

/**
 *
 * @author  Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class FluxCapacitorRunner {

    /**
     * Default common files
     */
    public static String DEFAULT_PARAMETER_FILE = "params.par";
    public static String DEFAULT_OUTPUT_FILE = "output"+ File.separator+"results.gtf";

    /**
     * Run the capacitor with the specified parameter file.
     *
     * @param parFile the parameter file
     * @return output the output of the capacitor
     */
    public static MappingStats runCapacitor(File parFile) throws Exception {
        FluxCapacitor capacitor= new FluxCapacitor();
        capacitor.setFile(parFile);
        Future<MappingStats> captain= Execute.getExecutor().submit(capacitor);
        MappingStats stats = captain.get();
        return stats;
    }

    /**
     * Creates a temporary directory structure for the current run of the Capacitor
     *
     * @param cwd current working directory
     * @param parameters Map of Capacitor parameters
     * @throws RequiredParameterMissingException if param name is not a known parameter name
     * @return the current paramter file
     */
    public static File createTestDir(File cwd, Map<String,Object> parameters) throws RequiredParameterMissingException,ParameterException,IOException {

        //check for mandatory parameters
        if (!parameters.containsKey("ANNOTATION_FILE"))
            throw new RequiredParameterMissingException("The parameter for annotation file is missing");
        if (!parameters.containsKey("MAPPING_FILE"))
            throw new RequiredParameterMissingException("The parameter for mapping file is missing");

        //get instance for the read descriptor
        if (parameters.containsKey("READ_DESCRIPTOR")) {
            UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
            descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        }

        //check if sorted files should be kept and set up directory
        if (parameters.containsKey("KEEP_SORTED")) {
            String sortedPath = parameters.get("KEEP_SORTED").toString();
            if (!sortedPath.startsWith(File.separator)) {
                File sortDir = new File(cwd, sortedPath);
                sortDir.mkdir();
                parameters.put("KEEP_SORTED", sortDir);
            }
        }

        //set up the output file
        File outDir = new File(cwd, "output");
        outDir.mkdir();
        if (!parameters.containsKey("STDOUT_FILE")) {
            File outFile = new File(cwd, DEFAULT_OUTPUT_FILE);
            parameters.put("STDOUT_FILE", outFile);
            outFile.delete();
        }

        //write the parameter file
        File parFile= new File(cwd,DEFAULT_PARAMETER_FILE);
        FluxCapacitorSettings settings= new FluxCapacitorSettings();
        for (Object m : parameters.entrySet()) {
            Map.Entry<String,Object> entry = (Map.Entry)m;
            settings.set(entry.getKey(), entry.getValue().toString());
        }
        settings.validate();

        BufferedWriter buffy = null;
        try {
            buffy = new BufferedWriter(new FileWriter(parFile));
            for (Object m : parameters.entrySet()) {
                Map.Entry<String,Object> entry = (Map.Entry)m;
                buffy.write(entry.getKey()+"\t"+entry.getValue().toString()+OSChecker.NEW_LINE);
            }
        } catch (IOException e) {
            throw e;
        } finally {
            if (buffy!=null)
                buffy.close();
        }

        return parFile;
    }
}
