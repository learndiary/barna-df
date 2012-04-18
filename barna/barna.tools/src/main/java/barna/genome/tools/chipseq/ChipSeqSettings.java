package barna.genome.tools.chipseq;

import barna.commons.parameters.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * General settings for chipseq tools
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class ChipSeqSettings extends ParameterSchema {
	    
		/**
	     * Load the setting from a file
	     *
	     * @param f the file
	     * @return settings the loaded and validated settings
	     * @throws Exception in case of a parser or validation error
	     */
	    public static ChipSeqSettings createSettings(File f) throws Exception {
	        if (f == null) {
	            throw new NullPointerException("Null parameter file not permitted!");
	        }
	        if (!f.exists()) {
	            throw new IllegalArgumentException("Parameter file " + f.getAbsolutePath() + " can not be found!");
	        }
	        InputStream in = null;
	        try {

	            f = new File(f.getCanonicalPath());  
	            ChipSeqSettings settings = new ChipSeqSettings();
	            in = new FileInputStream(f);
	            settings.parse(in);
	            settings.validate();
	            return settings;
	        } finally {
	            if (in != null) {
	                try {
	                    in.close();
	                } catch (IOException e) {
	                }
	            }
	        }
	    }
	    
	    public static final Parameter<File> FILE_INPUT = Parameters.fileParameter("FILE_INPUT", "input file");
	    public static final Parameter<File> FILE_OUTPUT = Parameters.fileParameter("FILE_OUTPUT", "primary output file", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = schema.get(FILE_OUTPUT);
                if(file == null){
                    throw new ParameterException(parameter, null, "Please specify the FILE_OUTPUT parameter");
                }else if(file.isDirectory()){
                    throw new ParameterException(parameter, file.getAbsolutePath(), "FILE_OUTPUT parameter is a directory, please specify a file");
                }
            }
        });
	    public static final Parameter<File> FILE_OUTPUT2 = Parameters.fileParameter("FILE_OUTPUT2", "secondary output file");
	    public static final Parameter<File> FILE_OUTPUT3 = Parameters.fileParameter("FILE_OUTPUT3", "tertiary output file");

	    public static final Parameter<Boolean> OUTPUT2 = Parameters.booleanParameter("OUTPUT2", "toggle secondary output", false);
	    public static final Parameter<Boolean> OUTPUT3 = Parameters.booleanParameter("OUTPUT3", "toggle tertiary output", false);
	    
	    public static final Parameter<String> READ_DESCRIPTOR = Parameters.stringParameter("READ_DESCRIPTOR", "descriptor of read ID");

	    public static final Parameter<Integer> DISTO_PEAK_MAIN_LEFT = Parameters.intParameter("DISTO_PEAK_MAIN_LEFT", "left flank of main peak");
	    public static final Parameter<Integer> DISTO_PEAK_MAIN_RIGHT = Parameters.intParameter("DISTO_PEAK_MAIN_RIGHT", "right flank of main peak");
	    public static final Parameter<Integer> DISTO_PEAK_MAIN_MAX = Parameters.intParameter("DISTO_PEAK_MAIN_MAXPOS", "position with maximum in main peak");
	    public static final Parameter<Integer> DISTO_PEAK_SEC_LEFT = Parameters.intParameter("DISTO_PEAK_SEC_LEFT", "position with maximum secondary peak to the left of main peak");
	    public static final Parameter<Integer> DISTO_PEAK_SEC_RIGHT = Parameters.intParameter("DISTO_PEAK_SEC_RIGHT", "position with maximum secondary peak to the right of main peak");

	}

