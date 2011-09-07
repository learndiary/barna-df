package fbi.genome.sequencing.rnaseq.reconstruction;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import fbi.commons.StringUtils;
import fbi.commons.parameters.Parameter;
import fbi.commons.parameters.ParameterException;
import fbi.commons.parameters.ParameterSchema;
import fbi.commons.parameters.ParameterValidator;
import fbi.commons.parameters.Parameters;
import fbi.genome.io.rna.UniversalReadDescriptor;
import fbi.genome.model.constants.Constants;

/**
 * Container class for settings of the <code>FluxCapacitor</code>.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class FluxCapacitorSettings extends ParameterSchema {

	 /**
	  * 
	  * @author Micha Sammeth (gmicha@gmail.com)
	  *
	  */
	 public static enum AnnotationMapping {
		 PAIRED, STRANDED, SINGLE, COMBINED	    
	 }
	 
	 /**
	  * Descriptor with parsing info for the read IDs.
	  */
	 public static final Parameter<UniversalReadDescriptor> READ_DESCRIPTOR = new Parameter<UniversalReadDescriptor>(
			 "READ_DESCRIPTOR",
			 " provides an expression how to parse the read IDs, or shorthand name ("
			 	+ StringUtils.toString(UniversalReadDescriptor.getMapSimpleDescriptors().keySet(), ',')+ ")",
			 null,
			 UniversalReadDescriptor.class,
			 null) {
		
		ParameterException parseException;
		UniversalReadDescriptor descriptor;
		 
		@Override
		protected void set(UniversalReadDescriptor value) {
			descriptor= value;
		}
		
		@Override
		protected void parse(String value) throws ParameterException {
			
			descriptor= new UniversalReadDescriptor();
			
			try {
				descriptor.init(value);
			} catch (RuntimeException e) {
				parseException= new ParameterException(
					this, value, e
				);
				throw parseException;
			}
			
		}
		
		@Override
		protected void validate(ParameterSchema schema)
				throws ParameterException {
			if (parseException!= null)
				throw parseException;
		}

		@Override
		protected UniversalReadDescriptor get() {
			return descriptor;
		}
		
		@Override
		public Parameter copy() {
			// TODO Auto-generated method stub
			return null;
		}
	};
		 

	 /**
	  * Information used during annotation mapping
	  */
	 public static final Parameter<AnnotationMapping> ANNOTATION_MAPPING = Parameters.enumParameter(
			 "ANNOTATION_MAPPING", 
			 " specifying information of the read descriptor that will be used for annotation mapping", 
			 AnnotationMapping.SINGLE, 
			 new ParameterValidator() {
			        @Override
			        public void validate(final ParameterSchema schema, final Parameter parameter) throws ParameterException {
			        	
			        	UniversalReadDescriptor d= schema.get(READ_DESCRIPTOR);
			        	AnnotationMapping a= schema.get(ANNOTATION_MAPPING);
			        	
			        	if ((!d.isPaired())&& (a.equals(AnnotationMapping.PAIRED)|| a.equals(AnnotationMapping.COMBINED)))
			        		throw new ParameterException("Annotation mapping "+a + " requires a paired-end read descriptor!");
			        	if ((!d.isStranded())&& a.equals(AnnotationMapping.STRANDED))
			        		throw new ParameterException("Annotation mapping "+a + " requires a stranded read descriptor!");
			        	if (a.equals(AnnotationMapping.COMBINED)&&
			        			(!(d.toString().equals(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE)))
			        			|| d.toString().equals(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE2_SENSE))))
			        		throw new ParameterException("Annotation mapping "+a + " requires a read descriptor "
			        				+ UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE+ " or "+ UniversalReadDescriptor.DESCRIPTORID_MATE2_SENSE+ "!");
			        }
	    	});
	 
	 	/**
	 	 * The file containing the annotation.
	 	 */
	    public static final Parameter<File> ANNOTATION_FILE = Parameters.fileParameter("ANNOTATION_FILE", "The annotation file", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                if (file == null) {
                    throw new ParameterException("You have to specify an annotation file");
                }
                if (!file.exists()) {
                    throw new ParameterException("The annotation file " + file.getAbsolutePath() 
                    		+ " could not be found!");
                }

            }
        });

	    /**
	     * The file containing the mapped reads.
	     */
	    public static final Parameter<File> MAPPING_FILE = Parameters.fileParameter("MAPPING_FILE", "The mapping file", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                if (file == null) {
                    throw new ParameterException("You have to specify an mapping file");
                }
                if (!file.exists()) {
                    throw new ParameterException("The mapping file " + file.getAbsolutePath() 
                    		+ " could not be found!");
                }

            }
        });

	    /**
	     * The file for default output.
	     */
	    public static final Parameter<File> STDOUT_FILE = Parameters.fileParameter("STDOUT_FILE", "The file for default output", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                if (!file.getParentFile().exists()) {
                    throw new ParameterException("Folder for output file " + file.getAbsolutePath() 
                    		+ " could not be found!");
                }

            }
        });

	    /**
	     * The log file.
	     */
	    public static final Parameter<File> STDERR_FILE = Parameters.fileParameter("STDERR_FILE", "The file for log messages", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                if (!file.getParentFile().exists()) {
                    throw new ParameterException("Folder for log file " + file.getAbsolutePath() 
                    		+ " could not be found!");
                }

            }
        });

	    /**
	     * The file where profiles are stored in.
	     */
	    public static final Parameter<File> PROFILE_FILE = Parameters.fileParameter("PROFILE_FILE", "The file to which profiles are stored", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                if (file == null) {
                    throw new ParameterException("You have to specify an file for profiles");
                }
                if (!file.exists()) {
                    throw new ParameterException("The file for profiles " + file.getAbsolutePath() 
                    		+ " could not be found!");
                }

            }
        });

	    /**
	     * The temporary directory.
	     */
	    public static final Parameter<File> TMP_DIR = Parameters.fileParameter("TMP_DIR", "The temporary directory", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                if (file == null) {
                    schema.set(parameter, new File(System.getProperty(Constants.PROPERTY_TMPDIR)));
                }
                if (!file.exists()) {
                    throw new ParameterException("The temporary directory " + file.getAbsolutePath() 
                    		+ " could not be found!");
                }

            }
        });
	    
	    /**
	     * Load the setting from a file
	     *
	     * @param f the file
	     * @return settings the loaded and validated settings
	     * @throws Exception in case of a parser or validation error
	     */
	    public static FluxCapacitorSettings createSettings(File f) throws Exception {
	        if (f == null) {
	            throw new NullPointerException("Null parameter file not permitted!");
	        }
	        if (!f.exists()) {
	            throw new IllegalArgumentException("Parameter file " + f.getAbsolutePath() + " can not be found!");
	        }
	        InputStream in = null;
	        try {

	            FluxCapacitorSettings settings = new FluxCapacitorSettings();
	            settings.parameterFile = f;
	            in = new FileInputStream(f);
	            settings.parse(in);
	            settings.validate();
	            return settings;
	        } finally {
	            if (in != null) {
	                try {
	                    in.close();
	                } catch (IOException e) {
	                	throw new RuntimeException(e);
	                }
	            }
	        }
	    }

	    /**
	     * A <code>boolean</code> value specifying whether locus sorting of reads 
	     * is carried out in RAM-memory or on disk.
	     */
	    public static final Parameter<Boolean> SORT_IN_RAM = Parameters.booleanParameter("SORT_IN_RAM", "Sort reads in RAM memory, not on disk", false);
	    
	    /**
	     * The parameter file
	     */
	    private File parameterFile;

	    /**
	     * Get the parameter file or null
	     *
	     * @return parameterFile the parameter file or null
	     */
	    public File getParameterFile() {
	        return parameterFile;
	    }
	    
	    public int getMaxThreads() {
	        return 1;
	    }


}
