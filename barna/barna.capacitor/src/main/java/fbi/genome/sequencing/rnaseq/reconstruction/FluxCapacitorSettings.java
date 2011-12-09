/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package fbi.genome.sequencing.rnaseq.reconstruction;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import fbi.commons.StringUtils;
import fbi.commons.parameters.FileNameParser;
import fbi.commons.parameters.Parameter;
import fbi.commons.parameters.ParameterException;
import fbi.commons.parameters.ParameterSchema;
import fbi.commons.parameters.ParameterValidator;
import fbi.commons.parameters.Parameters;
import fbi.genome.io.FileHelper;
import fbi.genome.io.rna.UniversalReadDescriptor;
import fbi.genome.model.constants.Constants;

/**
 * Container class for settings of the <code>FluxCapacitor</code>.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class FluxCapacitorSettings extends ParameterSchema {

	 protected static class UniversalReadDescriptorParameter extends Parameter<UniversalReadDescriptor> {
		 	ParameterException parseException;
			UniversalReadDescriptor descriptor;
			
			public UniversalReadDescriptorParameter() {
				super("READ_DESCRIPTOR",
					  " Expression how to parse the read IDs, or one of the shorthand names ("
						+ StringUtils.toString(UniversalReadDescriptor.getMapSimpleDescriptors().keySet(), ',')+ ")",
					  null,
					  UniversalReadDescriptor.class,
					  null);
			}
			
			public UniversalReadDescriptorParameter(UniversalReadDescriptorParameter anotherURDP) {
				super(anotherURDP);
				this.parseException= anotherURDP.parseException;
				this.descriptor= anotherURDP.descriptor;
			}
			 
			protected void set(UniversalReadDescriptor value) {
				descriptor= value;
			}
			
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
			
			protected void validate(ParameterSchema schema)
					throws ParameterException {
				if (parseException!= null)
					throw parseException;
			}

			protected UniversalReadDescriptor get() {
				return descriptor;
			}
			
			public Parameter copy() {
				UniversalReadDescriptorParameter clone= 
					new UniversalReadDescriptorParameter(this);  
				
				return clone;
			}
	 }
	
	
	 /**
	  * 
	  * @author Micha Sammeth (gmicha@gmail.com)
	  *
	  */
	 public static enum AnnotationMapping {
		 PAIRED, STRANDED, SINGLE, COMBINED	    
	 }
	 
	 /**
	 * Helper class to allow us to parse file names relative to a parent directory
	 */
	static class RelativePathParser implements FileNameParser {
	    /**
	     * the Parent directory
	     */
	    File parentDir = null;
	
	    @Override
	    public File parse(String string) {
	        return FileHelper.fromRelative(string, parentDir);
	    }
	}
	
    /**
     * Helper to parse relative filenames
     */
    private static RelativePathParser relativePathParser = new RelativePathParser();

	/**
	  * Descriptor with parsing info for the read IDs.
	  */
	 public static final Parameter<UniversalReadDescriptor> READ_DESCRIPTOR =
		 	new UniversalReadDescriptorParameter();
		 

	 /**
	  * Information used during annotation mapping
	  */
	 public static final Parameter<AnnotationMapping> ANNOTATION_MAPPING = Parameters.enumParameter(
			 "ANNOTATION_MAPPING", 
			 " Information from the read descriptor that will be used for annotation mapping", 
			 AnnotationMapping.SINGLE, 
			 new ParameterValidator() {
			        @Override
			        public void validate(final ParameterSchema schema, final Parameter parameter) throws ParameterException {
			        	
			        	UniversalReadDescriptor d= schema.get(READ_DESCRIPTOR);
			        	AnnotationMapping a= schema.get(ANNOTATION_MAPPING);
			        	
			        	// paired read descriptor requires paired-end descriptor, not vice versa
			        	if ((!d.isPaired())&& (a.equals(AnnotationMapping.PAIRED)|| a.equals(AnnotationMapping.COMBINED)))
			        		throw new ParameterException("Annotation mapping "+a + " requires a paired-end read descriptor!");
			        	// stranded annotation mapping requires stranded descriptor, not vice versa
			        	if ((!d.isStranded())&& a.equals(AnnotationMapping.STRANDED))
			        		throw new ParameterException("Annotation mapping "+a + " requires a stranded read descriptor!");
			        	if (a.equals(AnnotationMapping.COMBINED)&&
			        			(!(d.toString().contains(UniversalReadDescriptor.TAG_MATE1SENSE)
			        			|| d.toString().contains(UniversalReadDescriptor.TAG_MATE2SENSE)))) {
			        		
			        		throw new ParameterException("Annotation mapping "+a + " requires a read descriptor "
			        				+ UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE+ " or "+ UniversalReadDescriptor.DESCRIPTORID_MATE2_SENSE+ "!");
			        	}
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
        }, relativePathParser);

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
        }, relativePathParser);

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
        }, relativePathParser);

	    /**
	     * The log file.
	     */
	    public static final Parameter<File> STDERR_FILE = Parameters.fileParameter("STDERR_FILE", "The file for log messages", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                if (file== null)
                	return;
                if (!file.getParentFile().exists()) {
                    throw new ParameterException("Folder for log file " + file.getAbsolutePath() 
                    		+ " could not be found!");
                }

            }
        }, relativePathParser);


	    /**
	     * Flag to output coverage statistic
	     */
	    public static final Parameter<Boolean> COVERAGE_STATS = Parameters.booleanParameter("COVERAGE_STATS", "Flag to output coverage statistics", false, new ParameterValidator() {
	    		 @Override
	             public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
	                 boolean set= (Boolean) schema.get(parameter);
	                 File file= (File) schema.get(COVERAGE_FILE);
	                 if (set&& file== null)
	                	 throw new ParameterException("Parameter "+ COVERAGE_FILE.getName()
	                			 + " has to be set to output coverage statistics.");
	             }
	    });

	    /**
	     * The file where profiles are stored in.
	     */
	    public static final Parameter<File> COVERAGE_FILE = Parameters.fileParameter("COVERAGE_FILE", "The file to which coverage profiles are stored", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                // if set for writing, check whether the parent directory is valid
                if ((file != null)&& (!file.exists())&& ((!file.getParentFile().exists())|| (!file.getParentFile().canWrite()))) {
                    throw new ParameterException("Parent folder " + file.getParentFile().getAbsolutePath() 
                    		+ " to write coverage file "+ file.getName()+ " cannot be found or is write-protected.");
                }
            }
        }, relativePathParser);

	    /**
	     * The temporary directory.
	     */
	    public static final Parameter<File> TMP_DIR = Parameters.fileParameter("TMP_DIR", "The temporary directory", new File(System.getProperty("java.io.tmpdir")), new ParameterValidator() {
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                if (file == null) {
                    schema.set(parameter, new File(System.getProperty(Constants.PROPERTY_TMPDIR)));
                }
                if (!file.exists()) {
                    throw new ParameterException("The temporary directory " + file.getAbsolutePath() 
                    		+ " could not be found!");
                }
                if (!file.canWrite()) {
                    throw new ParameterException("The temporary directory " + file.getAbsolutePath()
                    		+ " is not writable!");
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
	            relativePathParser.parentDir = f.getParentFile();
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
	     * Flag whether sorted input files (annotation, mappings) should be kept,
	     * <b>iff</b> they were unsorted. 
	     */
	    public static final Parameter<Boolean> KEEP_SORTED_FILES = Parameters.booleanParameter("KEEP_SORTED_FILES", "Keeps input files ("+ 
	    		ANNOTATION_FILE.getName()+ ", "+ MAPPING_FILE+ ")", false);
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
