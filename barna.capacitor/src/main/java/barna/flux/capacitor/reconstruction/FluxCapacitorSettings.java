/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.capacitor.reconstruction;

import barna.commons.parameters.*;
import barna.commons.utils.StringUtils;
import barna.io.RelativePathParser;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.constants.Constants;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.EnumSet;

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
     * Helper to parse relative filenames
     */
    static RelativePathParser relativePathParser = new RelativePathParser();

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
                if (file != null && !file.getParentFile().exists()) {
                    throw new ParameterException("Folder for output file " + file.getAbsolutePath() 
                    		+ " could not be found!");
                }

            }
        }, relativePathParser);


	    /**
	     * The file for outputting the learned profiles.
	     */
		public static final Parameter<File> PROFILE_FILE= Parameters.fileParameter("PROFILE_FILE", "The file for outputting profiles", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                if (file== null)
                	return;
                if (!file.getParentFile().exists()) {
                    throw new ParameterException("Folder for output file " + file.getAbsolutePath() 
                    		+ " could not be found!");
                }

            }
            
        }, relativePathParser);

	    /**
	     * The file for fragments of correctly paired reads (insert sizes).
	     */
	    public static final Parameter<File> INSERT_FILE = Parameters.fileParameter("INSERT_FILE", "The file for output of inserts", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                if (file == null)
                    return;
                if (!file.getParentFile().exists()) {
                    throw new ParameterException("Folder for output file " + file.getAbsolutePath()
                            + " could not be found!");
                }

                // check pre-conditions for read pairing
                UniversalReadDescriptor d = schema.get(READ_DESCRIPTOR);
                AnnotationMapping a = schema.get(ANNOTATION_MAPPING);
                if (!(d.isPaired() && (a.equals(AnnotationMapping.PAIRED) || a.equals(AnnotationMapping.COMBINED)))) {
                    throw new ParameterException("Read pairing required for annotating inserts: " +
                            (d.isPaired() ? ANNOTATION_MAPPING.getName() + " " + a.toString() : READ_DESCRIPTOR.getName() + " " + d.toString()));
                }

            }

        }, relativePathParser);

	    /**
	     * The number of reads mapping to the genomic reference.
	     */
	    public static final Parameter<Integer> NR_READS_MAPPED = Parameters.intParameter("NR_READS_MAPPED", "The number of reads mapped", -1, -1, Integer.MAX_VALUE, null);
	    
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
                if ((file != null) && (!file.exists()) && ((!file.getParentFile().exists()) || (!file.getParentFile().canWrite()))) {
                    throw new ParameterException("Parent folder " + file.getParentFile().getAbsolutePath()
                            + " to write coverage file " + file.getName() + " cannot be found or is write-protected.");
                }
            }
        }, relativePathParser);

	    /**
	     * The file where profiles are stored in.
	     */
	    public static final Parameter<File> STATS_FILE = Parameters.fileParameter("STATS_FILE", "The file to which the run characteristics are written", null, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File file = (File) schema.get(parameter);
                // if set for writing, check whether the parent directory is valid
                if ((file != null)&& (!file.exists())&& ((!file.getParentFile().exists())|| (!file.getParentFile().canWrite()))) {
                    throw new ParameterException("Parent folder " + file.getParentFile().getAbsolutePath() 
                    		+ " to write stats file "+ file.getName()+ " cannot be found or is write-protected.");
                }
            }
        }, relativePathParser);
        /**
         * If true, and a stats file is specified, the stats are added to the existing stats. This is useful, i.e, if you if you
         * run the chromosomes in single runs and want to sum up the stats in one file. NOTE that the stats
         * are actually added and the stats are summed up. You will 'loose' the old stats
         */
        public static final Parameter<Boolean> STATS_FILE_APPEND = Parameters.booleanParameter("STATS_FILE_APPEND", "Append to the stats file. \n" +
                "This adds results from this run to an existing capacitor stats file.", false);

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

        public static enum CountElements {SPLICE_JUNCTIONS,INTRONS};

        /**
         * Parameter for counting reads that falls into specific elements
         */
        public static final Parameter<EnumSet<CountElements>> COUNT_ELEMENTS = Parameters.enumSetParameter(
                "COUNT_ELEMENTS",
                " Count elements specified in the list",
                EnumSet.noneOf(CountElements.class),
                CountElements.class,
                null);

        /**
         * Parameter for skipping deconvolution
         */
        public static final Parameter<Boolean> NO_DECOMPOSE = Parameters.booleanParameter(
                "NO_DECOMPOSE",
                " Avoid running deconvolution step on the dataset",
                false,
                null);
	    
	    /**
	     * Load the setting from a file. NOTE that this does not validate the settings!
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
	            relativePathParser.setParentDir(f.getParentFile());
	            settings.parameterFile = f;
	            in = new FileInputStream(f);
	            settings.parse(in);
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
