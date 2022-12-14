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
import barna.io.FileHelper;
import barna.io.RelativePathParser;
import barna.model.Transcript;
import barna.model.constants.Constants;
import barna.model.rna.UniversalReadDescriptor;
import net.sf.samtools.SAMFileReader;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

/**
 * Container class for settings of the <code>FluxCapacitor</code>.
 *
 * @author Micha Sammeth (gmicha@gmail.com)
 */
public class FluxCapacitorSettings extends ParameterSchema {

    protected static class UniversalReadDescriptorParameter extends Parameter<UniversalReadDescriptor> {
        ParameterException parseException;
        UniversalReadDescriptor descriptor;

        public UniversalReadDescriptorParameter() {
            super("READ_DESCRIPTOR",
                    " Expression how to parse the read IDs, or one of the shorthand names ("
                            + StringUtils.toString(UniversalReadDescriptor.getMapSimpleDescriptors().keySet(), ',') + ")",
                    null,
                    UniversalReadDescriptor.class,
                    null);
        }

        public UniversalReadDescriptorParameter(UniversalReadDescriptorParameter anotherURDP) {
            super(anotherURDP);
            this.parseException = anotherURDP.parseException;
            this.descriptor = anotherURDP.descriptor;
        }

        protected void set(UniversalReadDescriptor value) {
            descriptor = value;
        }

        public UniversalReadDescriptor parse(String value) throws ParameterException {

            descriptor = UniversalReadDescriptor.createTestDescriptor();

            try {
                descriptor.init(value);
            } catch (RuntimeException e) {
                parseException = new ParameterException(
                        this, value, e
                );
                throw parseException;
            }
            return this.descriptor;
        }

        protected void validate(ParameterSchema schema)
                throws ParameterException {
            if (parseException != null)
                throw parseException;
            File mapping = schema.get(FluxCapacitorSettings.MAPPING_FILE);
            if (FileHelper.getExtension(mapping).toUpperCase().equals("BAM") && schema.get(READ_DESCRIPTOR) != null) {
                throw new ParameterException("You cannot specify a READ_DESCRIPTOR for BAM files");
            }
            if (FileHelper.getExtension(mapping).toUpperCase().contains("BED") && schema.get(READ_DESCRIPTOR) == null) {
                throw new ParameterException("You must specify a READ_DESCRIPTOR for BED files");
            }
        }

        protected UniversalReadDescriptor get() {
            if (descriptor != null)
                return descriptor;
            else
                return this.getDefault();
        }

        public Parameter copy() {
            UniversalReadDescriptorParameter clone =
                    new UniversalReadDescriptorParameter(this);
            clone.longOption(getLongOption()).shortOption(getShortOption());
            return clone;
        }

//            public String getValuesString() {
//                if (descriptor!=null)
//                    return descriptor.toString();
//                else
//                    return this.getDefault().toString();
//            }
    }


    /**
     *
     * @author Micha Sammeth (gmicha@gmail.com)
     *
     */

    /**
     * Enum for annotation mapping types
     */
    public static enum AnnotationMapping {
        AUTO, SINGLE, PAIRED, SINGLE_STRANDED, PAIRED_STRANDED;

        public boolean isStranded() {
            return (this.equals(SINGLE_STRANDED) || this.equals(PAIRED_STRANDED));
        }

        public boolean isPaired() {
            return (this.equals(PAIRED) || this.equals(PAIRED_STRANDED));
        }

        public boolean isSingle() {
            return (this.equals(SINGLE) || this.equals(SINGLE_STRANDED));
        }
    }

    /**
     * Sets annotation mapping according to the (non-)paired mapping reader and the read strandedness of the settings.
     * @param readerPaired <code>true</code> if read pairing was detected in the mapping file,
     *                     <code>false</code> otherwise
     * @return automatically detected annotation mapping value
     */
    public AnnotationMapping setAnnotationMappingAuto(boolean readerPaired) {
        FluxCapacitorSettings.ReadStrand readStrand = get(FluxCapacitorSettings.READ_STRAND);
        if (readerPaired) {
            if (readStrand.equals(FluxCapacitorSettings.ReadStrand.NONE)) {
                set(FluxCapacitorSettings.ANNOTATION_MAPPING, AnnotationMapping.PAIRED_STRANDED);
            } else {
                set(FluxCapacitorSettings.ANNOTATION_MAPPING, AnnotationMapping.PAIRED);
            }
        }
        else {
            if (readStrand.equals(FluxCapacitorSettings.ReadStrand.NONE)) {
                set(FluxCapacitorSettings.ANNOTATION_MAPPING, AnnotationMapping.SINGLE);
            } else {
                set(FluxCapacitorSettings.ANNOTATION_MAPPING, AnnotationMapping.SINGLE_STRANDED);
            }
        }

        return get(FluxCapacitorSettings.ANNOTATION_MAPPING);
    }

    /**
     * Lazy wrapper method
     * @return <code>true</code> if the annotation mapping is paired,
     * <code>false</code> otherwise
     */
    public boolean isPaired() {
        return get(ANNOTATION_MAPPING).isPaired();
    }

    /**
     * Lazy wrapper method
     * @return <code>true</code> if the annotation mapping is stranded,
     * <code>false</code> otherwise
     */
    public boolean isStranded() {
        return get(ANNOTATION_MAPPING).isStranded();
    }

    /**
     * Enum for strandedness types
     */
    public static enum ReadStrand {
        NONE, SENSE, ASENSE, MATE1_SENSE, MATE2_SENSE;

        public boolean isNone() {
            return this.equals(NONE);
        }

        public boolean isPaired() {
            return (this.equals(MATE1_SENSE) || this.equals(MATE2_SENSE));
        }

        public boolean isSingle() {
            return (this.equals(SENSE) || this.equals(ASENSE));
        }
    }

    /**
     * Helper to parse relative filenames
     */
    static RelativePathParser relativePathParser = new RelativePathParser();

    /**
     * Descriptor with parsing info for the read IDs.
     */
    public static final Parameter<UniversalReadDescriptor> READ_DESCRIPTOR =
            new UniversalReadDescriptorParameter().longOption("read-descriptor").shortOption('d');

    public static final Parameter<ReadStrand> READ_STRAND = Parameters.enumParameter(
            "READ_STRAND",
            " Information about read strandedness",
            ReadStrand.NONE,
            null);

    /**
     * Annotation sources that are disregarded during profiling.
     */
    public static final Parameter<List<String>> PROFILE_INCLUDE =
            Parameters.listParameter(
                    "PROFILE_INCLUDE",
                    " Annotation sources (GTF field 2) that are considered during profiling",
                    null,
                    new ParameterValidator() {
                        @Override
                        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                            FluxCapacitorSettings.checkMutuallyExclusive(schema.get(PROFILE_INCLUDE), schema.get(PROFILE_EXCLUDE));
                        }
                    });


    static protected void checkMutuallyExclusive(List<String> one, List<String> two) throws ParameterException {
        for (int i = 0; one!= null&& i < one.size(); i++) {
            for (int j = 0; two!= null&& j < two.size(); j++) {
                if (one.get(i).equals(two.get(j)))
                    throw new ParameterException("You cannot include AND exclude source "+ one.get(i)+ " for profiling!");
            }
        }

    }

    /**
     * Annotation sources that are disregarded during profiling.
     */
    public static final Parameter<List<String>> PROFILE_EXCLUDE =
            Parameters.listParameter(
                    "PROFILE_EXCLUDE",
                    " Annotation sources (GTF field 2) that are disregarded during profiling",
                    null,
                    new ParameterValidator() {
                        @Override
                        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                            FluxCapacitorSettings.checkMutuallyExclusive(schema.get(PROFILE_INCLUDE), schema.get(PROFILE_EXCLUDE));
                        }
                    });
    /**
     * Information used during annotation mapping
     */
    public static final Parameter<AnnotationMapping> ANNOTATION_MAPPING = Parameters.enumParameter(
            "ANNOTATION_MAPPING",
            " Information from the read descriptor that will be used for annotation mapping",
            AnnotationMapping.AUTO,
            new ParameterValidator() {
                @Override
                public void validate(final ParameterSchema schema, final Parameter parameter) throws ParameterException {

                    UniversalReadDescriptor d = schema.get(READ_DESCRIPTOR);
                    AnnotationMapping a = schema.get(ANNOTATION_MAPPING);
                    ReadStrand r = schema.get(READ_STRAND);
                    if (d != null) {
                        // paired read descriptor requires paired-end descriptor, not vice versa
                        if ((!d.isPaired()) && a.isPaired())
                            throw new ParameterException("Annotation mapping " + a + " requires a paired-end read descriptor!");
                        // stranded annotation mapping requires stranded descriptor, not vice versa
                        if ((!d.isStranded()) && a.isStranded()) {
                            if (r.isNone())
                                throw new ParameterException("Annotation mapping " + a + " requires a stranded read descriptor or strand information!");
                        }
//                            if (a.equals(AnnotationMapping.PAIRED_STRANDED)&&
//                                    (!(d.toString().contains(UniversalReadDescriptor.TAG_MATE1SENSE)
//                                    || d.toString().contains(UniversalReadDescriptor.TAG_MATE2SENSE)))) {
//
//                                throw new ParameterException("Annotation mapping " + a + " requires a read descriptor "
//                                        + UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE+ " or "+ UniversalReadDescriptor.DESCRIPTORID_MATE2_SENSE+ "!");
//                            }
                    } else {
                        if (a.isStranded() && r.isNone()) {
                            throw new ParameterException("Annotation mapping " + a + " requires strand information.");
                        }
                    }
                    if (a.isPaired() && r.isSingle())
                        throw new ParameterException("Annotation mapping " + a + " requires paired reads.");
                }
            }).longOption("annotation-mapping").shortOption('m');

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
    }, relativePathParser).longOption("annotation").shortOption('a');

    /**
     * The file containing the mapped reads.
     */
    public static final Parameter<File> MAPPING_FILE = Parameters.fileParameter("MAPPING_FILE", "The mapping file", null, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File file = (File) schema.get(parameter);
            if (file == null) {
                throw new ParameterException("You have to specify a mapping file");
            }
            if (!file.exists()) {
                throw new ParameterException("The mapping file " + file.getAbsolutePath()
                        + " could not be found!");
            }

        }
    }, relativePathParser).longOption("input").shortOption('i');

    /**
     * The file containing the read bias profile.
     */
    public static final Parameter<File> PROFILE_FILE = Parameters.fileParameter("PROFILE_FILE", "The profile file", null, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File file = (File) schema.get(parameter);
                /*if (file == null) {
                    throw new ParameterException("You have to specify a profile file");
                }
                if (!file.exists()) {
                    throw new ParameterException("The profile file " + file.getAbsolutePath()
                            + " could not be found!");
                }*/

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
    }, relativePathParser).longOption("output").shortOption('o');


    /**
     * The file for outputting the learned profiles.
     */
    public static final Parameter<File> PROFILE_OUTPUT = Parameters.fileParameter("PROFILE_OUTPUT", "The file for outputting profiles", null, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File file = (File) schema.get(parameter);
            if (file == null)
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
            File inputFile = schema.get(MAPPING_FILE);
            if (d != null && (inputFile != null && !inputFile.getName().toLowerCase().endsWith(".bam"))) {
                AnnotationMapping a = schema.get(ANNOTATION_MAPPING);
                if (!(d.isPaired() && a.isPaired())) {
                    throw new ParameterException("Read pairing required for annotating inserts: " +
                            (d.isPaired() ? ANNOTATION_MAPPING.getName() + " " + a.toString() : READ_DESCRIPTOR.getName() + " " + d.toString()));
                }
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
            if (file == null)
                return;
            if (!file.getParentFile().exists()) {
                throw new ParameterException("Folder for log file " + file.getAbsolutePath()
                        + " could not be found!");
            }

        }
    }, relativePathParser);


    /**
     * The file where profiles are stored in.
     */
    public static final Parameter<File> COVERAGE_FILE = Parameters.fileParameter("COVERAGE_FILE", "Calculate coverage profile write it to the specified file", null, new ParameterValidator() {
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
            if ((file != null) && (!file.exists()) && ((!file.getParentFile().exists()) || (!file.getParentFile().canWrite()))) {
                throw new ParameterException("Parent folder " + file.getParentFile().getAbsolutePath()
                        + " to write stats file " + file.getName() + " cannot be found or is write-protected.");
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

    public static enum CountElements {SPLICE_JUNCTIONS, INTRONS}

    ;

    /**
     * Parameter for counting reads that falls into specific elements
     */
    public static final Parameter<EnumSet<CountElements>> COUNT_ELEMENTS = Parameters.enumSetParameter(
            "COUNT_ELEMENTS",
            " Count specified elements. Possible elements are : " + Arrays.toString(CountElements.values()),
            EnumSet.noneOf(CountElements.class),
            CountElements.class,
            null);

    /**
     * Parameter for skipping deconvolution
     */
    public static final Parameter<Boolean> DISABLE_DECONVOLUTION = Parameters.booleanParameter(
            "DISABLE_DECONVOLUTION",
            "Disable the deconvolution step",
            false,
            null).longOption("disable-deconvolution");

    /**
     * Parameter for settting SAMtools validation stringency
     */
    public static final Parameter<SAMFileReader.ValidationStringency> SAM_VALIDATION_STRINGENCY = Parameters.enumParameter(
            "SAM_VALIDATION_STRINGENCY",
            " Set SAMtools validation stringency for validating records. One of STRICT|LENIENT|SILENT",
            SAMFileReader.ValidationStringency.SILENT,
            null).longOption("sam-validation-stringency");

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
     * A <code>boolean</code> value specifying whether the profiling
     * is carried out.
     */
        /*public static final Parameter<Boolean> PROFILE = Parameters.booleanParameter(
                "PROFILE",
                "Get the read bias profile",
                false,
                new ParameterValidator() {
                    @Override
                    public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                        boolean set = (Boolean)schema.get(parameter);
                        if (!set) {
                            File file = schema.get(PROFILE_FILE);
                            if (file == null) {
                                throw new ParameterException("You have to specify a profile file");
                            }
                        }
                    }
                }).longOption("profile").shortOption('p');*/

    /**
     * A <code>boolean</code> value specifying whether locus sorting of reads
     * is carried out in RAM-memory.
     */
    public static final Parameter<Boolean> SORT_IN_RAM = Parameters.booleanParameter("SORT_IN_RAM",
            "Sort reads in RAM memory, not on disk. This is the default behaviour. You can force " +
            "sorting on disk with the --sort-on-disk option",
            true).longOption("sort-in-ram").shortOption('r');

    /**
     * A <code>boolean</code> value specifying whether locus sorting of reads
     * is carried out on disk.
     */
    public static final Parameter<Boolean> SORT_ON_DISK = Parameters.booleanParameter("SORT_ON_DISK",
            "Sort reads on disk",
            false).longOption("sort-on-disk");

    /**
     * An <code>int</code> value specifying the minimum mapping score to use a mapping
     * for quantification
     */
    public static final Parameter<Integer> MIN_SCORE = Parameters.intParameter("MIN_SCORE",
            "Minimum mapping score. Mappings with score < min_score are discarded (mapq for BAM, score for BED)",
            -1).longOption("min-score").shortOption('q');

    /**
     * Minimum length of introns that are considered to be functional and not gaps/indels in genomic alignments of cDNA.
     * Neighboring exons in the same transcript with a distance &lt; min_ilen are joined.
     */
    public static final Parameter<Integer> MIN_ILEN = Parameters.intParameter("MIN_ILEN",
            "Minimum length of introns of the annotation that are considered real and not indels/gaps, " +
                    "\"introns\" with a length < MIN_ILEN are removed by joining the flanking exons",
            25, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            int val = (Integer) schema.get(parameter);
            if (val<= 0|| val> 255) {
                throw new ParameterException(parameter.getName()+ " has to be positive and <= 255");
            }
            Transcript.maxLengthIntronIsGap= (byte) (val- 1);
        }
    }).longOption("minilen");

    /**
     * Minimum number of reads to be left in a locus by the linear solver. If it's &gt;0 the lp-solver will be forced to
     * deconvolute the locus even if it would me more expensive that removing all the reads from it (BARNA-375).
     *
     * It can be specified as either:
     * - absolute number of reads (&gt;1)
     * - fraction of the observed reads (0,1)
     */
    public static final Parameter<Double> MIN_OBS = Parameters.doubleParameter("MIN_OBS",
            "Minimum number of reads to be left in a locus by the linear solver.\n" +
                    " It can be expressed as an absolute number of reads (>1) or as fraction of the observed read count [0..1].\n" +
                    "With MIN_OBS > 0 the lp-solver will be forced to provide a deconvolution in all loci with > 0 " +
                    "annotation-mapped reads or pairs. The derived quantifications may be variable in " +
                    "loci with less observations than changes that are necessary to perform the deconvolution.",
            0, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            double val = (Double) schema.get(parameter);
            if (val< 0) {
                throw new ParameterException(parameter.getName()+ " has to be positive");
            }
            Transcript.maxLengthIntronIsGap= (byte) (val- 1);
        }
    }).longOption("min-obs");

    /**
     * A <code>boolean</code> value specifying if the SAM flags have to be used to scan a BAM file
     * for quantification
     */
    public static final Parameter<Boolean> IGNORE_SAM_FLAGS = Parameters.booleanParameter("IGNORE_SAM_FLAGS",
            "Ignore SAM flags when scanning the mapping file",
            false).longOption("ignore-sam-flags");

    /**
     * A <code>boolean</code> value specifying if only primary alignments should be considered
     * for quantification
     */
    public static final Parameter<Boolean> SAM_PRIMARY_ONLY = Parameters.booleanParameter("SAM_PRIMARY_ONLY",
            "Only use primary alignments for quantification",
            false, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            boolean sam_primary_only = (Boolean)schema.get(parameter);
            boolean weighted = !schema.get(DISABLE_MULTIMAP_WEIGHTING);
            if (sam_primary_only && weighted)
                schema.set(FluxCapacitorSettings.DISABLE_MULTIMAP_WEIGHTING, true);
        }
    }).longOption("sam-primary-only");

    /**
     * A <code>boolean</code> value specifying if pairing information from the SAM file should be used
     * for quantification
     */
    public static final Parameter<Boolean> IGNORE_SAM_PAIRING_INFORMATION = Parameters.booleanParameter("IGNORE_SAM_PAIRING_INFORMATION",
            "Ignore SAM pairing information in the quantification",
            false, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            boolean sam_primary_only = schema.get(SAM_PRIMARY_ONLY);
            if (sam_primary_only)
                schema.set(parameter, true);
        }
    }).longOption("ignore-sam-pairing-information");

    /**
     * A <code>boolean</code> value specifying if only unique alignments should be considered
     * for quantification
     */
    public static final Parameter<Boolean> SAM_UNIQUE_ONLY = Parameters.booleanParameter("SAM_UNIQUE_ONLY",
            "Only use unique alignments for quantification",
            false).longOption("sam-unique-only");

    /**
     * A <code>boolean</code> value specifying to exclude file checking before the run
     */
    public static final Parameter<Boolean> DISABLE_FILE_CHECK = Parameters.booleanParameter("DISABLE_FILE_CHECK",
            "Disable scanning of input files before the run",
            false).longOption("disable-file-check");

    /**
     * A <code>boolean</code> value specifying to weight mapping counts by the number of multi-maps
     */
    public static final Parameter<Boolean> DISABLE_MULTIMAP_WEIGHTING = Parameters.booleanParameter("DISABLE_MULTIMAP_WEIGHTING",
            "Disable weighted counts for multi-maps",
            false, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            boolean sam_primary_only = schema.get(SAM_PRIMARY_ONLY);
            if (sam_primary_only)
                schema.set(parameter, true);
        }
    }).longOption("disable-multimap-weighting");

    /**
     * Flag whether sorted input files (annotation, mappings) should be kept,
     * <b>iff</b> they were unsorted.
     */
    public static final Parameter<File> KEEP_SORTED = Parameters.fileParameter("KEEP_SORTED", "Keeps input files (" +
            ANNOTATION_FILE.getName() + ", " + MAPPING_FILE + ")", null, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File file = (File) schema.get(parameter);
            if (file != null && !file.getParentFile().exists()) {
                throw new ParameterException("Folder for keeping sorted files " + file.getAbsolutePath()
                        + " could not be found!");
            }
        }
    }, relativePathParser);


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

    @Override
    public void validate() throws ParameterException {
        super.validate();
        if(this.get(FluxCapacitorSettings.SORT_ON_DISK)){
            this.set(FluxCapacitorSettings.SORT_IN_RAM, false);
        }
    }


    /**
     * Instance describing the attributes of reads and the readID encoding.
     */
    UniversalReadDescriptor descriptor;

    /**
     * Provides a non-null descriptor of the read attributes, constructs a new instance if not already initialized.
     * @return descriptor encapsulating the attributes of reads
     */
    public UniversalReadDescriptor getReadDescriptor() {

        if (descriptor== null) {

            // create a basic descriptor
            if (get(FluxCapacitorSettings.READ_DESCRIPTOR)== null)
                descriptor= new UniversalReadDescriptor((isPaired()? UniversalReadDescriptor.DESCRIPTORID_PAIRED:
                                                                    UniversalReadDescriptor.DESCRIPTORID_SIMPLE));
            else
                descriptor=  get(FluxCapacitorSettings.READ_DESCRIPTOR);

            // re-init for stranded protocols
            FluxCapacitorSettings.ReadStrand readStrand = get(FluxCapacitorSettings.READ_STRAND);
            if (!descriptor.isStranded() && !readStrand.equals(FluxCapacitorSettings.ReadStrand.NONE)) {
                if (descriptor.isPaired()) {
                    descriptor.init(readStrand.equals(FluxCapacitorSettings.ReadStrand.MATE1_SENSE) ?
                            UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE :
                            UniversalReadDescriptor.DESCRIPTORID_MATE2_SENSE);
                } else {
                    descriptor.init(readStrand.equals(FluxCapacitorSettings.ReadStrand.SENSE) ?
                            UniversalReadDescriptor.DESCRIPTORID_SENSE :
                            UniversalReadDescriptor.DESCRIPTORID_ANTISENSE);
                }
            }

        }

        return descriptor;
    }


}
