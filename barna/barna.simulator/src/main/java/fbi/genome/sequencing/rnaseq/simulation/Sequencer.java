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

package fbi.genome.sequencing.rnaseq.simulation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.StringUtils;
import fbi.commons.io.IOHandler;
import fbi.commons.io.IOHandlerFactory;
import fbi.genome.io.FileHelper;
import fbi.genome.io.gtf.GTFwrapper;
import fbi.genome.io.rna.FMRD;
import fbi.genome.model.Exon;
import fbi.genome.model.Gene;
import fbi.genome.model.Graph;
import fbi.genome.model.Transcript;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.sequencing.rnaseq.simulation.error.MarkovErrorModel;
import fbi.genome.sequencing.rnaseq.simulation.error.ModelPool;
import fbi.genome.sequencing.rnaseq.simulation.error.QualityErrorModel;

/**
 * Sequence the library
 */
public class Sequencer implements Callable<Void> {
    public static final byte BYTE_TAB = '\t', BYTE_DOT = '.', BYTE_COMMA = ',', BYTE_0 = '0', BYTE_1 = '1', BYTE_PLUS = 43, BYTE_MINUS = '-', BYTE_GT = 62, BYTE_AT = 64, BYTE_NL = '\n';
    public final static String NAME_SEQ = "Sequencing";
    static final String SFX_FASTA = "fasta", SFX_FASTQ = "fastq";
    private static final byte BYTE_a = 97, BYTE_t = 116;
    
    private static final ByteArrayCharSequence CHR_POLYA = new ByteArrayCharSequence("polyA");
    /**
     * Helper to create only one instance for the initial profiler map
     */
    private static final Long long0 = 1l;

    /**
     * flux molecule identifier
     */
    private FluxSimulatorSettings settings;
    /**
     * The profiler
     */
    private Profiler profiler;
    /**
     * Error model
     */
    private ModelPool babes = null;

    /**
     * Map passed to the profiler mapping from global id to
     * other attributes (read count, coverage attributes...)
     */
    private Hashtable<CharSequence, Number[]> map;
    /**
     * Read probability
     */
    private double p = -1;
    /**
     * Store total number of reads after a run
     */
    private int totalReads = -1;

    public Sequencer(FluxSimulatorSettings settings, Profiler profiler) {
        this.settings = settings;
        this.profiler = profiler;
    }


    public Void call() throws Exception {
        Log.info("SEQUENCING", "getting the reads");
        File inFile = settings.get(FluxSimulatorSettings.LIB_FILE);

        File zipFile = FileHelper.createTempFile("sim", "master.gz");
        zipFile.deleteOnExit();

        /*
        Write initial zip file and collect the line number. This represents
        the number of fragments written and is used to compute the probability
        for a read to be sequenced
         */
        long noOfFragments = writeInitialFile(inFile, zipFile);
        p = settings.get(FluxSimulatorSettings.READ_NUMBER) / (double) noOfFragments;


        File referenceFile = settings.get(FluxSimulatorSettings.REF_FILE);

        sequence(zipFile, referenceFile);
        return null;
    }

    /**
     * Write and return the initial zip file
     *
     * @param libraryFile the library file
     * @param zipFile     the target file
     * @return lines lines zipped
     * @throws IOException in case of errors
     */
    long writeInitialFile(File libraryFile, File zipFile) throws IOException {
        if (libraryFile == null) {
            throw new NullPointerException("NULL library file not permitted");
        }
        if (zipFile == null) {
            throw new NullPointerException("NULL target file not permitted");
        }
        ByteArrayCharSequence cs = new ByteArrayCharSequence(100);
        IOHandler io = IOHandlerFactory.createDefaultHandler();
        ZipOutputStream zipOut = null;

        long nrOfFrags = 0;
        long lines = 0;
        try {
            Log.message("\tinitialize");
            Log.progressStart("zipping");

            // read the sorted file and put it in a zip form
            InputStream sortedIn = new FileInputStream(libraryFile);
            io.addStream(sortedIn);


            // the target stream
            zipOut = new ZipOutputStream(new FileOutputStream(zipFile));


            long totBytes = libraryFile.length(), currBytes = 0;
            ByteArrayCharSequence lastID = null;

            while (io.readLine(sortedIn, cs) != -1) {
                currBytes += cs.length() + 1;
                //nrOfFrags++;
                lines++;
                Log.progress(currBytes, totBytes);

                cs.resetFind();
                ByteArrayCharSequence id = cs.getToken(2);
                int dups = cs.getTokenInt(3);
                nrOfFrags+=dups;
                if (!id.equals(lastID)) {
                    zipOut.putNextEntry(new ZipEntry(id.toString()));
                    lastID = id.cloneCurrentSeq();
                }
                zipOut.write(cs.chars, cs.start, cs.length());
                zipOut.write(BYTE_NL);
            }
            Log.progressFinish(lines + " lines zipped ("+nrOfFrags +" fragments)", true);
        } finally {
            io.close();
            if (zipOut != null) {
                try {
                    zipOut.close();
                } catch (IOException ignore) {
                	throw new RuntimeException("Unable to close zip file "+ zipFile.getAbsolutePath(), ignore);
                }
            }
        }
        return nrOfFrags;
    }

    /**
     * Load the error model
     *
     * @return success true if loaded
     */
    public boolean loadErrors() {
        // load model
        String errorFile = settings.get(FluxSimulatorSettings.ERR_FILE);
        boolean fastOutput = settings.get(FluxSimulatorSettings.FASTA);
        if (fastOutput && errorFile != null && errorFile.length() > 0) {
            QualityErrorModel errorModel;
            try {
                InputStream input = null;
                String name = null;
                if(errorFile.equals("35")){
                    input = getClass().getResource("/35_error.model").openStream();
                    name = "35 bases model";
                }else if(errorFile.equals("76")){
                    input = getClass().getResource("/76_error.model").openStream();
                    name = "76 bases model";
                }else {
                    File file = new File(errorFile);
                    if(!file.canRead()){
                        throw new RuntimeException("unable to read error model from file " + errorFile);
                    }
                    input = new FileInputStream(file);
                    name = file.getAbsolutePath();
                }

                errorModel = MarkovErrorModel.loadErrorModel(name, input);
                babes = new ModelPool(settings.get(FluxSimulatorSettings.FASTA), errorModel);
            } catch (Exception e) {
                Log.error("Unable to load error model : " + e.getMessage(), e);
                throw new RuntimeException("Unable to load error model : " + e.getMessage(), e);
            }

            // check length
            int readLength = settings.get(FluxSimulatorSettings.READ_LENGTH);
            int modelLength = errorModel.getReadLength();
            if(readLength > modelLength){
                throw new RuntimeException("The error model supports a read length of " + modelLength + " but\n" +
                        "you are trying to create reads of length "+ readLength +"! This is not supported. Please \n" +
                        "use a different error model or reduce your read length!");
            }
            if(readLength < modelLength){
                Log.warn("The error model supports a read length of " + modelLength + " but\n" +
                        "you are trying to create reads of length " + readLength + "!");
            }
            return true;
        }else if( fastOutput){
            // just plain fasta
            babes = new ModelPool(true, null);
        }

        return false;
    }


    boolean sequence(File zipFile, File referenceFile) {
        if (zipFile == null) {
            throw new NullPointerException("NULL initial file not permitted!");
        }

        try {
            Log.progressStart("sequencing");

            GTFwrapper reader = new GTFwrapper(referenceFile.getAbsolutePath());
            reader.setReadAheadLimit(500);
            reader.setSilent(true);
            reader.setStars(true);

            // vars
            Gene[] g;

            File tmpFile = FileHelper.createTempFile("flux", NAME_SEQ);
            File tmpFasta = null;

            File genomeDir = settings.get(FluxSimulatorSettings.GEN_DIR);
            String errorModelFile = settings.get(FluxSimulatorSettings.ERR_FILE);
            boolean fasta = settings.get(FluxSimulatorSettings.FASTA);
            boolean hasErrorModel = errorModelFile != null && errorModelFile.trim().length() > 0;

            if (genomeDir != null && fasta) {
                tmpFasta = FileHelper.createTempFile("flux", NAME_SEQ);
                Graph.overrideSequenceDirPath = genomeDir.getAbsolutePath();
            }

            map = new Hashtable<CharSequence, Number[]>(profiler.size(), 1f);

            // init IO
            int readLength = settings.get(FluxSimulatorSettings.READ_LENGTH);
            boolean pairs = settings.get(FluxSimulatorSettings.PAIRED_END);

            // init the hash
            // create ZIP
            // hash entries
            Hashtable<CharSequence, ZipEntry> zipHash = new Hashtable<CharSequence, ZipEntry>(profiler.size());
            ZipFile zFile = new ZipFile(zipFile);
            Enumeration e = zFile.entries();
            ZipEntry ze;
            while (e.hasMoreElements()) {
                ze = (ZipEntry) e.nextElement();
                zipHash.put(ze.getName(), ze);
            }
            zFile.close();


            SequenceWriter writer = new SequenceWriter(tmpFile, tmpFasta, readLength);
            Processor processor = new Processor(writer, pairs, zipFile, zipHash);
            long fileLen= referenceFile.length();
            for (reader.read(); (g = reader.getGenes()) != null; reader.read()) {
                for (Gene aG : g) {
                    processor.process(aG);
                    Log.progress(reader.getBytesRead(), fileLen);
                }
            }
            // stats
            processor.close();
            writer.close();

            Log.progressFinish();
            Log.message("");
            Log.message("\t" + processor.fragments + " fragments found");
            Log.message("\t" + writer.totalReads + " reads sequenced");
            Log.message("\t" + writer.countPolyAReads + " reads fall in poly-A tail");
            Log.message("\t" + writer.countTruncatedReads + " truncated reads");
            if(tmpFasta != null && babes != null){
                Log.message("");
                Log.message("\tQuality stats: ");
                Log.message("\t" + StringUtils.fprint(babes.getAverageMutations(), 2)+" % average mutations per sequence");
                Log.message("\t" + StringUtils.fprint(babes.getAverageQuality(), 2)+" average quality ");
            }

            // store reads
            totalReads = writer.totalReads;

            Log.message("\n\tMoving temporary BED file");
            FileHelper.move(tmpFile, settings.get(FluxSimulatorSettings.SEQ_FILE));
            Log.progressFinish();

            if (tmpFasta != null) {
                Log.message("\n\tCopying qFasta file");
                File fileFASTA = getFASTAfile();
                FileHelper.move(tmpFasta, fileFASTA);
                Log.progressFinish();
            }
            
            // write profile
            HashMap<CharSequence, Number> map2= null;
            for (int i = 0; i < 3; i++) {
            	Iterator<CharSequence> iter= map.keySet().iterator();
            	if (map2== null)
            		map2= new HashMap<CharSequence, Number>(map.size(), 1f);
            	else
            		map2.clear();
            	while(iter.hasNext()) {
            		CharSequence id= iter.next();
            		Number[] n= map.get(id);
            		map2.put(id, n[i]);
            	}
            	int colNr= ProfilerFile.PRO_COL_NR_SEQ+ i;
            	if (i> 0)
            		++colNr; // add. col fraction
                ProfilerFile.appendProfile(settings.get(
                		FluxSimulatorSettings.PRO_FILE), 
                		colNr, 
                		map2, i== 0);
			}
             
            return true;
        } catch (Exception e) {
            Log.error("Error while sequencing : " + e.getMessage(), e);
            return false;
        } finally {
            if (zipFile != null) {
                if (!zipFile.delete()) {
                    Log.error("Unable to delete zip file " + zipFile.getAbsolutePath());
                }
            }
        }
    }

    private boolean hasQualities() {
        return babes != null && babes.hasErrorModel();
    }

    private File getFASTAfile() {
        return FileHelper.replaceSfx(settings.get(FluxSimulatorSettings.SEQ_FILE), "." + (hasQualities() ? SFX_FASTQ : SFX_FASTA));
    }

    public static void createQname(BEDobject2 obj2, ByteArrayCharSequence cs, ModelPool babes) {

        byte[] a = obj2.chars;
        cs.clear();
        int p1 = obj2.getNameP1(), p2 = obj2.getNameP2();
        cs.ensureLength(0, 1 + (p2 - p1));
        byte[] b = cs.chars;
        b[0] = (babes == null) ? BYTE_GT : BYTE_AT;
        ++cs.end;
        assert (p1 > 0 && p2 > 0);
        System.arraycopy(a, p1, b, 1, (p2 - p1));
        cs.end += (p2 - p1);
        cs.append(BYTE_NL);
    }


    /**
     * Create polya read
     *
     * @param obj       the bed object to fill
     * @param start     the read start
     * @param end       the read end
     * @param t         the transcript
     * @param molNr     the molecule number
     * @param absDir    the global direction
     * @param fragStart fragment start
     * @param fragEnd   fragment end
     * @param left      read direction
     * @param pairedEndSide paired end side, either 1 or 2 or 0 for no paired ends
     * @return bed the filled bed object
     */
    public static BEDobject2 createReadPolyA(BEDobject2 obj, int start, int end, Transcript t, long molNr, byte absDir, int fragStart, int fragEnd, boolean left, int pairedEndSide) {
        obj.clear();
        obj.append(CHR_POLYA);
        obj.append(BYTE_TAB);
        obj.append(0);
        obj.append(BYTE_TAB);
        obj.append(end - start+ 1);	// (+1) for end being excluded in BED, included in tx coordinates
        obj.append(BYTE_TAB);
        FMRD.appendReadName(obj, t, molNr, absDir, fragStart, fragEnd, start, end, left, pairedEndSide);
        obj.append(BYTE_TAB);
        obj.append(BYTE_0);
        obj.append(BYTE_TAB);
        obj.append((absDir >= 0 ? BYTE_PLUS : BYTE_MINUS));
        obj.append(BYTE_TAB);
        obj.append(BYTE_DOT);
        obj.append(BYTE_TAB);
        obj.append(BYTE_DOT);
        obj.append(BYTE_TAB);
        obj.append(BYTE_0);
        obj.append(BYTE_COMMA);
        obj.append(BYTE_0);
        obj.append(BYTE_COMMA);
        obj.append(BYTE_0);
        obj.append(BYTE_TAB);
        obj.append(BYTE_1);
        obj.append(BYTE_TAB);
        obj.append(end - start + 1);
        obj.append(BYTE_TAB);
        obj.append(BYTE_0);

        return obj;
    }

    public final static byte[] BYTE_ARRAY_FROM_STRAND_TO_BLOCKS = new byte[]{'\t', '.', '\t', '.', '\t', '0', ',', '0', ',', '0'};


    /**
     * Prepare a read and fill the given bed object
     *
     * @param obj       the bed object to fill
     * @param start     the start read start position in transcript coordinates (0-based)
     * @param end       the end read end position in transcript coordinates (0-based)
     * @param t         the transcript with genomic positions (1-based)
     * @param molNr     molecule number of the transcript
     * @param absDir    absolute direction of the transcript
     * @param fragStart fragment start start position of the fragment in transcript coordinates (0-based)
     * @param fragEnd   fragment end end position of the fragment in transcript coordinates (0-based)
     * @param left      flag indicating whether the read is the left end of the fragment
     * @param pairedEndSide either 1 or 2 or anything else if no pairedend reads are produced
     * @return polyA the number of nucleotides in the poly-A tail
     */
    public static int createRead(BEDobject2 obj, int start, int end, Transcript t, long molNr, byte absDir, int fragStart, int fragEnd, boolean left, int pairedEndSide) {


        int originalStart = start;
        int originalEnd = end;
        int tlen = t.getExonicLength();

        if (start > tlen- 1) {
            createReadPolyA(obj, start, end, t, molNr, absDir, fragStart, fragEnd, left, pairedEndSide);    // read in polyA tail
            return (end- start+ 1);
        }
        int offsStart = 0, offsEnd = 0;
        if (start < 0) {
            offsStart = start; // t-coords 0-based, offs negative
            start = 0;
        }
        if (end > tlen- 1) {
            offsEnd = end - tlen+ 1;    // positive, pos's after annotated end
            end = tlen- 1;
        } else if (end < 1) {
            offsEnd = end;    // negative, pos's before annotated start
            end = 0;
        }


        // bed boundaries
        byte strand = t.getStrand();
        // 0-based transcript coordinates
        int bedStart = Math.abs(t.getGenomicPosition(start)),
                bedEnd = Math.abs(t.getGenomicPosition(end)); 
        int idxExA = t.getExonIdx(strand * bedStart),
                idxExB = t.getExonIdx(strand * bedEnd);    // exon indices
        if (idxExA == -1 || idxExB == -1) {
            throw new RuntimeException("Invalid read (out of exons): " +
            		"tx strand/length ["+ strand+ ","+ tlen+ "], index first/second exon ["+ idxExA + "," + idxExB+ 
            		"], bedStartEnd ["+ bedStart+ ","+ bedEnd+ "], originalStartEnd ["+ originalStart+ ","+ originalEnd+
            		"], offsetStartEnd ["+ offsStart+ ","+ offsEnd+ "]");
        }
        bedEnd = offsEnd >= 0 ? bedEnd + (offsEnd * strand)
                : bedStart + (offsEnd * strand);    // use original bedstart, before!
        bedStart += offsStart * strand;            // correct out of range
        if (bedStart > bedEnd) {
            if (t.getStrand() >= 0) {
                throw new RuntimeException("Invalid read (end before start): " +
                		"tx strand/length ["+ strand+ ","+ tlen+ "], index first/second exon ["+ idxExA + "," + idxExB+ 
                		"], bedStartEnd ["+ bedStart+ ","+ bedEnd+ "], originalStartEnd ["+ originalStart+ ","+ originalEnd+
                		"], offsetStartEnd ["+ offsStart+ ","+ offsEnd+ "]");
            }
            int h = bedStart;
            bedStart = bedEnd;
            bedEnd = h;    // swap for neg strand
        }
        --bedStart; // lower the lower pos, BED:0-based

        // build object
        obj.clear();
        obj.append(t.getChromosome());
        obj.append(BYTE_TAB);
        obj.append(bedStart);
        obj.append(BYTE_TAB);
        obj.append(bedEnd);
        obj.append(BYTE_TAB);
        FMRD.appendReadName(obj, t,
                molNr, absDir, fragStart, fragEnd,
                originalStart, originalEnd, left, pairedEndSide);
        obj.append(BYTE_TAB);
        obj.append(BYTE_0);
        obj.append(BYTE_TAB);
        obj.append(absDir >= 0 ? BYTE_PLUS : BYTE_MINUS);
        obj.append(BYTE_ARRAY_FROM_STRAND_TO_BLOCKS, 0,
                BYTE_ARRAY_FROM_STRAND_TO_BLOCKS.length);    // spare if there are blocks?
        if (idxExA == idxExB) {    // no blocks, spare?
            obj.append(BYTE_TAB);
            obj.append(BYTE_1);    // 1 block
            obj.append(BYTE_TAB);
            obj.append(bedEnd - bedStart);
            obj.append(BYTE_TAB);
            obj.append(BYTE_0);
            return offsEnd;
        }

        if (strand < 0) {
            int h = idxExA;
            idxExA = idxExB;
            idxExB = h;    // swap for iterating oder left->right
        }

        int nr = Math.abs(idxExB - idxExA) + 1;
        obj.append(BYTE_TAB);
        obj.append(nr);
        Exon[] ee = t.getExons();    // idx in trpt dir
        for (int i = idxExA; ; i += (strand >= 0 ? 1 : -1)) {
            int gp = ((strand >= 0 && i == 0 && offsStart < 0) || (strand < 0 && i == (ee.length - 1) && offsEnd > 0)) ? bedStart + 1
                    : Math.max(Math.abs(ee[i].getStart()), bedStart + 1);
            int x = gp - (bedStart + 1);
            obj.setNextBlockStart(x);
            x = ((strand >= 0 && i == (ee.length - 1) && offsEnd > 0) || (strand < 0 && i == 0 && offsStart < 0) ? bedEnd
                    : Math.min(Math.abs(ee[i].getEnd()), bedEnd)) - gp + 1;
            obj.setNextBlockSize(x);
            if (i == idxExB) {
                break;
            }
        }
        return offsEnd;
    }


    /**
     * Returns an error message if something is broken or missing and null if everything is fine
     *
     * @return message error message or null
     */
    public String isReady() {
        if (settings == null) {
            return "No Setting specified!";
        }
        File file = settings.get(FluxSimulatorSettings.LIB_FILE);
        if (file == null || !file.exists()) {
            if (file == null) {
                return "No Fragmentation file specified! Check your parameters file!";
            }
            if (!file.exists()) {
                return "Fragmentation file " + file.getAbsolutePath() + " not found. Make sure fragmentation was done !";
            } else if (file.length() == 0) {
                return "Fragmentation file " + file.getAbsolutePath() + " is empty. Make sure fragmentation was done !";
            }
        }
        return null;
    }

    public static ByteArrayCharSequence createQSeq(ByteArrayCharSequence cs, BEDobject2 obj, byte absDir, byte tDir, int len, int flen, ModelPool babes) {

        //int flen= fend- fstart+ 1;
        cs.ensureLength(cs.end, len);
        int x;
        for (x = 0; x < CHR_POLYA.length(); ++x) {
            if (obj.charAt(x) != CHR_POLYA.charAt(x)) {
                break;
            }
        }
        int fStart = cs.start, seqStart = cs.end;    // start of fasta obj
        if (x != CHR_POLYA.length()) {        // not in pA, at least partially
            obj.readSequence(cs);    // not completely in poly-A
            cs.toUpperCase(seqStart, cs.end);
        }

        // Issue 36 (20100516): sequencing poly-A in fragments < readLen
        // change from len (readLength) to Math.min(flen, len)
        int diff = Math.min(flen, len) - (cs.end - seqStart);    // fill trailing As
        if (diff > 0) {


            // prevent polyA in the middle of a transcript
            if (absDir == tDir) {
                Arrays.fill(cs.chars, cs.end, cs.end + diff, BYTE_a);
            } else {
                System.arraycopy(cs.chars, fStart, cs.chars, fStart + diff, diff);
                Arrays.fill(cs.chars, fStart, fStart + diff, BYTE_t);
            }
            cs.end += diff;

            //++cntPolyA; // count only reads that fully fall into pA, see createRead()
        }

        // create qual seq
        if (babes != null) {
            //int mark= cs.end;
            babes.apply(cs, seqStart);
        }

        return cs;
    }

    /**
     * Get the total number of reads created in the last run
     *
     * @return reads total number of reads
     */
    public int getTotalReads() {
        return totalReads;
    }

    private Number[] getNumbersInit() {
		Number[] n= new Number[3];
		n[0]= new Long(0);
		n[1]= new Double(0);
		n[1]= new Double(0);
	
		return n;
	}

	/**
     * Process reads and pass them to the writer
     */
    class Processor {
        /**
         * Reader cache
         */
        private ByteArrayCharSequence cs;
        /**
         * The sequence writer
         */
        private SequenceWriter writer;
        /**
         * Paired end reads
         */
        private boolean pairedEnd;
        /**
         * The zip hash
         */
        private Hashtable<CharSequence, ZipEntry> zipHash;

        /**
         * Random sampler pick reads
         */
        private Random rnd = new Random();
        /**
         * sens or anti-sense sampler
         */
        private Random rndFiftyFifty = new Random();
        /**
         * Count fragments processed
         */
        private int fragments = 0;
        /**
         * Count sens reads written
         */
        private int cntPlus = 0;
        /**
         * Count antisens reads written
         */
        private int cntMinus = 0;
        /**
         * The zip file
         */
        private ZipFile zip;
        
        /**
         * Temporary array for transcript coverage
         */
        private int[] tmpCoverage= null; 



        public Processor(SequenceWriter writer, boolean pairedEnd, final File zipFile, final Hashtable<CharSequence, ZipEntry> zipHash) throws IOException {
            this.writer = writer;
            this.pairedEnd = pairedEnd;
            this.zipHash = zipHash;

            // init caches
            cs = new ByteArrayCharSequence(128);
            zip = new ZipFile(zipFile);

        }

        public void process(Gene gene) {
            
        	if (gene == null) {
                throw new NullPointerException("Null gene not permitted");
            }
            
            // process every transcript in the gene
            String baseID = gene.getGeneID() + FluxSimulatorSettings.SEP_LOC_TID;
            for (int j = 0; j < gene.getTranscripts().length; j++) {
                Transcript t = gene.getTranscripts()[j];
                int elen= t.getExonicLength();
                if (tmpCoverage== null|| tmpCoverage.length< elen) 
                	tmpCoverage= new int[elen];
                else
                	Arrays.fill(tmpCoverage, 0);
                String compID = baseID + t.getTranscriptID();


                ZipEntry ze = zipHash.remove(compID);
                if (ze == null) {
                    continue;    // not in frg file
                }


                InputStream is = null;
                BufferedReader buffy = null;
                try {
                    is = zip.getInputStream(ze);
                    buffy = new BufferedReader(new InputStreamReader(is));

                    int k = 0;
                    String s = null;
                    while ((s = buffy.readLine()) != null) {
                        cs.set(s);
                        int fstart = cs.getTokenInt(0);
                        int fend = cs.getTokenInt(1);
                        int dups = cs.getTokenInt(3);
                        // TODO ++dups ?
                        dups = Math.max(dups, 1);	// file provides nr. of duplicates, not molecules


                        double q = p * dups;
                        int frags = (int) q;	
                        double rest = q-frags;

                        // write fragments
                        for(int dd = 0; dd< frags;dd++){
                            fragments++;
                            ++k;
                            if(pairedEnd){
                                int dir = rndFiftyFifty.nextDouble() <= 0.5 ? 1:2;
                                writer.writeRead(true, dir,t, tmpCoverage, fstart, fend, k);
                                ++cntPlus;

                                writer.writeRead(false, dir==1?2:1, t, tmpCoverage, fstart, fend, k);
                                ++cntMinus;
                            }else{
                                if (rndFiftyFifty.nextDouble() < 0.5) {
                                    writer.writeRead(true, 0 ,t, tmpCoverage, fstart, fend, k);
                                    ++cntPlus;
                                }else{
                                    writer.writeRead(false, 0, t, tmpCoverage, fstart, fend, k);
                                    ++cntMinus;
                                }
                            }

                        }

                        // try for the rest
                        double r = rnd.nextDouble();
                        if (r < rest) {
                            if(pairedEnd){
                                int dir = rndFiftyFifty.nextDouble() <= 0.5 ? 1:2;
                                writer.writeRead(true, dir, t, tmpCoverage, fstart, fend, k);
                                ++cntPlus;

                                writer.writeRead(false,dir==1?2:1, t, tmpCoverage, fstart, fend, k);
                                ++cntMinus;
                            }else{
                                if (rndFiftyFifty.nextDouble() < 0.5) {
                                    writer.writeRead(true,0, t, tmpCoverage, fstart, fend, k);
                                    ++cntPlus;
                                }else{
                                    writer.writeRead(false,0, t, tmpCoverage, fstart, fend, k);
                                    ++cntMinus;
                                }
                            }
                            ++k;
                            fragments++;
                        }
                    } // end all fragments
                
                // catch I/O errors
                } catch (IOException e) {
                    Log.error("Error while reading zip entry in sequencer: " + e.getMessage(), e);
                } finally {
                    if (is != null) {
                        try {
                            is.close();                            
                        } catch (IOException ignore) {
                        	throw new RuntimeException("Could not close zip inputstream from "+ zip.getName(), ignore);
                        }
                    }
                    if (buffy != null) {
                        try {
                            buffy.close();
                        } catch (IOException ignore) {
                        	throw new RuntimeException("Could not close reader of zip inputstream from "+ zip.getName(), ignore);
                        }
                    }
                }

                Number[] n= null;
                if (map.containsKey(compID))
                	n= map.get(compID);
                else {
                	n= Sequencer.this.getNumbersInit();
                	map.put(compID, n);
                }

                // chi-square, exclude 0-positions
                double avgCov= 0d;
                for (int i = 0; i < elen; i++) {
                	if (tmpCoverage[i]== 0)
                		continue;
					avgCov+= tmpCoverage[i];
                }
                avgCov/= elen;
                double x2= 0d;
                for (int i = 0; i < elen; i++) { 
                	if (tmpCoverage[i]== 0)
                		continue;
					x2+= (tmpCoverage[i]- avgCov)* (tmpCoverage[i]- avgCov);
                }
				x2/= avgCov;
				n[1]= new Long(Math.round(x2));
                
				// CV, exclude 0-positions
				double mean= 0, min= Double.MAX_VALUE;
				int cnt= 0;
				for (int i = 0; i < elen; i++) {
					if (tmpCoverage[i]== 0)
						continue;
					++cnt;
					double a= tmpCoverage[i];
					// Anscombe residuals [Hansen et al. 2010]
					a= (3d/ 2d)* (Math.pow(tmpCoverage[i], 2d/3d)- Math.pow(avgCov, 2d/3d))/ Math.pow(avgCov, 1d/6d);
					mean+= a;
					if (a< min)
						min= a;
				}
				mean/= cnt;
				mean+= 2* Math.abs(min);
				double cv= 0;
				for (int i = 0; i < elen; i++) {
					if (tmpCoverage[i]== 0)
						continue;
					double a= tmpCoverage[i];
					a= (3d/ 2d)* (Math.pow(tmpCoverage[i], 2d/3d)- Math.pow(avgCov, 2d/3d))/ Math.pow(avgCov, 1d/6d);
					a+= 2* Math.abs(min);
					cv+= (a- mean)* (a- mean);
				}
				cv/= cnt;
				cv= Math.sqrt(cv);	// sdev
				cv/= mean;
				n[2]= new Double(cv);
                
            }	// end all transcripts
        }

        /**
         * Closes the zip stream and de-allocates temporary arrays.
         */
        public void close() {
        	tmpCoverage= null;
            if (zip != null) {
                try {
                    zip.close();
                } catch (IOException ignore) {
                }
            }
        }

    }

    /**
     * Write reads
     */
    class SequenceWriter {
        /**
         * the bed writer
         */
        private BufferedWriter bedOut;
        /**
         * the fastq writer
         */
        private BufferedWriter qFastaOut;
        /**
         * Cache BED object
         */
        private BEDobject2 obj;
        /**
         * Reader cache
         */
        private ByteArrayCharSequence cs;

        /**
         * The read length
         */
        private int rLen;
        /**
         * Number of polya reads
         */
        private int countPolyAReads;
        /**
         * Nunmber of trucated reads
         */
        private int countTruncatedReads;
        /**
         * Number of total readss
         */
        private int totalReads;


        public SequenceWriter(File bedFile, File qFasta, int rLen) throws IOException {
            this.rLen = rLen;
            if (bedFile != null) {
                bedOut = new BufferedWriter(new FileWriter(bedFile));
            }
            if (qFasta != null) {
                qFastaOut = new BufferedWriter(new FileWriter(qFasta));
            }
            // init caches
            cs = new ByteArrayCharSequence(128);
            obj = new BEDobject2(128);
        }

        /**
         * Process a transcript and write the read
         *
         * @param left   read direction
         * @param pairedEndSide either 1 or 2 (or anything if no pairedend reads are produced)
         * @param t      the transcript
         * @param fstart the fragment start
         * @param fend   the fragment end
         * @param k      the molecule number
         * @throws IOException in case of any errors
         */
        public void writeRead(boolean left,int pairedEndSide, Transcript t, int[] cov, int fstart, int fend, int k) throws IOException {
            byte absDir = (byte) (t.getStrand() >= 0 ? 1 : -1);
            byte antiDir = (byte) (t.getStrand() >= 0 ? -1 : 1);

            // check if the read is truncated
            int flen = fend - fstart + 1;
            if (flen < rLen) {
                ++countTruncatedReads;
            }

            // update profiler counts
            // TODO use the profiler to get the global ID ?
            ByteArrayCharSequence id = new ByteArrayCharSequence(t.getGene().getGeneID() + FluxSimulatorSettings.SEP_LOC_TID + t.getTranscriptID());
            Number[] n= null;
            if (map.containsKey(id)) 
            	n= map.get(id);
            else {
            	n= getNumbersInit();
            	map.put(id, n);
            }
            n[0]= new Long(n[0].longValue()+ 1);
            

            totalReads++;

            // coverage stats
            if (left) {
            	for (int i = Math.max(fstart, 0); i < Math.max(Math.min(fstart + rLen - 1, fend), 0); i++) {
					++cov[i];
				}
            } else {
            	for (int i = Math.max(Math.max(fend - rLen + 1, fstart), 0); i< Math.max(0, fend); i++) {
					++cov[i];
				}
            }

            // bed object
            if (bedOut != null) {
            	
            	int polyA= 0;
                if (left) {
                    polyA = createRead(obj,
                            fstart, Math.min(fstart + rLen - 1, fend),     // start, end
                            t, k, absDir,
                            fstart, fend, left, pairedEndSide);
                } else {
                    polyA = createRead(obj,
                            Math.max(fend - rLen + 1, fstart), fend,     // start, end
                            t, k, antiDir,
                            fstart, fend, left, pairedEndSide);
                }
                if (polyA> 0) {
                    countPolyAReads++;
                }

                bedOut.write(obj.toString());
                bedOut.write("\n");
            }


            // fasta seq
            if (qFastaOut != null) {
                createQname(obj, cs, babes);
                createQSeq(cs, obj, absDir, t.getStrand(), rLen, flen, babes);
                qFastaOut.write(cs.toString());
                qFastaOut.write("\n");
            }
        }


        /**
         * Close the writer
         */
        public void close() {
            if (bedOut != null) {
                try {
                    bedOut.close();
                } catch (IOException ignore) {
                	throw new RuntimeException(ignore);
                }
            }
            if (qFastaOut != null) {
                try {
                    qFastaOut.close();
                } catch (IOException ignore) {
                	throw new RuntimeException(ignore);
                }
            }
        }

    }

}
