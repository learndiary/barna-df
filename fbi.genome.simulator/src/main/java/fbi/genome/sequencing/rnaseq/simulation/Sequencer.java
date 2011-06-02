package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.file.FileHelper;
import fbi.commons.io.IOHandler;
import fbi.commons.io.IOHandlerFactory;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.io.rna.FMRD;
import fbi.genome.model.Exon;
import fbi.genome.model.Gene;
import fbi.genome.model.Graph;
import fbi.genome.model.Transcript;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.sequencing.rnaseq.simulation.error.ModelPool;

import java.io.*;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

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
     * read count
     */
    private Hashtable<ByteArrayCharSequence, Long> map;
    /**
     * Read propability
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

        File zipFile = File.createTempFile("sim", "master.gz");
        zipFile.deleteOnExit();

        /*
        Write initial zip file and collect the line number. This represents
        the number of fragments written and is used to compute the priopability
        for a read to be sequenced
         */
        int noOfFragments = writeInitialFile(inFile, zipFile);
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
    int writeInitialFile(File libraryFile, File zipFile) throws IOException {
        if (libraryFile == null) {
            throw new NullPointerException("NULL library file not permitted");
        }
        if (zipFile == null) {
            throw new NullPointerException("NULL target file not permitted");
        }
        ByteArrayCharSequence cs = new ByteArrayCharSequence(100);
        IOHandler io = IOHandlerFactory.createDefaultHandler();
        ZipOutputStream zipOut = null;

        int nrOfFrags = 0;
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
                nrOfFrags++;
                Log.progress(currBytes, totBytes);

                cs.resetFind();
                ByteArrayCharSequence id = cs.getToken(2);
                if (!id.equals(lastID)) {
                    zipOut.putNextEntry(new ZipEntry(id.toString()));
                    lastID = id.cloneCurrentSeq();
                }
                zipOut.write(cs.chars, cs.start, cs.length());
                zipOut.write(BYTE_NL);
            }
            zipOut.close();
            Log.progressFinish(nrOfFrags + " lines zipped", true);
        } finally {
            io.close();
            if (zipOut != null) {
                try {
                    zipOut.close();
                } catch (IOException ignore) {
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
    public boolean loadErrors() { // todo : this should happen internally and protected or private
        // load model
        if (settings.get(FluxSimulatorSettings.ERR_FILE) != null) {
            babes = ModelPool.read(settings.get(FluxSimulatorSettings.ERR_FILE), settings);
            if (babes == null) {
                return false;
            }
            // check qualities for issued #48
            if (!babes.hasQualities() && settings.get(FluxSimulatorSettings.FASTQ)) {
                Log.warn("FastQ output requested, but the model does not support qualities. Disabled FastQ output!");
                settings.set(FluxSimulatorSettings.FASTQ, false);
            }
            return true;
        }
        return false;
    }


    boolean sequence(File zipFile, File referenceFile) {
        if (zipFile == null) {
            throw new NullPointerException("NULL initial file not permitted!");
        }

        try {
            Log.progressStart("sequencing");

            GFFReader reader = new GFFReader(referenceFile.getAbsolutePath());
            reader.setReadAheadLimit(500);
            reader.setSilent(true);
            reader.setStars(true);

            // vars
            Gene[] g;

            File tmpFile = File.createTempFile("flux", NAME_SEQ, settings.get(FluxSimulatorSettings.TMP_DIR));
            File tmpFasta = null;
            if (settings.get(FluxSimulatorSettings.FASTQ) && settings.get(FluxSimulatorSettings.GEN_DIR) != null) {
                tmpFasta = File.createTempFile("flux", NAME_SEQ, settings.get(FluxSimulatorSettings.TMP_DIR));
                Graph.overrideSequenceDirPath = settings.get(FluxSimulatorSettings.GEN_DIR).getAbsolutePath();
            }

            map = new Hashtable<ByteArrayCharSequence, Long>(profiler.size());

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
            for (reader.read(); (g = reader.getGenes()) != null; reader.read()) {
                for (Gene aG : g) {
                    processor.process(aG);
                }
            }
            // stats
            writer.close();

            Log.progressFinish();
            Log.message("");
            Log.message("\t" + processor.fragments + " fragments found");
            Log.message("\t" + writer.totalReads + " reads sequenced");
            Log.message("\t" + writer.countPolyAReads + " reads fall in poly-A tail");
            Log.message("\t" + writer.countTruncatedReads + " truncated reads");

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
            ProfilerFile.appendProfile(settings.get(FluxSimulatorSettings.PRO_FILE), ProfilerFile.PRO_COL_NR_SEQ, map);
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
        return !(babes == null || !babes.hasQualities());
    }

    private File getFASTAfile() {
        return FileHelper.replaceSfx(settings.get(FluxSimulatorSettings.SEQ_FILE), "." + (hasQualities() ? SFX_FASTQ : SFX_FASTA));
    }

    private void createQname(BEDobject2 obj2, ByteArrayCharSequence cs) {

        byte[] a = obj2.chars;
        cs.clear();
        int p1 = obj2.getNameP1(), p2 = obj2.getNameP2();
        cs.ensureLength(0, 1 + (p2 - p1));
        byte[] b = cs.chars;
        b[0] = (babes == null || !babes.hasQualities()) ? BYTE_GT : BYTE_AT;
        ++cs.end;
        assert (p1 > 0 && p2 > 0);
        System.arraycopy(a, p1, b, 1, p2 - p1);
        cs.end += p2 - p1;
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
     * @return bed the filled bed object
     */
    private BEDobject2 createReadPolyA(BEDobject2 obj, int start, int end, Transcript t, long molNr, byte absDir, int fragStart, int fragEnd, boolean left) {
        obj.clear();
        obj.append(CHR_POLYA);
        obj.append(BYTE_TAB);
        obj.append(0);
        obj.append(BYTE_TAB);
        obj.append(end - start);
        obj.append(BYTE_TAB);
        FMRD.appendReadName(obj, t, molNr, absDir, fragStart, fragEnd, start, end, settings.get(FluxSimulatorSettings.PAIRED_END), left);
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
     * @param start     the start
     * @param end       the end
     * @param t         the transcript
     * @param molNr     molecule number
     * @param absDir    direction
     * @param fragStart fragment start
     * @param fragEnd   fragment end
     * @param left      sens/antisense
     * @return polyA true if polya read
     */
    private boolean createRead(BEDobject2 obj, int start, int end, Transcript t, long molNr, byte absDir, int fragStart, int fragEnd, boolean left) {


        int originalStart = start;
        int originalEnd = end;
        int tlen = t.getExonicLength();
        boolean isPolyA = false;
        if (start > tlen) {
            createReadPolyA(obj, start, end, t, molNr, absDir, fragStart, fragEnd, left);    // read in polyA tail
            return true;
        }
        int offsStart = 0, offsEnd = 0;
        if (start < 1) {
            offsStart = start - 1; // t-coords 1-based, offs negative
            start = 1;
        }
        if (end > tlen) {
            offsEnd = end - tlen;    // positive, pos's after annotated end
            end = tlen;
            if (!left) {
                isPolyA = true;
            }
        } else if (end < 1) {
            offsEnd = end - 1;    // negative, pos's before annotated start
            end = 1;
        }


        // bed boundaries
        byte strand = t.getStrand();
        int bedStart = Math.abs(t.getGenomicPosition(start - 1)),
                bedEnd = Math.abs(t.getGenomicPosition(end - 1));    // uncorrected positions, 1-based
        int idxExA = t.getExonIdx(strand * bedStart),
                idxExB = t.getExonIdx(strand * bedEnd);    // exon indices
        if (idxExA == -1 || idxExB == -1) {
            Log.error("[INCONSISTENT] strand " + strand + " idx " + idxExA + ", " + idxExB);
        }
        bedEnd = offsEnd >= 0 ? bedEnd + (offsEnd * strand)
                : bedStart + (offsEnd * strand);    // use original bedstart, before!
        bedStart += offsStart * strand;            // correct out of range
        if (bedStart > bedEnd) {
            if (t.getStrand() >= 0) {
                Log.error("[INCONSISTENT] start " + bedStart + ", end " + bedEnd + ", strand " + t.getStrand());
            }
            int h = bedStart;
            bedStart = bedEnd;
            bedEnd = h;    // swap for neg strand
        }
        --bedStart; // lower the lower pos, 0-based

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
                originalStart, originalEnd, settings.get(FluxSimulatorSettings.PAIRED_END), left);
        obj.append(BYTE_TAB);
        obj.append(BYTE_0);
        obj.append(BYTE_TAB);
        obj.append(absDir >= 0 ? BYTE_PLUS : BYTE_MINUS);
        obj.append(BYTE_ARRAY_FROM_STRAND_TO_BLOCKS, 0,
                BYTE_ARRAY_FROM_STRAND_TO_BLOCKS.length);    // spare if there are blocks?
        if (idxExA == idxExB) {    // no blocks, spare?
            // spare?
            obj.append(BYTE_TAB);
            obj.append(BYTE_1);    // 1 block
            obj.append(BYTE_TAB);
            obj.append(bedEnd - bedStart);
            obj.append(BYTE_TAB);
            obj.append(BYTE_0);
            return isPolyA;
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
        return isPolyA;
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

    private ByteArrayCharSequence createQSeq(ByteArrayCharSequence cs, BEDobject2 obj, byte absDir, byte tDir, int len, int flen/*int fstart, int fend, Transcript t*/) {

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
            String baseID = gene.getGeneID() + FluxSimulatorSettings.SEP_LOC_TID;
            for (int j = 0; j < gene.getTranscripts().length; j++) {
                Transcript t = gene.getTranscripts()[j];
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
                        fragments++;
                        cs.set(s);
                        int fstart = cs.getTokenInt(0);
                        int fend = cs.getTokenInt(1);

                        double r = rnd.nextDouble();
                        if (r < p) {
                            double q = rndFiftyFifty.nextDouble();
                            if (pairedEnd || q < 0.5) {
                                writer.writeRead(true, t, fstart, fend, k);
                                ++cntPlus;
                            }
                            if (pairedEnd || q > 0.5) {
                                writer.writeRead(false, t, fstart, fend, k);
                                ++cntMinus;
                            }
                        }
                        ++k;
                    }
                } catch (IOException e) {
                    Log.error("Error while reading zip entry in sequencer: " + e.getMessage(), e);
                } finally {
                    if (buffy != null) {
                        try {
                            buffy.close();
                        } catch (IOException ignore) {
                        }
                    }
                    if (is != null) {
                        try {
                            is.close();
                        } catch (IOException ignore) {
                        }
                    }
                }
            }
        }

        public void close() {
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
         * @param t      the transcript
         * @param fstart the fragment start
         * @param fend   the fragment end
         * @param k      the molecule number
         * @throws IOException in case of any errors
         */
        public void writeRead(boolean left, Transcript t, int fstart, int fend, int k) throws IOException {
            byte absDir = (byte) (t.getStrand() >= 0 ? 1 : -1);
            byte antiDir = (byte) (t.getStrand() >= 0 ? -1 : 1);

            // check if the read is truncated
            int flen = fend - fstart + 1;
            if (flen < rLen) {
                ++countTruncatedReads;
            }

            // update profiler counts
            // todo: use the profiler to get the global ID ?
            ByteArrayCharSequence id = new ByteArrayCharSequence(t.getGene().getGeneID() + FluxSimulatorSettings.SEP_LOC_TID + t.getTranscriptID());
            if (map.containsKey(id)) {
                map.put(id, map.get(id) + 1);
            } else {
                map.put(id, long0);
            }

            totalReads++;


            // bed object
            if (bedOut != null) {
                if (left) {
                    boolean polyA = createRead(obj,
                            fstart, Math.min(fstart + settings.get(FluxSimulatorSettings.READ_LENGTH) - 1, fend),     // start, end
                            t, k, absDir,
                            fstart, fend, left);
                    if (polyA) {
                        countPolyAReads++;
                    }
                } else {
                    boolean polyA = createRead(obj,
                            Math.max(fend - settings.get(FluxSimulatorSettings.READ_LENGTH) + 1, fstart), fend,     // start, end
                            t, k, antiDir,
                            fstart, fend, left);
                    if (polyA) {
                        countPolyAReads++;
                    }
                }

                bedOut.write(obj.toString());
                bedOut.write("\n");
            }


            // fasta seq
            if (qFastaOut != null) {
                createQname(obj, cs);
                createQSeq(cs, obj, absDir, t.getStrand(), rLen, flen);
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
                }
            }
            if (qFastaOut != null) {
                try {
                    qFastaOut.close();
                } catch (IOException ignore) {
                }
            }
        }

    }

}
