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

package barna.flux.simulator;

import barna.commons.ByteArrayCharSequence;
import barna.commons.log.Log;
import barna.commons.system.OSChecker;
import barna.commons.utils.StringUtils;
import barna.flux.simulator.error.MarkovErrorModel;
import barna.flux.simulator.error.ModelPool;
import barna.flux.simulator.error.QualityErrorModel;
import barna.flux.simulator.fragmentation.FragmentDB;
import barna.io.FileHelper;
import barna.io.gtf.GTFwrapper;
import barna.model.Exon;
import barna.model.Gene;
import barna.model.Graph;
import barna.model.Transcript;
import barna.model.bed.BEDobject2;
import barna.model.commons.Coverage;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

/**
 * Sequence the library
 */
public class Sequencer implements Callable<Void> {
    public static final byte BYTE_TAB = '\t', BYTE_DOT = '.', BYTE_COMMA = ',', BYTE_0 = '0', BYTE_1 = '1', BYTE_PLUS = 43, BYTE_MINUS = '-', BYTE_GT = 62, BYTE_AT = 64, BYTE_NL = '\n';
    public final static String NAME_SEQ = "Sequencing";
    static final String SFX_FASTA = "fasta", SFX_FASTQ = "fastq";
    private static final byte BYTE_a = 97, BYTE_t = 116;
    public static final byte BYTE_DELIM_FMOLI = ':';
    public static final byte BYTE_DELIM_BARNA = '/';
    public final static byte[] BYTE_ARRAY_FROM_STRAND_TO_BLOCKS = new byte[]{'\t', '.', '\t', '.', '\t', '0', ',', '0', ',', '0'};


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

    /**
     * Number of fields written to .PRO file.
     */
    private int proFields = 4;

    /**
     * Create unique read ids for paired reads and skip the sense/antisense information
     */
    private boolean uniqueIds = false;

    public Sequencer(FluxSimulatorSettings settings, Profiler profiler) {
        this.settings = settings;
        this.profiler = profiler;
        if(settings != null){
            setUniqueIds(settings.get(FluxSimulatorSettings.UNIQUE_IDS));
        }
    }

    /**
     * Returns true if the sequencer generates unique ids for paired reads and
     * skips the A/S sense/anti-sense information from the read ids.
     *
     * @return uniqueIds true if unique paiered read ids are generated
     */
    public boolean isUniqueIds() {
        return uniqueIds;
    }

    /**
     * Set to true if the sequencer should generates unique ids for paired reads and
     * skips the A/S sense/anti-sense information from the read ids.
     *
     * @param uniqueIds generate unique reads
     */
    public void setUniqueIds(boolean uniqueIds) {
        this.uniqueIds = uniqueIds;
    }

    public Void call() throws Exception {
        Log.info("SEQUENCING", "getting the reads");
        File inFile = settings.get(FluxSimulatorSettings.LIB_FILE);

        /*
        Write initial zip file and collect the line number. This represents
        the number of fragments written and is used to compute the probability
        for a read to be sequenced
         */
        FragmentDB fragmentIndex = createFragmentIndex(inFile);
        p = Math.min(1,
                settings.get(FluxSimulatorSettings.READ_NUMBER) / (double) fragmentIndex.getNumberOfFragments());
        if (settings.get(FluxSimulatorSettings.PAIRED_END))
            p /= 2;  // half the probability for a fragment to be sequenced

        File referenceFile = settings.get(FluxSimulatorSettings.REF_FILE);

        sequence(fragmentIndex, referenceFile);
        return null;
    }

    /**
     * Write and return the initial zip file
     *
     * @param libraryFile the library file
     * @return db fragment index
     * @throws IOException in case of errors
     */
    FragmentDB createFragmentIndex(File libraryFile) throws IOException {
        if (libraryFile == null) {
            throw new NullPointerException("NULL library file not permitted");
        }
        FragmentDB index = new FragmentDB(libraryFile);
        index.setPrintStatus(true);
        index.createIndex();
        return index;
    }

    /**
     * Load the error model
     *
     * @return success true if loaded
     */
    public boolean loadErrors() {
        // load model
        String errorFile = settings.get(FluxSimulatorSettings.ERR_FILE);
        if (!errorFile.startsWith("/") || !new File(errorFile).exists()) {
            errorFile = new File(settings.getParameterFile().getParentFile(), errorFile).getAbsolutePath();
        }
        boolean fastOutput = settings.get(FluxSimulatorSettings.FASTA);
        if (errorFile != null && errorFile.length() > 0) {
            QualityErrorModel errorModel;
            try {
                InputStream input = null;
                String name = null;
                /*
                Fix BARNA-166 and make sure we only use the name of the file
                 */
                String fileName = new File(errorFile).getName();
                if (fileName.equals("35")) {
                    input = getClass().getResource("/35_error.model").openStream();
                    name = "35 bases model";
                } else if (fileName.equals("76")) {
                    input = getClass().getResource("/76_error.model").openStream();
                    name = "76 bases model";
                } else {
                    File file = new File(errorFile);
                    if (!file.canRead()) {
                        throw new RuntimeException("unable to read error model from file " + errorFile);
                    }
                    input = new FileInputStream(file);
                    name = file.getAbsolutePath();
                }

                errorModel = MarkovErrorModel.loadErrorModel(name, input);
                babes = new ModelPool(settings.get(FluxSimulatorSettings.FASTA), errorModel, settings.get(FluxSimulatorSettings.READ_LENGTH));
            } catch (Exception e) {
                Log.error("Unable to load error model : " + e.getMessage(), e);
                throw new RuntimeException("Unable to load error model : " + e.getMessage(), e);
            }

            // check length
            int readLength = settings.get(FluxSimulatorSettings.READ_LENGTH);
            int modelLength = errorModel.getReadLength();
            if (readLength != modelLength) {
                Log.warn("The error model supports a read length of " + modelLength + " but\n" +
                        "you are trying to create reads of length " + readLength +
                        ". We are scaling.");
            }

            return true;
        } else if (fastOutput) {
            // just plain fasta
            babes = new ModelPool(true, null, settings.get(FluxSimulatorSettings.READ_LENGTH));
        }

        return false;
    }


    boolean sequence(FragmentDB index, File referenceFile) {
        if (index == null) {
            throw new NullPointerException("NULL index not permitted!");
        }

        if (referenceFile == null) {
            throw new NullPointerException("NULL reference not permitted!");
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

            SequenceWriter writer = new SequenceWriter(tmpFile, tmpFasta, readLength);
            Processor processor = new Processor(writer, pairs, index);
            long fileLen = referenceFile.length();
            for (reader.read(); (g = reader.getGenes()) != null; reader.read()) {
                for (Gene aG : g) {
                    processor.process(aG);
                    Log.progress(reader.getBytesRead(), fileLen);
                }
            }
            // stats
            processor.close();
            writer.close();

            Log.progressFinish(StringUtils.OK, true);
            Log.message("");
            Log.message("\t" + index.getNumberOfFragments() + " fragments found (" + index.getNumberOfLines() + " without PCR duplicates)");
            Log.message("\t" + writer.totalReads + " reads sequenced");
            Log.message("\t" + writer.countPolyAReads + " reads fall in poly-A tail");
            Log.message("\t" + writer.countTruncatedReads + " truncated reads");
            if (tmpFasta != null && babes != null) {
                Log.message("");
                Log.message("\tQuality stats: ");
                String avg = StringUtils.fprint(babes.getAverageMutations(), 8);
                if (avg.equals("0.00000000") && babes.getAverageMutations() > 0.0) {
                    avg = Double.toString(babes.getAverageMutations());
                }
                Log.message("\t" + avg + " % average mutations per sequence");
                Log.message("\t" + StringUtils.fprint(babes.getAverageQuality(), 2) + " average quality ");
            }

            // store reads
            totalReads = writer.totalReads;

            Log.message("\n\tMoving temporary BED file");
            FileHelper.move(tmpFile, settings.get(FluxSimulatorSettings.SEQ_FILE));
            Log.progressFinish();

            if (tmpFasta != null) {
                Log.message("\n\tCopying Fasta file");
                File fileFASTA = getFASTAfile();
                FileHelper.move(tmpFasta, fileFASTA);
                Log.progressFinish();
            }

            // write profile
            HashMap<CharSequence, Number> map2 = null;
            for (int i = 0; i < proFields; i++) {
                Iterator<CharSequence> iter = map.keySet().iterator();
                if (map2 == null)
                    map2 = new HashMap<CharSequence, Number>(map.size(), 1f);
                else
                    map2.clear();
                while (iter.hasNext()) {
                    CharSequence id = iter.next();
                    Number[] n = map.get(id);
                    map2.put(id, n[i]);
                }
                int colNr = ProfilerFile.PRO_COL_NR_SEQ + i;
                if (i > 0)
                    ++colNr; // add. col fraction
                ProfilerFile.appendProfile(settings.get(
                        FluxSimulatorSettings.PRO_FILE),
                        colNr,
                        map2, i == 0);
            }

            return true;
        } catch (Exception e) {
            Log.error("Error while sequencing : " + e.getMessage(), e);
            throw new RuntimeException("Error while sequencing : " + e.getMessage(), e);
        } finally {
            index.close();
        }
    }

    private boolean hasQualities() {
        return babes != null && babes.hasErrorModel();
    }

    private File getFASTAfile() {
        return FileHelper.replaceSfx(settings.get(FluxSimulatorSettings.SEQ_FILE), "." + (hasQualities() ? SFX_FASTQ : SFX_FASTA));
    }

    public void createQname(BEDobject2 obj2, ByteArrayCharSequence cs, ModelPool babes) {

        byte[] a = obj2.chars;
        cs.clear();
        int p1 = obj2.getNameP1(), p2 = obj2.getNameP2();
        cs.ensureLength(0, 1 + (p2 - p1));
        byte[] b = cs.chars;
        b[0] = (babes == null) ? BYTE_GT : BYTE_AT;
        ++cs.end;
        assert (p1 > 0 && p2 > 0);
        try {
            System.arraycopy(a, p1, b, 1, (p2 - p1));
        } catch (ArrayIndexOutOfBoundsException e) {
            Log.error("Problem when generating read ID: " + p1 + ", " + p2 + ", " + a.length + ", " + b.length);
            throw new RuntimeException(e);
        }
        cs.end += (p2 - p1);
        //cs.append(BYTE_NL);
        cs.append(OSChecker.NEW_LINE);
    }

    public void appendReadName(BEDobject2 obj, Transcript t, int molNr, int fragStart, int fragEnd, boolean sense, int pairedEndSide) {
        // FURI
        obj.append(t.getGene().getLocusID());
        obj.append(BYTE_DELIM_FMOLI);
        obj.append(t.getTranscriptID());
        obj.append(BYTE_DELIM_FMOLI);
        obj.append((int) (molNr + 1));
        obj.append(BYTE_DELIM_FMOLI);
        obj.append(t.getExonicLength());
        obj.append(BYTE_DELIM_FMOLI);
        obj.append(fragStart);
        obj.append(BYTE_DELIM_FMOLI);
        obj.append(fragEnd);

        if (pairedEndSide == 1 || pairedEndSide == 2) {
            if(!isUniqueIds()){
                // always append sense antisens to single end reads
                obj.append(BYTE_DELIM_FMOLI);
                obj.append(sense ? "S" : "A");
            }

            obj.append(BYTE_DELIM_BARNA);
            if(isUniqueIds()){
                // make sure we use sense/antisense information fot the pairing
                // in case we generate unique reads
                obj.append(sense ? 1 : 2);
            }else{
                obj.append(pairedEndSide);
            }
        }else{
            // always append sense antisens to single end reads
            obj.append(BYTE_DELIM_FMOLI);
            obj.append(sense ? "S" : "A");
        }
    }

    /**
     * Create polya read
     *
     * @param obj           the bed object to fill
     * @param start         the read start
     * @param end           the read end
     * @param t             the transcript
     * @param molNr         the molecule number
     * @param absDir        the global direction
     * @param fragStart     fragment start
     * @param fragEnd       fragment end
     * @param left          read direction
     * @param pairedEndSide paired end side, either 1 or 2 or 0 for no paired ends
     * @return bed the filled bed object
     */
    public BEDobject2 createReadPolyA(BEDobject2 obj, int start, int end, Transcript t, int molNr, byte absDir, int fragStart, int fragEnd, boolean left, int pairedEndSide) {
        obj.clear();
        obj.append(CHR_POLYA);
        obj.append(BYTE_TAB);
        obj.append(0);
        obj.append(BYTE_TAB);
        obj.append(end - start + 1);    // (+1) for end being excluded in BED, included in tx coordinates
        obj.append(BYTE_TAB);
        appendReadName(obj, t, molNr, fragStart, fragEnd, left, pairedEndSide);
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




    /**
     * Prepare a read and fill the given bed object
     *
     * @param obj           the bed object to fill
     * @param start         the start read start position in transcript coordinates (0-based)
     * @param end           the end read end position in transcript coordinates (0-based)
     * @param t             the transcript with genomic positions (1-based)
     * @param molNr         molecule number of the transcript
     * @param aDir          absolute direction of the transcript
     * @param fragStart     fragment start start position of the fragment in transcript coordinates (0-based)
     * @param fragEnd       fragment end end position of the fragment in transcript coordinates (0-based)
     * @param left          flag indicating whether the read is the left end of the fragment
     * @param pairedEndSide either 1 or 2 or anything else if no pairedend reads are produced
     * @return polyA the number of nucleotides in the poly-A tail
     */
    public int createRead(BEDobject2 obj, int start, int end, Transcript t, int molNr, byte aDir, int fragStart, int fragEnd, boolean left, int pairedEndSide) {

        byte tDir = t.getStrand();
        int originalStart = start;
        int originalEnd = end;
        int tlen = t.getExonicLength();

        if (start > tlen - 1) {
            createReadPolyA(obj, start, end, t, molNr, aDir, fragStart, fragEnd, left, pairedEndSide);    // read in polyA tail
            return (end - start + 1);
        }
        int offsStart = 0, offsEnd = 0;
        if (start < 0) {
            offsStart = start; // t-coords 0-based, offs negative
            start = 0;
        }
        if (end > tlen - 1) {
            offsEnd = end - tlen + 1;    // positive, pos's after annotated end
            end = tlen - 1;
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
                    "tx strand/length [" + strand + "," + tlen + "], index first/second exon [" + idxExA + "," + idxExB +
                    "], bedStartEnd [" + bedStart + "," + bedEnd + "], originalStartEnd [" + originalStart + "," + originalEnd +
                    "], offsetStartEnd [" + offsStart + "," + offsEnd + "]");
        }
        bedEnd = offsEnd >= 0 ? bedEnd + (offsEnd * strand)
                : bedStart + (offsEnd * strand);    // use original bedstart, before!
/*        if (offsEnd> 0) {
            if (tDir> 0)
                bedEnd+= offsEnd;
            else
                bedStart-= offsEnd;
        } else {    // <0, before start
            if (tDir> 0)
                bedStart+= offsEnd;
            else
                bedEnd-= offsEnd;
        }*/
        bedStart += offsStart * strand;            // correct out of range
/*        if (tDir > 0)
            bedStart -= offsStart;
        else
            bedEnd += offsStart;*/
        if (bedStart > bedEnd) {
            if (t.getStrand() >= 0) {
                throw new RuntimeException("Invalid read (end before start): " +
                        "tx strand/length [" + strand + "," + tlen + "], index first/second exon [" + idxExA + "," + idxExB +
                        "], bedStartEnd [" + bedStart + "," + bedEnd + "], originalStartEnd [" + originalStart + "," + originalEnd +
                        "], offsetStartEnd [" + offsStart + "," + offsEnd + "]");
            }
            int h = bedStart;
            bedStart = bedEnd;
            bedEnd = h;    // swap for neg strand
        }
        --bedStart; // lower the lower pos, BED:0-based

        //bedStart= Math.max(0, bedStart);    // prevent underflow

        // build object
        obj.clear();
        obj.append(t.getChromosome());
        obj.append(BYTE_TAB);
        obj.append(bedStart);
        obj.append(BYTE_TAB);
        obj.append(bedEnd);
        obj.append(BYTE_TAB);
        appendReadName(obj, t,
                molNr, fragStart, fragEnd,
                left, pairedEndSide);
        obj.append(BYTE_TAB);
        obj.append(BYTE_0);
        obj.append(BYTE_TAB);
        obj.append(aDir >= 0 ? BYTE_PLUS : BYTE_MINUS);
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

    public ByteArrayCharSequence createQSeq(ByteArrayCharSequence cs, BEDobject2 obj, int t3p, byte tDir, int len, int flen, ModelPool babes) {

        //int flen= fend- fstart+ 1;
        cs.ensureLength(cs.end, len);
        // check whether polyA read,
        // compare against chr identifier
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
            int diff = (tDir > 0 ? obj.getEnd() : -(obj.getStart() + 1)) - t3p;
            // remove trailing N's outside of chr sequence if in polyA-tail
            if (diff > 0) {
                diff = Math.min(diff, cs.end - seqStart);
                if (obj.getStrand() < 0)
                    try {
                        System.arraycopy(cs.chars, seqStart + diff, cs.chars, seqStart, cs.end - (seqStart + diff));
                    } catch (ArrayIndexOutOfBoundsException e) {
                        System.err.println(obj);
                        System.err.println(t3p + "; " + diff + ", " + seqStart + ": " + cs.chars.length);
                        throw (e);
                    }
                cs.end -= diff;
            }
        }

        // Issue 36 (20100516): sequencing poly-A in fragments < readLen
        // change from len (readLength) to Math.min(flen, len)
        // to prevent polyA in the middle of a transcript
        int diff = Math.min(flen, len) - (cs.end - seqStart);    // fill trailing As
        if (diff > 0) {

            if (obj.getStrand() * tDir > 0) { // sense reads
                Arrays.fill(cs.chars, cs.end, cs.end + diff, BYTE_a);
            } else {
                //System.err.println("before:"+ cs.subSequence(seqStart, cs.end));
                try {
                    System.arraycopy(cs.chars, seqStart, cs.chars, seqStart + diff, cs.end - seqStart);
                } catch (ArrayIndexOutOfBoundsException e) {
                    System.err.println(diff + ", " + seqStart + ": " + cs.chars.length);
                }
                Arrays.fill(cs.chars, seqStart, seqStart + diff, BYTE_t);
                //System.err.println("after:"+ cs.subSequence(seqStart, cs.end));
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
         * The sequence writer
         */
        private SequenceWriter writer;
        /**
         * Paired end reads
         */
        private boolean pairedEnd;
        /**
         * the fragment index
         */
        private FragmentDB index;

        /**
         * Random sampler pick reads
         */
        private Random rnd = new Random();
        /**
         * sens or anti-sense sampler
         */
        private Random rndFiftyFifty = new Random();
        /**
         * Count sens reads written
         */
        private int cntPlus = 0;
        /**
         * Count antisens reads written
         */
        private int cntMinus = 0;

        /**
         * Coverage profile of a sequenced transcript.
         */
        private Coverage coverage = null;

        public Processor(SequenceWriter writer, boolean pairedEnd, FragmentDB index) throws IOException {
            this.writer = writer;
            this.pairedEnd = pairedEnd;
            this.index = index;
        }

        boolean noAmpWarn= false;
        public void process(Gene gene) {

            if (gene == null) {
                throw new NullPointerException("Null gene not permitted");
            }

            // process every transcript in the gene
            String baseID = gene.getLocusID() + FluxSimulatorSettings.SEP_LOC_TID;
            for (int j = 0; j < gene.getTranscripts().length; j++) {
                Transcript t = gene.getTranscripts()[j];
                int elen = t.getExonicLength();
                if (coverage == null)
                    coverage = new Coverage(elen);
                else
                    coverage.reset(elen);
                String compID = baseID + t.getTranscriptID();


                // get the iterator
                int readsSequenced = 0;
                try {

                    Iterable<ByteArrayCharSequence> entries = index.getEntries(compID);

                    int k = 0;
                    for (ByteArrayCharSequence cs : entries) {
                        int fstart = cs.getTokenInt(0);
                        int fend = cs.getTokenInt(1);
                        int dups = 0;
                        try {
                            dups= cs.getTokenInt(3);
                        } catch (IllegalArgumentException e) {
                            if (!noAmpWarn) {
                                Log.warn("No 3rd field (amplified molecules) found in library, assuming \'0\'.");
                                noAmpWarn= true;
                            }
                        }
                        dups = Math.max(dups, 1);    // file provides nr. of duplicates, not molecules


                        double q = p * dups;
                        int frags = (int) q;
                        double rest = q - frags;

                        // write fragments
                        for (int dd = 0; dd < frags; dd++) {
                            readsSequenced++;
                            ++k;
                            if (pairedEnd) {
                                int dir = rndFiftyFifty.nextDouble() <= 0.5 ? 1 : 2;
                                writer.writeRead(true, dir, t, coverage, fstart, fend, k);
                                ++cntPlus;

                                writer.writeRead(false, dir == 1 ? 2 : 1, t, coverage, fstart, fend, k);
                                ++cntMinus;
                            } else {
                                if (rndFiftyFifty.nextDouble() < 0.5) {
                                    writer.writeRead(true, 0, t, coverage, fstart, fend, k);
                                    ++cntPlus;
                                } else {
                                    writer.writeRead(false, 0, t, coverage, fstart, fend, k);
                                    ++cntMinus;
                                }
                            }

                        }

                        // try for the rest
                        double r = rnd.nextDouble();
                        if (r < rest) {
                            ++k;
                            if (pairedEnd) {
                                int dir = rndFiftyFifty.nextDouble() <= 0.5 ? 1 : 2;
                                writer.writeRead(true, dir, t, coverage, fstart, fend, k);
                                ++cntPlus;

                                writer.writeRead(false, dir == 1 ? 2 : 1, t, coverage, fstart, fend, k);
                                ++cntMinus;
                            } else {
                                if (rndFiftyFifty.nextDouble() < 0.5) {
                                    writer.writeRead(true, 0, t, coverage, fstart, fend, k);
                                    ++cntPlus;
                                } else {
                                    writer.writeRead(false, 0, t, coverage, fstart, fend, k);
                                    ++cntMinus;
                                }
                            }
                            readsSequenced++;
                        }
                    } // end all fragments

                    // catch I/O errors
                } catch (IOException e) {
                    throw new RuntimeException("Error while reading from fragment index sequencer: " + e.getMessage(), e);
                }

                Number[] n = null;
                if (map.containsKey(compID))
                    n = map.get(compID);
                else {
                    n = new Number[proFields];
                    map.put(compID, n);
                }
                // reads produced
                n[0] = new Integer(pairedEnd ? 2 * readsSequenced : readsSequenced);
                // covered positions
                n[1] = new Double(coverage.getFractionCovered());
                // chi-square, exclude 0-positions
                n[2] = new Long(coverage.getChiSquare(true, false));
                // CV, exclude 0-positions
                n[3] = new Double(coverage.getCV(true, true));

            }    // end all transcripts
        }

        /**
         * Closes the zip stream and de-allocates no longer needed objects.
         */
        public void close() {
            coverage = null;
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
         * @param left          read direction
         * @param pairedEndSide either 1 or 2 (or anything if no pairedend reads are produced)
         * @param t             the transcript
         * @param fstart        the fragment start
         * @param fend          the fragment end
         * @param k             the molecule number
         * @throws IOException in case of any errors
         */
        public void writeRead(boolean left, int pairedEndSide, Transcript t, Coverage cov, int fstart, int fend, int k) throws IOException {
            byte absDir = (byte) (t.getStrand() >= 0 ? 1 : -1);
            byte antiDir = (byte) (t.getStrand() >= 0 ? -1 : 1);

            // check if the read is truncated
            int flen = fend - fstart + 1;
            if (flen < rLen) {
                ++countTruncatedReads;
            }

            // update profiler counts
            // TODO use the profiler to get the global ID ?
            ByteArrayCharSequence id = new ByteArrayCharSequence(t.getGene().getLocusID() + FluxSimulatorSettings.SEP_LOC_TID + t.getTranscriptID());
            totalReads++;

            // coverage stats
            if (left) {
                for (int i = Math.max(fstart, 0); i < Math.max(Math.min(fstart + rLen - 1, fend), 0); i++) {
                    cov.increment(i);
                }
            } else {
                for (int i = Math.max(Math.max(fend - rLen + 1, fstart), 0); i < Math.max(0, fend); i++) {
                    cov.increment(i);
                }
            }

            // bed object
            if (bedOut != null) {

                int polyA = 0;
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
                if (polyA > 0) {
                    countPolyAReads++;
                }

                bedOut.write(obj.toString());
                bedOut.write(barna.commons.system.OSChecker.NEW_LINE);
            }


            // fasta seq
            if (qFastaOut != null) {
                createQname(obj, cs, babes);
                createQSeq(cs, obj, t.get3PrimeEdge(), t.getStrand(), rLen, flen, babes);
                qFastaOut.write(cs.toString());
                qFastaOut.write(barna.commons.system.OSChecker.NEW_LINE);
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
