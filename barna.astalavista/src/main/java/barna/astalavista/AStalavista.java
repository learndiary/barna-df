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

package barna.astalavista;

import barna.commons.Execute;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.FluxTool;
import barna.commons.log.Log;
import barna.commons.parameters.ParameterException;
import barna.geneid.*;
import barna.io.GeneAheadReaderThread;
import barna.io.gtf.GTFwrapper;
import barna.model.*;
import barna.model.commons.MyFile;
import barna.model.constants.Constants;
import barna.model.splicegraph.*;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;


public class AStalavista implements FluxTool<Void>{

    AStalavistaSettings settings= null;

    public static void main(String[] args) {

        Execute.initialize(2);
        JSAP jsap = new JSAP();
        AStalavista myAsta= new AStalavista();
        List<Parameter> parameter = myAsta.getParameter();
        if(parameter != null){
            for (Parameter p : parameter) {
                try {
                    jsap.registerParameter(p);
                } catch (JSAPException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }
        try{
            JSAPResult toolParameter = jsap.parse(args);
            if (!myAsta.validateParameter(toolParameter)){
                System.exit(-1);
            }
        } catch (Exception e) {
            Log.error("Parameter error : " + e.getMessage(), e);
            e.printStackTrace();
            System.exit(-1);
        }

        Future<Void> captain= Execute.getExecutor().submit(myAsta);
        try {
            captain.get();
        } catch (InterruptedException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (ExecutionException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        Execute.shutdown();
    }

    public static Species species = new Species("human");

    static {
        species.setGenomeVersion("hg18");
    }

    public static void setSpecies(Species newSpecies) {
        species = newSpecies;
    }

    private String outputFname= null;


    /**
     * Input file containing with the transcriptome annotation in GTF.
     */
    private static File inputFile;


    // Genes read ahead
    private static int readAheadLimit= -1;  // TODO not static

    /**
     * Flag to only detect events that are not at the transcript edges.
     */
    static boolean onlyInternal= true;


    public static boolean DEBUG= false;

    /**
     * Counter for introns with invalid lengths/sequence attributes.
     */
    static long invalidIntrons= 0;

    /**
     * Number of introns analyzed.
     */
    static long totalIntrons= 0;

    /**
     * Flag to switch on/off intron checks.
     */
    static boolean acceptableIntrons= false;   // TODO not static, remove--it is always true

    /**
     * Writer for splice site scores.
     */
    BufferedWriter siteScoreWriter= null;

    /**
     * Parameters for geneID.
     */
    GParam geneidParam= null;

    /**
     * Variant hash, maps location to SNP base.
     */
    HashMap<String,String> variants= null;



    @Override
    public Void call() throws Exception {

        // TODO parse arguments
        //inputFile= new MyFile(parseArguments(SplicingGraph.writerThread, args).getAbsolutePath());

        //
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr11.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr6.gtf
        //
        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807.gtf
        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807_mRNAs_fromUCSC070919.gtf
        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807_mRNAs_fromUCSC070919_splicedESTs_fromUCSC070919.gtf
        boolean output= false, output2= true;
//		if (rusc)
//			outputFname= "delme.asta";

        EventExtractor.writerThread.start();

        // init and start threads
        long t0= System.currentTimeMillis();
        if (output2) {
            // writerThread
            // outputStats(SplicingGraph.writerThread, new OutputStreamWriter(System.err));
            // TODO write settings
            //Date ti= new Date(t0);
            //System.out.println("["+ti+"]  started, k= "+EventExtractor.n+" species "+EventExtractor.species+", input file "+inputFile.getAbsolutePath()+", output file= "+outputFname);
        }
        //GTFChrReader reader= new GTFChrReader(file.getAbsolutePath());
        //ChromosomeReaderThread readerThread= new ChromosomeReaderThread(reader);
        GTFwrapper reader= new GTFwrapper(inputFile.getAbsolutePath());
        if (!reader.isApplicable()) {
            inputFile= reader.sort();
            reader= new GTFwrapper(inputFile.getAbsolutePath());
        }
        if (readAheadLimit> 0)
            reader.setReadAheadLimit(readAheadLimit);
        reader.setNoIDs(null);
        //reader.sweepToChromosome("chr17");
        GeneAheadReaderThread readerThread= new GeneAheadReaderThread(reader);
        readerThread.setOutput(output);
        readerThread.setOutput2(output2);
        readerThread.start();

        while (true) {
            Gene[] g= readerThread.getGenes();
            if (g== null)
                break;

            String chromo = g[0].getChromosome();
            int evBefore = EventExtractor.counter;
            for (int i = 0; i < g.length; i++) {

                // score splice sites
                if (siteScoreWriter!= null) {
                    try {
                        g[i].markAlternativeSpliceSites();
                        scoreSites(g[i].getSpliceSites());
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                    continue;
                }

                if (g[i].getTranscriptCount() == 1) {
                    g[i] = null;
                    continue;
                }
//						if (g[i].getTranscriptCount()> 5000)  {
//							if (output2) {
//								Date ti= new Date(System.currentTimeMillis());
//								System.err.println("["+ti+"] "+chromo+" skipped locus "+g[i].getTranscripts()[0].getTranscriptID());
//							}
//							continue;
//						}

                g[i].setSpecies(species);

                EventExtractor extractor= new EventExtractor(g[i], settings);
                Thread extractorThread= new Thread(extractor);
                extractorThread.start();
                extractorThread.join();

            }

            System.err.println(chromo);

            // prepare for next batch
            if (output2) {
                Date ti = new Date(System.currentTimeMillis());
                int div = (int) (EventExtractor.cumulGC +
                        EventExtractor.cumulGF +
                        EventExtractor.cumulEV) / 1000;
                int frac = 0;
                if (div > 0)
                    frac = (EventExtractor.counter - evBefore) / div;
                else
                    frac = (EventExtractor.counter - evBefore);

                System.err.println("[" + ti + "] " + chromo +
                        " graph construct " + (EventExtractor.cumulGC / 1000) +
                        //" sec, fuzzy flanks "+(cumulGF/1000) +
                        " sec, contraction " + (EventExtractor.cumulGT / 1000) +
                        " sec, extract events " + (EventExtractor.cumulEV / 1000) +
                        " sec, found " + (EventExtractor.counter - evBefore) +
                        " events, " + frac + " ev/sec.");
            }
            System.gc();
            EventExtractor.writerThread.interrupt();

        }

        System.err.println("took "+((System.currentTimeMillis()- t0)/1000)+" sec.");

        if (siteScoreWriter!= null)
            siteScoreWriter.close();

        try {
            EventExtractor.writerThread.setKill(true);
            EventExtractor.writerThread.interrupt();
            EventExtractor.writerThread.join();
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
        }
        System.err.println("found "+ EventExtractor.counter+" events.");
        if (acceptableIntrons) {
            DecimalFormat df = new DecimalFormat("#.##");
            System.err.println("discarded " + invalidIntrons + " introns, " +
                    "found " + (totalIntrons - invalidIntrons) + " valid ones when checking splice sites: " +
                    "ratio (invalid/total) = " + df.format(((double) invalidIntrons) / totalIntrons));
        }


        return null;
    }

    /**
     * Provided with a splice site of a certain type and GeneID models,
     * the method branches to the right routine to compute the splice
     * site score.
     * @param spliceSite splice site to be scored
     * @param seq sequence of the splice site that is to be scored
     * @return the score according to the appropriate GeneID model,
     * or <code>NaN</code> if no the site corresponding GeneID model
     * is available.
     */
    protected float scoreSite(SpliceSite spliceSite, String seq) {

        seq= seq.toUpperCase();

        if (spliceSite.isDonor())
            return GeneID.scoreDonor(seq, geneidParam.getDonorProfile());
        if (spliceSite.isAcceptor())
            return GeneID.scoreAcceptor(seq, geneidParam.getAcceptorProfile(), null, null);

        return Float.NaN; // throw exception?
    }

    /**
     * Retrieves the score of splice sites as obtained from the genomic sequence,
     * and also of variants of those, if annotated.
     * @param spliceSites
     */
    protected void scoreSites(SpliceSite[] spliceSites) {

        String seq= null;
        float score;
        int flank5, flank3;
        for (int i = 0; i < spliceSites.length; i++) {

            if (spliceSites[i].isDonor()) {
                // prefix_donor = (offset+ order)
                // suffix_donor = (dimension- offset- order- 2)
                // DonorProfile: order= 1, offset= 1, dimension= 9
                // ==> prefix 2, suffix 5
                // "GC GT ACCCC"
                flank5= geneidParam.getDonorProfile().getOffset()+
                        geneidParam.getDonorProfile().getOrder();
                flank3= geneidParam.getDonorProfile().getDimension()-
                        geneidParam.getDonorProfile().getOffset()-
                        geneidParam.getDonorProfile().getOrder()- 2;

            } else if (spliceSites[i].isAcceptor()) {

                // prefix_acceptor = offset- 2
                // prefix_acceptor = (dimension- offset- 2)
                // AcceptorProfile: order= 1, offset= 24, dimension= 27
                // ==> prefix 24, suffix 3
                // "CTCTCTCTCTCTCTCTCTCTCT AG CGC"
                flank5= geneidParam.getAcceptorProfile().getOffset()- 2;
                flank3= geneidParam.getAcceptorProfile().getDimension()-
                        geneidParam.getAcceptorProfile().getOffset();

            } else
                continue;

            // get variants, if annotated
            Vector<String> vvec= getVariants(spliceSites[i], flank5, flank3);
            int nrCombinations= (int) Math.pow(2, vvec== null? 0: vvec.size());

            // arrays to store reference results and
            String[] seqs= new String[nrCombinations];
            float[] scores= new float[nrCombinations];
            String[] varTuples= new String[nrCombinations];

            // genomic sequence
            seq= Graph.readSequence(spliceSites[i], flank5, flank3);
            score= scoreSite(spliceSites[i], seq);
            // outputSite(spliceSites[i], null, seq, score);
            seqs[0]= seq;
            scores[0]= score;
            varTuples[0]= "reference";

            // recursion to do all tuple combinations
            if (vvec!= null) {

                if (vvec.size()> 1)
                    Log.warn("Site "+ spliceSites[i]+ " has "+ vvec.size()+ " variants.");

                // map to last included position, relative to exon boundary (splice site pos)
                if (spliceSites[i].isAcceptor()) {
                    flank5+= 2;
                } else {// donor
                    flank5-= 1;
                }
                int nr= scoreVariants(vvec, "", 0, 0, 0, spliceSites[i], flank5, flank3, seq.toLowerCase(), seqs, scores, varTuples);
                assert((nr+ 1)== nrCombinations);
            }

            outputSitesVCF(spliceSites[i], vvec, seqs, scores, varTuples);
        }

    }



    @Override
    public String getName() {
        return "astalavista";
    }

    @Override
    public String getDescription() {
        return "The AStalavista event retriever";
    }

    @Override
    public String getLongDescription() {
        return null;
    }


    @Override
    public List<Parameter> getParameter() {

        // converts parameter file parameters to CLI parameters

        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        AStalavistaSettings settings= new AStalavistaSettings();
        Collection<barna.commons.parameters.Parameter> pars=
                settings.getParameters().values();
        for (barna.commons.parameters.Parameter parameter : pars) {

            Class c= parameter.getType();
            Parameter p= null;
            if (c.equals(Boolean.class)) {
                p= JSAPParameters.switchParameter(
                        parameter.getLongOption(),
                        parameter.getShortOption())
                        .defaultValue(parameter.getDefault().toString())
                        .type(c)
                        .help(parameter.getDescription())
                        .get();
            } else {
                p= JSAPParameters.flaggedParameter(
                    parameter.getLongOption(),
                    parameter.getShortOption())
                    .type(c)
                    .help(parameter.getDescription())
                    .valueName(parameter.getName())
                    .get();
            }
            // TODO required() not implemented
            if (parameter.getLongOption()!= null|| parameter.getShortOption()!= 0)
                parameters.add(p);
        }

       return parameters;
    }

    /**
     * Scans for variants in the range of the splice site sequence. If SNPs are found,
     * the correspondingly mutated sequences are scored.
     *
     * @param ss the splice site
     * @param flank5 upstream sequence flanking the dinucleotide
     * @param flank3 downstream sequence flanking the dinucleotide
     */
    protected Vector<String> getVariants(SpliceSite ss, int flank5, int flank3) {

        if (variants== null)
            return null;

        String chr= ss.getGene().getChromosome();
        chr= chr.substring(3);  // curiosity with this VCF file omitting "chr" prefixes

        // boundaries, take into account the di-nucleotides
        int from= Math.abs(ss.getPos()- flank5- (ss.isAcceptor()? 2: -1));
        int to= Math.abs(ss.getPos()+ flank3+ (ss.isDonor()? 2: -1));

        Vector<String> vvar= new Vector<String>();
        for (int i = Math.min(from, to), j= 0; i <= Math.max(from, to); ++i, ++j) {

            // key: chrNr + @ + position
            String key= chr+ "@"+ Integer.toString(i);
            if (!(variants.containsKey(key)))
                continue;

            // val: snpID + @ + ref string + @ + variant string
            String val= variants.get(key);

            // VCF already provides strand-specific bases
            vvar.add(key+ "@"+ val);
        }

        if (vvar.size()== 0)
            return null;
        return vvar;

    }

    protected int scoreVariants(Vector<String> vvec, String varID, int rec, int nr, int idx, SpliceSite ss, int flank5, int flank3, String seq,
                                String[] seqs, float[] scores, String[] varTuples) {

        // break
        if (idx>= vvec.size())
            return nr;

        // recursion
        for (int i = idx; i < vvec.size(); ++i) {

            ++nr;

            // 0:chrNr, 1:position, 2:snpID, 3:ref string, 4:variant string
            String[] vv= vvec.elementAt(i).split("@");
            int snPos= Integer.parseInt(vv[1]);
            boolean deletion= vv[3].length()> vv[4].length();
            boolean insertion= vv[3].length()< vv[4].length();
            boolean substitution= vv[3].length()== vv[4].length();

            int del= vv[3].length()- 1;
            if (ss.getPos()< 0) {
                if (vv[4].length()> 1)  // || vv[4].length()> 0)
                    System.currentTimeMillis();
                snPos+= del;
                vv[3]= Graph.reverseSequence(Graph.complementarySequence(vv[3]));
                vv[4]= Graph.reverseSequence(Graph.complementarySequence(vv[4]));
            }
            int j= Math.abs(Math.abs(ss.getPos()- flank5)- snPos);

            // VCF already provides strand-specific bases
            String varSeq= null;
            // check
            if (rec== 0) {
                String sub= seq.substring(j, Math.min(j+ del+ 1, seq.length()));
                if (!(vv[3].substring(0, sub.length()).equalsIgnoreCase(sub)))
                    System.currentTimeMillis();
                //assert(vv[3].substring(0, sub.length()).equalsIgnoreCase(sub));
            }

            if (j+ vv[4].length()- del<= seq.length()) {

                varSeq= seq.substring(0, j)+ vv[4];
                if (j+ del< seq.length())
                    varSeq+= seq.substring(j+ del+ 1);

                // trim start
                int start= 0;
                if (insertion&& j< flank5)
                    start= vv[4].length()- vv[3].length();
                // trim end
                if (varSeq.length()> seq.length())  // dirty
                    varSeq= varSeq.substring(0, seq.length());
            } else
                varSeq= seq.substring(0, j)+ vv[4].substring(0, seq.length()- j);
            if (varSeq.length()< seq.length()) {
                int missing= seq.length()- varSeq.length();
                // fill with upstream seq
                if (j+ missing< flank5) {
                    String s= Graph.readSequence(ss, flank5+ missing- (ss.isAcceptor()? 2: 0), 0);
                    varSeq= s.substring(0, missing).toLowerCase()+ varSeq;
                } else {    // fill with downstream seq
                    int skip= j+ del+ 1- seq.length(); // to be skipped ds
                    String s= Graph.readSequence(ss, 0, flank3+ skip+ missing);
                    varSeq+= s.substring(2+ flank3+ skip).toLowerCase();
                }
            }
            float score= scoreSite(ss, varSeq);
            String vvarID= varID+ (varID.length()> 0? ",": "")+ vv[2];
            //outputSite(ss, vvarID, varSeq, score);
            seqs[nr]= varSeq;
            scores[nr]= score;
            varTuples[nr]= vvarID;

            // start recursion
            nr= scoreVariants(vvec, vvarID, rec+ 1, nr, i + 1, ss, flank5, flank3, varSeq, seqs, scores, varTuples);
        }

        return nr;

    }

    /**
     *
     * @param ss
     * @param variants
     * @param sequences
     * @param scores
     * @param varTuples
     */
    private void outputSitesVCF(SpliceSite ss, Vector<String> variants, String[] sequences, float[] scores, String[] varTuples) {

        // CHROM: number/letter without "chr"
        String chr= ss.getGene().getChromosome();
        if (chr.startsWith("chr"))
            chr= chr.substring(3);
        StringBuilder sb= new StringBuilder(chr);
        sb.append("\t");

        // POS: last/first exonic position flanking splice site di-nucleotide
        sb.append(Integer.toString(Math.abs(ss.getPos())));
        sb.append("\t");

        // ID: strand, coord, chromosome
        sb.append(ss.toString());
        sb.append(chr);
        sb.append("\t");

        // REF: genomic splice site sequence
        sb.append(sequences[0]);
        sb.append("\t");

        // VAR: missing value ".", or comma-separated list of variant sequences
        if (sequences.length== 1)
            sb.append(".\t");
        else {
            for (int i = 1; i < sequences.length; i++) {
                sb.append(sequences[i]);
                sb.append(",");
            }
            sb.replace(sb.length() - 1, sb.length(), "\t");
        }

        // SCORE: score of the reference splice site
        sb.append(Float.toString(scores[0]));
        sb.append("\t");

        // FILLTER: "PASS" if the site is possible, otherwise a semicolon-separated list of codes for filters that fail.
        // e.g. “q10;s50” might indicate that at this site the quality is below 10
        sb.append(scores[0]< (-1000)? "q-1000":"PASS");
        sb.append("\t");

        // INFO
        // modality: alternative/constitutive
        sb.append("MOD=");
        sb.append(ss.isAlternative()? "ALT;": "CON;");
        // variant list
        for (int i = 1; i < varTuples.length; ++i) {
            sb.append("ALT");
            sb.append(Integer.toString(i));
            sb.append("=");
            sb.append(varTuples[i]);
            sb.append(";");
        }
        if (scores.length> 1)
            sb.append("VAR_SCORES=");
        for (int i = 1; i < scores.length; ++i) {
            sb.append(Float.toString(scores[i]));
            sb.append(",");
        }
        if (scores.length> 1)
            sb.replace(sb.length()- 1, sb.length(), ";");
            // non-redundant snp list
        if (variants!= null&& variants.size()> 0) {
            sb.append("SNPS=");
            for (int i = 0; i < variants.size(); i++) {
                String[] a= variants.elementAt(i).split("@");
                sb.append(a[2]);
                sb.append(",");
            }
            sb.deleteCharAt(sb.length()- 1);
        }

        sb.append("\n");

        try {
            siteScoreWriter.write(sb.toString());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void outputSite(SpliceSite ss, String varID, String seq, float score) {

        String id= (varID== null? "genomic_sequence    ": varID);
        if (id.length()< 20) {
            char[] ext= new char[20- id.length()];
            Arrays.fill(ext, ' ');
            id= id+ String.copyValueOf(ext);
        }

        String scoreStr= Float.toString(score);
        if (scoreStr.length()< 15) {
            char[] ext= new char[15- scoreStr.length()];
            Arrays.fill(ext, ' ');
            scoreStr= scoreStr+ String.copyValueOf(ext);
        }

        String line= ss.getGene().getChromosome()+ "\t"+
                ss.toString()+ "\t"+
                (ss.isAlternative()? "ALT": "CON")+ "\t"+
                id+ "\t"+
                scoreStr+ "\t"+
                seq+ "\n";

        try {
            siteScoreWriter.write(line);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    @Override
    public boolean validateParameter(JSAPResult args) {

        // output help
        if (args.userSpecified(AStalavistaSettings.PRINT_PARAMETERS.getName())) {
            AStalavistaSettings settings= new AStalavistaSettings();
            settings.write(System.out);
            return false;
        }

        // create or merge settings
        try {
            settings= createSettings(settings, args);
        } catch (ParameterException e) {
            Log.error(e.getMessage(), e);
            return false;
        }

        // TODO from now on do all checks on the settings stub
        // TODO move "logic" tests of parameters to the settings ParameterSchema.validateParameters() and call that in the end
        if (settings.get(AStalavistaSettings.REF_FILE)== null) {
            Log.error("Hey, you forgot to specify a valid input file! ");
            return false;
        } else {
            inputFile= settings.get(AStalavistaSettings.REF_FILE);
            EventExtractor.writerThread= new WriterThread();
            if (EventExtractor.writerThread.outputFname== null&& (!SplicingGraph.writeStdOut)) {
                EventExtractor.writerThread.outputFname= inputFile.getAbsolutePath()+"_astalavista.gtf.gz";
            }
            if (EventExtractor.writerThread.outputFname!= null&& new MyFile(EventExtractor.writerThread.outputFname).exists()) {
                // Confirm o..+"\n by typing \'yes\':"
                Log.error("Overwriting output file "+ EventExtractor.writerThread.outputFname+".");
            }
        }

        // genome dir
        if (settings.get(AStalavistaSettings.GEN_DIR)== null) {
            // TODO acceptableIntrons= true;
        } else {
            Graph.overrideSequenceDirPath= settings.get(AStalavistaSettings.GEN_DIR).getAbsolutePath();
        }

        // init variants
        if (settings.get(AStalavistaSettings.VARIANTS)!= null) {
            variants= getVariants(settings.get(AStalavistaSettings.VARIANTS));
        }

        // events
        if (!settings.get(AStalavistaSettings.EVENT_TYPES).isEmpty()) {
            if (settings.get(AStalavistaSettings.EVENT_TYPES).contains(AStalavistaSettings.EventTypes.ASExt)
                    || settings.get(AStalavistaSettings.EVENT_TYPES).contains(AStalavistaSettings.EventTypes.DS)
                    || settings.get(AStalavistaSettings.EVENT_TYPES).contains(AStalavistaSettings.EventTypes.VS))

                onlyInternal= false;
            else
                onlyInternal= true;

            if (settings.get(AStalavistaSettings.EVENT_TYPES).contains(AStalavistaSettings.EventTypes.ASExt)
                    || settings.get(AStalavistaSettings.EVENT_TYPES).contains(AStalavistaSettings.EventTypes.ASInt))
                EventExtractor.retrieveASEvents= true;
            else
                EventExtractor.retrieveASEvents= false;

            if (settings.get(AStalavistaSettings.EVENT_TYPES).contains(AStalavistaSettings.EventTypes.DS))
                EventExtractor.retrieveDSEvents= true;
            else
                EventExtractor.retrieveDSEvents= false;

            if (settings.get(AStalavistaSettings.EVENT_TYPES).contains(AStalavistaSettings.EventTypes.VS))
                EventExtractor.retrieveVSEvents= true;
            else
                EventExtractor.retrieveVSEvents= false;
        }
        System.exit(0);


        // intron filtering
        if ((SplicingGraph.canonicalSS|| acceptableIntrons)&& Graph.overrideSequenceDirPath== null) {
            Log.error("You want me to check introns for valid/canonical splice sites, but you did not provide a valid sequence directory");
            return false;
        }

        // check splice site scoring stuff
        if (settings.get(AStalavistaSettings.SCORE_SITES)!= null) {

            try {
                siteScoreWriter= new BufferedWriter(new FileWriter(args.getFile(AStalavistaSettings.SCORE_SITES.getName())));
            } catch (IOException e) {
                Log.error(e.getMessage(), e);
            }

            if(settings.get(AStalavistaSettings.GEN_DIR)== null) {
                Log.error("Splice site scoring requires the genomic sequence, set parameter \'GEN_DIR\'");
                return false;
            }

            try {
                if(settings.get(AStalavistaSettings.GENEID_PARAM)== null) {
                    Log.warn("No GeneID parameter file for scoring models, using default");
                    geneidParam= Profile.readParam(GeneIDconstants.PARAMETERFILE, new GeneIDsettings())[0];
                } else
                    geneidParam= Profile.readParam(settings.get(AStalavistaSettings.GENEID_PARAM).getAbsolutePath(),
                            new GeneIDsettings())[0];
            } catch (Exception e) {
                Log.error(e.getMessage(), e);
                return false;
            }
        }


        return true;
    }

    /**
     * Method to convert command line arguments to parameter file stub,
     * possibly overwriting values set via additional parameter file.
     *
     * @param settings pre-defined settings or <code>null</code>
     * @param args parsed command line arguments
     * @return newly created or extended settings
     * @throws ParameterException in case a parameter does not get what
     * it expects
     */
    AStalavistaSettings createSettings(AStalavistaSettings settings, JSAPResult args) throws ParameterException {

        // lazily create
        if (settings== null)
            settings = new AStalavistaSettings();

        // input file
        if (args.userSpecified(AStalavistaSettings.REF_FILE.getLongOption()))
            settings.set(AStalavistaSettings.REF_FILE, args.getFile(AStalavistaSettings.REF_FILE.getLongOption()));

        // genome directory
        if (args.userSpecified(AStalavistaSettings.GEN_DIR.getLongOption())) {
            settings.set(AStalavistaSettings.GEN_DIR, args.getFile(AStalavistaSettings.GEN_DIR.getLongOption()));
        }

        // variants, VCF file
        if (args.userSpecified(AStalavistaSettings.VARIANTS.getLongOption())) {
            settings.set(AStalavistaSettings.VARIANTS, args.getFile(AStalavistaSettings.VARIANTS.getLongOption()));
        }

        // events
        if (args.userSpecified(AStalavistaSettings.EVENT_TYPES.getLongOption())) {
            String s= args.getString(AStalavistaSettings.EVENT_TYPES.getLongOption());
            settings.set(AStalavistaSettings.EVENT_TYPES.getName(),
                    args.getString(AStalavistaSettings.EVENT_TYPES.getLongOption()));
        }

        // splice site scoring stuff
        if (args.userSpecified(AStalavistaSettings.SCORE_SITES.getLongOption()))
            settings.set(AStalavistaSettings.SCORE_SITES, args.getFile(AStalavistaSettings.SCORE_SITES.getLongOption()));

        if (args.userSpecified(AStalavistaSettings.GENEID_PARAM.getLongOption()))
            settings.set(AStalavistaSettings.GENEID_PARAM, args.getFile(AStalavistaSettings.GENEID_PARAM.getLongOption()));

        // output options
//        if (args.userSpecified(AStalavistaSettings.OUTPUT_LOCUS.getName()))
//                // parses enum string and sets EnumSet correspondingly
//                settings.set(AStalavistaSettings.OUTPUT_LOCUS.getName(), args.getString(AStalavistaSettings.OUTPUT_LOCUS.getName()));

        return settings;
    }

    /**
     * Reads vcf file and fills a hash with position x snp information.
     * @param vcf file with the variants in vcf
     * @return hash representing the information of the provided file
     */
    protected HashMap<String, String> getVariants(File vcf) {

        try {
            HashMap<String, String> map= new HashMap<String, String>((int) (vcf.length()/ 1000));
            BufferedReader buffy= new BufferedReader(new FileReader(vcf));
            StringTokenizer t;
            for (String s= null; (s= buffy.readLine())!= null; ) {
                t= new StringTokenizer(s, "\t");
                String loc= t.nextToken()+ "@"+ t.nextToken();  // chrNr + @ + position
                String snpID= t.nextToken();
                String ref= t.nextToken();
                String var= t.nextToken();
                map.put(loc, snpID+ "@"+ ref+ "@"+ var); // snpID + @ + ref String + @ + var string
            }

            return map;
        } catch (Exception e) {
            Log.error("Error reading VCF file");
            throw new RuntimeException(e);
        }
    }

     public void outputStats(Writer writer) {
        BufferedWriter buffy= new BufferedWriter(writer);
        try {
            buffy.write("# started\t"+new Date(System.currentTimeMillis())+barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# input\t"+inputFile.getAbsolutePath()+barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# output");
            if (!SplicingGraph.writeStdOut)
                buffy.write("\t"+outputFname+barna.commons.system.OSChecker.NEW_LINE);
            else
                buffy.write("\tstdout\n");
            if (barna.model.Graph.overrideSequenceDirPath== null) {
                if (DEBUG)
                    buffy.write("# genome\t"+ species+ barna.commons.system.OSChecker.NEW_LINE);
            } else
                buffy.write("# genome\t"+ barna.model.Graph.overrideSequenceDirPath+barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# dimension\t"+ EventExtractor.n+barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# internalOnly\t"+ SplicingGraph.onlyInternal +barna.commons.system.OSChecker.NEW_LINE);
            //buffy.write("# canonicalSS "+canonicalSS+barna.commons.system.OSChecker.NEW_LINE);
            //buffy.write("# acceptableIntrons "+acceptableIntrons+barna.commons.system.OSChecker.NEW_LINE);
            if (SplicingGraph.acceptableIntrons)
                buffy.write("# intronConfidenceLevel "+SplicingGraph.intronConfidenceLevel+barna.commons.system.OSChecker.NEW_LINE);
            if (!SplicingGraph.onlyInternal)
                buffy.write("# edgeConfidenceLevel "+ Transcript.getEdgeConfidenceLevel()+barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# as_events\t");
            if (SplicingGraph.retrieveASEvents)
                buffy.write("true");
            else
                buffy.write("false");
            buffy.write(barna.commons.system.OSChecker.NEW_LINE);
            if (SplicingGraph.retrieveDSEvents) {
                buffy.write("# ds_events\t");
                if (SplicingGraph.retrieveDSEvents)
                    buffy.write("true");
                else
                    buffy.write("false");
            }
            buffy.write(barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# vs_events\t");
            if (SplicingGraph.retrieveVSEvents)
                buffy.write("true");
            else
                buffy.write("false");
            buffy.write(barna.commons.system.OSChecker.NEW_LINE);

            buffy.write(barna.commons.system.OSChecker.NEW_LINE);
            buffy.flush();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public java.io.File parseArguments(String[] args) {

        System.err.println("\nThis is ASta"
                +", graph-based AS event retriever of the AStalavista package.");

        boolean helpRequested= false;
        for (int i = 0; args!= null&& i < args.length; i++) {
            if (args[i].equals("--help"))
                helpRequested= true;
        }

        if (helpRequested|| args== null|| args.length< 2) {	// -i input file
            System.err.println("Here is a list of the options I understand:\n");
            System.err.println("-i, --input <input file>");
            System.err.println("This is a bit important, I cannot work without an input annotation. I want a GTF file with " +
                    "transcript annotations (exon features, with a mandatory optional attribute named \'transcript_id\') " +
                    "IN THE SAME COLUMN (i.e., if the transcript identifier of the 1st line is in column #10, it has to be in " +
                    "all lines of the file in column #10. The rest of the file should comply with the standard as specified " +
                    "at http://mblab.wustl.edu/GTF2.html.\n"+
                    "There may also be CDS features, but they become only interesting when checking for additional things " +
                    "as NMD probability etc.."+
                    barna.commons.system.OSChecker.NEW_LINE);
            System.err.println("-o, --output <output file|\'stdout\'>");
            System.err.println("Optional, the name of the output file (fully qualified path) OR the keyword \'stdout\' for " +
                    "writing everything to the standard output stream. " +
                    "If nothing is specified, the output will be written to a file \'<input file>_astalavista.gtf.gz\'. " +
                    barna.commons.system.OSChecker.NEW_LINE);
            System.err.println("-g, --genome <path to directory>");
            System.err.println("Path to the directory containing sequence files corresponding to the <seqname> field " +
                    "in the input GTF. A genome directory is required if a intron confidence value is specified." +
                    barna.commons.system.OSChecker.NEW_LINE);
            System.err.println("-k, --dimension <int value>");
            System.err.println("Dimension >1 of the events to be extracted. Default is 2 (i.e., \'pairwise events\'). " +
                    barna.commons.system.OSChecker.NEW_LINE);
            System.err.println("-tmp");
            System.err.println("Set temporary directory" +
                    barna.commons.system.OSChecker.NEW_LINE);
            System.err.println("-ext, +ext");
            System.err.println("(De-)activate external events, i.e. events that include the transcript start or the " +
                    "poly-adenylation site" +
                    barna.commons.system.OSChecker.NEW_LINE);
            System.err.println("-ic, --intronConfidence [int value]");
            System.err.println("Level of intron confidence. The default is to trust all introns. Introns are assigned " +
                    "a confidency class:\n" +
                    "\t 0 if 'RefSeq' appears in the source field of the annotation\n" +
                    "\t 1 if 'mRNA' appears in the source field of the annotation\n" +
                    "\t 2 if 'EST' appears in the source field of the annotation\n" +
                    "\t 3 if if none of the above applies\n" +
                    "all introns of confidency level > intronConfidence will be checked for proper splice sites when extracting events." +
                    barna.commons.system.OSChecker.NEW_LINE);
            System.err.println("-as, +as");
            System.err.println("Deactivate (\'-as\') or activate (\'+as\') the retrieval of Alternative Splicing events. See documentation " +
                    "for the definition of events that suffice an alternative splicing event." +
                    barna.commons.system.OSChecker.NEW_LINE);
            System.err.println("-ds, +ds");
            System.err.println("Deactivate (\'-ds\') or activate (\'+ds\') the retrieval of aDditional splicing events. See documentation " +
                    "for the definition of events that suffice an additional splicing event." +
                    barna.commons.system.OSChecker.NEW_LINE);
            System.err.println("-s, --seqsite");
            System.err.println("Output splice site sequences with events. Requires a reference genome."+
                    barna.commons.system.OSChecker.NEW_LINE);
            System.err.println("--flankType");
            System.err.println("Output the type of the event flanks, i.e., \'constitutive\' or \'alternative\'."+
                    barna.commons.system.OSChecker.NEW_LINE);

            // reactivated on 20100112
            System.err.println("-ec, --edgeConfidence [int value]");
            System.err.println("Level of confidence for edges (i.e., annotated transcription starts/poly-adenylation sites). " +
                    "The default is to trust no annotated edge and to extend overlapping first/last exons of a transcript to " +
                    "their most extreme position. :\n" +
                    "\t 0 if 'RefSeq' appears in the source field of the annotation\n" +
                    "\t 1 if 'mRNA' appears in the source field of the annotation\n" +
                    "\t 2 if 'EST' appears in the source field of the annotation\n" +
                    "\t 3 if if none of the above applies\n" +
                    "all transcript edges of confidency level > edgeConfidence will be extended in case the annotation shows " +
                    "another exon with the same adjacent splice site and an earlier/later start/end." +
                    barna.commons.system.OSChecker.NEW_LINE);
            System.err.println("AStalavista.");
            System.exit(0);
        }

        java.io.File file= null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-i")|| args[i].equals("--input")) {
                if (i+1>= args.length) {
                    System.err.println("Hey, you forgot to specify the input file!");
                    System.exit(-1);
                }
                file= new java.io.File(args[++i]);
                continue;
            }
            if (args[i].equals("-s")|| args[i].equals("--seqsite")) {
                WriterThread.outputSeq= true;
                continue;
            }
            if (args[i].equals("-o")|| args[i].equals("--output")) {
                String s= args[++i];
                if (s.equalsIgnoreCase("stdout"))
                    SplicingGraph.writeStdOut= true;
                else
                    outputFname= s+".gz";
                continue;
            }
            if (args[i].equals("-g")|| args[i].equals("--genome")) {
                if (i+1>= args.length) {
                    System.err.println("You forgot to specify the genome!");
                    System.exit(-1);
                }
                MyFile checkFile= new MyFile(args[i+1]);
                if (checkFile.exists()&& checkFile.isDirectory())
                    barna.model.Graph.overrideSequenceDirPath= checkFile.getAbsolutePath();
                else {
                    String[] s= args[++i].split("_");
                    if (s.length!= 2) {
                        System.err.println("Invalid genome directory: "+args[i+1]);
                        System.exit(-1);
                    }
                    Species spe= new Species(s[0]);
                    spe.setGenomeVersion(s[1]);
                    setSpecies(spe);
                }
                SplicingGraph.acceptableIntrons= true;
                continue;
            }

            if (args[i].equals("-k")|| args[i].equals("--dimension")) {
                try {
                    int x= Integer.parseInt(args[++i]);
                    if (x< -1|| ((x> -1)&& (x< 2)))
                        System.err.println(args[i]+" is not a valid dimension, ignored");
                    EventExtractor.n= x;
                } catch (NumberFormatException e) {
                    System.err.println(args[i]+" is not a valid dimension, ignored"); // :)
                }
                continue;
            }
            // -c is now RESERVED for "command" multicaster
//				if (args[i].equals("-c")|| args[i].equals("--canonical")) {
//					canonicalSS= true;
//					continue;
//				}
            //			if (args[i].equals("-c")|| args[i].equals("--cluster")) {
            //				readAheadLimit= Integer.parseInt(args[++i]);
            //				continue;
            //			}

            if (args[i].equalsIgnoreCase("-3pc")) {
                ASEvent.check3Pcomplete= true;
                continue;
            }
            if (args[i].equals("-tmp")) {
                System.setProperty(Constants.PROPERTY_TMPDIR, args[++i]);
                continue;
            }
            if (args[i].equals("-ds")) {
                SplicingGraph.retrieveDSEvents= false;
                continue;
            }
            if (args[i].equals("+ds")) {
                SplicingGraph.retrieveDSEvents= true;
                continue;
            }
            if (args[i].equals("-as")) {
                SplicingGraph.retrieveASEvents= false;
                continue;
            }
            if (args[i].equals("+as")) {
                SplicingGraph.retrieveASEvents= true;
                continue;
            }
            if (args[i].equalsIgnoreCase("-nmd")) {
                ASEvent.checkNMD= true;
                continue;
            }

            if (args[i].equals("+vs")) {
                SplicingGraph.retrieveVSEvents= true;
                continue;
            }
            if (args[i].equals("-vs")) {
                SplicingGraph.retrieveVSEvents= false;
                continue;
            }
            if (args[i].equals("-ext")) {
                SplicingGraph.onlyInternal= true;	// net false
                continue;
            }
            if (args[i].equals("+ext")) {
                SplicingGraph.onlyInternal= false;
                continue;
            }
            if (args[i].equalsIgnoreCase("-ic")|| args[i].equalsIgnoreCase("--intronConfidence")) {
                if (i+1== args.length)
                    System.err.println("You did not provide an intron confidence.");
                try {
                    SplicingGraph.intronConfidenceLevel= Byte.parseByte(args[i+1]);
                    ++i;	// ignore if missing
                    SplicingGraph.acceptableIntrons= true;
                } catch (NumberFormatException e) {
                    System.err.println("Intron confidence must be an integer value, you gave me "+args[i+1]);
                }
                continue;
            }

            // reactivated 20100112
            if (args[i].equalsIgnoreCase("-ec")|| args[i].equalsIgnoreCase("--edgeConfidence")) {
                if (i+1== args.length)
                    System.err.println("You did not provide an edge confidence.");
                try {
                    Transcript.setEdgeConfidenceLevel(Byte.parseByte(args[i+1]));
                    ++i;	// ignore if missing
                } catch (NumberFormatException e) {
                    System.err.println("Exon confidence must be an integer value, you gave me "+args[i+1]);
                }
                continue;
            }

            if (args[i].equals("--flankType")) {
                ASEvent.setOutputFlankMode(true);
                continue;
            }


        }

        if (ASEvent.isOutputFlankMode()&& barna.model.Graph.overrideSequenceDirPath== null) {
            System.err.println("[OOOPS] Parameters require genomic sequence. Use the option -g.");
            System.exit(-1);
        }

        // check for necessary components
        if (file== null|| !file.exists()) {
            System.err.print("Hey, you forgot to specify a valid input file! ");
            if (file!= null)
                System.err.println("Cannot find: "+file.getAbsolutePath());
            else
                System.err.println();
            System.exit(-1);
        }

        if ((SplicingGraph.canonicalSS|| SplicingGraph.acceptableIntrons)&& barna.model.Graph.overrideSequenceDirPath== null) {
            System.err.println("You want me to check introns for valid/canonical splice sites, but you did not provide a valid sequence directory");
            System.exit(-1);
        }

        if (outputFname== null&& (!SplicingGraph.writeStdOut)) {
            outputFname= file.getAbsolutePath()+"_astalavista.gtf.gz";
        }
        if (outputFname!= null&& new MyFile(outputFname).exists()) {
            // Confirm o..+"\n by typing \'yes\':"
            System.err.println("Overwriting output file "+outputFname+".");
//				try {
//					StringBuffer sb= new StringBuffer(3);
//					for (int i = System.in.read(); i != '\n';i= System.in.read()) {
//						sb.append((char) i);
//					}
//					sb.deleteCharAt(sb.length()-1);
//					if (sb.toString().trim().equalsIgnoreCase("yes")) {
//						System.err.println("AStalavista.");
//						System.exit(0);
//					}
//				} catch (Exception e) {
//					e.printStackTrace();
//					System.err.println("Output file exists, I give up.\nAStalavista.");
//					System.exit(0);
//				}
            while (new MyFile(outputFname).exists())	// TODO strange access probs under win32
                try {
                    new MyFile(outputFname).delete();
                    Thread.sleep(100);
                } catch (InterruptedException e) {
                    ; // e.printStackTrace();
                }
            System.err.println(outputFname+ " deleted.");
            //System.out.println("File "+outputFname+" exists, check - I give up now.");
            //System.exit(0);
        }


        GTFwrapper checkReader= new GTFwrapper(file.getAbsolutePath());
        //System.err.println("DEBUG -- temporarily deactivated file check");
        if (!checkReader.isApplicable()) {
            System.err.println("sorting input file, temporary directory "+System.getProperty(Constants.PROPERTY_TMPDIR));
            file= checkReader.sort();
            System.err.println("Here is a sorted version of your file: "+file.getAbsolutePath());
        }


        return file;
    }


}
