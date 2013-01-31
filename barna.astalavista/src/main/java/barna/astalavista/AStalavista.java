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
import barna.io.FileHelper;
import barna.io.GeneAheadReaderThread;
import barna.io.gtf.GTFwrapper;
import barna.model.*;
import barna.model.commons.MyFile;
import barna.model.splicegraph.*;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.*;
import java.lang.reflect.Field;
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
    static boolean acceptableIntrons= false;   // TODO not static

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


        // INIT
        long t0= System.currentTimeMillis();
        Log.message("# started\t" + new Date(t0));

        // print parameters
        ArrayList<barna.commons.parameters.Parameter> pars= new ArrayList<barna.commons.parameters.Parameter>(
                settings.getParameters().values());
        Collections.sort(pars, new barna.commons.parameters.Parameter.ParameterByNameComparator());
        for (barna.commons.parameters.Parameter p : pars) {
            Log.message("# "+ p.getName()+ "\t"+ settings.get(p));
        }

        // genome dir
        if (settings.get(AStalavistaSettings.CHR_SEQ)!= null) {
            Graph.overrideSequenceDirPath= settings.get(AStalavistaSettings.CHR_SEQ).getAbsolutePath();
        }

        // init variants
        if (settings.get(AStalavistaSettings.VARIANT_FILE)!= null) {
            variants= getVariants(settings.get(AStalavistaSettings.VARIANT_FILE));
        }

        // init input / reader thread
        GTFwrapper reader= new GTFwrapper(settings.get(AStalavistaSettings.IN_FILE).getAbsolutePath());
        if (!reader.isApplicable()) {
            File f= reader.sort();
            settings.set(AStalavistaSettings.IN_FILE, f);
            reader= new GTFwrapper(settings.get(AStalavistaSettings.IN_FILE).getAbsolutePath());
        }
        if (readAheadLimit> 0)
            reader.setReadAheadLimit(readAheadLimit);
        reader.setNoIDs(null);
        //reader.sweepToChromosome("chr17");
        GeneAheadReaderThread readerThread= new GeneAheadReaderThread(reader);
        readerThread.setOutput(false);
        readerThread.setOutput2(true);
        readerThread.start();



        // init site writer
        if (!settings.get(AStalavistaSettings.SITES).isEmpty()) {
            try {
                File f= null;
                if (settings.get(AStalavistaSettings.SITES_FILE)!= null)
                    f= settings.get(AStalavistaSettings.SITES_FILE);
                else
                    f= new File(FileHelper.append(settings.get(AStalavistaSettings.IN_FILE).getAbsolutePath(), "_sites", true, "vcf"));
                siteScoreWriter= new BufferedWriter(new FileWriter(f));
            } catch (IOException e) {
                Log.error(e.getMessage(), e);
            }
        }

        // init event writer thread
        EventExtractor.writerThread= new WriterThread();
        if (settings.get(AStalavistaSettings.EVENTS_FILE)!= null)
            EventExtractor.writerThread.outputFname= settings.get(AStalavistaSettings.EVENTS_FILE).getAbsolutePath();
        else
            EventExtractor.writerThread.outputFname= settings.get(AStalavistaSettings.IN_FILE).getAbsolutePath()+"_astalavista.gtf.gz";

        if (EventExtractor.writerThread.outputFname!= null&& new MyFile(EventExtractor.writerThread.outputFname).exists()) {
            // Confirm o..+"\n by typing \'yes\':"
            Log.warn("Overwriting output file " + EventExtractor.writerThread.outputFname + ".");
        }
        EventExtractor.writerThread.start();


        // MAIN LOOP
        Log.progressStart("Iterating Annotation");
        long inBytes= settings.get(AStalavistaSettings.IN_FILE).length();
        Log.progress(0, inBytes);
        while (true) {
            Gene[] g= readerThread.getGenes();
            Log.progress(readerThread.getBytesRead(),inBytes);
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
                g[i].setSpecies(species);

                // sets types of events to be extracted, k, etc..
                EventExtractor extractor= new EventExtractor(g[i], settings);
                Thread extractorThread= new Thread(extractor);
                extractorThread.start();
                extractorThread.join();

            }

            // prepare for next batch
            Date ti = new Date(System.currentTimeMillis());
            int div = (int) (EventExtractor.cumulGC +
                    EventExtractor.cumulGF +
                    EventExtractor.cumulEV) / 1000;
            int frac = 0;
            if (div > 0)
                frac = (EventExtractor.counter - evBefore) / div;
            else
                frac = (EventExtractor.counter - evBefore);

            Log.debug("[" + ti + "] " + chromo +
                    " graph construct " + (EventExtractor.cumulGC / 1000) +
                    //" sec, fuzzy flanks "+(cumulGF/1000) +
                    " sec, contraction " + (EventExtractor.cumulGT / 1000) +
                    " sec, extract events " + (EventExtractor.cumulEV / 1000) +
                    " sec, found " + (EventExtractor.counter - evBefore) +
                    " events, " + frac + " ev/sec.");
            System.gc();
            EventExtractor.writerThread.interrupt();

        }
       Log.progressFinish("done.", true);
       Log.info("took "+((System.currentTimeMillis()- t0)/1000)+" sec.");

        if (siteScoreWriter!= null)
            siteScoreWriter.close();

        try {
            EventExtractor.writerThread.setKill(true);
            EventExtractor.writerThread.interrupt();
            EventExtractor.writerThread.join();
        } catch (InterruptedException e) {
            ; // :)
        }

        Log.info("found "+ EventExtractor.counter+" events.");

        if (settings.get(AStalavistaSettings.EVENTS_OPT).contains(AStalavistaSettings.EventOptions.IOK)) {
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
            if (c.toString().toLowerCase().contains("enum"))
                System.currentTimeMillis();
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
        if (args.userSpecified(AStalavistaSettings.HELP.getName())) {
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

        if (settings.get(AStalavistaSettings.IN_FILE)== null) {
            Log.error("Hey, you forgot to specify a valid input file! ");
            return false;
        }

        if (settings.requiresGenomicSequence()
            && settings.get(AStalavistaSettings.CHR_SEQ)== null) {

            Log.error("You want me to check introns for valid/canonical splice sites, " +
                    "but you did not provide a valid sequence directory");
        }

        if (settings.get(AStalavistaSettings.CHR_SEQ)!= null) {

            File f= settings.get(AStalavistaSettings.CHR_SEQ);
            if (f.exists()&& f.isDirectory()) {
                Graph.overrideSequenceDirPath= f.getAbsolutePath();
            } else {
                Log.message("Trying to parse species_version pair");
                String[] s= f.getName().split("_");
                if (s.length!= 2) {
                    Log.error("The genome " + f.getAbsolutePath() + " could not be found!");
                    return false;
                }
                Species spe= new Species(s[0]);
                spe.setGenomeVersion(s[1]);
                AStalavista.setSpecies(spe);

            }
        }

        // check splice site scoring stuff
        if (settings.get(AStalavistaSettings.SITES_OPT).contains(AStalavistaSettings.SiteOptions.SSS)) {

            if(settings.get(AStalavistaSettings.CHR_SEQ)== null) {
                Log.error("Splice site scoring requires the genomic sequence, set parameter \'CHR_SEQ\'");
                return false;
            }

            try {
                if(settings.get(AStalavistaSettings.GENE_ID)== null) {
                    Log.warn("No GeneID parameter file for scoring models, using default");
                    geneidParam= Profile.readParam(GeneIDconstants.PARAMETERFILE, new GeneIDsettings())[0];
                } else
                    geneidParam= Profile.readParam(settings.get(AStalavistaSettings.GENE_ID).getAbsolutePath(),
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
        if (settings== null) {
            settings = new AStalavistaSettings();
        }

        // copy
        Collection<barna.commons.parameters.Parameter> pars=
                settings.getParameters().values();
        for (barna.commons.parameters.Parameter p : pars) {
            if (args.userSpecified(p.getLongOption())) {
                p.parse(args.getObject(p.getLongOption()).toString());
            }
        }

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

    /**
     * Summarizes the settings before the run
     */
    protected void printSettings() {

        Log.message("# started\t" + new Date(System.currentTimeMillis()));
        Log.message("# INPUT #");
        Log.message("# annotation\t" + settings.get(AStalavistaSettings.IN_FILE).getAbsolutePath());
        Log.message("# chromosomes\t"+ barna.model.Graph.overrideSequenceDirPath);
        if (settings.get(AStalavistaSettings.GENE_ID)!= null)
            Log.message("# chromosomes\t"+ settings.get(AStalavistaSettings.GENE_ID));
        if (settings.get(AStalavistaSettings.VARIANT_FILE)!= null)
            Log.message("# variants\t"+ settings.get(AStalavistaSettings.VARIANT_FILE));

        if (!settings.get(AStalavistaSettings.EVENTS).isEmpty()) {
            Log.message("# EVENTS #");
            Log.message("# output\t"+ outputFname);
            Log.message("# dimension\t" + settings.get(AStalavistaSettings.EVENTS_DIMENSION));
            Log.message("# as_events\t"+ (SplicingGraph.retrieveASEvents? "true": "false"));
            Log.message("# ds_events\t"+ (SplicingGraph.retrieveDSEvents? "true": "false"));
            Log.message("# vs_events\t"+ (SplicingGraph.retrieveVSEvents? "true": "false"));
            Log.message("# internalOnly\t" + SplicingGraph.onlyInternal);
            if (!SplicingGraph.onlyInternal)
                Log.message("# edgeConfidenceLevel " + Transcript.getEdgeConfidenceLevel());
            Log.message("# canonicalSS\t"+ SplicingGraph.canonicalSS);
            Log.message("# acceptableIntrons\t"+ settings.get(AStalavistaSettings.EVENTS_OPT).contains(AStalavistaSettings.EventOptions.IOK));
            if (SplicingGraph.acceptableIntrons)
                Log.message("# intronConfidenceLevel " + SplicingGraph.intronConfidenceLevel);
        }

        if (!settings.get(AStalavistaSettings.SITES).isEmpty()) {
            Log.message("# SITES #");
            Log.message("# scores\t"+ settings.get(AStalavistaSettings.SITES_OPT).contains(AStalavistaSettings.SiteOptions.SSS));
        }

        System.exit(0);
    }

}
