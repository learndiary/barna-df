package barna.astalavista;

import barna.commons.Execute;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.log.Log;
import barna.geneid.GParam;
import barna.geneid.GeneID;
import barna.geneid.GeneIDconstants;
import barna.geneid.GeneIDsettings;
import barna.geneid.Profile;
import barna.io.FileHelper;
import barna.model.Gene;
import barna.model.Graph;
import barna.model.SpliceSite;
import barna.model.splicegraph.*;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 2/5/13
 * Time: 2:48 PM
 */
public class Scorer extends AStalavista {

    public static final int TYPE_UNKNOWN= -1;
    public static final int TYPE_U2_GTAG= 0;
    public static final int TYPE_U2_GCAG= 1;
    public static final int TYPE_U12_GTAG= 2;
    public static final int TYPE_U12_ATAC= 3;
    public static final String[] TYPE_DESCRIPTION= {"U2_GTAG", "U2_GCAG", "U12_GTAG", "U12_ATAC"};

    /**
     * Maximum length of string flanking a splice site required
     * to compute all given donor / acceptor site models.
     * Usually a 2x2 matrix:<br>
     * [0][0] donor upstream 5'-flank<br>
     * [0][1] donor downstream 3'-flank<br>
     * [1][0] acceptor upstream 5'-flank<br>
     * [1][1] acceptor downstream 3'-flank<br>
     */
    protected int[][] maxFlanks= null;

    /**
     * Map intron types to their corresponding donor model
     */
    protected Profile[] donorProfiles= new Profile[4];

    /**
     * Map intron types to their corresponding acceptor model
     */
    protected Profile[] acceptorProfiles= new Profile[4];

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

    public static void main(String[] args) {

        Execute.initialize(2);
        Scorer myScorer= new Scorer();

        // construct to register parameters in JSAP
        List<Parameter> parameter = JSAPParameters.getJSAPParameter(new ScorerSettings());
        JSAP jsap = JSAPParameters.registerParameters(parameter);

        // parse
        try{
            JSAPResult toolParameter = jsap.parse(args);
            if (!myScorer.validateParameter(toolParameter)){
                System.exit(-1);
            }
        } catch (Exception e) {
            Log.error("Parameter error : " + e.getMessage(), e);
            e.printStackTrace();
            System.exit(-1);
        }

        Future<Void> captain= Execute.getExecutor().submit(myScorer);
        try {
            captain.get();
        } catch (InterruptedException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (ExecutionException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        Execute.shutdown();
    }


    @Override
    public String getName() {
        return "scorer";
    }

    @Override
    public String getDescription() {
        return "Splice site scorer";
    }

    @Override
    public String getLongDescription() {
        return "The splice site scorer obtains geneID scores for all splice sites in the annotation, " +
                "possibly including individual variants.";
    }


    @Override
    public Void call() throws Exception {

        // init variants
        if (settings.get(ScorerSettings.VARIANT_FILE)!= null) {
            variants= getVariants(settings.get(ScorerSettings.VARIANT_FILE));
        }

        return super.call();
    }

    @Override
    protected void callLoop(Gene g) throws Exception {
        // score splice sites
        g.markAlternativeSpliceSites();

        // deprecated, difficult to distinguish U12 acceptors
        SplicingGraph gr= new SplicingGraph(g);
        gr.constructGraph();
        scoreSites(gr, g.getSpliceSites());

        // difficult to avoid site redundancy
        //scoreSites(gr.getNodesInGenomicOrder());
    }

    @Override
    protected void callBegin() throws Exception {
        super.callBegin();
        siteScoreWriter= getSiteScoreWriter();
    }

    @Override
    protected void callFinish() throws Exception {
        super.callFinish();
        if (siteScoreWriter!= null)
            siteScoreWriter.close();
    }

    @Override
    public boolean validateParameter(JSAPResult args) {

        if (!super.validateParameter(new ScorerSettings(), args))
            return false;

        // check splice site scoring stuff
        if (settings.get(ScorerSettings.SITES_OPT).contains(ScorerSettings.SiteOptions.SSS)) {

            if(settings.get(AStalavistaSettings.CHR_SEQ)== null) {
                Log.error("Splice site scoring requires the genomic sequence, provide a value for parameter \'"+
                    AStalavistaSettings.CHR_SEQ.getName()+ "\' in the parameter file, or via " +
                        "the command line flags -"+ (AStalavistaSettings.CHR_SEQ.getShortOption())+
                        " or --"+ (AStalavistaSettings.CHR_SEQ.getLongOption())+ "!");
                return false;
            }

            try {
                if(settings.get(ScorerSettings.GENE_ID)== null) {
                    Log.warn("No GeneID parameter file for scoring models, using default (for human)");
                    geneidParam= Profile.readParam(GeneIDconstants.PARAMETERFILE, new GeneIDsettings())[0];
                } else
                    geneidParam= Profile.readParam(settings.get(ScorerSettings.GENE_ID).getAbsolutePath(),
                            new GeneIDsettings())[0];

                // determine max flanks to read genomic sequence only once later-on
                maxFlanks= new int[2][];
                maxFlanks[0]= new int[2];
                maxFlanks[1]= new int[2];

                assert(geneidParam.getDonorProfile()!= null&& geneidParam.getAcceptorProfile()!= null);

                // U2 GTAG
                donorProfiles[TYPE_U2_GTAG]= geneidParam.getDonorProfile();
                acceptorProfiles[TYPE_U2_GTAG]= geneidParam.getAcceptorProfile();
                maxFlanks[0][0]= GParam.getDonorFlank5(geneidParam.getDonorProfile());
                maxFlanks[0][1]= GParam.getDonorFlank3(geneidParam.getDonorProfile());
                maxFlanks[1][0]= GParam.getAcceptorFlank5(geneidParam.getAcceptorProfile());
                maxFlanks[1][1]= GParam.getAcceptorFlank3(geneidParam.getAcceptorProfile());

                // U2 GCAG
                if (geneidParam.getU2gcagDonorProfile()!= null) {
                    donorProfiles[TYPE_U2_GCAG]= geneidParam.getU2gcagDonorProfile();
                    acceptorProfiles[TYPE_U2_GCAG]= geneidParam.getAcceptorProfile();
                    maxFlanks[0][0]= Math.max(maxFlanks[0][0],
                            GParam.getDonorFlank5(geneidParam.getU2gcagDonorProfile()));
                    maxFlanks[0][1]= Math.max(maxFlanks[0][1],
                            GParam.getDonorFlank3(geneidParam.getU2gcagDonorProfile()));
                }

                // U12 GTAG
                if (geneidParam.getU12gtagDonorProfile()!= null && geneidParam.getU12gtagAcceptorProfile()!= null) {
                    donorProfiles[TYPE_U12_GTAG]= geneidParam.getU12gtagDonorProfile();
                    acceptorProfiles[TYPE_U12_GTAG]= geneidParam.getU12gtagAcceptorProfile();
                    maxFlanks[0][0]= Math.max(maxFlanks[0][0],
                            GParam.getDonorFlank5(geneidParam.getU12gtagDonorProfile()));
                    maxFlanks[0][1]= Math.max(maxFlanks[0][1],
                            GParam.getDonorFlank3(geneidParam.getU12gtagDonorProfile()));
                    maxFlanks[1][0]= Math.max(maxFlanks[1][0],
                            GParam.getAcceptorFlank5(geneidParam.getU12gtagAcceptorProfile()));
                    maxFlanks[1][1]= Math.max(maxFlanks[1][1],
                            GParam.getAcceptorFlank3(geneidParam.getU12gtagAcceptorProfile()));
                }

                // U12 ATAC
                if (geneidParam.getU12atacDonorProfile()!= null && geneidParam.getU12atacDonorProfile()!= null) {
                    donorProfiles[TYPE_U12_ATAC]= geneidParam.getU12atacDonorProfile();
                    acceptorProfiles[TYPE_U12_ATAC]= geneidParam.getU12atacAcceptorProfile();
                    maxFlanks[0][0]= Math.max(maxFlanks[0][0],
                            GParam.getDonorFlank5(geneidParam.getU12atacDonorProfile()));
                    maxFlanks[0][1]= Math.max(maxFlanks[0][1],
                            GParam.getDonorFlank3(geneidParam.getU12atacDonorProfile()));
                    maxFlanks[1][0]= Math.max(maxFlanks[1][0],
                            GParam.getAcceptorFlank5(geneidParam.getU12atacAcceptorProfile()));
                    maxFlanks[1][1]= Math.max(maxFlanks[1][1],
                            GParam.getAcceptorFlank3(geneidParam.getU12atacAcceptorProfile()));
                }


            } catch (Exception e) {
                Log.error(e.getMessage(), e);
                return false;
            }
        }

        return true;
    }


    /**
     * Summarizes the settings before the run
     */
    @Deprecated
    protected void printSettings(AStalavistaSettings settings) {
        super.printSettings(settings);
        if (settings.get(ScorerSettings.GENE_ID)!= null)
            Log.message("# chromosomes\t"+ settings.get(ScorerSettings.GENE_ID));
        if (settings.get(ScorerSettings.VARIANT_FILE)!= null)
            Log.message("# variants\t"+ settings.get(ScorerSettings.VARIANT_FILE));
        if (!settings.get(ScorerSettings.SITES).isEmpty()) {
            Log.message("# SITES #");
            Log.message("# scores\t"+ settings.get(ScorerSettings.SITES_OPT).contains(ScorerSettings.SiteOptions.SSS));
        }
    }

    private BufferedWriter getSiteScoreWriter() {
        if (siteScoreWriter== null) {
            // init site writer
            try {
                File f= null;
                if (settings.get(ScorerSettings.SITES_FILE)!= null)
                    f= settings.get(ScorerSettings.SITES_FILE);
                else
                    f= new File(FileHelper.append(settings.get(AStalavistaSettings.IN_FILE).getAbsolutePath(), "_sites", true, "vcf"));
                siteScoreWriter= new BufferedWriter(new FileWriter(f));
            } catch (IOException e) {
                Log.error(e.getMessage(), e);
            }
        }
        return siteScoreWriter;
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


    /**
     *
     * @param ss
     * @param variants
     * @param sequences
     * @param scores
     * @param varTuples
     */
    private void outputSitesVCF(SpliceSite ss, Vector<String> variants, String[] sequences, float[] scores, int[] types, String[] varTuples) {

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
        sb.append("TYPE=");
        sb.append(types[0]< 0? "UNKNOWN": TYPE_DESCRIPTION[types[0]]);
        sb.append(";");
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

        if (types.length> 1)
            sb.append("VAR_TYPES=");
        for (int i = 1; i < scores.length; ++i) {
            sb.append(types[0]< 0? "UNKNOWN": TYPE_DESCRIPTION[types[0]]);
            sb.append(",");
        }
        if (types.length> 1)
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

    @Override
    protected AStalavistaSettings getSettings() {
        if (settings== null)
            return new ScorerSettings();
        return settings;
    }

    protected int scoreVariants(SplicingGraph graph, int[] type, Vector<String> vvec, String varID, int rec, int nr, int idx, SpliceSite ss, int flank5, int flank3, String seq,
                                String[] seqs, float[] scores, int[] types, String[] varTuples) {

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
            float score= scoreSite(graph, ss, varSeq, type);
            String vvarID= varID+ (varID.length()> 0? ",": "")+ vv[2];
            //outputSite(ss, vvarID, varSeq, score);
            seqs[nr]= varSeq;
            scores[nr]= score;
            types[nr]= type[0];
            varTuples[nr]= vvarID;

            // start recursion
            nr= scoreVariants(graph, type,  vvec, vvarID, rec+ 1, nr, i + 1, ss, flank5, flank3, varSeq, seqs, scores, types, varTuples);
        }

        return nr;

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
    protected float scoreSite(SplicingGraph graph, SpliceSite spliceSite, String seq, int[] type) {

        seq= seq.toUpperCase();

        Node n= graph.getNode(spliceSite);
        if (spliceSite.isDonor()) {

            String nt= seq.substring(maxFlanks[0][0], maxFlanks[0][0]+ 2);

            float sscore= Float.NEGATIVE_INFINITY;
            int stype= TYPE_UNKNOWN;
            for (int i = 0; i < n.getOutEdges().size(); i++) {

                if (!n.getOutEdges().elementAt(i).isIntronic())
                    continue;

                SpliceSite cosite= n.getOutEdges().elementAt(i).getHead().getSite();
                String seq2= nt+ Graph.readSequence(cosite, 0, 0).toUpperCase();

                if (seq2.equals("GTAG")) {
                    // decide based on donor only
                    float scoreU2= scoreDonor(seq, geneidParam.getDonorProfile());
                    float scoreU12= scoreDonor(seq, geneidParam.getU12gtagDonorProfile());
                    if (scoreU2< -100|| scoreU12> 5) {
                        // no, too many U12
                        //getNormScore(scoreU2, geneidParam.getDonorProfile())< getNormScore(scoreU12, geneidParam.getU12gtagDonorProfile())) {
                        type[0]= TYPE_U12_GTAG;
                        return scoreU12;
                    } else {
                        type[0]= TYPE_U2_GTAG;
                        return scoreU2;
                    }
                } else if (seq2.equals("GCAG")) {
                    float score= scoreDonor(seq, geneidParam.getU2gcagDonorProfile());
                    type[0]= TYPE_U2_GCAG;
                    return score;
                } else if (seq2.equals("ATAC")) {
                    float score= scoreDonor(seq, geneidParam.getU12atacDonorProfile());
                    type[0]= TYPE_U12_ATAC;
                    return score;
                }
            }

        } else if (spliceSite.isAcceptor()) {

            String nt= seq.substring(maxFlanks[0][0], maxFlanks[0][0]+ 2);

            for (int i = 0; i < n.getInEdges().size(); i++) {

                if (!n.getInEdges().elementAt(i).isIntronic())
                    continue;

                SpliceSite cosite= n.getInEdges().elementAt(i).getTail().getSite();
                String coseq= Graph.readSequence(cosite, maxFlanks[0][0], maxFlanks[0][1]).toUpperCase();
                String seq2= coseq.substring(maxFlanks[0][0], maxFlanks[0][0]+ 2)+ nt;

                if (seq2.equals("GTAG")) {
                    // decide based on donor only
                    float scoreU2= scoreDonor(coseq, geneidParam.getDonorProfile());
                    float scoreU12= scoreDonor(coseq, geneidParam.getU12gtagDonorProfile());
                    if (scoreU2< -100|| scoreU12> 5) {
                        //getNormScore(scoreU2, geneidParam.getDonorProfile())< getNormScore(scoreU12, geneidParam.getU12gtagDonorProfile())) {
                        type[0]= TYPE_U12_GTAG;
                        float score= scoreAcceptor(seq, geneidParam.getU12gtagAcceptorProfile(), null, null);
                        return score;
                    } else {
                        type[0]= TYPE_U2_GTAG;
                        float score= scoreAcceptor(seq, geneidParam.getAcceptorProfile(), null, null);
                        return scoreU2;
                    }
                } else if (seq2.equals("GCAG")) {
                    float score= scoreAcceptor(seq, geneidParam.getAcceptorProfile(), null, null);
                    type[0]= TYPE_U2_GCAG;
                    return score;
                } else if (seq2.equals("ATAC")) {
                    float score= scoreAcceptor(seq, geneidParam.getU12atacAcceptorProfile(), null, null);
                    type[0]= TYPE_U12_ATAC;
                    return score;
                }

            }
        }

        type[0]= TYPE_UNKNOWN;
        return Float.NEGATIVE_INFINITY;
    }


    protected float scoreSiteDELME(SplicingGraph graph, SpliceSite spliceSite, String seq, int[] type) {

        seq= seq.toUpperCase();

        if (spliceSite.isDonor()) {

            String nt= seq.substring(maxFlanks[0][0], maxFlanks[0][0]+ 2);

            if (nt.equals("GC")) {
                type[0]= TYPE_U2_GCAG;
                float score= scoreDonor(seq, geneidParam.getU2gcagDonorProfile());
                return score;
                // difficult to score U12 wo whole intron
//            } else if (nt.equals("AT")) {
//                type[0]= TYPE_U12_ATAC;
//                float score= scoreDonor(seq, geneidParam.getU12atacDonorProfile());
//                return score;
            } else if (nt.equals("GT")|| nt.equals("AT")) {
                // score first U12, quite a scrict sequence ..GTATCCTtt
                float scoreU2= scoreDonor(seq, geneidParam.getDonorProfile());
                float scoreU12= Float.NEGATIVE_INFINITY;
                float donorScore= Float.NEGATIVE_INFINITY;
                if (nt.equals("GT"))
                    donorScore= scoreDonor(seq, geneidParam.getU12gtagDonorProfile());
                else
                    donorScore= scoreDonor(seq, geneidParam.getU12atacDonorProfile());
                Node n= graph.getNode(spliceSite);
                for (int i = 0; i < n.getOutEdges().size(); i++) {

                    if (!n.getOutEdges().elementAt(i).isIntronic())
                        continue;

                    SpliceSite acceptor= n.getOutEdges().elementAt(i).getHead().getSite();
                    String acceptorSeq= Graph.readSequence(acceptor, maxFlanks[1][0], maxFlanks[1][1]).toUpperCase();
                    float acceptorScore= Float.NEGATIVE_INFINITY;
                    if (nt.equals("GT"))
                        acceptorScore= scoreAcceptor(acceptorSeq, geneidParam.getU12gtagAcceptorProfile(), null, null);
                    else
                        acceptorScore= scoreAcceptor(acceptorSeq, geneidParam.getU12atacAcceptorProfile(), null, null);

                    if (donorScore+ acceptorScore> scoreU12)
                        scoreU12= donorScore+ acceptorScore;
                }

                // decision
                if (scoreU12> geneidParam.getU12SpliceScoreThresh()|| (scoreU2< -1000&& scoreU12> -10)) {
                    if (nt.equals("GT"))
                        type[0]= TYPE_U12_GTAG;
                    else
                        type[0]= TYPE_U12_ATAC;
                    return donorScore;
                } else {
                    type[0]= TYPE_U2_GTAG;
                    return scoreU2;
                }

            }
        }

        // else
        if (spliceSite.isAcceptor()) {
            String nt= seq.substring(maxFlanks[1][0], maxFlanks[1][0]+ 2);
            // score first u2 with ppt
            float scoreU2= scoreAcceptor(seq, geneidParam.getAcceptorProfile(), null, null); // TODO use u2 bp, ppt
            float scoreU12= Float.NEGATIVE_INFINITY;
            float acceptorScore= Float.NEGATIVE_INFINITY;
            if (nt.equals("AC")|| (nt.equals("AG"))) {

                if (nt.equals("AC"))
                    acceptorScore= scoreAcceptor(seq, geneidParam.getU12atacAcceptorProfile(), null, null);
                else
                    acceptorScore= scoreAcceptor(seq, geneidParam.getU12gtagAcceptorProfile(), null, null);

                Node n= graph.getNode(spliceSite);
                for (int i = 0; i < n.getInEdges().size(); i++) {

                    if (!n.getInEdges().elementAt(i).isIntronic())
                        continue;

                    SpliceSite donor= n.getInEdges().elementAt(i).getTail().getSite();
                    String donorSeq= Graph.readSequence(donor, maxFlanks[0][0], maxFlanks[0][1]).toUpperCase();
                    String donorNt= donorSeq.substring(maxFlanks[0][0], maxFlanks[0][0]+ 2);
                    float donorScore= Float.NEGATIVE_INFINITY;
                    if (nt.equals("AC")) {
                        donorScore= scoreDonor(donorSeq, geneidParam.getU12atacDonorProfile());
                    } else
                        donorScore= scoreDonor(donorSeq, geneidParam.getU12gtagDonorProfile());

                    if (donorScore+ acceptorScore> scoreU12)
                        scoreU12= donorScore+ acceptorScore;
                }
            }
            if (scoreU12> geneidParam.getU12SpliceScoreThresh()|| (scoreU2< -1000&& scoreU12> -10)) {
                if (nt.equals("AC")) {
                    type[0]= TYPE_U12_ATAC;
                    return acceptorScore;
                } else {
                    type[0]= TYPE_U12_GTAG;
                    return acceptorScore;
                }
            } else {
                type[0]= TYPE_U2_GTAG;
                return scoreU2;
            }


        }

        // else: corrupt splice site
        type[0]= TYPE_UNKNOWN;
        return Float.NaN; // throw exception?
    }

    protected float scoreSiteDELME(SpliceSite spliceSite, String seq) {

        seq= seq.toUpperCase();

        if (spliceSite.isDonor())
            return GeneID.scoreDonor(seq, geneidParam.getDonorProfile());
        if (spliceSite.isAcceptor())
            return GeneID.scoreAcceptor(seq, geneidParam.getAcceptorProfile(), null, null);

        return Float.NaN; // throw exception?
    }

    /*
     * @deprecated not used
     */
    protected void scoreIntron(SimpleEdge intron) {

        // get donor/acceptor up- and downstream
        assert(maxFlanks.length== 2&& maxFlanks[0].length== 2&& maxFlanks[1].length== 2);
        String donorSeq= Graph.readSequence(intron.getTail().getSite(), maxFlanks[0][0], maxFlanks[0][1]).toUpperCase();
        String acceptorSeq= Graph.readSequence(intron.getHead().getSite(), maxFlanks[1][0], maxFlanks[1][1]).toUpperCase();

        // pre-decide type by dinucleotide
        String donNt= donorSeq.substring(maxFlanks[0][0], maxFlanks[0][0]+ 2).toUpperCase();
        String accNt= acceptorSeq.substring(maxFlanks[1][0],maxFlanks[1][0]+ 2).toUpperCase();
        String diNt= donNt+ accNt;


        // calculate all scores, for letting variants take their own decision
        float[] donorScores= new float[4], acceptorScores= new float[4];
        float dMax= Float.NEGATIVE_INFINITY, aMax= Float.NEGATIVE_INFINITY;
        int dMaxPos= -1, aMaxPos= -1;
        for (int i = 0; i < donorScores.length; i++) {
            if (diNt.equals("GCAG")&& i== 1)
                System.currentTimeMillis();
            donorScores[i]= scoreDonor(donorSeq, donorProfiles[i]);
            if (donorScores[i]!= Float.NaN&& donorScores[i]> dMax) {
                dMax= donorScores[i];
                dMaxPos= i;
            }
            acceptorScores[i]= scoreAcceptor(acceptorSeq, acceptorProfiles[i], null, null); // TODO no PPT and BP used
            if (acceptorScores[i]!= Float.NaN&& acceptorScores[i]> aMax) {
                aMax= acceptorScores[i];
                aMaxPos= i;
            }
        }



        if (dMaxPos> 0|| aMaxPos> 0)
            System.currentTimeMillis();

        if (diNt.equals("GTAG"))
            System.currentTimeMillis();


        //scoreSiteVariants(intron.getTail().getSite(), acceptorScores, maxFlanks[0][0], maxFlanks[0][1], donorSeq);
        //scoreSiteVariants(intron.getHead().getSite(), donorScores, maxFlanks[1][0], maxFlanks[1][1], acceptorSeq);

    }


    /**
     * Score an entire Intron
     * @param donorSeq sequence at donor site
     * @param acceptorSeq seuqence at acceptor site
     * @param profileDon donor site profile
     * @param profileAcc acceptor site profile
     * @param profilePPT poly-pyrimidine tract profile, can be <code>null</code>
     * @param profileBP branch-point profile, can be <code>null</code>
     * @return score of the intron as a sum of both scores of its splice sites
     */
    protected float scoreIntron(String donorSeq, String acceptorSeq, Profile profileDon, Profile profileAcc, Profile profilePPT, Profile profileBP) {

        if(profileDon== null|| profileAcc== null)
            return Float.NaN;

        float scoreDon= scoreDonor(donorSeq, profileDon);
        float scoreAcc= scoreAcceptor(acceptorSeq, profileAcc, profilePPT, profileBP);

        return (scoreDon+ scoreAcc);
    }


    /**
     * Adapt sequence for profile and delegate to HMM score computation
     * @param baseSeq (super-)sequence of the splice site area
     * @param profile donor profile for computation
     * @return score of the donor site according to model
     */
    protected float scoreDonor(String baseSeq, Profile profile) {
        if (profile== null)
            return Float.NaN;
        String seq= baseSeq;
        if (GParam.getDonorFlank5(profile)< maxFlanks[0][0])
            seq= baseSeq.substring(maxFlanks[0][0]- GParam.getDonorFlank5(profile));
        if (GParam.getDonorFlank3(profile)< maxFlanks[1][1])
            seq= baseSeq.substring(0, seq.length()- (maxFlanks[1][1]- GParam.getDonorFlank3(profile)));
        return GeneID.scoreDonor(seq, profile);
    }

    /**
     * Computes a relative score for the donor
     * @param score the real score to be normalized
     * @param profile donor profile for computation
     * @return score of the site normalized according to model
     */
    protected float getNormScore(float score, Profile profile) {

        float min= getLimitScore(profile, true);
        float max= getLimitScore(profile, false);

        float nval= ((score- min)/ (max- min));
        return nval;
    }

    protected HashMap<Profile, Float> mapMin= null, mapMax= null;

    /**
     * Returns the minimum/maximum sum of a profile, ignoring
     * artificially low numbers that represent (-Infinity).
     * @param profile the profile to be assessed
     * @param min flag to trigger the search for the minimum
     *            respectively maximum
     * @return limiting value of the sum
     */
    protected float getLimitScore(Profile profile, boolean min) {
        if (min) {
            if (mapMin== null)
                mapMin= new HashMap<Profile, Float>();
            if (mapMin.containsKey(profile))
                return mapMin.get(profile);
            else {
                float val= calcLimitScore(profile, min);
                mapMin.put(profile, val);
                return val;
            }

        } else {    /* max */
            if (mapMax== null)
                mapMax= new HashMap<Profile, Float>();
            if (mapMax.containsKey(profile))
                return mapMax.get(profile);
            else {
                float val= calcLimitScore(profile, min);
                mapMax.put(profile, val);
                return val;
            }

        }
    }

    /**
     * Computes minimum or maximum score of a markov chain,
     * disregarding constants that express -Inf (-9999).
     * @param profile the Markov chain
     * @param min flag to indicate whether the min or max is to be computed
     * @return the corresponding minimum/maximum sum of the Markov chain
     */
    protected float calcLimitScore(Profile profile, boolean min) {

        float[][] t= profile.getTransitionValues();
        float sum= 0f;
        for (int i = 0; i < profile.getDimension(); i++) {

            float x= (min? Float.MAX_VALUE: Float.MIN_VALUE);
            for (int j = 0; j < t[i].length; j++) {
                if (t[i][j]<= -1000)    // empiric value from parameter file, no constant available
                    continue;
                if ((min&& t[i][j]< x)|| ((!min)&& t[i][j]> x))
                    x= t[i][j];
            }
            sum+= x;

        }

        return sum;
    }


    /**
     * Adapt sequence for profile and delegate to HMM score computation
     * @param baseSeq (super-)sequence of the splice site area
     * @param profileAcc acceptor profile for computation
     * @param profilePPT poly-pyrimidine tract profile, can be <code>null</code>
     * @param profilePPT branch-point profile, can be <code>null</code>
     * @return score of the acceptor site according to model
     */
    protected float scoreAcceptor(String baseSeq, Profile profileAcc, Profile profilePPT, Profile profileBP) {
        if (profileAcc== null)
            return Float.NaN;
        String seq= baseSeq;
        if (GParam.getAcceptorFlank5(profileAcc)< maxFlanks[1][0])
            seq= baseSeq.substring(maxFlanks[1][0]- GParam.getAcceptorFlank5(profileAcc));
        if (GParam.getAcceptorFlank3(profileAcc)< maxFlanks[1][1])
            seq= baseSeq.substring(0, seq.length()- (maxFlanks[1][1]- GParam.getAcceptorFlank3(profileAcc)));

        return GeneID.scoreAcceptor(seq, profileAcc, profilePPT, profileBP);
    }


    protected void scoreSites(Node[] sites) {

        for (int i = 0; i < sites.length; i++) {
            for (int j= 0; sites[i].getSite().isDonor()&& j< sites[i].getOutEdges().size(); ++j) {

                SimpleEdge intron= sites[i].getOutEdges().elementAt(j);
                if (!intron.getHead().getSite().isAcceptor())
                    continue;

                // score every intron a splice site participates in independently
                // who says that a site cannot interact with different spliceosomes
                // in the context of different introns?
                scoreIntron(intron);
            }

        }
    }

    /**
     * Retrieves the score of splice sites as obtained from the genomic sequence,
     * and also of variants of those, if annotated.
     * @param spliceSites
     */
    protected void scoreSites(SplicingGraph graph, SpliceSite[] spliceSites) {

        // get donor/acceptor up- and downstream
        assert(maxFlanks.length== 2&& maxFlanks[0].length== 2&& maxFlanks[1].length== 2);

        String seq= null;
        float score;
        int flank5, flank3;
        for (int i = 0; i < spliceSites.length; i++) {

            if (spliceSites[i].isDonor()) {
                flank5= maxFlanks[0][0];
                flank3= maxFlanks[0][1];

            } else if (spliceSites[i].isAcceptor()) {
                flank5= maxFlanks[1][0];
                flank3= maxFlanks[1][1];

            } else
                continue;

            seq= Graph.readSequence(spliceSites[i], flank5, flank3);
            if (seq== null)
                continue;

            // get variants, if annotated
            Vector<String> vvec= getVariants(spliceSites[i], flank5, flank3);
            int nrCombinations= (int) Math.pow(2, vvec== null? 0: vvec.size());

            // arrays to store reference results and
            String[] seqs= new String[nrCombinations];
            float[] scores= new float[nrCombinations];
            int[] types= new int[nrCombinations];
            String[] varTuples= new String[nrCombinations];

            // type site
            int[] type= new int[1];
            score= scoreSite(graph, spliceSites[i], seq, type);
            // outputSite(spliceSites[i], null, seq, score);
            seqs[0]= seq;
            scores[0]= score;
            types[0]= type[0];
            varTuples[0]= "reference";

            // recursion to do all tuple combinations
            if (vvec!= null) {


                // TODO count instances for histogram (x sites with 1,2,3,... variants)f
//                if (vvec.size()> 1)
//                    Log.warn("Site "+ spliceSites[i]+ " has "+ vvec.size()+ " variants.");

                // map to last included position, relative to exon boundary (splice site pos)
                if (spliceSites[i].isAcceptor()) {
                    flank5+= 2;
                } else {// donor
                    flank5-= 1;
                }
                int nr= scoreVariants(graph, type, vvec, "", 0, 0, 0, spliceSites[i], flank5, flank3, seq.toLowerCase(), seqs, scores, types, varTuples);
                assert((nr+ 1)== nrCombinations);
            }

            outputSitesVCF(spliceSites[i], vvec, seqs, scores, types, varTuples);
        }

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
                if (s.startsWith("#"))
                    continue;
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



}
