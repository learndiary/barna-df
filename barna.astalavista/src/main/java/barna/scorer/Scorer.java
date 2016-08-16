package barna.scorer;

import barna.astalavista.AStalavista;
import barna.astalavista.AStalavistaSettings;
import barna.commons.Execute;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.log.Log;
import barna.geneid.*;
import barna.io.FileHelper;
import barna.model.Gene;
import barna.model.Graph;
import barna.model.SpliceSite;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.*;
import java.net.URL;
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


    /**
     * Run a Scorer with the specified parameters.
     *
     * @param settings the parameters in a settings object
     * @param params alternatively, a vector of command line arguments
     */
    public static void runIt(ScorerSettings settings, String[] params) {
        Scorer aScorer= new Scorer();
        if (params != null) {
            JSAP jsap = new JSAP();
            for (Parameter p : aScorer.getParameter()) {
                try {
                    jsap.registerParameter(p);
                } catch (JSAPException e) {
                    e.printStackTrace();
                }
            }
            aScorer.validateParameter(jsap.parse(params));
        } else {
            aScorer.setSettings(settings);
            if (!aScorer.validateSettings(settings))
                return;
        }

        Future<Void> captain= Execute.getExecutor().submit(aScorer);
        try {
            Object o = captain.get();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }
     }


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
        scoreSites(g.getSpliceSites());
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

        return validateSettings(settings);
    }

    public boolean validateSettings(AStalavistaSettings settings) {

        if (!super.validateSettings(settings))
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
                String pFileName= (settings.get(ScorerSettings.GENE_ID)== null)?
                        null: settings.get(ScorerSettings.GENE_ID).getAbsolutePath();

                // load it fro file provided, or take default matrices (human)
                geneidParam= Profile.readParam(pFileName, new GeneIDsettings())[0];
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

    /**
     * Scans for variants in the range of the splice site sequence. If SNPs are found,
     * the correspondingly mutated sequences are scored.
     *
     * @param ss the splice site
     * @param flank5 upstream sequence flanking the dinucleotide
     * @param flank3 downstream sequence flanking the dinucleotide
     * @return vector of strings
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

            // correct for neg.strand (VCF file is in genomic positions, indels are to be shift)
            int del= vv[3].length()- 1;
            if (ss.getPos()< 0) {
                snPos+= del;
                vv[3]= Graph.reverseSequence(Graph.complementarySequence(vv[3]));
            }
            int p= snPos- (ss.getPos()- flank5);
            if (p< 0|| p>= seq.length())
                continue;   // after correction (pe neg.strand) out of site area


            // BARNA-317 split for multiple substitution possibilities
            String[] vvv= vv[4].split(",");
            for (int j = 0; j < vvv.length; j++) {
                boolean deletion= vv[3].length()> vvv[j].length();
                boolean insertion= vv[3].length()< vvv[j].length();
                boolean substitution= vv[3].length()== vvv[j].length();

                if (ss.getPos()< 0) {
                    vvv[j]= Graph.reverseSequence(Graph.complementarySequence(vvv[j]));
                }

                // VCF already provides strand-specific bases
                String varSeq= null;
                if (p+ vvv[j].length()- del<= seq.length()) {

                    varSeq= seq.substring(0, p)+ vvv[j];
                    if (p+ del< seq.length())
                        varSeq+= seq.substring(p+ del+ 1);

                    // trim start
                    int start= 0;
                    if (insertion&& p< flank5)
                        start= vvv[j].length()- vv[3].length();
                    // trim end
                    if (varSeq.length()> seq.length())  // dirty
                        varSeq= varSeq.substring(0, seq.length());
                } else
                    varSeq= seq.substring(0, p)+ vvv[j].substring(0, seq.length()- p);
                if (varSeq.length()< seq.length()) {
                    int missing= seq.length()- varSeq.length();
                    // fill with upstream seq
                    if (p+ missing< flank5) {
                        String s= Graph.readSequence(ss, flank5+ missing- (ss.isAcceptor()? 2: 0), 0);
                        varSeq= s.substring(0, missing).toLowerCase()+ varSeq;
                    } else {    // fill with downstream seq
                        int skip= p+ del+ 1- seq.length(); // to be skipped ds
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
     * @param spliceSites vector of (splice) sites
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

            // int nrCombinations= (int) Math.pow(2, vvec== null? 0: vvec.size());
            // BARNA-317: unfortunately not that easy, one line may contain multiple ALT values (comma-separated)
            //
            int nrCombinations= 1;
            for(int j=0; (vvec!= null) && j< vvec.size(); ++j) {
                String[] vv= vvec.elementAt(j).split("@");
                String[] vvv= vv[4].split(",");
                nrCombinations*= (vvv.length+ 1); // total nr.= nr. variants + ref
            }

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


                // TODO count instances for histogram (x sites with 1,2,3,... variants)f
//                if (vvec.size()> 1)
//                    Log.warn("Site "+ spliceSites[i]+ " has "+ vvec.size()+ " variants.");

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
