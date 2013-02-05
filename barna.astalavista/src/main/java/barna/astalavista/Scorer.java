package barna.astalavista;

import barna.commons.log.Log;
import barna.geneid.*;
import barna.io.FileHelper;
import barna.model.Gene;
import barna.model.Graph;
import barna.model.SpliceSite;
import com.martiansoftware.jsap.JSAPResult;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Vector;

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
    protected void callFinish() throws Exception {
        super.callFinish();
        if (siteScoreWriter!= null)
            siteScoreWriter.close();
    }

    @Override
    public boolean validateParameter(JSAPResult args) {

        // check splice site scoring stuff
        if (settings.get(ScorerSettings.SITES_OPT).contains(ScorerSettings.SiteOptions.SSS)) {

            if(settings.get(AStalavistaSettings.CHR_SEQ)== null) {
                Log.error("Splice site scoring requires the genomic sequence, set parameter \'CHR_SEQ\'");
                return false;
            }

            try {
                if(settings.get(ScorerSettings.GENE_ID)== null) {
                    Log.warn("No GeneID parameter file for scoring models, using default");
                    geneidParam= Profile.readParam(GeneIDconstants.PARAMETERFILE, new GeneIDsettings())[0];
                } else
                    geneidParam= Profile.readParam(settings.get(ScorerSettings.GENE_ID).getAbsolutePath(),
                            new GeneIDsettings())[0];
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

    private void getSiteScoreWriter() {
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



}
