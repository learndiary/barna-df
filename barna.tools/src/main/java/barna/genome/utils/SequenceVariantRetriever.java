package barna.genome.utils;

import barna.io.vcf.VCFwrapper;
import barna.model.DirectedRegion;
import barna.model.Graph;

import java.io.*;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.zip.GZIPInputStream;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 1/9/13
 * Time: 3:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class SequenceVariantRetriever implements Callable {

    public static void main(String[] args){

        args= new String[] {
            "/Volumes/DB/Users/micha/Desktop/Geuvadis/novel_sites/polyA/polyA_signals_chr22",
            "/Volumes/DB/Users/micha/Desktop/Geuvadis/external/GEUVADIS.chr22.PH1PH2_465.IMPFRQFILT_BIALLELIC.genotypes.vcf",
            "/Volumes/Raptor/genomes/hg19"};

        if (args== null|| args.length!= 3) {
            System.err.println("Usage: "+ SequenceVariantRetriever.class.getSimpleName()+ " coordinates vcf genomeDir");
            System.exit(-1);
        }

        for (int i = 8; i >= 8; --i) {

            args= new String[] {
                    "/Volumes/DB/Users/micha/Desktop/Geuvadis/novel_sites/polyA/polyA_signals_chr"+ i,
                    "/Volumes/DB/Users/micha/Desktop/Geuvadis/external/GEUVADIS.chr"+i+".PH1PH2_465.IMPFRQFILT_BIALLELIC.genotypes.vcf.gz",
                    "/Volumes/Raptor/genomes/hg19"};


            SequenceVariantRetriever myRetriever= new SequenceVariantRetriever(args);
            myRetriever.call();
        }

    }

    String[] args= null;

    public SequenceVariantRetriever(String[] args) {
        this.args= args;
    }


    public Object call() {

        System.err.println("Start "+ args[0]);
        long t0= System.currentTimeMillis();

        // init genome dir
        Graph.overrideSequenceDirPath= args[2];

        // init variants
        System.err.print("Reading variants..");
        System.err.flush();
        File vcf= new File(args[1]);
        VCFwrapper vcfReader= new VCFwrapper(vcf);
        HashMap<String, String> variantMap= vcfReader.getVariants();
        System.err.println("done.");

        try {
            BufferedReader buffy= new BufferedReader(new FileReader(args[0]));
            BufferedWriter writer= new BufferedWriter(new FileWriter(args[0]+"_variants.txt"));
            for (String s; (s= buffy.readLine())!= null;) {
                String[] ss= s.split("\\s");
                DirectedRegion reg= new DirectedRegion(
                        Integer.parseInt(ss[1]), Integer.parseInt(ss[2]), (ss[3]=="+"?1:-1));
                reg.setChromosome(ss[0]);

                if(Math.abs(reg.getStart()- reg.getEnd())+ 1!= 6) {
                    System.err.println("skipped "+ s);
                    continue;
                }

                String[][] vv= Graph.readSequence(reg, variantMap);
                if (vv== null|| vv.length== 0)
                    continue;
                StringBuffer sb= new StringBuffer(100);
                for (int i = 0; vv!= null&& i < vv[0].length; i++) {
                   sb.append(vv[0][i]+ "="+ vv[1][i]+";");
                }
                writer.write(ss[0] + "\t" + ss[1] + "\t" + ss[2] + "\t" + sb.substring(0, sb.length() - 1) + "\t" + ss[3] + "\n");

            }
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.err.println(((System.currentTimeMillis()- t0)/ 1000) + " sec");
        return null;
    }
}
