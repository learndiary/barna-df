package barna.io.vcf;

import barna.commons.log.Log;
import barna.io.FileHelper;

import java.io.*;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 1/9/13
 * Time: 3:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class VCFwrapper {

    /**
     * The file this wrapper is based on.
     */
    File vcf= null;

    /**
     * Flag to read SNPs
     */
    boolean readSNP= true;

    /**
     * Flag to read indels
     */
    boolean readIndel= false;

    /**
     * A hash representing the variant information.
     */
    HashMap<String, String> map= null;

    public VCFwrapper(File vcf) {
        this.vcf= vcf;
    }

    /**
     * Reads vcf file and fills a hash with position x snp information.
     * @return hash representing the information of the provided file
     */
    public HashMap<String, String> getVariants() {

        if (map== null) {

            try {
                map= new HashMap<String, String>((int) (vcf.length()/ 1000));
                BufferedReader buffy= null;
                byte comp= FileHelper.getCompression(vcf);
                if (comp== FileHelper.COMPRESSION_NONE)
                    buffy= new BufferedReader(new FileReader(vcf));
                else if (comp== FileHelper.COMPRESSION_GZIP)
                    buffy= new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(vcf))));

                StringTokenizer t;
                int cntIgnore= 0;
                for (String s= null; (s= buffy.readLine())!= null; ) {
                    if (s.startsWith("#"))
                        continue;
                    t= new StringTokenizer(s, "\t");
                    String loc= t.nextToken()+ "@"+ t.nextToken();  // chrNr + @ + position
                    String snpID= t.nextToken();
                    if (snpID.startsWith("snp")&& (!readSNP)) {
                        ++cntIgnore;
                        continue;
                    }
                    if (snpID.startsWith("indel")&& (!readIndel)) {
                        ++cntIgnore;
                        continue;
                    }

                    String ref= t.nextToken();
                    String var= t.nextToken();
                    map.put(loc, snpID+ "@"+ ref+ "@"+ var); // snpID + @ + ref String + @ + var string
//                    if (map.size()% 1000== 0)
//                        System.err.println(map.size());
                }

                System.err.println("Read "+map.size()+ " variants, ignored "+ cntIgnore+ (!readSNP?" SNPs":"")+ (!readIndel?" indels":""));
                return map;
            } catch (Exception e) {
                Log.error("Error reading VCF file");
                throw new RuntimeException(e);
            }

        }

        return map;
    }
}
