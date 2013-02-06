package barna.genome.tools.overlap;

import barna.io.FileHelper;
import barna.model.commons.IntVector;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 12/12/12
 * Time: 12:07 AM
 * To change this template use File | Settings | File Templates.
 */
public class OverlapVcf {

    public static void main(String[] args) {

        int minInd= 0, minReads= 0;
        for (int i = 22; i >=1; --i) {
            System.err.println("chr"+ i);
            getMutRate(Integer.toString(i));
            //getMutRateDS(Integer.toString(i), true, minInd, minReads);
            //getMutRateDS(Integer.toString(i), false, minInd, minReads);
        }
//        System.err.println("chrX");
//        getMutRateDS("X");
//        System.err.println("chrY");
//        getMutRateDS("Y");

    }


    static void getMutRate(String chrID) {
        String fpath= "/Volumes/DB/Users/micha/Desktop/Geuvadis/external/GEUVADIS.chr"+chrID+".PH1PH2_465.IMPFRQFILT_BIALLELIC.genotypes.vcf.gz";
        // String rpath= "/Volumes/DB/Users/micha/Desktop/Geuvadis/novel_sites/chr"+chrID+"_spld_reads_TANDEM.bed";
        String rpath= "/Volumes/DB/Users/micha/Desktop/Geuvadis/novel_sites/annotv2.acceptors.pos.bed_"+ chrID;

        long tstart= System.currentTimeMillis();
        try {
            IntVector i1= new IntVector(), i2= new IntVector(), i3= new IntVector(), i4= new IntVector(),
                    i5= new IntVector(), i6= new IntVector(), i7= new IntVector(), i8= new IntVector();
            BufferedReader buffy= new BufferedReader(new FileReader(rpath));
            String s= null;
            int offset5=30, max3= 90;
            int c= 0;
            while ((s= buffy.readLine())!= null) {
                ++c;
                String[] ss= s.split("\\s");
                int start= Integer.parseInt(ss[1]);
                int end= Integer.parseInt(ss[2]);
                int len= (end- start);

                if(ss[ss.length- 1].equals("+")) {
                    if (ss[3].contains("^")) {
                        i1.add(start- offset5);
                        i2.add(end+ (max3- len- offset5));
                    } else {
                        i3.add(start- (max3- len- offset5));
                        i4.add(end+ offset5);
                    }
                } else {
                    if (ss[3].contains("^")) {
                        i5.add(start- (max3- len- offset5));
                        i6.add(end+ offset5);
                    } else {
                        i7.add(start- offset5);
                        i8.add(end+ (max3- len- offset5));
                    }
                }

                if (c%100== 0)
                    ;//System.err.println("\t"+ c+ " lines read.");;
            }

            System.err.println("\tregions read "+ i1.size()+", "+i2.size()+", "+i3.size()+", "+i4.size());

            int[] donPstart= i1.toIntArray(), donPend= i2.toIntArray(),
                    accPstart= i3.toIntArray(), accPend= i4.toIntArray(),
                    donNstart= i5.toIntArray(), donNend= i6.toIntArray(),
                    accNstart= i7.toIntArray(), accNend= i8.toIntArray();
            System.err.println("\tconverted");

            i1= null; i2= null; i3= null; i4= null;
            i5= null; i6= null; i7= null; i8= null;
            System.gc();

            // TODO sanity check >=

            int[] don= new int[91], acc= new int[91];

            try {
                buffy= new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fpath))));
                while((s= buffy.readLine())!= null) {

                    // System.err.println(s);
                    if (s.startsWith("#"))
                        continue;
                    int i;
                    String chr= null;
                    for (i = 0; i < s.length(); i++) {
                        if (s.charAt(i)== ' '|| s.charAt(i)=='\t') {
                            chr= s.substring(0, i);
                            break;
                        }
                    }
                    while(s.charAt(i)== ' '|| s.charAt(i)=='\t')
                        ++i;
                    int p= i;
                    while(s.charAt(i)!= ' '&& s.charAt(i)!='\t')
                        ++i;
                    int pos= Integer.parseInt(s.substring(p, i));

                    // donors
                    int pp= Arrays.binarySearch(donPstart, pos);
                    if (pp < 0)
                        pp= -(pp+ 1)- 1;
                    if (pp>= 0 && pos<= donPend[pp]) {
                        int delta= pos- donPstart[pp];
                        ++don[delta];
                    }
                    pp= Arrays.binarySearch(donNstart, pos);
                    if (pp < 0)
                        pp= -(pp+ 1)- 1;
                    if (pp>= 0 && pos<= donNend[pp]) {
                        int delta= donNend[pp]- pos;
                        ++don[delta];
                    }

                    pp= Arrays.binarySearch(accPstart, pos);
                    if (pp < 0)
                        pp= -(pp+ 1)- 1;
                    if (pp>= 0 && pos<= accPend[pp]) {
                        int delta= accPend[pp]- pos;
                        ++acc[delta];
                    }
                    pp= Arrays.binarySearch(accNstart, pos);
                    if (pp < 0)
                        pp= -(pp+ 1)- 1;
                    if (pp>= 0 && pos<= accNend[pp]) {
                        int delta= pos- accNstart[pp];
                        ++acc[delta];
                    }

                }
            } catch (Throwable t) {
                t.printStackTrace();
            }

            System.err.println("\twriting results");
            BufferedWriter bufDon= new BufferedWriter(new FileWriter(FileHelper.append(rpath,"_don",false, "")));
            BufferedWriter bufAcc= new BufferedWriter(new FileWriter(FileHelper.append(rpath,"_acc",false, "")));
            for (int i = 0; i < don.length; i++) {
                bufDon.write(Integer.toString(i+1)+ " "+ Integer.toString(don[i])+ "\n");
                bufAcc.write(Integer.toString(i+1)+ " "+ Integer.toString(acc[i])+ "\n");
            }
            bufDon.close();
            bufAcc.close();

        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        System.err.println((System.currentTimeMillis()- tstart)/ 1000);
    }

    static void getMutRateDS(String chrID, boolean ext, int minInd, int minReads) {
        String fpath= "/Volumes/DB/Users/micha/Desktop/Geuvadis/external/GEUVADIS.chr"+chrID+".PH1PH2_465.IMPFRQFILT_BIALLELIC.genotypes.vcf.gz";
        String rpath= "/Volumes/DB/Users/micha/Desktop/Geuvadis/novel_sites/chr"+chrID+"_spld_reads_TANDEM.ref-nref.bed";

        long tstart= System.currentTimeMillis();
        try {
            IntVector i1= new IntVector(), i2= new IntVector(), i3= new IntVector(), i4= new IntVector(),
                    i5= new IntVector(), i6= new IntVector(), i7= new IntVector(), i8= new IntVector();
            BufferedReader buffy= new BufferedReader(new FileReader(rpath));
            String s= null;
            int offset5=30, max3= 90;
            int c= 0;
            while ((s= buffy.readLine())!= null) {
                ++c;
                String[] ss= s.split("\\s");
                int start= Integer.parseInt(ss[1]);
                int end= Integer.parseInt(ss[2]);
                int len= (end- start);
                String[] id= ss[3].split("@");
                int indiv= Integer.parseInt(id[1]);
                int score= Integer.parseInt(ss[4]);
                if (indiv< minInd || score< minReads)
                    continue;

                if(ss[ss.length- 1].equals("+")) {
                    if (ss[3].contains("^")) {
                        i1.add((ext? end: start)- offset5);
                        i2.add((ext? end: start)+ offset5+ 1);
                    } else {
                        i3.add((ext? start: end)- offset5- 1);
                        i4.add((ext? start: end)+ offset5);
                    }
                } else {
                    if (ss[3].contains("^")) {
                        i5.add((ext? start:end)- offset5- 1);
                        i6.add((ext? start:end)+ offset5);
                    } else {
                        i7.add((ext? end:start)- offset5);
                        i8.add((ext? end:start)+ offset5+ 1);
                    }
                }

                if (c%100== 0)
                    ;//System.err.println("\t"+ c+ " lines read.");;
            }

            System.err.println("\tregions read");

            int[] donPstart= i1.toIntArray(), donPend= i2.toIntArray(),
                    accPstart= i3.toIntArray(), accPend= i4.toIntArray(),
                    donNstart= i5.toIntArray(), donNend= i6.toIntArray(),
                    accNstart= i7.toIntArray(), accNend= i8.toIntArray();
            System.err.println("\tconverted");

            i1= null; i2= null; i3= null; i4= null;
            i5= null; i6= null; i7= null; i8= null;
            System.gc();

            // TODO sanity check >=

            int[] don= new int[62], acc= new int[62];

            try {
                buffy= new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fpath))));
                while((s= buffy.readLine())!= null) {

                    // System.err.println(s);
                    if (s.startsWith("#"))
                        continue;
                    int i;
                    String chr= null;
                    for (i = 0; i < s.length(); i++) {
                        if (s.charAt(i)== ' '|| s.charAt(i)=='\t') {
                            chr= s.substring(0, i);
                            break;
                        }
                    }
                    while(s.charAt(i)== ' '|| s.charAt(i)=='\t')
                        ++i;
                    int p= i;
                    while(s.charAt(i)!= ' '&& s.charAt(i)!='\t')
                        ++i;
                    int pos= Integer.parseInt(s.substring(p, i));

                    // donors
                    int pp= Arrays.binarySearch(donPstart, pos);
                    if (pp < 0)
                        pp= -(pp+ 1)- 1;
                    if (pp>= 0 && pos<= donPend[pp]) {
                        int delta= pos- donPstart[pp];
                        ++don[delta];
                    }
                    pp= Arrays.binarySearch(donNstart, pos);
                    if (pp < 0)
                        pp= -(pp+ 1)- 1;
                    if (pp>= 0 && pos<= donNend[pp]) {
                        int delta= donNend[pp]- pos;
                        ++don[delta];
                    }

                    pp= Arrays.binarySearch(accPstart, pos);
                    if (pp < 0)
                        pp= -(pp+ 1)- 1;
                    if (pp>= 0 && pos<= accPend[pp]) {
                        int delta= accPend[pp]- pos;
                        ++acc[delta];
                    }
                    pp= Arrays.binarySearch(accNstart, pos);
                    if (pp < 0)
                        pp= -(pp+ 1)- 1;
                    if (pp>= 0 && pos<= accNend[pp]) {
                        int delta= pos- accNstart[pp];
                        ++acc[delta];
                    }

                }
            } catch (Throwable t) {
                t.printStackTrace();
            }

            System.err.println("\twriting results");
            String sfx= (ext? "_extd": "_real")+ "_"+ minInd+ "_"+ minReads;
            BufferedWriter bufDon= new BufferedWriter(new FileWriter(FileHelper.append(rpath,"_don"+sfx,false, "tab")));
            BufferedWriter bufAcc= new BufferedWriter(new FileWriter(FileHelper.append(rpath,"_acc"+sfx,false, "tab")));
            for (int i = 0; i < don.length; i++) {
                bufDon.write(Integer.toString(i+1)+ " "+ Integer.toString(don[i])+ "\n");
                bufAcc.write(Integer.toString(i+1)+ " "+ Integer.toString(acc[i])+ "\n");
            }
            bufDon.close();
            bufAcc.close();

        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        System.err.println((System.currentTimeMillis()- tstart)/ 1000);
    }
}
