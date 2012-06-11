package barna.flux.capacitor.graph;

import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.io.BufferedIterator;
import barna.io.bed.BEDwrapper;
import barna.io.gtf.GTFwrapper;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Gene;
import groovy.ui.SystemOutputInterceptor;
import junit.framework.Assert;
import junit.framework.TestCase;
import org.junit.Test;

import java.io.*;
import java.util.*;
import java.util.zip.ZipOutputStream;

/**
 * Created with IntelliJ IDEA.
 * User: Emilio Palumbo
 * Date: 5/31/12
 * Time: 12:39 PM
 */
public class AnnotationMapperTest extends TestCase {

    private final File gtfFile = new File(getClass().getResource("/mm9_chr1_chrX.gtf").getFile());
    private final File bedFile = new File(getClass().getResource("/chr1_chrX_uniq.bed").getFile());//("/home/emilio/tmp/chr1_chrX.bed");
    private FluxCapacitorSettings settings;
    private boolean paired =false;

    @Override
    public void setUp() throws Exception {
        super.setUp();
        setUpSettings();
    }

    private void setUpSettings() {
        UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMULATOR));
        settings = new FluxCapacitorSettings();
        settings.set(FluxCapacitorSettings.ANNOTATION_FILE,
                new File(gtfFile.getAbsolutePath()));
        settings.set(FluxCapacitorSettings.MAPPING_FILE,
                new File(bedFile.getAbsolutePath()));
        settings.set(FluxCapacitorSettings.READ_DESCRIPTOR,
                descriptor);
        settings.set(FluxCapacitorSettings.SORT_IN_RAM,
                false);
        settings.set(FluxCapacitorSettings.KEEP_SORTED_FILES,
                false);
        settings.set(FluxCapacitorSettings.ANNOTATION_MAPPING,
                FluxCapacitorSettings.AnnotationMapping.PAIRED);
        settings.set(FluxCapacitorSettings.STDOUT_FILE,
                null);
        settings.set(FluxCapacitorSettings.STATS_FILE,
                null);
        paired = settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).name().equals("PAIRED")?true:false;
    }

    private String[] getId(String [] split, int nBlocks) {
        String[] sjs = new String[nBlocks-1];
        for (int i =0;i<nBlocks-1;i++) {
            int sjStart = Integer.parseInt(split[1])+Integer.parseInt(split[11].split(",")[i])+Integer.parseInt(split[10].split(",")[i]);
            int sjEnd = Integer.parseInt(split[1])+Integer.parseInt(split[11].split(",")[i+1])+1;
            sjs[i] = sjStart + "-" + sjEnd ;
        }
        return sjs;
    }

    private String getAltStrand(String s) {
        if (s.equals("S"))
            return "A";
        if (s.equals("A"))
            return "S";
        return "";
    }


    private Map<String,Integer> getSJReads(Gene g, boolean paired) throws Exception{
        BufferedReader bedReader = new BufferedReader(new InputStreamReader(new FileInputStream(bedFile)));
        /*BufferedReader gtfReader = new BufferedReader(new InputStreamReader(new FileInputStream(gtfFile)));
        ArrayList<String[]> nodes = new ArrayList<String[]>();*/
        Map<String,Integer> reads = new TreeMap<String, Integer>();
        HashMap<String,ArrayList<String[]>> p1hash = new HashMap<String, ArrayList<String[]>>();
        HashMap<String,ArrayList<String[]>> p2hash = new HashMap<String, ArrayList<String[]>>();
        int nr = 0;
        /*String tx = null;
        for ( String line; (line = gtfReader.readLine())!= null;) {
            String[] gLine = line.split("\t");
            if (tx==null)
                    tx = gLine[8].split(";")[0].split("\\s")[1];
            if (gLine[2].equals("exon")) {
                if (tx.equals(gLine[8].split(";")[0].split("\\s")[1])) {
                    if (nodes.size() != 0) {
                        if (!nodes.contains(new String[]{nodes.get(nodes.size()-1)[1],gLine[3]}))
                            nodes.add(new String[]{nodes.get(nodes.size()-1)[1],gLine[3]});
                    }
                    if (!nodes.contains(new String[]{gLine[3],gLine[4]}))
                        nodes.add(new String[]{gLine[3],gLine[4]});
                } else {
                    tx = gLine[8].split(";")[0].split("\\s")[1];
                    if (!nodes.contains(new String[]{gLine[3],gLine[4]}))
                        nodes.add(new String[]{gLine[3],gLine[4]});
                }
            }
        }*/
        for ( String line; (line = bedReader.readLine())!= null;) {
            String[] bLine = line.split("\t");
            int nBlocks = Integer.parseInt(bLine[9]);
            if (g.getGeneID().equals(bLine[3].split(":")[0]+":"+ bLine[3].split(":")[1]))
            {
                if (paired) {
                    if (nBlocks<3) {
                        String readId = bLine[3].substring(0,bLine[3].length()-4);
                        if (bLine[3].endsWith("1")) {
                            if (!p1hash.containsKey(readId)) {
                                String sj = nBlocks==2?getId(bLine,nBlocks)[0]:null;
                                ArrayList<String[]> list = new ArrayList<String[]>();
                                list.add(new String[]{bLine[3].substring(bLine[3].length()-3,bLine[3].length()-2),sj});
                                p1hash.put(readId, list);
                            } else {
                                String sj = nBlocks==2?getId(bLine,nBlocks)[0]:null;
                                ArrayList<String[]> list = p1hash.get(readId);
                                list.add(new String[]{bLine[3].substring(bLine[3].length()-3,bLine[3].length()-2), sj});
                                p1hash.put(readId, list);
                            }
                        }
                        if (bLine[3].endsWith("2")) {
                            if (!p2hash.containsKey(readId)) {
                                String sj = nBlocks==2?getId(bLine,nBlocks)[0]:null;
                                ArrayList<String[]> list = new ArrayList<String[]>();
                                list.add(new String[]{bLine[3].substring(bLine[3].length()-3,bLine[3].length()-2),sj});
                                p2hash.put(readId, list);
                            } else {
                                String sj = nBlocks==2?getId(bLine,nBlocks)[0]:null;
                                ArrayList<String[]> list = p2hash.get(readId);
                                list.add(new String[]{bLine[3].substring(bLine[3].length()-3,bLine[3].length()-2), sj});
                                p2hash.put(readId, list);
                            }
                        }
                    }
                } else {
                    if (nBlocks==2) {
                        String[] ids = getId(bLine,nBlocks);
                        for (String id : ids) {
                            nr = 1;
                            if (reads.containsKey(id)) {
                                nr += reads.get(id);
                            }
                            reads.put(id,nr);
                        }
                    }
                }
            }
        }
        if (paired) {
            for (String id : p1hash.keySet()) {
                if (p2hash.containsKey(id)) {
                    ArrayList<String[]> p1 = p1hash.get(id);
                    ArrayList<String[]> p2 = p2hash.get(id);
                    for (String[] s1 : p1) {
                        for (String[] s2 : p2) {
                            if (s2[0].equals(getAltStrand(s1[0]))) {
                                if ((s1[1]!=null || s2[1] != null)) {
                                    String sj = s1[1]!=null?s1[1]:s2[1];
                                    nr = (s1[1]!=null?1:0) + (s2[1]!=null?1:0);
                                    if (reads.containsKey(sj)) {
                                        nr += reads.get(sj);
                                    }
                                    reads.put(sj,nr);
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
        return reads;
    }


    @Test
    public void testWriteSJReads() throws Exception {
        GTFwrapper gtf = new GTFwrapper(gtfFile);
        BEDwrapper bed = new BEDwrapper(bedFile);
        try {
            gtf = new GTFwrapper((gtf.sort()));
            gtf.setReadAll(true);
            gtf.setNoIDs(null);
            gtf.setReadFeatures(new String[]{"exon","CDS"});
            gtf.read();
            Gene g = gtf.getGenes()[0];
            int start = 0,end=0,tol=0;
            start = g.getStart();
            end=g.getEnd();
            if (g.getStrand()< 0) {
                start= -start;
                end= -end;
            }
            tol= 0;
            start= Math.max(1, start- tol);
            end= end+ tol;
            BufferedIterator iter = bed.readBedFile(g, start, end, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR),null);
            AnnotationMapper a = new AnnotationMapper(g);
            a.map(iter, settings);
            ZipOutputStream out = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream("/home/emilio/ann_mapper.zip")));
            a.writeSJReads(out, paired);
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
            Assert.fail();
        }
    }

    @Test
    public void testCompareSJReads() throws Exception {
        GTFwrapper gtf = new GTFwrapper(gtfFile);
        BEDwrapper bed = new BEDwrapper(bedFile);
        gtf = new GTFwrapper((gtf.sort()));
        gtf.setReadAll(true);
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon","CDS"});
        gtf.read();
        Gene g = gtf.getGenes()[0];
        int start = 0,end=0,tol=0;
        start = g.getStart();
        end=g.getEnd();
        if (g.getStrand()< 0) {
            start= -start;
            end= -end;
        }
        tol= 0;
        start= Math.max(1, start- tol);
        end= end+ tol;
        BufferedIterator iter = bed.readBedFile(g, start, end, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR),null);
        AnnotationMapper a = new AnnotationMapper(g);
        a.map(iter, settings);
        Map<String,Integer> m = a.getSJReads(paired);
        int count[] = new int[]{0,0};
        for (String e : m.keySet()) {
            System.err.println(e+ "\t"+ m.get(e));
            count[0]+=m.get(e);
        }
        System.err.println();
        Map<String,Integer> m1 = getSJReads(g, paired);
        for (String e : m1.keySet()) {
            System.err.println(e+ "\t"+ m1.get(e));
            count[1]+=m1.get(e);
        }
        assertEquals(count[1],count[0]);
    }

    @Test
    public void testCompareIntronReads() throws Exception {
        GTFwrapper gtf = new GTFwrapper(gtfFile);
        BEDwrapper bed = new BEDwrapper(bedFile);
        gtf = new GTFwrapper((gtf.sort()));
        gtf.setReadAll(true);
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon","CDS"});
        gtf.read();
        Gene g = gtf.getGenes()[0];
        int start = 0,end=0,tol=0;
        start = g.getStart();
        end=g.getEnd();
        if (g.getStrand()< 0) {
            start= -start;
            end= -end;
        }
        tol= 0;
        start= Math.max(1, start- tol);
        end= end+ tol;
        BufferedIterator iter = bed.readBedFile(g, start, end, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR),null);
        AnnotationMapper a = new AnnotationMapper(g);
        a.map(iter, settings);
        Map<String,Integer> m = a.getIntronReads(paired);
        int count[] = new int[]{0,0};
        for (String e : m.keySet()) {
            count[0]+=m.get(e);
        }
        //((SimpleEdgeIntronMappings)a.getNodesInGenomicOrder()[2].getOutEdges().elementAt(0)).getReadDist(new int[]{2,2,2,2,2,2,2});
        assertTrue(count[0]>0);
    }

}
