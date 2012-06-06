package barna.flux.capacitor.graph;

import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.io.BufferedIterator;
import barna.io.BufferedIteratorRAM;
import barna.io.bed.BEDwrapper;
import barna.io.gtf.GTFwrapper;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Gene;
import barna.model.Transcript;
import barna.model.bed.BEDMapping;
import junit.framework.Assert;
import junit.framework.TestCase;
import org.junit.Test;

import java.awt.datatransfer.StringSelection;
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
    private final File bedFile = new File(getClass().getResource("/chr1_chrX.bed").getFile());//("/home/emilio/tmp/chr1_chrX.bed");
    private FluxCapacitorSettings settings;
    private boolean paired =false;

    @Override
    public void setUp() throws Exception {
        super.setUp();
        setSettings();
    }

    private void setSettings() {
        UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMPLE));
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
                FluxCapacitorSettings.AnnotationMapping.SINGLE);
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


    private Map<String,Integer> getSJReads(Gene g, boolean paired) throws Exception{
        BufferedReader bedReader = new BufferedReader(new InputStreamReader(new FileInputStream(bedFile)));
        Map<String,Integer> reads = new TreeMap<String, Integer>();
        HashMap<String,String[]> hash = new HashMap<String, String[]>();
        HashMap<String,String[]> p2Hash = new HashMap<String, String[]>();
        int nr = 0;
        for ( String line; (line = bedReader.readLine())!= null;) {
            String[] bLine = line.split("\t");
            int nBlocks = Integer.parseInt(bLine[9]);
            if (g.getGeneID().equals(bLine[3].split(":")[0]+":"+ bLine[3].split(":")[1]))
            {
                if (paired) {
                    String readId = bLine[3].substring(0,bLine[3].length()-4);

                    if (!hash.containsKey(readId)) {
                        String[] sj = nBlocks==2?getId(bLine,nBlocks):null;
                        hash.put(readId,sj);
                    } else {
                        if (nBlocks<=2) {
                            String[] ids = hash.get(readId);
                            if (ids != null) {
                                nr = ids.length+nBlocks-1;
                            } else {
                                ids = getId(bLine,nBlocks);
                                nr = ids.length;
                            }
                            for (String id : ids) {
                                if (reads.containsKey(id)) {
                                    nr += reads.get(id);
                                }
                                reads.put(id,nr);
                            }
                            hash.remove(readId);
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
        /*if (paired) {
            for (String k : hash.keySet()) {
                if (hash.get(k)!=null) {
                    String[] ids = hash.get(k);
                    for (String id : ids) {
                        nr = 1;
                        if (reads.containsKey(id)) {
                            nr += reads.get(id);
                        }
                        reads.put(id,nr);
                    }
                }
            }
        }*/
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
            count[0]+=m.get(e);
        }
        Map<String,Integer> m1 = getSJReads(g, paired);
        for (String e : m1.keySet()) {
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
        assertTrue(count[0]>0);
    }

}
