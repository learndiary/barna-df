package barna.flux.capacitor.graph;

import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.io.BufferedIterator;
import barna.io.bed.BEDwrapper;
import barna.io.gtf.GTFwrapper;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Gene;
import com.sun.org.apache.bcel.internal.generic.RETURN;
import com.sun.org.apache.xpath.internal.NodeSet;
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

    private String path = "/home/emilio/fromMicha";
    private final File gtfFile = new File(path + "/hg19_ref_ucsc120203_sorted.gtf");//(getClass().getResource("/mm9_chr1_chrX.gtf").getFile());//(path+"/hg19_ref_ucsc120203_sorted.gtf");//
    private final File bedFile = new File(path + "/NA12546_NA12546.1.M_120209_gem_2_76-76-50-30_120313170321-1689404293_chr22.bed");//(getClass().getResource("/chr1_chrX.bed").getFile());//
    private FluxCapacitorSettings settings;
    Map<String, ArrayList<String[]>> nodes = new HashMap<String, ArrayList<String[]>>();

    @Override
    public void setUp() throws Exception {
        super.setUp();
        initSettings(UniversalReadDescriptor.DESCRIPTORID_SIMULATOR, FluxCapacitorSettings.AnnotationMapping.PAIRED);
    }

    private void initSettings(String descriptorStr, FluxCapacitorSettings.AnnotationMapping mapping) {
        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor(descriptorStr));
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
                mapping);
        settings.set(FluxCapacitorSettings.STDOUT_FILE,
                null);
        settings.set(FluxCapacitorSettings.STATS_FILE,
                null);
    }

    private String[] getId(String[] split, int nBlocks) {
        String[] sjs = new String[nBlocks - 1];
        for (int i = 0; i < nBlocks - 1; i++) {
            int sjStart = Integer.parseInt(split[1]) + Integer.parseInt(split[11].split(",")[i]) + Integer.parseInt(split[10].split(",")[i]);
            int sjEnd = Integer.parseInt(split[1]) + Integer.parseInt(split[11].split(",")[i + 1]) + 1;
            sjs[i] = sjStart + "-" + sjEnd;
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

    private void readGtf(Gene g) throws IOException {
        BufferedReader gtfReader = new BufferedReader(new InputStreamReader(new FileInputStream(gtfFile)));
        String tx = null;
        nodes.clear();
        int i = 0;
        int start = 0, end = 0, tol = 0;
        start = g.getStart();
        end = g.getEnd();
        if (g.getStrand() < 0) {
            start = -start;
            end = -end;
        }
        tol = 0;
        start = Math.max(1, start - tol);
        end = end + tol;
        //System.err.println("\n[TEST] Reading GTF for gene "+g.getGeneID()+" ...");
        for (String line; (line = gtfReader.readLine()) != null; ) {
            String[] gLine = line.split("\t");
            if (gLine[0].equals(g.getChromosome()) && gLine[6].equals(g.getStrand() < 0 ? "-" : "+") && Integer.parseInt(gLine[3]) >= start && Integer.parseInt(gLine[4]) <= end) {
                if (tx == null)
                    tx = gLine[8].split(";")[0].split("\\s")[1].replace("\"", "");
                if (gLine[2].equals("exon")) {
                    String[] strings = new String[]{gLine[3], gLine[4]};
                    if (tx.equals(gLine[8].split(";")[0].split("\\s")[1])) {
                        ArrayList<String[]> list = null;
                        if (!nodes.containsKey(tx)) {
                            list = new ArrayList<String[]>();
                            list.add(strings);
                        } else {
                            list = nodes.get(tx);
                            if (!list.contains(strings))
                                list.add(strings);
                        }
                        nodes.put(tx, list);
                    } else {
                        tx = gLine[8].split(";")[0].split("\\s")[1];
                        ArrayList<String[]> list = new ArrayList<String[]>();
                        list.add(strings);
                        nodes.put(tx, list);
                    }
                }
            }
        }
        //System.err.println("[TEST] GTF read.");
    }

    private boolean containsExon(Map<String, String[]> transcripts, String[] exon) {
        for (String t : transcripts.keySet()) {
            if (exon.length != transcripts.get(t).length)
                return false;
            if (exon.length == 2) {
                if (transcripts.get(t)[0].equals(exon[0]) && transcripts.get(t)[1].equals(exon[1])) {
                    return true;
                }
            }
            if (exon.length == 4) {
                if (transcripts.get(t)[0].equals(exon[0]) && transcripts.get(t)[1].equals(exon[1]) &&
                        transcripts.get(t)[2].equals(exon[2]) && transcripts.get(t)[3].equals(exon[3])) {
                    return true;
                }
            }
        }
        return false;
    }

    private Map<String, Integer> getSJReads(Gene g, boolean paired) throws Exception {
        BufferedReader bedReader = new BufferedReader(new InputStreamReader(new FileInputStream(bedFile)));
        Map<String, Integer> reads = new TreeMap<String, Integer>();
        HashMap<String, ArrayList<String[]>> p1hash = new HashMap<String, ArrayList<String[]>>();
        HashMap<String, ArrayList<String[]>> p2hash = new HashMap<String, ArrayList<String[]>>();
        UniversalReadDescriptor rd = new UniversalReadDescriptor();
        rd.init(UniversalReadDescriptor.DESCRIPTORID_CASAVA18);
        int nr = 0;
        int start = 0, end = 0, tol = 0;
        start = g.getStart();
        end = g.getEnd();
        if (g.getStrand() < 0) {
            start = -start;
            end = -end;
        }
        tol = 0;
        start = Math.max(1, start - tol);
        end = end + tol;
        for (String line; (line = bedReader.readLine()) != null; ) {
            String[] bLine = line.split("\t");
            boolean mapped = false;
            Map<String, String[]> transcripts = new HashMap<String, String[]>();
            int nBlocks = Integer.parseInt(bLine[9]);
            if (bLine[0].equals(g.getChromosome()) && Integer.parseInt(bLine[1]) + 1 >= start && Integer.parseInt(bLine[2]) <= end) {
                for (String tx : nodes.keySet()) {
                    for (int i = 0; i < nodes.get(tx).size(); i++) {
                        //for (String[] pos : nodes.get(tx)) {
                        String[] exon = nodes.get(tx).get(i);
                        if (nBlocks == 1) {
                            if (Integer.parseInt(exon[0]) <= Integer.parseInt(bLine[1]) + 1 && Integer.parseInt(exon[1]) >= Integer.parseInt(bLine[2])) {
                                mapped = true;
                                if (!containsExon(transcripts, exon))
                                    transcripts.put(tx, exon);
                                else
                                    transcripts.put("(" + tx + ")", exon);
                                break;
                            }
                        }
                        if (nBlocks == 2 && i < nodes.get(tx).size() - 1) {
                            String[] nextExon = nodes.get(tx).get(i + 1);
                            int[] block1 = {Integer.parseInt(bLine[1]) + Integer.parseInt(bLine[11].split(",")[0]) + 1, Integer.parseInt(bLine[1]) + Integer.parseInt(bLine[11].split(",")[0]) + Integer.parseInt(bLine[10].split(",")[0])};
                            int[] block2 = {Integer.parseInt(bLine[1]) + Integer.parseInt(bLine[11].split(",")[1]) + 1, Integer.parseInt(bLine[1]) + Integer.parseInt(bLine[11].split(",")[1]) + Integer.parseInt(bLine[10].split(",")[1])};
                            if (Integer.parseInt(exon[0]) <= block1[0] && Integer.parseInt(exon[1]) == block1[1] && Integer.parseInt(nextExon[0]) == block2[0] && Integer.parseInt(nextExon[1]) >= block2[1]) {
                                mapped = true;
                                String[] junction = {exon[0], exon[1], nextExon[0], nextExon[1]};
                                if (!containsExon(transcripts, junction))
                                    transcripts.put(tx, junction);
                                else
                                    transcripts.put("(" + tx + ")", junction);
                                break;
                            }
                        }
                    }
                    /*if (mapped)
                        break;*/
                }
                if (mapped)// && g.getGeneID().equals(bLine[3].split(":")[0]+":"+ bLine[3].split(":")[1]))
                {
                    if (paired) {
                        if (nBlocks < 3) {
                            String readId = null;
                            char mate = 0;
                            if (settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).toString().equals(rd.toString())) {
                                String[] readDescTokens = bLine[3].split(" ");
                                readId = readDescTokens[0];
                                mate = readDescTokens[1].charAt(0);
                            } else {
                                readId = bLine[3].substring(0, bLine[3].length() - 4);
                                mate = bLine[3].charAt(bLine[3].length() - 1);
                            }
                            String transcript = null;
                            for (String s : transcripts.keySet()) {
                                transcript = transcript == null ? s : transcript + "," + s;
                            }

                            if (mate == '1') {
                                if (!p1hash.containsKey(readId)) {
                                    String sj = nBlocks == 2 ? getId(bLine, nBlocks)[0] : null;
                                    ArrayList<String[]> list = new ArrayList<String[]>();
                                    list.add(new String[]{bLine[1], bLine[5], sj, transcript});
                                    p1hash.put(readId, list);
                                } else {
                                    String sj = nBlocks == 2 ? getId(bLine, nBlocks)[0] : null;
                                    ArrayList<String[]> list = p1hash.get(readId);
                                    list.add(new String[]{bLine[1], bLine[5], sj, transcript});
                                    p1hash.put(readId, list);
                                }
                            }
                            if (mate == '2') {
                                if (!p2hash.containsKey(readId)) {
                                    String sj = nBlocks == 2 ? getId(bLine, nBlocks)[0] : null;
                                    ArrayList<String[]> list = new ArrayList<String[]>();
                                    list.add(new String[]{bLine[1], bLine[5], sj, transcript});
                                    p2hash.put(readId, list);
                                } else {
                                    String sj = nBlocks == 2 ? getId(bLine, nBlocks)[0] : null;
                                    ArrayList<String[]> list = p2hash.get(readId);
                                    list.add(new String[]{bLine[1], bLine[5], sj, transcript});
                                    p2hash.put(readId, list);
                                }
                            }
                        }
                    } else {
                        if (nBlocks == 2) {
                            String[] ids = getId(bLine, nBlocks);
                            for (String id : ids) {
                                nr = 1;
                                if (reads.containsKey(id)) {
                                    nr += reads.get(id);
                                }
                                reads.put(id, nr);
                            }
                        }
                    }
                    mapped = false;
                }
            }
        }
        if (paired) {
            for (String id : p1hash.keySet()) {
                if (p2hash.containsKey(id)) {
                    ArrayList<String[]> p1 = p1hash.get(id);
                    ArrayList<String[]> p2 = p2hash.get(id);
                    int p1Start, p2Start;
                    byte p1Strand, p2Strand;
                    String[] tx1, tx2;
                    for (String[] s1 : p1) {
                        for (String[] s2 : p2) {
                            tx1 = s1[3].split(",");
                            tx2 = s2[3].split(",");
                            Arrays.sort(tx1);
                            Arrays.sort(tx2);
                            p1Start = Integer.parseInt(s1[0]);
                            p2Start = Integer.parseInt(s2[0]);
                            p1Strand = s1[1].equals("+") ? (byte) 1 : (byte) -1;
                            p2Strand = s2[1].equals("+") ? (byte) 1 : (byte) -1;
                            if (p1Strand != p2Strand &&
                                    ((p1Start < p2Start && p1Strand == g.getStrand()) ||
                                            (p2Start < p1Start && p2Strand == g.getStrand()))) {

                                if (!Arrays.equals(tx1, tx2)) {
                                    for (String s : tx1) {
                                        for (String ss : tx2) {
                                            if (!(s.startsWith("(") && ss.startsWith("("))) {
                                                if (s.replaceAll("[()]", "").equals(ss.replaceAll("[()]", "")) && !s.startsWith("(") && s1[2] != null) {
                                                    String sj = s1[2];
                                                    nr = 1;
                                                    if (reads.containsKey(sj)) {
                                                        nr += reads.get(sj);
                                                    }
                                                    reads.put(sj, nr);
                                                }

                                                if (s.replaceAll("[()]", "").equals(ss.replaceAll("[()]", "")) && !ss.startsWith("(") && s2[2] != null) {
                                                    String sj = s2[2];
                                                    nr = 1;
                                                    if (reads.containsKey(sj)) {
                                                        nr += reads.get(sj);
                                                    }
                                                    reads.put(sj, nr);
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    if (s1[2] != null) {
                                        String sj = s1[2];
                                        nr = 1;
                                        if (reads.containsKey(sj)) {
                                            nr += reads.get(sj);
                                        }
                                        reads.put(sj, nr);
                                    }

                                    if (s2[2] != null) {
                                        String sj = s2[2];
                                        nr = 1;
                                        if (reads.containsKey(sj)) {
                                            nr += reads.get(sj);
                                        }
                                        reads.put(sj, nr);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return reads;
    }

    private Map<String, Integer> getAllIntronicReads(Gene g, boolean paired) throws Exception {
        BufferedReader bedReader = new BufferedReader(new InputStreamReader(new FileInputStream(bedFile)));
        Map<String, Integer> reads = new TreeMap<String, Integer>();
        HashMap<String, ArrayList<String[]>> p1hash = new HashMap<String, ArrayList<String[]>>();
        HashMap<String, ArrayList<String[]>> p2hash = new HashMap<String, ArrayList<String[]>>();
        Map<String, Integer> introns = new TreeMap<String, Integer>();
        int nr = 0;
        int start = 0, end = 0, tol = 0;
        start = g.getStart();
        end = g.getEnd();
        if (g.getStrand() < 0) {
            start = -start;
            end = -end;
        }
        tol = 0;
        start = Math.max(1, start - tol);
        end = end + tol;

        for (String tx : nodes.keySet()) {
            for (int i = 0; i < nodes.get(tx).size() - 1; i++) {
                String[] exon = nodes.get(tx).get(i);
                String[] nextExon = nodes.get(tx).get(i + 1);
                String intron = exon[1]+"-"+nextExon[0];
                if (!introns.containsKey(intron)) {
                    introns.put(intron, 1);
                } else {
                    introns.put(intron, introns.get(intron) + 1);
                }
            }
        }
        for (String intron : introns.keySet()) {
            if (introns.get(intron)==g.getTranscriptCount())
                continue;
            else {
                Iterator<String> iter = introns.keySet().iterator();
                String nextIntron = null;
                while (iter.hasNext()) {
                    nextIntron = iter.next();
                    if (Integer.parseInt(intron.split("-")[0])==Integer.parseInt(nextIntron.split("-")[0])) {
                        if (Integer.parseInt(intron.split("-")[1])<Integer.parseInt(nextIntron.split("-")[1])) {
                            introns.put(intron, introns.get(intron)+introns.get(nextIntron));
                        }
                    }
                    if (Integer.parseInt(intron.split("-")[1])==Integer.parseInt(nextIntron.split("-")[1])) {
                        if (Integer.parseInt(intron.split("-")[0])<Integer.parseInt(nextIntron.split("-")[0])) {
                            introns.put(nextIntron, introns.get(nextIntron)+introns.get(intron));
                        }
                    }
                }
            }
        }

        for (String line; (line = bedReader.readLine()) != null; ) {
            String[] bLine = line.split("\t");
            String[] mIntron = new String[2];
            boolean mapped = false;
            int nBlocks = Integer.parseInt(bLine[9]);
            if (bLine[0].equals(g.getChromosome()) && Integer.parseInt(bLine[1]) + 1 >= start && Integer.parseInt(bLine[2]) <= end) {
                for (String intronString : introns.keySet()) {
                    if (introns.get(intronString) == g.getTranscriptCount()) {
                        if (nBlocks == 1) {
                            String[] intron = intronString.split("-");
                            if (Integer.parseInt(intron[0]) < Integer.parseInt(bLine[1]) + 1 && Integer.parseInt(intron[1]) > Integer.parseInt(bLine[2])) {
                                mapped = true;
                                mIntron[0] = intron[0];
                                mIntron[1] = intron[1];
                                break;
                            }
                        }
                    }
                }
                if (mapped)// && g.getGeneID().equals(bLine[3].split(":")[0]+":"+ bLine[3].split(":")[1]))
                {
                    if (paired) {
                        if (nBlocks < 3) {
                            String readId = bLine[3].substring(0, bLine[3].length() - 4);
                            if (bLine[3].endsWith("1")) {
                                if (!p1hash.containsKey(readId)) {
                                    String sj = nBlocks == 2 ? getId(bLine, nBlocks)[0] : null;
                                    ArrayList<String[]> list = new ArrayList<String[]>();
                                    list.add(new String[]{bLine[3].substring(bLine[3].length() - 3, bLine[3].length() - 2), sj});
                                    p1hash.put(readId, list);
                                } else {
                                    String sj = nBlocks == 2 ? getId(bLine, nBlocks)[0] : null;
                                    ArrayList<String[]> list = p1hash.get(readId);
                                    list.add(new String[]{bLine[3].substring(bLine[3].length() - 3, bLine[3].length() - 2), sj});
                                    p1hash.put(readId, list);
                                }
                            }
                            if (bLine[3].endsWith("2")) {
                                if (!p2hash.containsKey(readId)) {
                                    String sj = nBlocks == 2 ? getId(bLine, nBlocks)[0] : null;
                                    ArrayList<String[]> list = new ArrayList<String[]>();
                                    list.add(new String[]{bLine[3].substring(bLine[3].length() - 3, bLine[3].length() - 2), sj});
                                    p2hash.put(readId, list);
                                } else {
                                    String sj = nBlocks == 2 ? getId(bLine, nBlocks)[0] : null;
                                    ArrayList<String[]> list = p2hash.get(readId);
                                    list.add(new String[]{bLine[3].substring(bLine[3].length() - 3, bLine[3].length() - 2), sj});
                                    p2hash.put(readId, list);
                                }
                            }
                        }
                    } else {
                        if (nBlocks == 1) {
                            String[] ids = {mIntron[0] + "-" + mIntron[1]};
                            for (String id : ids) {
                                nr = 1;
                                if (reads.containsKey(id)) {
                                    nr += reads.get(id);
                                }
                                reads.put(id, nr);
                            }
                        }
                    }
                    mapped = false;
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
                                if (s1[1] != null) {
                                    String sj = s1[1];
                                    nr = 1;
                                    if (reads.containsKey(sj)) {
                                        nr += reads.get(sj);
                                    }
                                    reads.put(sj, nr);
                                }
                                if (s2[1] != null) {
                                    String sj = s2[1];
                                    nr = 1;
                                    if (reads.containsKey(sj)) {
                                        nr += reads.get(sj);
                                    }
                                    reads.put(sj, nr);
                                }
                                //break;
                            }
                        }
                    }
                }
            }
        }
        return reads;
    }

    private void writeSJReads() throws Exception {
        GTFwrapper gtf = new GTFwrapper(gtfFile);
        BEDwrapper bed = new BEDwrapper(bedFile);
        try {
            gtf = new GTFwrapper((gtf.sort()));
            gtf.setReadAll(true);
            gtf.setNoIDs(null);
            gtf.setReadFeatures(new String[]{"exon", "CDS"});
            gtf.read();
            Gene g = gtf.getGenes()[0];
            int start = 0, end = 0, tol = 0;
            start = g.getStart();
            end = g.getEnd();
            if (g.getStrand() < 0) {
                start = -start;
                end = -end;
            }
            tol = 0;
            start = Math.max(1, start - tol);
            end = end + tol;
            BufferedIterator iter = bed.readBedFile(g, start, end, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
            AnnotationMapper a = new AnnotationMapper(g);
            a.map(iter, settings);
            ZipOutputStream out = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream("/home/emilio/ann_mapper.zip")));
            a.writeSJReads(out, true);
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
            Assert.fail();
        }
    }

    @Test
    public void testCompareSJReadsSingle() throws Exception {
        GTFwrapper gtf = new GTFwrapper(gtfFile);
        BEDwrapper bed = new BEDwrapper(bedFile);
        initSettings(UniversalReadDescriptor.DESCRIPTORID_SIMPLE, FluxCapacitorSettings.AnnotationMapping.SINGLE);
        //gtf = new GTFwrapper((gtf.sort()));
        //gtf.setReadAll(true);
        gtf.setChromosomeWise(true);
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon", "CDS"});
        gtf.sweepToChromosome("chr22");
        gtf.read();
        //Gene g = gtf.getGenes()[0];        ;
        for (Gene g : gtf.getGenes()) {
            //if (g.getGeneID().contains("18893736-18899601")) {
                int start = 0, end = 0, tol = 0;
                start = g.getStart();
                end = g.getEnd();
                if (g.getStrand() < 0) {
                    start = -start;
                    end = -end;
                }
                tol = 0;
                start = Math.max(1, start - tol);
                end = end + tol;
                BufferedIterator iter = bed.readBedFile(g, start, end, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
                AnnotationMapper a = new AnnotationMapper(g);
                a.map(iter, settings);
                Map<String, Integer> m = a.getSJReads(false);
                int count[] = new int[]{0, 0};
                for (String e : m.keySet()) {
                    count[0] += m.get(e);
                }
                readGtf(g);
                Map<String, Integer> m1 = getSJReads(g, false);
                for (String e : m1.keySet()) {
                    count[1] += m1.get(e);
                }
                if (count[0]!=count[1])
                    System.err.println("Gene : " + g.getGeneID() + "\tAnnotation Mapper: " + count[0] + "\tTest: " + count[1]);
                //assertEquals(count[1],count[0]);
            //}
        }
    }

    @Test
    public void testCompareSJReadsPaired() throws Exception {
        GTFwrapper gtf = new GTFwrapper(gtfFile);
        BEDwrapper bed = new BEDwrapper(bedFile);
        initSettings(UniversalReadDescriptor.DESCRIPTORID_CASAVA18, FluxCapacitorSettings.AnnotationMapping.PAIRED);
        //gtf = new GTFwrapper((gtf.sort()));
        //gtf.setReadAll(true);
        gtf.setChromosomeWise(true);
        gtf.sweepToChromosome("chr22");
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon", "CDS"});
        gtf.read();
        for (Gene g : gtf.getGenes()) {
            //if (g.getGeneID().contains("50946645-50963209")) {
            int start = 0, end = 0, tol = 0;
            start = g.getStart();
            end = g.getEnd();
            if (g.getStrand() < 0) {
                start = -start;
                end = -end;
            }
            tol = 0;
            start = Math.max(1, start - tol);
            end = end + tol;
            BufferedIterator iter = bed.readBedFile(g, start, end, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
            AnnotationMapper a = new AnnotationMapper(g);
            a.map(iter, settings);
            Map<String, Integer> m = a.getSJReads(true);
            int count[] = new int[]{0, 0};
            for (String e : m.keySet()) {
                count[0] += m.get(e);
            }
            readGtf(g);
            Map<String, Integer> m1 = getSJReads(g, true);
            for (String e : m1.keySet()) {
                count[1] += m1.get(e);
            }
            if (count[0] != count[1])
                System.err.println("Gene : " + g.getGeneID() + "\tAnnotation Mapper: " + count[0] + "\tTest: " + count[1]);
            //assertEquals(count[1], count[0]);
            //}
        }
    }

    @Test
    public void testCompareIntronReadsSingle() throws Exception {
        GTFwrapper gtf = new GTFwrapper(gtfFile);
        BEDwrapper bed = new BEDwrapper(bedFile);
        initSettings(UniversalReadDescriptor.DESCRIPTORID_SIMPLE, FluxCapacitorSettings.AnnotationMapping.SINGLE);
        //gtf = new GTFwrapper((gtf.sort()));
        //gtf.setReadAll(true);
        gtf.setChromosomeWise(true);
        gtf.sweepToChromosome("chr22");
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon", "CDS"});
        gtf.read();
        for (Gene g : gtf.getGenes()) {
            if (g.getGeneID().contains("19744226-19771112")) {//32870707-32894818
            int start = 0, end = 0, tol = 0;
            start = g.getStart();
            end = g.getEnd();
            if (g.getStrand() < 0) {
                start = -start;
                end = -end;
            }
            tol = 0;
            start = Math.max(1, start - tol);
            end = end + tol;
            BufferedIterator iter = bed.readBedFile(g, start, end, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
            AnnotationMapper a = new AnnotationMapper(g);
            a.map(iter, settings);
            Map<String, Integer[]> m = a.getAllIntronicReads();
            int count[] = new int[]{0, 0};
            for (String e : m.keySet()) {
                count[0] += m.get(e)[0];
            }
            readGtf(g);
            Map<String, Integer> m1 = getAllIntronicReads(g, false);
            for (String e : m1.keySet()) {
                count[1] += m1.get(e);
            }
            if (count[0]!=count[1])
                System.err.println("Gene : " + g.getGeneID() + "\tAnnotationMapper: " + count[0] + "\tTest: " + count[1]);
            //assertTrue(count[0] > 0);
            }
        }
    }

}
