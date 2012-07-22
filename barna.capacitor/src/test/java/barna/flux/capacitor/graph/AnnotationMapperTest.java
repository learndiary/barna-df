package barna.flux.capacitor.graph;

import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.io.MSIterator;
import barna.io.bed.BEDReader;
import barna.io.gtf.GTFwrapper;
import barna.io.rna.UniversalReadDescriptor;
import barna.io.sam.SAMReader;
import barna.model.Gene;
import barna.model.Mapping;
import barna.model.bed.BEDMapping;
import barna.model.sam.SAMMapping;
import junit.framework.TestCase;
import org.junit.Test;

import java.io.*;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Emilio Palumbo
 * Date: 5/31/12
 * Time: 12:39 PM
 */
public class AnnotationMapperTest extends TestCase {

    private final File hgGtfFile = new File(getClass().getResource("/gencode_v12_hg_chr22_24030323-24041363.gtf").getFile());
    //private final File hgBedFile = new File(getClass().getResource("/test_hg_chr22_24030323-24041363.bed").getFile());
    private final File hgBedFile = new File("/home/emilio/fromMicha/test-chr22-24030323-24041363_new.bed");
    private final File hgBamIndexFile = new File("/home/emilio/fromMicha/test.bam");
    private final File hgBamFile = new File("/home/emilio/fromMicha/test1.bam");
    private final File mm9GtfFile = new File(getClass().getResource("/mm9_chr1_chrX.gtf").getFile());
    private final File mm9BedFile = new File(getClass().getResource("/chr1_chrX.bed").getFile());
    private FluxCapacitorSettings settings;
    Map<String, ArrayList<String[]>> nodes = new HashMap<String, ArrayList<String[]>>();

    @Override
    public void setUp() throws Exception {
        super.setUp();
        initSettings(UniversalReadDescriptor.DESCRIPTORID_CASAVA18, FluxCapacitorSettings.AnnotationMapping.PAIRED);
    }

    private void initSettings(String descriptorStr, FluxCapacitorSettings.AnnotationMapping mapping) {
        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor(descriptorStr));
        settings = new FluxCapacitorSettings();
        settings.set(FluxCapacitorSettings.ANNOTATION_FILE,
                new File(hgGtfFile.getAbsolutePath()));
        settings.set(FluxCapacitorSettings.MAPPING_FILE,
                new File(hgBedFile.getAbsolutePath()));
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

    private void readGtf(Gene g, File gtfFile) throws IOException {
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

    private Map<String, Integer> getSJReads(Gene g, boolean paired, File bedFile) throws Exception {
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
                        if (nBlocks > 1 && i < nodes.get(tx).size() - 1) {
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
                }
                if (mapped) {
                    if (paired) {
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
                    } else {
                        if (nBlocks > 1) {
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

    private Map<String, Integer> getAllIntronicReads(Gene g, boolean paired, File bedFile) throws Exception {
        BufferedReader bedReader = new BufferedReader(new InputStreamReader(new FileInputStream(bedFile)));
        Map<String, Integer> reads = new TreeMap<String, Integer>();
        HashMap<String, ArrayList<String[]>> p1hash = new HashMap<String, ArrayList<String[]>>();
        HashMap<String, ArrayList<String[]>> p2hash = new HashMap<String, ArrayList<String[]>>();
        Map<String, Integer> introns = new TreeMap<String, Integer>();
        TreeSet<String> exons = new TreeSet<String>();
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
            for (int i = 0; i < nodes.get(tx).size(); i++) {
                if (i < nodes.get(tx).size() - 1) {
                    String[] exon1 = nodes.get(tx).get(i);
                    String[] exon2 = nodes.get(tx).get(i + 1);
                    String intron = exon1[1] + "-" + exon2[0];

                    if (!introns.containsKey(intron)) {
                        introns.put(intron, 1);
                    } else {
                        introns.put(intron, introns.get(intron) + 1);
                    }
                }

                String exon = nodes.get(tx).get(i)[0] + "-" + nodes.get(tx).get(i)[1];
                if (!exons.contains(exon)) {
                    exons.add(exon);
                }


            }
        }

        for (String line; (line = bedReader.readLine()) != null; ) {
            String[] bLine = line.split("\t");
            String[] mIntron = new String[2];
            boolean mapped = true;
            int nBlocks = Integer.parseInt(bLine[9]);
            if (bLine[0].equals(g.getChromosome()) && Integer.parseInt(bLine[1]) + 1 >= start && Integer.parseInt(bLine[2]) <= end && nBlocks == 1) {
                for (String exonString : exons) {
                    String[] exon = exonString.split("-");
                    if (Integer.parseInt(exon[0]) < Integer.parseInt(bLine[2]) && Integer.parseInt(exon[1]) > Integer.parseInt(bLine[1]) + 1) {
                        mapped = false;
                        break;
                    }
                }
                if (mapped) {
                    for (String intronString : introns.keySet()) {
                        String[] intron = intronString.split("-");
                        if (Integer.parseInt(intron[0]) < Integer.parseInt(bLine[1]) + 1 && Integer.parseInt(intron[1]) > Integer.parseInt(bLine[2])) {
                            mIntron[0] = intron[0];
                            mIntron[1] = intron[1];
                        }
                    }
                    if (mIntron[0] != null && mIntron[1] != null) {
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
                }
            }
        }

        return reads;
    }

    @Test
    public void testCompareSJReadsSingle() throws Exception {
        GTFwrapper gtf = new GTFwrapper(hgGtfFile);
        BEDReader bed = new BEDReader(hgBedFile, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
        byte lastStr = 0;
        initSettings(UniversalReadDescriptor.DESCRIPTORID_SIMPLE, FluxCapacitorSettings.AnnotationMapping.SINGLE);
        //gtf = new GTFwrapper((gtf.sort()));
        gtf.setReadAll(true);
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon", "CDS"});
        gtf.read();
        for (Gene g : gtf.getGenes()) {
            if (lastStr!=0&&lastStr!=g.getStrand()) {
                bed.reset(g.getChromosome());
                lastStr = g.getStrand();
            }
            if (lastStr == 0)
                lastStr = g.getStrand();
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
            MSIterator<Mapping> iter = bed.read(g.getChromosome(), start, end);
            AnnotationMapper a = new AnnotationMapper(g);
            a.map(iter, settings);
            Map<String, Integer> m = a.getSJReads(false);
            int count[] = new int[]{0, 0};
            for (String e : m.keySet()) {
                count[0] += m.get(e);
            }
            readGtf(g, hgGtfFile);
            Map<String, Integer> m1 = getSJReads(g, false, hgBedFile);
            for (String e : m1.keySet()) {
                count[1] += m1.get(e);
            }
            assertEquals(count[0], count[1]);
        }
    }

    @Test
    public void testCompareSJReadsPaired() throws Exception {
        initSettings(UniversalReadDescriptor.DESCRIPTORID_PAIRED, FluxCapacitorSettings.AnnotationMapping.PAIRED);
        GTFwrapper gtf = new GTFwrapper(hgGtfFile);
        BEDReader bed = new BEDReader(hgBedFile, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
        byte lastStr = 0;
        //gtf = new GTFwrapper((gtf.sort()));
        gtf.setReadAll(true);
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon", "CDS"});
        gtf.read();
        for (Gene g : gtf.getGenes()) {
            if (lastStr!=0&&lastStr!=g.getStrand()) {
                bed.reset(g.getChromosome());
                lastStr = g.getStrand();
            }
            if (lastStr == 0)
                lastStr = g.getStrand();
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
            MSIterator<Mapping> iter = bed.read(g.getChromosome(),start,end);
//            MSIterator<Mapping> iter = bed.readBedFile(g, start, end, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
            AnnotationMapper a = new AnnotationMapper(g);
            a.map(iter, settings);
            Map<String, Integer> m = a.getSJReads(true);
            int count[] = new int[]{0, 0};
            for (String e : m.keySet()) {
                count[0] += m.get(e);
            }
//            readGtf(g, hgGtfFile);
//            Map<String, Integer> m1 = getSJReads(g, true, hgBedFile);
//            for (String e : m1.keySet()) {
//                count[1] += m1.get(e);
//            }
            assertEquals(count[0], 4);
        }
    }

//    @Test
//    public void testSJReadsPaired() throws Exception {
//        GTFwrapper gtf = new GTFwrapper(hgGtfChrFile);
//        BEDwrapper bed = new BEDwrapper(hgBedChrFile);
//        Map<String,String> sjs = new TreeMap<String, String>();
//        int count=0;
//        byte lastStr = 0;
//        initSettings(UniversalReadDescriptor.DESCRIPTORID_CASAVA18, FluxCapacitorSettings.AnnotationMapping.PAIRED);
//        //gtf = new GTFwrapper((gtf.sort()));
//        gtf.setReadAll(true);
//        gtf.setNoIDs(null);
//        gtf.setReadFeatures(new String[]{"exon", "CDS"});
//        gtf.read();
//        for (Gene g : gtf.getGenes()) {
//            if (lastStr!=0&&lastStr!=g.getStrand()) {
//                bed.reset(g.getChromosome());
//                lastStr = g.getStrand();
//            }
//            if (lastStr == 0)
//                lastStr = g.getStrand();
//            int start = 0, end = 0, tol = 0;
//            start = g.getStart();
//            end = g.getEnd();
//            if (g.getStrand() < 0) {
//                start = -start;
//                end = -end;
//            }
//            tol = 0;
//            start = Math.max(1, start - tol);
//            end = end + tol;
//            BufferedIterator iter = bed.readBedFile(g, start, end, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
//            AnnotationMapper a = new AnnotationMapper(g);
//            a.map(iter, settings);
//            Map<String, Integer> m = a.getSJReads(true);
//            for (String e : m.keySet()) {
//                sjs.put(e,"OK");
//            }
//        }
//
//        System.err.println("Splice junctions found: " + sjs.keySet().size());
//    }

    @Test
    public void testCompareMultiSJReadsSingle() throws Exception {
        GTFwrapper gtf = new GTFwrapper(mm9GtfFile);
        BEDReader bed = new BEDReader(mm9BedFile, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
        byte lastStr = 0;
        initSettings(UniversalReadDescriptor.DESCRIPTORID_SIMPLE, FluxCapacitorSettings.AnnotationMapping.SINGLE);
        //gtf = new GTFwrapper((gtf.sort()));
        gtf.setReadAll(true);
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon", "CDS"});
        gtf.read();
        for (Gene g : gtf.getGenes()) {
            if (lastStr!=0&&lastStr!=g.getStrand()) {
                bed.reset(g.getChromosome());
                lastStr = g.getStrand();
            }
            if (lastStr == 0)
                lastStr = g.getStrand();
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
            MSIterator<Mapping> iter = bed.read(g.getChromosome(), start, end);
            AnnotationMapper a = new AnnotationMapper(g);
            a.map(iter, settings);
            Map<String, Integer> m = a.getSJReads(false);
            int count[] = new int[]{0, 0};
            for (String e : m.keySet()) {
                count[0] += m.get(e);
            }
            readGtf(g, mm9GtfFile);
            Map<String, Integer> m1 = getSJReads(g, false, mm9BedFile);
            for (String e : m1.keySet()) {
                count[1] += m1.get(e);
            }
            //assertEquals(count[0], count[1]);
            System.err.println(count[0] + "\t" + count[1]);
        }
    }

    @Test
    public void testCompareIntronReadsSingle() throws Exception {
        GTFwrapper gtf = new GTFwrapper(hgGtfFile);
        BEDReader bed = new BEDReader(hgBedFile, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
        byte lastStr = 0;
        initSettings(UniversalReadDescriptor.DESCRIPTORID_SIMPLE, FluxCapacitorSettings.AnnotationMapping.SINGLE);
        //gtf = new GTFwrapper((gtf.sort()));
        gtf.setReadAll(true);
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon", "CDS"});
        gtf.read();
        for (Gene g : gtf.getGenes()) {
            if (lastStr!=0&&lastStr!=g.getStrand()) {
                bed.reset(g.getChromosome());
                lastStr = g.getStrand();
            }
            if (lastStr == 0)
                lastStr = g.getStrand();
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
            MSIterator<Mapping> iter = bed.read(g.getChromosome(), start, end);
            AnnotationMapper a = new AnnotationMapper(g);
            a.map(iter, settings);
            Map<String, Float[]> m = a.getAllIntronicReads(false);
            int count[] = new int[]{0, 0};
            for (String e : m.keySet()) {
                count[0] += m.get(e)[0];
            }
            readGtf(g, hgGtfFile);
            Map<String, Integer> m1 = getAllIntronicReads(g, false,hgBedFile);
            for (String e : m1.keySet()) {
                count[1] += m1.get(e);
            }
            assertEquals(count[0], count[1]);
        }
    }

    @Test
    public void testCompareIntronReadsPaired() throws Exception { //TODO update getAllIntronic reads to work for paired end reads
        initSettings(UniversalReadDescriptor.DESCRIPTORID_PAIRED, FluxCapacitorSettings.AnnotationMapping.PAIRED);
        GTFwrapper gtf = new GTFwrapper(hgGtfFile);
        BEDReader bed = new BEDReader(hgBedFile, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
        byte lastStr = 0;
        //gtf = new GTFwrapper((gtf.sort()));
        gtf.setReadAll(true);
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon", "CDS"});
        gtf.read();
        for (Gene g : gtf.getGenes()) {
            if (lastStr!=0&&lastStr!=g.getStrand()) {
                bed.reset(g.getChromosome());
                lastStr = g.getStrand();
            }
            if (lastStr == 0)
                lastStr = g.getStrand();
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
            MSIterator<Mapping> iter = bed.read(g.getChromosome(), start, end);
            AnnotationMapper a = new AnnotationMapper(g);
            a.map(iter, settings);
            Map<String, Float[]> m = a.getAllIntronicReads(true);
            int count[] = new int[]{0, 0};
            for (String e : m.keySet()) {
                count[0] += m.get(e)[0];
            }
            readGtf(g,hgGtfFile);
            Map<String, Integer> m1 = getAllIntronicReads(g, false, hgBedFile);
            for (String e : m1.keySet()) {
                count[1] += m1.get(e);
            }
            assertEquals(76, count[0]);
        }
    }

    @Test
    public void testCompareBEDtoBAM() throws Exception {
        GTFwrapper gtf = new GTFwrapper(hgGtfFile);
        SAMReader sam = new SAMReader(hgBamFile);
        BEDReader bed = new BEDReader(hgBedFile, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR),null);
        initSettings(UniversalReadDescriptor.DESCRIPTORID_SIMPLE, FluxCapacitorSettings.AnnotationMapping.SINGLE);
        gtf.setReadAll(true);
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon", "CDS"});
        gtf.read();
        for (Gene g : gtf.getGenes()) {
            int start = 0, end = 0;
            start = g.getStart();
            end = g.getEnd();
            if (g.getStrand() < 0) {
                start = -start;
                end = -end;
            }
            start = Math.max(1, start);

            MSIterator<Mapping> iter1 = bed.read(g.getChromosome(),start,end);
            MSIterator<Mapping> iter2 = sam.read(g.getChromosome(),start,end);

            int[] count = {0,0};

            while (iter1.hasNext()) {
                iter1.next();
                count[0]++;
            }

            while (iter2.hasNext()) {
                iter2.next();
                count[1]++;
            }

//            Mapping m1,m2;
//            while (iter1.hasNext()) {
//                m1=iter1.next();
//                m2=iter2.next();
//                assertEquals(m1.getChromosome(),m2.getChromosome());
//                assertEquals(m1.getName(),m2.getName());
//                assertEquals(m1.getStart(),m2.getStart());
//                assertEquals(m1.getEnd(),m2.getEnd());
//                assertEquals(m1.getStrand(),m2.getStrand());
//            }

            assertEquals(count[0],count[1]);
        }

    }

    @Test
    public void testCompareBAMSJReadsSingle() throws Exception {
        initSettings(UniversalReadDescriptor.DESCRIPTORID_PAIRED, FluxCapacitorSettings.AnnotationMapping.PAIRED);
        GTFwrapper gtf = new GTFwrapper(hgGtfFile);
        SAMReader sam = new SAMReader(hgBamFile);
        SAMReader samIndex = new SAMReader(hgBamIndexFile);
        BEDReader bed = new BEDReader(hgBedFile, true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR), null);
        byte lastStr = 0;
        //gtf = new GTFwrapper((gtf.sort()));
        gtf.setReadAll(true);
        gtf.setNoIDs(null);
        gtf.setReadFeatures(new String[]{"exon", "CDS"});
        gtf.read();
        for (Gene g : gtf.getGenes()) {
            if (lastStr!=0&&lastStr!=g.getStrand()) {
                bed.reset(g.getChromosome());
                lastStr = g.getStrand();
            }
            if (lastStr == 0)
                lastStr = g.getStrand();
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
            MSIterator<Mapping> iter1 = bed.read(g.getChromosome(), start, end);
            MSIterator<Mapping> iter2 = sam.read(g.getChromosome(), start, end);
            MSIterator<Mapping> iter3 = samIndex.read(g.getChromosome(), start, end);

            AnnotationMapper a = new AnnotationMapper(g);
            AnnotationMapper b = new AnnotationMapper(g);
            a.map(iter1, settings);
            b.map(iter3, settings);

//            assertEquals(a.nrMappingsLocus,b.nrMappingsLocus);
            assertEquals(a.getNrMappingsMapped(),b.getNrMappingsMapped());
//            assertEquals(a.nrMappingsNotMapped,b.nrMappingsNotMapped);

            Map<String, Integer> m = a.getSJReads(true);
            Map<String, Integer> m1 = b.getSJReads(true);
            int count[] = new int[]{0, 0};
            for (String e : m.keySet()) {
                count[0] += m.get(e);
            }
//            readGtf(g, hgGtfFile);
            //Map<String, Integer> m1 = getSJReads(g, false, hgBedFile);
            for (String e : m1.keySet()) {
                count[1] += m1.get(e);
            }
            assertEquals(count[0], count[1]);
        }
    }

}
