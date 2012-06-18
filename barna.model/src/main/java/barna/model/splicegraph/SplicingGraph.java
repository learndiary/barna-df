/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.model.splicegraph;

import barna.commons.utils.ArrayUtils;
import barna.model.*;
import barna.model.commons.IntVector;
import barna.model.commons.MyFile;
import barna.model.commons.MyMath;
import barna.model.commons.MyTime;
import barna.model.constants.Constants;
import barna.model.gff.GFFObject;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.zip.GZIPOutputStream;

/**
 * <tt>07/11/14</tt>: handles soft starts
 * <tt>07/11/22</tt>: writes out stats as comments in GTF
 * <tt>08/03/08</tt>: edgeConfidenceLevel delegated to <code>model.Transcript</code>
 * constructGraph() now takes <= for trusting transcript end
 * <tt>08/03/09</tt>: changed strategy how to connect soft starts/ends to root/leaf
 * constructGraph() saved in constructGraph_x3()
 * generateTuples adapted to correct for ss extraction acc. to graph and not in annotation
 * #events increased, investigate
 * <tt>08/08/18</tt>: major bugfix in createEdge(), concerning checking of valid introns
 * (acceptableIntrons, intronConfidenceLevel)
 * <tt>08/08/26</tt>: major bugfix in constructGraph(), single exon genes are no longer extended
 * due to inconsistencies with s.getTranscripts() later-on,
 * some bugfixes in the ASEvent attributes and output
 */

public class SplicingGraph {
    public final static String TRANSCRIPT_ID_TAG = "transcript_id";

    class MyLittleInt implements Comparable {
        int val;

        public MyLittleInt(int newVal) {
            this.val = newVal;
        }

        @Override
        public int hashCode() {
            return val;
        }

        @Override
        public boolean equals(Object obj) {
            return (compareTo(obj) == 0);
        }

        public int compareTo(Object o) {
            MyLittleInt other = (MyLittleInt) o;
            return (val - other.val);
        }

    }

    public class EdgeCoordComparator implements Comparator<SimpleEdge> {
        int readLen = 36;

        public EdgeCoordComparator(int readLen) {
            this.readLen = readLen;
        }

        public int compare(SimpleEdge o1, SimpleEdge o2) {
            Transcript[] t = decodeTset(intersect(o1.getTranscripts(), o2.getTranscripts()));
//			if (o1.toString().equals("28905362-28905391^28905391^28905479-"))
//				System.currentTimeMillis();
            int[] coord1 = o1.getFrac(t[0], readLen), coord2 = o2.getFrac(t[0], readLen);
            if (coord1[0] < coord2[0]) {
//				if (coord1[1]> coord2[1])
//					System.currentTimeMillis();
                return -1;
            }
            if (coord1[0] > coord2[0]) {
//				if (coord1[1]< coord2[1])
//					System.currentTimeMillis();
                return 1;
            }
            return 0;
        }
    }

    public EdgeCoordComparator defaultEdgeCoordComparator = null;

    public void createDefaultCoordComparator(int readLen) {
        defaultEdgeCoordComparator = new EdgeCoordComparator(readLen);
    }


    public static boolean outputSeq = false;
    public static final String version = "2.1", build = "080826";

    HashMap<SpliceSite, SpliceSite> softStartHash = new HashMap<SpliceSite, SpliceSite>(), softEndHash = new HashMap<SpliceSite, SpliceSite>();
    public static int counter = 0;
    static int readAheadLimit = -1;
    /**
     *
     */
    public static WriterThread writerThread = null;    // let it null here, see addEvent()
    public static boolean writeStdOut = false;
    public static boolean retrieveASEvents = true;
    public static boolean retrieveDSEvents = false;
    public static boolean retrieveVSEvents = false;

    public static final byte ETYPE_AL = 0, ETYPE_EX = 1, ETYPE_IN = 2, ETYPE_SJ = 3, ETYPE_XJ = 4, ETYPE_PE = 5;

    public static boolean canonicalSS = false;
    public static boolean acceptableIntrons = false; // schmu-buh, immer true sobald intronConfidence gesetzt
    public static boolean DEBUG = false;
    public static byte intronConfidenceLevel = Transcript.ID_SRC_MOST_INCONFIDENT;
    public static boolean onlyInternal = true;
    static long invalidIntrons = 0, totalIntrons = 0;

    public static class EventExtractorThread extends Thread {

        Gene[] g;
        Thread upstreamThread;
        public static int n = 2;
        boolean output = false, output2 = true;
        public static Species species = new Species("human");

        static {
            species.setGenomeVersion("hg18");
        }

        public static void setSpecies(Species newSpecies) {
            species = newSpecies;
        }

        public EventExtractorThread(Thread newUpstreamThread) {
            super();
            setName("event_extraction_thread");
            upstreamThread = newUpstreamThread;
        }

        @Override
        public void run() {

            String chromo = null;
            int evBefore = SplicingGraph.counter;
            long cumulGC = 0l, cumulGF = 0l, cumulGT = 0l, cumulEV = 0l;
            for (int i = 0; i < g.length; i++) {

                if (chromo == null)
                    chromo = g[0].getChromosome();

                if (g[i].getTranscriptCount() == 1) {
                    g[i] = null;
                    continue;
                }
//						if (g[i].getTranscriptCount()> 5000)  {
//							if (output2) {
//								Date ti= new Date(System.currentTimeMillis());
//								System.err.println("["+ti+"] "+chromo+" skipped locus "+g[i].getTranscripts()[0].getTranscriptID());
//							}
//							continue;
//						}

                g[i].setSpecies(species);

                String tmp = null;
                try {
                    long t1 = System.currentTimeMillis();

                    // DEBUG
//							Transcript[] trpts= g[i].getTranscripts();
//							Arrays.sort(trpts, defaultTranscriptByNameComparator);
//							if (!trpts[0].getTranscriptID().equals("AA002101"))
//								continue;

                    SplicingGraph gr = new SplicingGraph(g[i]);
                    if (output) {
                        System.err.println(gr.trpts[0].getTranscriptID() + "  transcripts " + g[i].getTranscriptCount());
                    }
                    //gr.init(n);
                    gr.constructGraph();
                    Node[] no = gr.getNodesInGenomicOrder();
                    long t2 = System.currentTimeMillis();
                    long dA = (t2 - t1);
                    cumulGC += dA;
                    if (output) {
                        tmp = gr.trpts[0].getTranscriptID() + " construct " + (dA / 1000) + " n " + gr.nodeHash.size() + " e " + gr.edgeHash.size();
                        System.err.println(tmp);
                        System.err.flush();
                    }

                    // in again
                    // DEBUG: take out fuzzy flanks
                    gr.collapseFuzzyFlanks(); // formerly called in init()
                    //gr.cleanESTborderEdges();	// for Claudia, kill events linked to EST edge variations

                    long t22 = System.currentTimeMillis();
                    long dF = (t22 - t2);
                    cumulGF += dF;
                    if (output) {
                        tmp = gr.trpts[0].getTranscriptID() + " flanks " + (dF / 1000) + " n " + gr.nodeHash.size() + " e " + gr.edgeHash.size();
                        System.err.println(tmp);
                        System.err.flush();
                    }

                    //gr.contractGraph(2);	// not n !!!
                    no = gr.getNodesInGenomicOrder();
                    long t3 = System.currentTimeMillis();
                    long dT = (t3 - t22);
                    cumulGT += dT;
                    if (output) {
                        tmp = gr.trpts[0].getTranscriptID() + " contract " + (dT / 1000) + " n " + gr.nodeHash.size() + " e " + gr.edgeHash.size();
                        System.err.println(tmp);
                        System.err.flush();
                    }
                    //int oldSize= evMap.size();

                    gr.getEventsByPartitions(n);

                    //cnt+= evMap.size()- oldSize;
                    long dB = (System.currentTimeMillis() - t3);
                    cumulEV += dB;
//							if (output) 
//								//System.err.println("events "+(evMap.size()-oldSize));
//								System.err.println(" extract "+(dB/1000));
                } catch (OutOfMemoryError e) {
                    MyTime ti = new MyTime(System.currentTimeMillis());
                    System.err.println("[" + ti + "]");
                    e.printStackTrace();
                    System.err.println(g[i].getTranscripts()[0].getTranscriptID() + " " + g[i].getTranscriptCount() + " " + g[i].getSpliceSites().length);
                    if (tmp != null)
                        System.err.println(tmp);
                } finally {
                    g[i] = null;
                    writerThread.interrupt();
                }


            }    // for

            if (output2) {
                Date ti = new Date(System.currentTimeMillis());
                int div = (int) (cumulGC + cumulGF + cumulEV) / 1000;
                int frac = 0;
                if (div > 0)
                    frac = (SplicingGraph.counter - evBefore) / div;
                else
                    frac = (SplicingGraph.counter - evBefore);

                System.err.println("[" + ti + "] " + chromo +
                        " graph construct " + (cumulGC / 1000) +
                        //" sec, fuzzy flanks "+(cumulGF/1000) +
                        " sec, contraction " + (cumulGT / 1000) +
                        " sec, extract events " + (cumulEV / 1000) +
                        " sec, found " + (SplicingGraph.counter - evBefore) +
                        " events, " + frac + " ev/sec.");
            }
            System.gc();
            writerThread.interrupt();

        }

        public boolean isOutput() {
            return output;
        }

        public void setOutput(boolean output) {
            this.output = output;
        }

        public boolean isOutput2() {
            return output2;
        }

        public void setOutput2(boolean output2) {
            this.output2 = output2;
        }

        public Gene[] getG() {
            return g;
        }

        public void setG(Gene[] g) {
            this.g = g;
        }
    }

    public static class WriterThread extends Thread {
        Queue<ASEvent> queue = new ConcurrentLinkedQueue<ASEvent>();
        Thread writingThread;
        boolean kill = false, outputASTA = false, outputGTF = true;
        int maxQevents = 10000, minQevents = 1000;
        static boolean writeHeader = true;
        public String outputFname;


        public WriterThread() {
            super();
            setName("event_writing_thread");
        }

        public void addEvent(ASEvent ev) {
            queue.add(ev);
            writingThread.interrupt();
        }

        public void setKill(boolean newKill) {
            kill = newKill;
        }


        public void run() {
            writingThread = Thread.currentThread();

            //BufferedWriter buffy= new BufferedWriter(new OutputStreamWriter(
            OutputStreamWriter globalWriter = null;

            try {
                if (outputFname != null) {
                    MyFile f = new MyFile(outputFname);
                    if (!f.exists())
                        writeHeader = true;
                    GZIPOutputStream zipperStream = new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(f, false)));    // WAS true, for appending
                    globalWriter = new OutputStreamWriter(zipperStream);
                } else
                    globalWriter = new OutputStreamWriter(System.out);

                //				if (writeHeader) {
                //					//Writer writer= new OutputStreamWriter(zipperStream);
                //					//outputStats(writer);
                //					outputStats(globalWriter);
                //					writeHeader= false;
                //					//writer.close();
                //				}
            } catch (IOException e) {
                e.printStackTrace();
            }

            while (true) {
                try {
                    while (!queue.isEmpty()) {
                        ASEvent ev = queue.poll();
                        ev.init();
                        String s = null;
                        if (outputGTF)
                            s = ev.toStringGTF(outputSeq);    // true
                        if (outputASTA)
                            s = ev.toStringASTA();
                        //for (int i = 0; i < s.length(); i++) 
                        //zipperStream.write((byte) s.charAt(i));
                        globalWriter.write(s.toCharArray(), 0, s.length());
                        //zipperStream.write('\n');
                        globalWriter.write('\n');
                    }
                    //zipperStream.flush(); zipperStream.close();

                    // better compression if not flushed?
//					globalWriter.flush();

                    // close only at the end
//					if (outputFname!= null)
//						globalWriter.close();

                    //buffy.flush(); buffy.close();

                } catch (IOException e) {
                    e.printStackTrace();
                }

                if (kill) {
                    try {
                        globalWriter.flush();
                        if (outputFname != null)
                            globalWriter.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    return;
                }

                try {
                    writingThread.sleep(0);
                } catch (InterruptedException e) {
                    ; // :)
                }
            }
        }

        public int getMaxQevents() {
            return maxQevents;
        }

        public void setMaxQevents(int maxQevents) {
            this.maxQevents = maxQevents;
        }

        public int getMinQevents() {
            return minQevents;
        }

        public void setMinQevents(int minQevents) {
            this.minQevents = minQevents;
        }
    }


    public static void main(String[] args) {
        //_070808_test();
        //gphase.Constants.DATA_DIR= "/home/msammeth";
        //System.out.println(gphase.Constants.DATA_DIR);


        //_240808_test_multithread(args);

        /*IntronModel iModel= new IntronModel();
		iModel.read(new File("test.model"));
		genome.model.Graph.overrideSequenceDirPath= "c:\\genomes\\H.sapiens\\golden_path_200603\\chromFa";
		extractSpliceJunctions(30, 30, iModel, 
				new File("M:\\annotations\\hg18_all-mrna_UCSC090507.gtf"), 
				new File("C:\\testJunctions.fasta"));
	    */
    }


    public boolean isRoot(Node v) {
        if (v.equals(root))
            return true;
        return false;
    }

    public boolean isLeaf(Node v) {
        if (v.equals(leaf))
            return true;
        return false;
    }

    public static class TranscriptByNameComparator implements Comparator<Transcript> {
        public int compare(Transcript t1, Transcript t2) {
            return t1.getTranscriptID().compareTo(t2.getTranscriptID());
        }
    }

    public static boolean isNull(long[] sig) {
        for (int i = 0; i < sig.length; i++)
            if (sig[i] != 0l)
                return false;
        return true;
    }

    public static long[] intersect(long[] sig1, long[] sig2) {
        assert (sig1.length == sig2.length);
        long[] reSig = new long[sig1.length];
        for (int i = 0; i < sig2.length; i++)
            reSig[i] = sig1[i] & sig2[i];
        return reSig;
    }

    public static boolean intersects(long[] a, long[] b) {
        for (int i = 0; i < b.length; i++) {
            if ((a[i] & b[i]) != 0l)
                return true;
        }
        return false;
    }

    public static long[] unite(long[] sig1, long[] sig2) {
        assert (sig1.length == sig2.length);
        long[] reSig = new long[sig1.length];
        for (int i = 0; i < sig2.length; i++)
            reSig[i] = sig1[i] | sig2[i];
        return reSig;
    }

    public static long[] xor(long[] a, long[] b) {
        long[] c = new long[a.length];
        for (int i = 0; i < c.length; i++)
            c[i] = a[i] ^ b[i];
        return c;
    }

    public static long[] without(long[] a, long[] b) {
        long[] inter = intersect(a, b);
        long[] c = new long[a.length];
        for (int i = 0; i < c.length; i++) {
            c[i] = a[i] ^ inter[i];
        }
        return c;
    }


    public static long[] unite(Vector<long[]> a) {
        if (a.size() == 0)
            return null;
        if (a.size() == 1)
            return a.elementAt(0);
        long[] inter = unite(a.elementAt(0), a.elementAt(1));
        for (int i = 2; i < a.size(); i++)
            inter = unite(inter, a.elementAt(i));
        return inter;
    }

    public static long[] createNullArray(int size) {
        long[] a = new long[size];
        for (int i = 0; i < a.length; i++)
            a[i] = 0l;
        return a;
    }

    public long[] createAllArray() {
        long[] a = encodeTset(trpts);
        return a;
    }

    public static int intersect(Vector<long[]> ps1, Vector<long[]> ps2, Vector<long[]> interSet, boolean onlyFullIntersect) {
        int cntFullIntersect = 0;
        for (int k = 0; k < ps1.size(); k++) {
            for (int h = 0; h < ps2.size(); h++) {
                long[] sect = intersect(ps1.elementAt(k), ps2.elementAt(h));
                if (isNull(sect))
                    continue;
                int fis = 0;
                if (equalSet(sect, ps1.elementAt(k)))
                    fis = 1;
                if (interSet != null) {
                    if (!onlyFullIntersect)
                        interSet.add(sect);
                    else if (fis == 1)
                        interSet.add(ps1.elementAt(k));
                }
                cntFullIntersect += fis;
                if (onlyFullIntersect && fis == 1)
                    break;
            }
        }
        return cntFullIntersect;
    }

    public static boolean equalSet(long[] a, long[] b) {
        for (int i = 0; i < b.length; i++) {
            if (a[i] != b[i])
                return false;
        }
        return true;
    }


    boolean checkValidPrimers(int[] idx, int[] idx2, Partition[] playground) {
        // ensure there is no partition containing all elements
        // => intersection of parents is empty
        HashMap<PartitionSet, PartitionSet> partitions = (HashMap<PartitionSet, PartitionSet>)
                playground[idx[0]].parents.clone();
        for (int i = 1; i < idx.length; i++) {
            Object[] o = partitions.values().toArray();
            for (int j = 0; j < o.length; j++) {
                if (playground[idx[i]].parents.get(o[j]) == null)
                    partitions.remove(o[j]);
                if (partitions.size() == 0)
                    return true;
            }
        }

        for (int i = 0; i < idx2.length; i++) {
            Object[] o = partitions.values().toArray();
            for (int j = 0; j < o.length; j++) {
                if (playground[idx2[i]].parents.get(o[j]) == null)
                    partitions.remove(o[j]);
                if (partitions.size() == 0)
                    return true;
            }
        }

        return false;
    }

    boolean checkValid_new(int pivot, int[] idx, Partition[] playground) {
        // ensure there is no partition containing all elements
        // => intersection of parents is empty
        HashMap<PartitionSet, PartitionSet> partitions = (HashMap<PartitionSet, PartitionSet>)
                playground[pivot].parents.clone();
        Object[] o;
        for (int i = 0; i < idx.length; i++) {
            o = partitions.values().toArray();
            for (int j = 0; j < o.length; j++) {
                if (playground[idx[i]].parents.get(o[j]) == null)
                    partitions.remove(o[j]);
                if (partitions.size() == 0)
                    return true;
            }
        }

        return false;
    }

    boolean outputCDS = true;

    /**
     * @param k
     * @param src
     * @param snk
     * @param buckets
     */
    void generateTuples(int k, Node src, Node snk, Vector<Vector<Partition>> buckets) {

        int s = 0;
        for (int i = 0; i < buckets.size(); i++)
            s += buckets.elementAt(i).size();
        if (k <= 1)
            k = s;
        Partition[] playground = new Partition[s];
        int[] borders = new int[buckets.size()];
        s = 0;
        int t = 1;
        for (int i = 0; i < buckets.size(); i++) {
            for (int j = 0; j < buckets.elementAt(i).size(); j++)
                playground[s++] = buckets.elementAt(i).elementAt(j);
            if (i < buckets.size() - 1)
                borders[t++] = s;
        }

        int[] idx, idx2;
        // for every pivot region
        for (int i = 1; i < borders.length; i++) {    // min 1 other 'region'
            int left = borders[i - 1];
            if (playground.length - left < k)
                break;

            // now combine all sep combinations in the pivot partition with all (k-sep) combination from the rest
            for (int sep = 1; sep < k; sep++) {    // min 1 has to be separated

                if (sep > (borders[i] - left) || (k - sep) > (playground.length - borders[i]))    // no single combi possible in one of the parted sets
                    continue;    // no break, no?!

                idx = new int[sep];    // init pivot partition combination
                int x = left;
                for (int j = 0; j < idx.length; ++j)
                    idx[j] = x++;


                idx2 = new int[k - sep];
                //for (int j = left; j < borders[i]- sep; j++) {
                while (true) {    // iterate combis on pivot partition, actually idx[0]< borders[i]- sep, but performed by negPtr	

                    x = borders[i];
                    for (int j = 0; j < idx2.length; j++)    // init rest partition 
                        idx2[j] = x++;

                    while (true) {    // iterate combis on rest partition

                        // TODO DEBUG only, kill
//						Transcript[][] ttx= new Transcript[k][];
//						for (int h = 0; h < idx.length; h++) 
//							ttx[h]= decodeTset(playground[idx[h]].transcripts);
//						for (int h = 0; h < idx2.length; h++) 
//							ttx[idx.length+h]= decodeTset(playground[idx2[h]].transcripts);

                        // accept everything for complete events<>outputCDS
                        boolean valid = true;
                        if (!outputCDS)
                            valid = checkValid(idx, idx2, playground);
                        if (valid) {

                            // get transcripts
                            Transcript[][] tt = new Transcript[k][];
                            for (int h = 0; h < idx.length; h++)
                                tt[h] = decodeTset(playground[idx[h]].transcripts);
                            for (int h = 0; h < idx2.length; h++)
                                tt[idx.length + h] = decodeTset(playground[idx2[h]].transcripts);

                            // omit boundaries of ESTs
                            SpliceSite srcSite = src.getSite(), snkSite = snk.getSite();

                            // extract event
                            SpliceSite[][] ss = new SpliceSite[k][];
                            for (int h = 0; h < idx.length; h++)
                                ss[h] = tt[h][0].getSpliceSitesBetween(srcSite, snkSite);
                            for (int h = 0; h < idx2.length; h++)
                                ss[idx.length + h] = tt[idx.length + h][0].getSpliceSitesBetween(srcSite, snkSite);

                            // take min/max position for soft starts/ends
                            // NO: take the one recorded in the hashes
                            for (int h = 0; h < ss.length; h++) {
                                if (ss[h].length < 1)    // no <2, splicechain can only have one site, 
                                    continue;            // but it can be the 1st/last of the transcript

                                SpliceSite refs = null;
                                if (ss[h][0] == tt[h][0].getSpliceSitesAll()[0]) {    // cannot rely on TSS, maybe its an acceptor
                                    if (ss[h].length == 1) {
                                        int p = Arrays.binarySearch(tt[h][0].getSpliceSitesAll(), ss[h][0], SpliceSite.getDefaultPositionTypeComparator());
                                        if (p + 1 < tt[h][0].getSpliceSitesAll().length)
                                            refs = tt[h][0].getSpliceSitesAll()[p + 1];
                                    } else
                                        refs = ss[h][1]; // next in chain
                                }

                                SpliceSite subs = softStartHash.get(refs);
                                if (subs != null && ss[h][0].getSourceType() > Transcript.getEdgeConfidenceLevel()) {    //ss[h][0].getType()== SpliceSite.TYPE_SOFT_START) {
//									for (int j = 1; j < tt[h].length; j++) {
//										if (tt[h][j].getSpliceSitesAll()[0].getPos()< ss[h][0].getPos())
//											ss[h][0]= tt[h][j].getSpliceSitesAll()[0]; 
//									}
//									if (softStartHash.get(ss[h][1])== null)
//										System.currentTimeMillis();
                                    if (subs == srcSite) {
                                        SpliceSite[] ss2 = new SpliceSite[ss[h].length - 1];    // 080822 patch for throwing out the 1st site if it has been extended until src
                                        System.arraycopy(ss[h], 1, ss2, 0, ss2.length);
                                        ss[h] = ss2;
                                    } else
                                        ss[h][0] = subs;
                                }

                                if (ss[h].length == 0)
                                    continue;    // 080822 can have run empty due to removal 

                                refs = null;
                                if (ss[h][ss[h].length - 1] == tt[h][0].getSpliceSitesAll()[tt[h][0].getSpliceSitesAll().length - 1]) {    // cannot rely that it is a TES, can be a SS! 
                                    if (ss[h].length == 1) {
                                        int p = Arrays.binarySearch(tt[h][0].getSpliceSitesAll(), ss[h][ss[h].length - 1], SpliceSite.getDefaultPositionTypeComparator());
                                        if (p - 1 >= 0)
                                            refs = tt[h][0].getSpliceSitesAll()[p - 1];
                                    } else
                                        refs = ss[h][ss[h].length - 2]; // prev in chain
                                }
                                subs = softEndHash.get(refs);
                                if (subs != null && ss[h][ss[h].length - 1].getSourceType() > Transcript.getEdgeConfidenceLevel()) {        //getType()== SpliceSite.TYPE_SOFT_END) {
//									for (int j = 1; j < tt[h].length; j++) {
//										if (tt[h][j].getSpliceSitesAll()[tt[h][j].getSpliceSitesAll().length- 1].getPos()> ss[h][ss[h].length- 1].getPos())
//											ss[h][ss[h].length- 1]= tt[h][j].getSpliceSitesAll()[tt[h][j].getSpliceSitesAll().length- 1]; 
//									}
//									if ( softEndHash.get(ss[h][ss[h].length- 2])== null)
//										System.currentTimeMillis();
                                    if (subs == snkSite) {
                                        SpliceSite[] ss2 = new SpliceSite[ss[h].length - 1];    // 080822 patch for throwing out the 1st site if it has been extended until src
                                        System.arraycopy(ss[h], 0, ss2, 0, ss2.length);
                                        ss[h] = ss2;
                                    } else
                                        ss[h][ss[h].length - 1] = subs;
                                }
                            }

                            if (src.getSite().getPos() == 1396240
                                    && snk.getSite().getPos() == 1396324)
                                System.currentTimeMillis();
                            ASEvent ev = new ASEvent(tt, ss);
                            ev.setAnchors(srcSite, snkSite);
                            ev.setDimension(playground.length);
//							if (ev.toStringStructure().equals("1[2^,1[2^5-7^8-9^,1[2^6-7^8-9^,3[7^8-9^,4[7^8-9^"))
//								System.currentTimeMillis();
//							if(snkSite.getPos()== 155209868/*&& snkSite.getPos()== 155214653*/)
//								System.currentTimeMillis();
                            if (outputCDS) {
                                appendCDS(ev, idx, idx2, playground);
                            }


//							if (isRoot(src)|| isLeaf(snk)) {
                            // filter no_events and other not desired ones
                            byte type = ASEvent.TYPE_UNDEFINED;
                            type = ev.getType();
                            if ((type == ASEvent.TYPE_AS_EVENT && retrieveASEvents)
                                    || (type == ASEvent.TYPE_DS_EVENT && retrieveDSEvents)
                                    || (type == ASEvent.TYPE_VS_EVENT && retrieveVSEvents)) {
                                outputEvent(ev);    // not bucket-size
                            }
//							} else {
//								outputEvent(ev);
//							}
                        }

                        // next combi in rest partition
                        int negPtr = idx2.length - 1;
                        while (negPtr >= 0) {
                            if (idx2[negPtr] < playground.length - ((k - sep) - negPtr)) {
                                int c = ++idx2[negPtr];
                                for (int m = negPtr + 1; m < idx2.length; m++)
                                    idx2[m] = ++c;
                                break;
                            }
                            // else
                            --negPtr;
                        }
                        if (negPtr < 0)
                            break;
                    }

                    // next pivot partition combination
                    int negPtr = idx.length - 1;
                    while (negPtr >= 0) {
                        if (idx[negPtr] < borders[i] - (sep - negPtr)) {
                            int c = ++idx[negPtr];
                            for (int m = negPtr + 1; m < idx.length; m++)
                                idx[m] = ++c;
                            break;
                        }
                        // else
                        --negPtr;
                    }
                    if (negPtr < 0)
                        break;
                }
            }
        }

    }

    Vector<long[][]> generateTuples(int k, long[][][] buckets) {    //Vector<Vector<Path>> buckets

        int[] idx = new int[k];
        for (int i = 0; i < idx.length; i++)
            idx[i] = i;

        int[] combi = new int[k];
//		int permutNb= 1;
//		for (int i = 0; i < combi.length; i++) 
//			;
        Vector<long[][]> tuples = new Vector<long[][]>();
        while (idx[0] < (buckets.length - k + 1)) {

            // now all combinations
            for (int j = 0; j < combi.length; j++)
                combi[j] = 0;
            while (true) {

                long[][] p = new long[k][];
                for (int j = 0; j < combi.length; j++)
                    p[j] = buckets[idx[j]][combi[j]];
                tuples.add(p);

                int negPtr = combi.length - 1;
                while (negPtr >= 0) {
                    if (combi[negPtr] < (buckets[idx[negPtr]].length - 1)) {
                        ++combi[negPtr];
                        break;
                    }
                    // else
                    combi[negPtr] = 0;
                    --negPtr;
                }
                if (negPtr < 0)
                    break;
            }

            //

            int negPtr = idx.length - 1;
            while (negPtr >= 0) {
                if (idx[negPtr] < (buckets.length - 1)) {
                    int c = ++idx[negPtr];
                    for (int i = negPtr + 1; i < idx.length; i++)
                        idx[negPtr] = ++c;
                    break;
                }
                // else
                --negPtr;
            }
        }

        return tuples;
    }

    static TranscriptByNameComparator defaultTranscriptByNameComparator = new TranscriptByNameComparator();
    static Node.PositionTypeComparator defaultNodeByPositionTypeComparator = new Node.PositionTypeComparator();
    protected Gene gene;
    public Transcript[] trpts;    // for flux capacitor
    int taSize;
    protected HashMap<SpliceSite, Node> nodeHash = new HashMap<SpliceSite, Node>();
    protected HashMap<String, SimpleEdge> edgeHash = new HashMap<String, SimpleEdge>();
    public Node leaf;
    public Node root;
    protected Node[] nodesInGenomicOrder = null;
    protected Vector<ASEvent> eventV = null;

    /**
     * All edges in the graph ordered by their genomic positions.
     */
    protected AbstractEdge[] exonicEdgesInGenomicOrder = null;

    /**
     * Converts gene model to splicegraph.
     *
     * @param g a gene
     */
    public SplicingGraph(Gene g) {
        this.gene = g;
        trpts = g.getTranscripts();
        Arrays.sort(trpts, defaultTranscriptByNameComparator);
        taSize = trpts.length / 64;
        if (trpts.length % 64 != 0)
            ++taSize;
    }

    public long[] encodeTset(Transcript[] t) {
        long[] taVector = new long[taSize];
        for (int i = 0; i < taVector.length; i++)
            taVector[i] = 0l;
        for (int i = 0; i < t.length; i++) {
            int p = Arrays.binarySearch(trpts, t[i], defaultTranscriptByNameComparator);
            assert (p >= 0);
            int cnt = 0;
            while (p >= 64) {
                p -= 64;
                ++cnt;
            }
            taVector[cnt] |= MyMath.pow(2, p);
        }
        return taVector;
    }

    public long[] encodeTset(Transcript t) {
        long[] taVector = new long[taSize];
        for (int i = 0; i < taVector.length; i++)
            taVector[i] = 0l;
        int p = Arrays.binarySearch(trpts, t, defaultTranscriptByNameComparator);
        assert (p >= 0);
        int cnt = 0;
        while (p >= 64) {
            p -= 64;
            ++cnt;
        }
        taVector[cnt] |= MyMath.pow(2, p);

        return taVector;
    }

    public long[] encodeTx(int idx) {
        assert (idx >= 0);
        long[] taVector = new long[taSize];
        for (int i = 0; i < taVector.length; i++)
            taVector[i] = 0l;
        addTx(idx, taVector);
        return taVector;
    }

    public boolean addTx(int idx, long[] taVector) {
        assert (idx >= 0);
        int cnt = 0;
        while (idx >= 64) {
            idx -= 64;
            ++cnt;
        }
        long val = MyMath.pow(2, idx);
        boolean wasIn = ((taVector[cnt] & val) == val);
        taVector[cnt] |= val;

        return wasIn;
    }

    public static int getTranscriptNb(long[] c) {
        int cnt = 0;
        for (int i = 0; i < c.length; i++) {
            int base = 0;
            long val = 1;
            int sum = 0;
            while (base < 64) {

                if ((c[i] & val) != 0l) {
                    ++cnt;
                    sum += val;
                    if (sum == c[i])
                        break;
                }
                ++base;
                val *= 2;
            }
        }
        return cnt;
    }

    public Node getNode(SpliceSite ss) {
        return nodeHash.get(ss);
    }

    public SimpleEdge getEdge(Node v, Node w) {
        return edgeHash.get(v.getSite().toString() + w.getSite().toString());
    }

    public Node createNode(SpliceSite ss) {
        Node n = getNode(ss);
        if (n == null) {
            if (canonicalSS && ss.isSpliceSite() && (!ss.isCanonical()))
                return null;
            n = new Node(ss, encodeTset(ss.getTranscripts()));
            nodeHash.put(ss, n);
        }
        return n;
    }

    public SimpleEdge createEdge(Node v, Node w, long[] newTset, byte type) {

        return createEdge(v, w, newTset, type, v.getSite().isLeftFlank() && w.getSite().isRightFlank());

    }

    /**
     * All edge creation delegated to here, to be customized by sub-classes
     *
     * @param v
     * @param w
     * @param newTset
     * @return
     */
    protected SimpleEdge createSimpleEdge(Node v, Node w, long[] newTset) {
        SimpleEdge e = new SimpleEdge(v, w);
        e.setTranscripts(newTset);
        return e;
    }

    /**
     * All super-edge creation delegated to here, to be customized by sub-classes
     *
     * @param newEdges
     * @param newTset
     * @param isPairedEnd
     * @return
     */
    protected SuperEdge createSuperEdge(AbstractEdge[] newEdges, long[] newTset, boolean isPairedEnd) {
        SuperEdge se = new SuperEdge(newEdges, newTset, isPairedEnd);
        return se;
    }


    public SimpleEdge addEdge(Node v, Node w, long[] newTset) {
        SimpleEdge e = createSimpleEdge(v, w, newTset);
        //Edge chk= edgeHash.get(v.getSite().toString()+w.getSite().toString());
        edgeHash.put(v.getSite().toString() + w.getSite().toString(), e);
        return e;
    }

    public void init(int n) {
        constructGraph();
        //cleanGraphByNodes();
        collapseFuzzyFlanks();
        contractGraph(n);
    }

    Node[] removeEstChain(Node n, boolean fromFront) {

        // stop
        int min = 2;
        Transcript[] t = decodeTset(n.getTranscripts());
        if (t.length >= min)    // stop if at least 2 support
            return new Node[]{n};

        // do
        Object[] edges;
        if (fromFront)
            edges = n.getOutEdges().toArray();
        else
            edges = n.getInEdges().toArray();

        Vector<Node> v = new Vector();
        for (int i = 0; i < edges.length; i++) {
            SimpleEdge e = (SimpleEdge) edges[i];
            if (e.type > Transcript.ID_SRC_REFSEQ) { // was: == ID_SRC_MRNA
                Node nextN = null;
                if (fromFront)
                    nextN = e.getHead();
                else
                    nextN = e.getTail();

//					int x;
//					if (fromFront)
//						x= nextN.getInEdges().size();
//					else
//						x= nextN.getOutEdges().size();

                // remove last node
                removeEdge(e);
                if (n.getInEdges().size() == 0 && n.getOutEdges().size() == 0)
                    removeNode(n);

//					if (x== 1) {	// the last one from before
                Node[] nn = removeEstChain(nextN, fromFront);
                for (int j = 0; j < nn.length; j++)
                    ArrayUtils.addUnique(v, nn[j]);
//					}
            }
        }

        Node[] nn = new Node[v.size()];
        for (int i = 0; i < nn.length; i++)
            nn[i] = v.elementAt(i);
        return nn;
    }

    // better on graph-level first/last two edges
    // (alignment errors)
    void collapseFuzzyFlanks(Transcript[] trpts) {
        // filter first exon
        for (int i = 0; i < trpts.length; i++) {
//			if (!trpts[i].getSource().contains("EST"))
//				continue;
            for (int j = (i + 1); j < trpts.length; j++) {
                if (trpts[i].getExons().length == trpts[j].getExons().length)
                    continue;
                int k = 0;
                for (; k < trpts[i].getExons().length; k++) {
                    if (i > 0 && trpts[i].getExons()[k].get5PrimeEdge() != trpts[j].getExons()[k].get3PrimeEdge())
                        break;
                    //if (i< trpts[i].getExons().length- 1&& trpts[i].getExons()[k].get5PrimeEdge()!= trpts[j].getExons()[k].get3PrimeEdge())
                }
            }
        }
    }

    void collapseFuzzyFlanks(boolean forRoot) {
        Iterator<SimpleEdge> iterEdge = null, iterEdge2;
        if (forRoot)
            iterEdge = root.getOutEdges().iterator();
        else
            iterEdge = leaf.getInEdges().iterator();
        SimpleEdge e, f;
        HashMap<Node, Vector<SimpleEdge>> map = new HashMap<Node, Vector<SimpleEdge>>();
        Vector<SimpleEdge> v;
        Node n, u;
        while (iterEdge.hasNext()) {
            e = iterEdge.next();
            // 20100112: replaced
/*			Transcript[] t= decodeTset(e.getTranscripts());
			int x= 0;
			for (; x < t.length; x++) 
				if (!t[x].getSource().contains("Est"))
					break;
			if (x< t.length)
				continue;
*/
            if (forRoot)
                n = e.getHead();    // transcription start
            else
                n = e.getTail();    // transcription end
            if (!n.getSite().isSoftEdge())
                continue;    // 20100112: replacement

            // 20100512
            // removed, all with same 1st/last ss get clustered together
            //if (n.getInEdges().size()== 1&& n.getOutEdges().size()== 1) {	// TODO <k ??
            if (forRoot)
                iterEdge2 = n.getOutEdges().iterator();
            else
                iterEdge2 = n.getInEdges().iterator();
            while (iterEdge2.hasNext()) {
                f = iterEdge2.next();    // first/last exonic edge
                if (forRoot)
                    u = f.getHead();
                else
                    u = f.getTail();
                v = map.get(u);    // pivot splice site (first, resp'y last)
                if (v == null) {
                    v = new Vector<SimpleEdge>();
                    map.put(u, v);
                }
                // only remove second one first, see later whether first edge still needed
                v.add(f);    // first last exonic edges
            }
            //}
        }

        // remove second edges
        Iterator<Node> iterNode = map.keySet().iterator();
        while (iterNode.hasNext()) {
            n = iterNode.next();    // pivot splice site (first, resp'y last)
            v = map.get(n);        // vector(edges) 

            // 20100512: find uttermost tx edge
            f = null;
            for (int i = 0; i < v.size(); i++) {
                u = forRoot ? v.elementAt(i).getTail() : v.elementAt(i).getHead();
                if (f == null || (forRoot && f.getTail().getSite().getPos() < u.getSite().getPos())
                        || (!forRoot) && f.getHead().getSite().getPos() > u.getSite().getPos())
                    f = v.elementAt(i);
            }

            // make new transcript group
            long[] unity = f.getTranscripts();
            if (v.size() > 1) {    // more than one chain of 2 edges
                for (int i = 0; i < v.size(); i++) {
                    e = v.elementAt(i);
                    if (e != f && (forRoot && (!isRoot(e.getTail()))) || ((!forRoot) && (!isLeaf(e.getHead())))) {        // .. because of deleted 2nd edges
                        unity = unite(unity, e.getTranscripts());
                        int delta = forRoot ? e.getTail().getSite().getPos() - f.getTail().getSite().getPos() :
                                e.getHead().getSite().getPos() - f.getHead().getSite().getPos();
                        Transcript[] tt = decodeTset(f.getTranscripts());
                        for (int j = 0; j < tt.length; j++) {
                            if (tt[j].getTranscriptID().equals("ENST00000351671"))
                                System.currentTimeMillis();
                            tt[j].setExonicLength(tt[j].getExonicLength() + delta);
                            if (forRoot)
                                tt[j].getExons()[0].set5PrimeEdge(f.getTail().getSite().getPos());
                            else
                                tt[j].getExons()[tt[j].getExons().length - 1].set3PrimeEdge(f.getHead().getSite().getPos());
                        }
                        removeEdge(e);    // ..nodes will be purged of the end of collapse
                    }
                }
                f.setTranscripts(unity);    // 20100512

                Vector<SimpleEdge> w = forRoot ? f.getTail().getInEdges() :
                        f.getHead().getOutEdges();
                assert (w.size() == 1);    // only connecting to root/leaf
                w.elementAt(0).setTranscripts(unity);
                if (forRoot)
                    w.elementAt(0).getHead().transcripts = unity;
                else
                    w.elementAt(0).getTail().transcripts = unity;

/*				if (forRoot) { 
					addEdge(root, n, unity);	// 20100512: changed from n to u					
				} else {
					addEdge(n, leaf, unity);
				}
*/
            }

        }

        // check first edges, remove no longer needed
        Object[] o = null;
        if (forRoot)
            o = root.getOutEdges().toArray();
        else
            o = leaf.getInEdges().toArray();
        for (int i = 0; i < o.length; i++) {
            e = (SimpleEdge) o[i];
            if ((forRoot && e.getHead().getOutEdges().size() == 0) ||
                    ((!forRoot) && e.getTail().getInEdges().size() == 0))    // no longer needed
                removeEdge(e);
        }
    }

    /**
     * extends evidence >= ec (edge confidence)
     * with the same 1st splice site
     * to the longest exonic evidence.
     *
     * @see #cleanESTborderEdges()
     * @deprecated now in constructGraph
     */
    public void collapseFuzzyFlanks() {

        // kill divergences only by start/end
        collapseFuzzyFlanks(true);
        collapseFuzzyFlanks(false);
    }

    /**
     * Kills all events that start have one variant
     * with exclusively EST evidence (hardcoded!!!)
     * <> you need a first/last common site for EST events
     *
     * @see #collapseFuzzyFlanks(boolean)
     */
    void cleanESTborderEdges() {
        cleanESTborderEdges(root.getOutEdges());
        cleanESTborderEdges(leaf.getInEdges());
    }

    /**
     */
    void cleanESTborderEdges(Vector<SimpleEdge> edges) {
        SimpleEdge[] ee = new SimpleEdge[edges.size()];
        for (int i = 0; i < ee.length; i++)
            ee[i] = edges.elementAt(i);

        // ConcurrentModificationException
//		Iterator<Edge> iter= edges.iterator();
//		while (iter.hasNext()) {
//			Edge e= iter.next();
        for (int x = 0; x < ee.length; ++x) {
            SimpleEdge e = ee[x];
//			try {
//				e= iter.next();
//			} catch (ConcurrentModificationException ex) {
//				System.currentTimeMillis();
//			}
            Transcript[] t = decodeTset(e.getTranscripts());
            int i = 0;
            for (; i < t.length; i++)
                if (t[i].getSourceType() <= Transcript.getEdgeConfidenceLevel())
                    break;
            if (i == t.length)
                removeEdge(e);
        }
    }

    public void contractGraph(int k) {

        // TODO DEBUG, kill
//		if (trpts[0].getTranscriptID().equals("AA002101")) {
//			Transcript[] tx, ty;
//			if (edgeHash.get("163893764-163894039^")!= null)
//				tx= decodeTset(edgeHash.get("163893764-163894039^").getTranscripts());
//			if (edgeHash.get("163893764-163900465)")!= null)
//				ty= decodeTset(edgeHash.get("163893764-163900465)").getTranscripts());
//			System.currentTimeMillis();	
//		}

//		if (trpts[0].getTranscriptID().equals("NM_001005410"))
//			System.currentTimeMillis();

        // contract the rest of the nodes
        Node[] nodes = getNodesInGenomicOrder();    // have to do all for "full breaks" by edge filtering (chr2:179,380,178-179,380,190)
//		if (this.trpts[0].getTranscriptID().equals("AA037639"))
//			System.currentTimeMillis();
        for (int i = 0; i < nodes.length; i++) {
            if (nodes[i].getInEdges().size() > 0)
                continue;
            Vector<SimpleEdge> outV = nodes[i].getOutEdges();
            for (int j = 0; j < outV.size(); ++j) {
                if (!outV.elementAt(j).isProcessed())    // was: isContracted()
                    contractGraph(k, outV.elementAt(j));
            }
        }
        nodesInGenomicOrder = null;

        // purge hashes
        Object[] keys = edgeHash.keySet().toArray();
        Object key;
        for (int i = 0; i < keys.length; i++) {
            SimpleEdge e = edgeHash.get(keys[i]);
            if (e.isContracted())
                continue;
            edgeHash.remove(keys[i]);
            e.getTail().removeOutEdge(e);
            e.getHead().removeInEdge(e);
        }

        keys = nodeHash.keySet().toArray();
        for (int i = 0; i < keys.length; i++) {
            Node v = nodeHash.get(keys[i]);
            if (v.getInEdges().size() > 0 || v.getOutEdges().size() > 0)
                continue;
            nodeHash.remove(keys[i]);
        }

        // TODO DEBUG, kill
//		if (trpts[0].getTranscriptID().equals("AA002101")) {
//			Transcript[] tx, ty;
//			if (edgeHash.get("163893764-163894039^")!= null)
//				tx= decodeTset(edgeHash.get("163893764-163894039^").getTranscripts());
//			if (edgeHash.get("163893764-2147483647?")!= null)
//				ty= decodeTset(edgeHash.get("163893764-2147483647?").getTranscripts());
//			System.currentTimeMillis();	
//		}
    }

    //DFS
    void contractGraph(int k, SimpleEdge e) {
        Node n = e.getHead();
        if (e.isProcessed())
            return;
        e.setProcessed(true);

        Vector<SimpleEdge> outV = n.getOutEdges();
        Vector<SimpleEdge> inV = n.getInEdges();
        SimpleEdge srcEdge = e, endEdge = null, f = null;
        long[] newTset = e.getTranscripts();    // has to be since edges are deleted
        boolean validity = srcEdge.valid;
        while (outV.size() < k && inV.size() < k && outV.size() > 0) {
            f = outV.iterator().next();
            newTset = intersect(newTset, f.getTranscripts());
            if (isNull(newTset))
                break;    // because of deleted edges
            endEdge = f;
            validity &= endEdge.valid;
            n = endEdge.getHead();
            outV = n.getOutEdges();
            inV = n.getInEdges();
        }

        // contract
        if (endEdge != null) {

            f = createSimpleEdge(srcEdge.getTail(), endEdge.getHead(), newTset);
            // srcEdge.getTranscripts() not, since edges are deleted; really TODO
            f.setContracted(true);
            f.setProcessed(true);
            f.valid = validity;
            f.type = srcEdge.type;
            edgeHash.put(f.getTail().getSite().toString() + f.getHead().getSite().toString(), f);    // cannot be in there

        } else {
            e.setContracted(true);
        }

        Object[] o = outV.toArray();
        for (int i = 0; i < o.length; i++) {
            contractGraph(k, (SimpleEdge) o[i]);
        }
    }

    public void removeEdge(SimpleEdge e) {
//		if (e.getTail().getSite().getPos()== -2147483648&& e.getHead().getSite().getPos()== 36710258)
//			System.currentTimeMillis();
        edgeHash.remove(new Integer(e.getTail().getSite().getPos() + e.getHead().getSite().getPos()));
        e.getTail().removeOutEdge(e);
        e.getHead().removeInEdge(e);
    }

    public void removeNode(Node n) {
        nodeHash.remove(n.getSite());
    }

    public Node[] getNodesInGenomicOrder() {
        if (nodesInGenomicOrder == null) {
            Iterator<Node> iter = nodeHash.values().iterator();
            int cnt = 0;
            nodesInGenomicOrder = new Node[nodeHash.size()];
            while (iter.hasNext())
                nodesInGenomicOrder[cnt++] = iter.next();
            Arrays.sort(nodesInGenomicOrder, defaultNodeByPositionTypeComparator);
        }

        return nodesInGenomicOrder;
    }

    private void outputEvent(final ASEvent event) {

        if (writerThread == null || (!writerThread.isAlive()))
            getEventV().add(event);
        else {
//			String s= event.toString();
//			if (s.startsWith("1-2^,1"))
//				System.currentTimeMillis();
//			SpliceSite[] ss= event.getSpliceUniverse(false);
//			for (int i = 0; i < ss.length; i++) {
//				int ctr= 0;
//				for (int j = 0; j < event.getSpliceChains().length; j++) {
//					int p= Arrays.binarySearch(event.getSpliceChains()[j], ss[i], SpliceSite.getDefaultPositionComparator());
//					if (p>= 0)
//						++ctr;
//				}
//				if (ctr== event.getSpliceChains().length)
//					System.currentTimeMillis();
//			}
            writerThread.addEvent(event);
            writerThread.interrupt();    // worked w/o ?!
            // block if too many events in queue
            while (writerThread.queue.size() > writerThread.getMaxQevents()) {
                try {
                    Thread.sleep(100);
                    writerThread.interrupt();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        ++counter;
    }

    public DirectedRegion[] getInternalConstitutiveExonsIntrons() {
        Node[] nodes = getNodesInGenomicOrder();
        long[] trptsAlive = createNullArray(taSize);
        int intronCtr = 0, exonCtr = 0;
        Vector<DirectedRegion> v = new Vector<DirectedRegion>();
        for (int i = 0; i < nodes.length; i++) {
            if (isLeaf(nodes[i]))
                continue;

            // first collect
            for (int j = 0; j < nodes[i].getInEdges().size(); j++) {
                SimpleEdge e = nodes[i].getInEdges().elementAt(j);
                if (isRoot(e.getTail())
                        || (!e.valid) || (!equalSet(trptsAlive, e.getTranscripts())))
                    continue;    // not valid or not constitutive
                SpliceSite tail = e.getTail().getSite(),
                        head = e.getHead().getSite();

                if (tail.isTSS() || head.isTES())
                    continue;    // not internal

                StringBuffer sb = new StringBuffer();
                Transcript[] tt = decodeTset(e.getTranscripts());
                for (int k = 0; k < tt.length; k++) {
                    sb.append(tt[k].getTranscriptID());
                    if (k < tt.length - 1)
                        sb.append("/");
                }

                if (tail.isDonor()) {
                    if (!head.isAcceptor())
                        System.err.println("Intron error in getInternalConstitutiveExonsIntrons()");
                    else {
                        DirectedRegion reg = new DirectedRegion(tail.getPos() + 1, head.getPos() - 1, trpts[0].getStrand());
                        reg.setID("intron");    // +(++intronCtr)
                        reg.setChromosome(trpts[0].getChromosome());
                        reg.addAttribute(Transcript.GTF_ATTRIBUTE_SOURCE, Transcript.getSourceStr(e.type));
                        reg.addAttribute(GFFObject.TRANSCRIPT_ID_TAG, sb.toString());
                        v.add(reg);
                    }
                } else {    // exons
                    if (!head.isDonor())
                        System.err.println("Exon error in getInternalConstitutiveExonsIntrons()");
                    else {
                        if (!acceptableIntrons || (acceptableIntrons && Exon.checkAcceptableIntron(head, tail))) {    // TODO invalid exons are taken into account when determing constitutive introns
                            DirectedRegion reg = new DirectedRegion(tail.getPos(), head.getPos(), trpts[0].getStrand());
                            reg.setID("exon");    // +(++exonCtr)
                            reg.setChromosome(trpts[0].getChromosome());
                            reg.addAttribute(Transcript.GTF_ATTRIBUTE_SOURCE, Transcript.getSourceStr(e.type));
                            reg.addAttribute(GFFObject.TRANSCRIPT_ID_TAG, sb.toString());
                            v.add(reg);
                        }
                    }
                }
            }

            if (!isRoot(nodes[i]))
                trptsAlive = unite(trptsAlive, encodeTset(nodes[i].getSite().getTranscripts()));    // add new born

            for (int j = 0; j < nodes[i].getOutEdges().size(); j++) {    // kill dead oens
                if ((!nodes[i].getOutEdges().elementAt(j).valid)
                        || isLeaf(nodes[i].getOutEdges().elementAt(j).getHead()))
                    trptsAlive = without(trptsAlive, nodes[i].getOutEdges().elementAt(j).getTranscripts());
            }

        }

        return (DirectedRegion[]) ArrayUtils.toField(v);
    }

    public Transcript[] decodeTset(long[] c) {
        int count = 0;
        Vector<Transcript> v = new Vector<Transcript>();
        for (int i = 0; i < c.length; i++) {
            int base = 0;
            long val = 1;
            long sum = 0;
            while (base < 64) {

                if ((c[i] & val) != 0l) {
                    v.add(trpts[count + base]);
                    sum += val;
                    if (sum == c[i])
                        break;
                }
                ++base;
                val *= 2;
            }
            count += 64;
        }
        Transcript[] t = new Transcript[v.size()];
        for (int i = 0; i < t.length; i++)
            t[i] = v.elementAt(i);
        return t;
    }

    public static int getNextTxIdx(long[] partition, int tIdx) {

        // get token
        int base = tIdx + 1;
        int ctr = 0;
        for (; base >= 64; base -= 64, ++ctr) ;
        if (ctr >= partition.length)
            return -1;

        // search this part
        long x = 1;
        for (int i = 0; i < base; ++i, x *= 2) ;
        for (; base < 64; ++base, x *= 2) {    // omit: partition[ctr]>= x&& , overflows			
            if ((partition[ctr] & x) == x)
                return base + (ctr * 64);
        }

        // find next
        base = 0;
        ++ctr;
        for (; ctr < partition.length && partition[ctr] == 0; ++ctr) ;    // skip empty
        if (ctr >= partition.length)
            return -1;

        // search next part
        for (x = 1; base < 64; ++base, x *= 2) {
            if ((partition[ctr] & x) == x)
                return base + (ctr * 64);
        }
        assert (true);    // there must be a next in this idx
        return -1;

    }

    public int getNextTxIdx_old(long[] partition, int tIdx) {

        int base = tIdx + 1;
        int ctr = 0;
        // TODO MAX_VALUE: A constant holding the 
        // maximum value a long can have, 2^{63}-1.
        for (; base > 64; base -= 64, ++ctr) ;
        if (ctr >= partition.length)
            return -1;

        // search this part
        for (; base < 64; ++base) {    // 101028: 2^4= 0
            long x = MyMath.pow(2, base);
            if ((partition[ctr] & x) == x)
                return base + (ctr * 64);
        }

        // find next
        base = 0;
        ++ctr;
        for (; ctr < partition.length && partition[ctr] == 0; ++ctr) ;    // skip empty
        if (ctr >= partition.length)
            return -1;

        // search next part
        for (; base < 64; ++base) {
            long x = MyMath.pow(2, base);
            if ((partition[ctr] & x) == 0)
                return base + (ctr * 64);
        }
        assert (true);    // there must be a next in this idx
        return -1;
    }

    public int decodeCount(long[] c) {
        int count = 0, cntForms = 0;
        for (int i = 0; i < c.length; i++) {
            int base = 0;
            long val = 1;
            int sum = 0;
            while (base < 64) {

                if ((c[i] & val) != 0l) {
                    ++cntForms;
                    sum += val;
                    if (sum == c[i])
                        break;
                }
                ++base;
                val *= 2;
            }
            count += 64;
        }

        return cntForms;
    }

    public int decodeFirstIdx(long[] c) {
        int count = 0;
        for (int i = 0; i < c.length; i++) {
            int base = 0;
            long val = 1;
            while (base < 64) {

                if ((c[i] & val) != 0l)
                    return count + base;
                ++base;
                val *= 2;
            }
            count += 64;
        }

        return -1;
    }

    void cleanGraphByNodes() {
        Object[] nodeO = nodeHash.keySet().toArray();
        for (int j = 0; j < nodeO.length; ++j) {
            Node v = (Node) nodeHash.get(nodeO[j]);
            if ((!v.getSite().isSpliceSite()) || v.getSite().isCanonical())
                continue;
            // delete
            Object[] edgeO = v.getInEdges().toArray();
            for (int i = 0; i < edgeO.length; i++)
                removeEdge((SimpleEdge) edgeO[i]);
            edgeO = v.getOutEdges().toArray();
            for (int i = 0; i < edgeO.length; i++)
                removeEdge((SimpleEdge) edgeO[i]);
            removeNode(v);
        }
    }

    public void writeSpliceJunctionSeqs(int eFlankDon, int eFlankAcc, IntronModel iModel, File f) {
        try {
            BufferedWriter buffy = null;

            if (f == null)
                buffy = new BufferedWriter(new OutputStreamWriter(System.out));
            else
                buffy = new BufferedWriter(new FileWriter(f, true));

            Node[] n = getNodesInGenomicOrder();
            for (int i = 0; i < n.length; i++) {
                if (!n[i].getSite().isDonor())
                    continue;
                for (int j = i + 1; j < n.length; j++) {
                    if (!n[j].getSite().isAcceptor())
                        continue;
                    Vector[] v = getSequence(n[i], eFlankDon, n[j], eFlankAcc, iModel);
                    if (v == null)
                        continue;
                    SimpleEdge e = getEdge(n[i], n[j], false);
                    String sfx = "";
                    if (e != null && e.isIntronic())
                        sfx = "*";
                    for (int k = 0; k < v[0].size(); k++) {
                        for (int m = 0; m < v[1].size(); m++) {
                            buffy.write(">");
                            buffy.write(trpts[0].getGene().getGeneID());
                            buffy.write("_");
                            int start = Math.abs(n[i].getSite().getPos()), end = Math.abs(n[j].getSite().getPos());
//							if (start> end) {
//								int h= start;n[j].getSite().getPos();start= end; end= h;
//							}
                            buffy.write(Integer.toString(start));
                            buffy.write("_");
                            buffy.write(Integer.toString(end));
                            buffy.write(sfx + "\n");
                            // seq
                            buffy.write(v[0].elementAt(k).toString());
                            buffy.write(v[1].elementAt(m).toString());
                            buffy.write("\n");
                        }
                    }
                }
            }

            buffy.flush();
            if (f != null)
                buffy.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public Vector<String>[] getSequence(Node donor, int eFlankDon, Node acceptor, int eFlankAcc, IntronModel iModel) {
        if (!iModel.isValid(acceptor.getSite().getPos() - donor.getSite().getPos() + 1))
            return null;
        String seqDon = barna.model.Graph.readSequence(null,
                donor.getSite().getGene().getChromosome(),
                donor.getSite().getGene().isForward(),
                donor.getSite().getPos() + 1,
                donor.getSite().getPos() + iModel.getMaxDonorLength()    // +1 -1
        );
        String seqAcc = barna.model.Graph.readSequence(null,
                donor.getSite().getGene().getChromosome(),
                donor.getSite().getGene().isForward(),
                acceptor.getSite().getPos() - iModel.getMaxAcceptorLength(),
                acceptor.getSite().getPos()
        );
        if (!iModel.isValid(seqDon, seqAcc))
            return null;

        Vector<String> vecDonSeqs = new Vector<String>(2, 2), vecAccSeqs = new Vector<String>(2, 2);
        char[] donPfx = new char[eFlankDon], accPfx = new char[eFlankAcc];
        walkSeq(donor, donPfx, 0, false, vecDonSeqs);
        walkSeq(acceptor, accPfx, 0, true, vecAccSeqs);

        return new Vector[]{vecDonSeqs, vecAccSeqs};
    }

    void walkSeq(Node n, char[] pfx, int len, boolean forward, Vector<String> v) {

        if ((forward && isLeaf(n)) || ((!forward) && isRoot(n))) {
            Arrays.fill(pfx, len, pfx.length - 1, '.');
            v.add(new String(pfx));
            return;
        }

        Vector<SimpleEdge> vEdge = (forward) ? n.getOutEdges() : n.getInEdges();
        for (int i = 0; i < vEdge.size(); i++) {
            SimpleEdge e = vEdge.elementAt(i);
            if (!e.isExonic())
                continue;
            int p1 = forward ? e.getTail().getSite().getPos() :
                    Math.max(e.getHead().getSite().getPos() - (pfx.length - len) + 1,
                            e.getTail().getSite().getPos());
            int p2 = forward ? Math.min(e.getTail().getSite().getPos() + (pfx.length - len) - 1,
                    e.getHead().getSite().getPos()) :
                    e.getHead().getSite().getPos();
            String seq = barna.model.Graph.readSequence(null,
                    n.getSite().getGene().getChromosome(),
                    n.getSite().getGene().isForward(),
                    p1, p2
            );
            char[] nuPfx = (char[]) pfx.clone();
            if (seq.length() != nuPfx.length - len)
                System.currentTimeMillis();
            System.arraycopy(seq.toCharArray(), 0, nuPfx, len, seq.length());
            if (len + seq.length() < nuPfx.length)
                walkSeq(forward ? e.getHead() : e.getTail(), nuPfx, len + seq.length(), forward, v);
            else {
                assert (seq.length() + len == nuPfx.length);
                v.add(new String(nuPfx));
            }
        }
    }

    /**
     * 080504: separate iterations for soft starts/ends to avoid pbs with single exon transcripts
     */
    public void constructGraph() {

//			for (int i = 0; i < trpts.length; i++) {
//				if (trpts[i].getTranscriptID().equals("NM_206855"))
//					System.currentTimeMillis();
//			}
//			if (trpts[0].getTranscriptID().equals("AA002101"))
//				System.currentTimeMillis();

//			if (trpts[0].getTranscriptID().equals("NM_001005410"))
//				System.currentTimeMillis();

        // create root/leaf
        SpliceSite s = new SpliceSite(Integer.MIN_VALUE, SpliceSite.TYPE_NOT_INITED, gene);
        s.setTranscripts(trpts);
        root = createNode(s, encodeTset(trpts));
        s = new SpliceSite(Integer.MAX_VALUE, SpliceSite.TYPE_NOT_INITED, gene);
        s.setTranscripts(trpts);
        leaf = createNode(s, encodeTset(trpts));

        // for the soft ends, extend them to max exon boundary seen..
        softStartHash = new HashMap<SpliceSite, SpliceSite>();

        // first the intermediate
        for (int i = 0; i < trpts.length; i++) {
            SpliceSite[] ss = trpts[i].getSpliceSitesAll(true);
            for (int j = 2; j <= ss.length - 2; j++) {    // at least one intron
                Node v = createNode(ss[j - 1], encodeTset(new Transcript[]{trpts[i]}));
                if (v == null)    // if non-canonical sites are skipped
                    continue;
                Node w = createNode(ss[j], encodeTset(new Transcript[]{trpts[i]}));
                if (w == null)
                    continue;
                createEdge(v, w, encodeTset(new Transcript[]{trpts[i]}), trpts[i].getSourceType());
            }
        }

        // collect soft starts
        for (int i = 0; i < trpts.length; i++) {
            SpliceSite[] ss = trpts[i].getSpliceSitesAll(true);
            if (ss[0].getSourceType() > Transcript.getEdgeConfidenceLevel()
                    && ss.length > 2) {    // 080826 bugfix
                s = softStartHash.get(ss[1]);
                if (s == null || ss[0].getPos() < s.getPos())
                    s = ss[0];
                softStartHash.put(ss[1], s);
            } else {
                Node v = createNode(ss[0], encodeTset(new Transcript[]{trpts[i]}));
                Node w = createNode(ss[1], encodeTset(new Transcript[]{trpts[i]}));
                createEdge(v, w, encodeTset(new Transcript[]{trpts[i]}), trpts[i].getSourceType());
                createEdge(root, v, encodeTset(new Transcript[]{trpts[i]}), trpts[i].getSourceType());
            }
        }
        // collect soft ends
        softEndHash = new HashMap<SpliceSite, SpliceSite>(softStartHash.size());
        for (int i = 0; i < trpts.length; i++) {
            SpliceSite[] ss = trpts[i].getSpliceSitesAll(true);
            if (ss[ss.length - 1].getSourceType() > Transcript.getEdgeConfidenceLevel()
                    && ss.length > 2) {
                s = softEndHash.get(ss[ss.length - 2]);
                if (s == null || ss[ss.length - 1].getPos() > s.getPos())
                    s = ss[ss.length - 1];
                softEndHash.put(ss[ss.length - 2], s);
            } else {     // hard-wire hard ends
                Node v = createNode(ss[ss.length - 2], encodeTset(new Transcript[]{trpts[i]}));
                Node w = createNode(ss[ss.length - 1], encodeTset(new Transcript[]{trpts[i]}));
                createEdge(v, w, encodeTset(new Transcript[]{trpts[i]}), trpts[i].getSourceType());
                createEdge(w, leaf, encodeTset(new Transcript[]{trpts[i]}), trpts[i].getSourceType());
            }
        }


        // connect soft starts
        Iterator<SpliceSite> iter = softStartHash.keySet().iterator();
        while (iter.hasNext()) {
            s = iter.next();    // get intermediate site to be connected by soft edges and its soft-edge support			
            SpliceSite x = softStartHash.get(s);    // start site
            assert (x.isLeftFlank());
            byte fide = Byte.MAX_VALUE;
            Vector<Transcript> tv = new Vector<Transcript>(s.getTranscripts().length);
            // get transcript support
            for (int i = 0; i < s.getTranscripts().length; i++) {
                SpliceSite[] tiss = s.getTranscripts()[i].getSpliceSitesAll(true);
                if (tiss[0].getSourceType() <= Transcript.getEdgeConfidenceLevel()    // not an EST start site
                        || (tiss[1] != s)
                        || tiss.length == 2)// not a candidate
                    continue;    // no soft edge
                if (s.getTranscripts()[i].getSourceType() < fide)
                    fide = s.getTranscripts()[i].getSourceType();
                tv.add(s.getTranscripts()[i]);
            }
            Transcript[] t = new Transcript[tv.size()];    // TODO
            for (int i = 0; i < t.length; i++) {
                t[i] = tv.elementAt(i);
                if (t[i].get5PrimeEdge() != x.getPos())
                    t[i].getExons()[0].set5PrimeEdge(x.getPos());    // 20100512: correct tx model
            }
//				s.setTranscripts(t);
//				for (int i = 0; i < t.length; i++) {
//					SpliceSite[] ss= t[i].getSpliceSitesAll();
//					ss[ss.length-1]= s;
//				}
            long[] tset = encodeTset(t);
            // connect
            Node v = createNode(x, tset);
            Node w = createNode(s);
            createEdge(v, w, tset, fide);
            createEdge(root, v, tset, fide);
        }


        // connect soft ends
        iter = softEndHash.keySet().iterator();
        while (iter.hasNext()) {
            s = iter.next();    // get intermediate site to be connected by soft edges and its soft-edge support			
            SpliceSite x = softEndHash.get(s);    // soft end site
            assert (x.isRightFlank());
            byte fide = Byte.MAX_VALUE;
            Vector<Transcript> tv = new Vector<Transcript>(s.getTranscripts().length);
            // get transcript support
            for (int i = 0; i < s.getTranscripts().length; i++) {
                SpliceSite[] tiss = s.getTranscripts()[i].getSpliceSitesAll(true);
                if (tiss[0].getSourceType() <= Transcript.getEdgeConfidenceLevel()    // not an EST start site
                        || tiss[tiss.length - 2] != s
                        || tiss.length == 2)    // not a candidate
                    continue;    // no soft edge
                if (s.getTranscripts()[i].getSourceType() < fide)
                    fide = s.getTranscripts()[i].getSourceType();
                tv.add(s.getTranscripts()[i]);
            }
            Transcript[] t = new Transcript[tv.size()];
            for (int i = 0; i < t.length; i++) {
                t[i] = tv.elementAt(i);
                if (t[i].get3PrimeEdge() != x.getPos())
                    t[i].getExons()[t[i].getExons().length - 1].set3PrimeEdge(x.getPos());    // 20100512: correct tx model
            }
//				s.setTranscripts(t);
//				for (int i = 0; i < t.length; i++) {
//					SpliceSite[] ss= t[i].getSpliceSitesAll();
//					ss[ss.length-1]= s;
//				}
            long[] tset = encodeTset(t);
            // connect
            Node v = createNode(s);
            Node w = createNode(x, tset);
            createEdge(v, w, tset, fide);
            createEdge(w, leaf, tset, fide);
        }

    }

    /**
     * create node with a specific transcript set and - in case -
     * add this tset to already existing node.
     *
     * @param ss
     * @param tset
     * @return
     */
    public Node createNode(SpliceSite ss, long[] tset) {
        Node n = getNode(ss);
        if (n == null) {
            if (canonicalSS && ss.isSpliceSite() && (!ss.isCanonical()))
                return null;
            n = new Node(ss, tset);
            nodeHash.put(ss, n);
        } else {
            n.transcripts = unite(n.transcripts, tset);
        }
        return n;
    }

    public HashMap<String, SimpleEdge> getEdgeHash() {
        return edgeHash;
    }

    public static WriterThread getWriterThread() {
        return writerThread;
    }

    public static void setWriterThread(WriterThread writerThread) {
        SplicingGraph.writerThread = writerThread;
    }

    public Vector<ASEvent> getEventV() {
        if (eventV == null) {
            eventV = new Vector<ASEvent>();

        }

        return eventV;
    }

    /**
     * sites stay the same, just reconnect:
     * - traverse the graph with nodes in genomic order
     * - maintain array with exonic status of transcripts
     * - preserve root->TSS, TES->leaf
     *
     * @return
     */
    public void transformToFragmentGraph() {

        Node[] nodes = getNodesInGenomicOrder();
        // set with tx that are exonic at the currently iterated position
        long[] exonic = encodeTset(new Transcript[0]);

        long[] active = encodeTset(new Transcript[0]);

        //long[] intronic= encodeTset(new Transcript[0]);

        for (int i = 1; i < nodes.length - 1; i++) {    // wo root/leaf

            // init current exonic, intronic
            // remove existing edges
            Object[] inEdges = nodes[i].getInEdges().toArray();
            long[] activeNext = active; // keeps track of tx that started and did not end yet
            for (int j = 0; j < inEdges.length; j++) {
                SimpleEdge e = (SimpleEdge) inEdges[j];
                // find new TSSs, actually only one per node can be found
                if (e.getTail() == root)
                    activeNext = unite(activeNext, e.getTranscripts());    // update active list
                    // clean-up exonic edges before segmenting
                else {
                    if (e.isExonic())
                        removeEdgeNew(e);
                    // DO NOT remove non-segmented intron edges,
                    // breaks probably the graph iterations
                }
            }

            // create new exonic in-edges with last segment
            Transcript[] t = decodeTset(exonic);
            byte minConf = Byte.MAX_VALUE;
            for (int j = 0; j < t.length; j++)
                minConf = (byte) Math.min(minConf, t[j].getSourceType());
            if (t.length > 0) {
                SimpleEdge eee = createEdge(nodes[i - 1], nodes[i], exonic, minConf, true);
            } else {
                if (!isNull(active)) {
                    SimpleEdge allintronic = createEdge(nodes[i - 1], nodes[i], without(active, exonic), SimpleEdge.ALL_INTRONIC, false);
                }
            }


            active = activeNext;
            Vector<SimpleEdge> outEdgeV = nodes[i].getOutEdges();
            for (int j = 0; j < outEdgeV.size(); j++) {
                // remove tx that ended from set of active tx
                if (outEdgeV.elementAt(j).getHead() == leaf) {
                    active = without(active, outEdgeV.elementAt(j).getTranscripts());    // update active list

                    // update set of tx that are in exonic state
                } else {
                    if (outEdgeV.elementAt(j).getTail().getSite().isLeftFlank()) {    // intron
                        exonic = unite(exonic, outEdgeV.elementAt(j).getTranscripts());
                    } else {    // exon
                        exonic = without(exonic, outEdgeV.elementAt(j).getTranscripts());
                    }
                }
            }

            // remove in-active tx from exonic set
            exonic = intersect(active, exonic);    // still active?
            //intronic= without(active, exonic);
        }
    }

    public void removeEdgeNew(SimpleEdge e) {
        edgeHash.remove(e.getTail().getSite().toString() + e.getHead().toString() + (e.isExonic() ? "1" : "0"));
        e.getTail().removeOutEdge(e);
        e.getHead().removeInEdge(e);
    }

    IntronModel iModel = new IntronModel();

    public SimpleEdge createEdge(Node v, Node w, long[] newTset, byte type, boolean exonic) {

        SimpleEdge e = getEdge(v, w, exonic);
        if (e == null) {// || !e.isAllIntronic()&&type == SimpleEdge.ALL_INTRONIC) {
            e = createSimpleEdge(v, w, newTset);
            e.type = type;
            if (exonic)
                e.exonic = true;
            edgeHash.put(v.getSite().toString() + w.getSite().toString() + (exonic ? "1" : "0"), e);
            // init validity
            if (acceptableIntrons) {
                if (v.getSite().isDonor() && w.getSite().isAcceptor()) {
                    if (w.getSite().getPos() - v.getSite().getPos() < Exon.MIN_INTRON_LENGTH_HUMAN)
                        e.valid = false;
//					String donSeq= genome.model.Graph.readSequence(v.getSite(), 0, 3).toUpperCase();
//					String accSeq= genome.model.Graph.readSequence(w.getSite(), 0, 0).toUpperCase();
                    if (e.valid //BUGFIX100111: && intronConfidenceLevel> 0 
                            && e.type > intronConfidenceLevel
                            //&& (iModel== null|| !iModel.isValid(donSeq, accSeq))	// GTAG filter
                            && (!Exon.checkAcceptableIntron(v.getSite(), w.getSite()))) {

                        e.valid = false;
                        ++invalidIntrons;
                    }
                    ++totalIntrons;
                }
            }
        } else {
            e.setTranscripts(unite(e.getTranscripts(), newTset));
            if (type < e.type) {
                // update validity
                if ((!e.valid) //BUGFIX100111: && intronConfidenceLevel> 0  
                        && e.type > intronConfidenceLevel
                        && (type <= intronConfidenceLevel
                        //BUGFIX100111: 
                        //&& 
                        || (w.getSite().getPos() - v.getSite().getPos() >= Exon.MIN_INTRON_LENGTH_HUMAN
                        && Exon.checkAcceptableIntron(v.getSite(), w.getSite())))) {

                    e.valid = true;
                    --invalidIntrons;
                }
                e.type = type;
            }
        }

        return e;
    }

    public SimpleEdge getEdge(Node v, Node w, boolean exonic) {
        return edgeHash.get(v.getSite().toString() + w.getSite().toString() + (exonic ? "1" : "0"));
    }

    public Transcript getAnyTranscript(long[] c) {
        for (int i = 0; i < c.length; i++) {
            int base = 0;
            long val = 1;
            while (base < 64) {

                if ((c[i] & val) != 0l)
                    return trpts[(i * 64) + base];
                ++base;
                val *= 2;
            }
        }
        return null;
    }

    public Vector<AbstractEdge> getEdges(long[] sig, boolean pend) {


        Node[] n = getNodesInGenomicOrder();
        Vector<AbstractEdge> v = new Vector<AbstractEdge>(2, 2);
        for (int i = 0; i < n.length; i++) {
            for (int j = 0; j < n[i].getOutEdges().size(); j++) {
                SimpleEdge e = n[i].getOutEdges().elementAt(j);
                if (!e.isExonic() || SplicingGraph.isNull(SplicingGraph.intersect(e.getTranscripts(), sig)))
                    continue;
                //double x= profile.getAreaFrac(e.getFrac(t, readLen), readLen, TProfile.DIR_FORWARD);
                //System.out.println(" "+x+" "+e);
                if (!pend) {
                    v.add(e);
                }
                for (int k = 0; e.getSuperEdges() != null &&
                        k < e.getSuperEdges().size(); k++) {
                    SuperEdge se = e.getSuperEdges().elementAt(k);
                    if (se.getEdges()[0] != e || SplicingGraph.isNull(SplicingGraph.intersect(se.getTranscripts(), sig)))
                        continue;
                    if ((se.isPend() && pend) || ((!se.isPend()) && (!pend))) {
                        v.add(se);
                    }
                    if (((!se.isPend()) && (pend)))
                        for (int m = 0; se.getSuperEdges() != null && m < se.getSuperEdges().size(); m++) {
                            if (SplicingGraph.isNull(SplicingGraph.intersect(se.getSuperEdges().elementAt(m).getTranscripts(), sig)))
                                continue;
                            v.add(se.getSuperEdges().elementAt(m));
                        }
                }
            }
        }


        return v;
    }

    /**
     * excludes:
     * - intronic edges
     * - superedges
     *
     * @return
     */
    public AbstractEdge[] getExonicEdgesInGenomicOrder() {
        if (exonicEdgesInGenomicOrder == null) {
            Vector<AbstractEdge> v = new Vector<AbstractEdge>();
            Iterator<SimpleEdge> iter = edgeHash.values().iterator();
            while (iter.hasNext()) {
                AbstractEdge e = iter.next();
                if (e.isIntronic()
                        || e instanceof SuperEdge
                        || e.getTail() == root
                        || e.getHead() == leaf)
                    continue;
                v.add(e);
            }

            exonicEdgesInGenomicOrder = new SimpleEdge[v.size()];
            for (int i = 0; i < v.size(); i++)
                exonicEdgesInGenomicOrder[i] = v.elementAt(i);
            Arrays.sort(exonicEdgesInGenomicOrder, SimpleEdge.getDefaultPositionComparator());
        }

        return exonicEdgesInGenomicOrder;
    }


    /**
     * returns the first (5'most) exonic edge
     *
     * @return
     */
    public AbstractEdge getFirstExonicEdge() {

        AbstractEdge f = null;
        Iterator<SimpleEdge> iter = edgeHash.values().iterator();
        while (iter.hasNext()) {
            AbstractEdge e = iter.next();
            if (e.isIntronic()
                    || e instanceof SuperEdge
                    || e.getTail() == root
                    || e.getHead() == leaf)
                continue;
            if ((f == null) || (SimpleEdge.getDefaultPositionComparator().compare(e, f) < 0))
                f = e;
        }

        assert (f != null);
        return f;
    }

    protected Vector<SimpleEdge> edgeVector = new Vector<SimpleEdge>();

    /**
     * Checks whether an atomary exonic stretch connects
     * <code>node</code> to <code>node2</code> by iterating
     * all out-edges of <code>node</code>.
     *
     * @param node  upstream site of exonic edge
     * @param node2 downnstream site of exonic edge
     * @return <code>true</code> iff there is an exon / exonic
     *         edge connecting <code>node</code> to <code>node2</code>,
     *         <code>false</code> otherwise.
     */
    public boolean checkExonicEdge(Node node, Node node2) {
        Vector<SimpleEdge> ev = node.getOutEdges();
        for (int i = 0; i < ev.size(); i++) {
            SimpleEdge g = ev.elementAt(i);
            if (g.isExonic() || g.isAllIntronic()) {
                if (g.getHead() == node2)
                    return true;
                // TODO optimize: can maximally be 1 for Segment-Graphs
                // break; 
            }
        }
        return false;
    }


    private int[] su = null;

    /**
     * Returns an array with position of all sites in genomic
     * order.<br>
     * Note: method optimized to initialize array of positions
     * exclusively once, i.e. further added nodes will not be
     * included.
     *
     * @return an array with position of sites in this locus
     */
    protected int[] getSpliceUniverse() {
        if (su == null) {
            Node[] nodes = getNodesInGenomicOrder();
            su = new int[nodes.length];
            for (int i = 0; i < su.length; i++) {
                su[i] = nodes[i].getSite().getPos();
            }
        }

        return su;
    }

    public ASEvent[] getEvents(int k) {
        if (eventV == null) {
            eventV = new Vector<ASEvent>();
            getEventsByPartitions(k);
        }

        ASEvent[] ev = new ASEvent[eventV.size()];
        for (int i = 0; i < ev.length; i++) {
            ev[i] = eventV.elementAt(i);
        }
        return ev;
    }

    public String getPathes(int start, int end, int lenMin, int lenMax) {

        // get exons with bounds
        long t0 = System.currentTimeMillis();
        HashMap<AbstractEdge, Integer> startE = new HashMap<AbstractEdge, Integer>(), endE = new HashMap<AbstractEdge, Integer>();
        //HashMap<Edge, Edge> edgesBetw= new HashMap<Edge, Edge>();
        int edgeSum = 0, max = 0;
        Iterator<SimpleEdge> iter = edgeHash.values().iterator();
        while (iter.hasNext()) {
            AbstractEdge e = iter.next();
            if (!isExon(e))
                continue;
            int offs = 0;
            if (start >= e.getTail().getSite().getPos()
                    && start <= e.getHead().getSite().getPos()) {
                offs = e.getHead().getSite().getPos() - start + 1;
                startE.put(e, offs);
            }
            if (end >= e.getTail().getSite().getPos()
                    && end <= e.getHead().getSite().getPos()) {
                if (offs > 0) {
                    offs = end - start + 1;
                    startE.put(e, offs);
                } else
                    offs = end - e.getTail().getSite().getPos() + 1;
                endE.put(e, offs);
                if (e.getHead().getSite().getPos() > max)
                    max = e.getHead().getSite().getPos();
            }
            if (e.getTail().getSite().getPos() >= start && e.getHead().getSite().getPos() <= end)
                edgeSum += (e.getHead().getSite().getPos()
                        - e.getTail().getSite().getPos() + 1); //edgesBetw.put(e, e);	// mayb too many
            else
                edgeSum += offs;    // half ones
        }
        if (startE.size() == 0)
            return "0 MISSING_START " + start + " " + ((System.currentTimeMillis() - t0) / 1000);
        if (endE.size() == 0)
            return "0 MISSING_END " + end + " " + ((System.currentTimeMillis() - t0) / 1000);

        // discard very bad cases
        if (edgeSum < lenMin)
            return "0 ITS_TOO_BIG " + edgeSum + " " + ((System.currentTimeMillis() - t0) / 1000);

        // check existing paths
        Iterator<AbstractEdge> startIter = startE.keySet().iterator();
        Path p = null;
        while (p == null && startIter.hasNext()) {
            AbstractEdge startN = startIter.next();
            Iterator<AbstractEdge> endIter = endE.keySet().iterator();
            while (p == null && endIter.hasNext()) {
                AbstractEdge endN = endIter.next();
                long[] validP = intersect(startN.transcripts, endN.transcripts);
                if (isNull(validP))
                    continue;
                HashMap<AbstractEdge, Integer> m = new HashMap<AbstractEdge, Integer>(1, 1f);
                m.put(endN, endE.get(endN));
                Path px = new Path();
                px.addEdge(startN, startE.get(startN));
                p = findPath(m, px, max, lenMin, lenMax, validP);
            }
        }
        if (p != null)
            return "1 REAL_PATH " + p.toString() + " " + ((System.currentTimeMillis() - t0) / 1000);

        // check all rest paths
        startIter = startE.keySet().iterator();
        p = null;
        while (p == null && startIter.hasNext()) {
            AbstractEdge e = startIter.next();
            Path px = new Path();
            px.addEdge(e, startE.get(e));
            p = findPath(endE, px, max, lenMin, lenMax, null);
        }
        if (p != null)
            return "1 HYPO_PATH " + p.toString() + " " + ((System.currentTimeMillis() - t0) / 1000);

        return "0 TRIED_ALL null " + ((System.currentTimeMillis() - t0) / 1000);
    }

    /**
     * iterates SJ and exons. NO PAIR-ENDS!
     *
     * @param e
     * @return
     */
    public AbstractEdge nextEdge(AbstractEdge e, long[] partition, boolean oneWayOnly) {

        SuperEdge f = null, se = null;
        Vector<SuperEdge> seV = null;
        if (e instanceof SuperEdge) {
            se = (SuperEdge) e;
            seV = se.edges[0].getSuperEdges();
        } else
            seV = e.getSuperEdges();

        // try to get shortest superEdge
        // starting with the same tail-edge
        // that is at least as long as se (if exists)
        for (int i = 0; i < seV.size(); i++) {
            SuperEdge seX = seV.elementAt(i);
            //			if (partition[0]== 17920&& e.toString().startsWith("54819681-54819741^")) {
            //				Transcript[] tt= decodeTset(seX.transcripts);
            //				System.currentTimeMillis();
            //			}
            if (se != null && (seX == se || seX.edges[0] != se.edges[0]
                    || seX.edges[seX.edges.length - 1].getHead().getSite().getPos() <= se.edges[se.edges.length - 1].getTail().getSite().getPos()))
                continue;
            if (se == null && seX.edges[0] != e)
                continue;
            if (f != null && (seX.edges[seX.edges.length - 1].getHead().getSite().getPos() > f.edges[f.edges.length - 1].getHead().getSite().getPos()))
                continue;
            if (seX.isPend()) // dont walk this way					
                continue;
            long[] p = intersect(seX.transcripts, partition);
            if ((oneWayOnly && !equalSet(p, partition)) ||
                    ((!oneWayOnly) && isNull(p)))
                continue;
            f = seX;
        }
        if (f != null)
            return f;    // found a next super-edge
        if (se != null)
            return se.edges[1];    // go back to next normal edge

        // 0length segments can happen, see chr1:54,132,696-54,132,704 
        if (e.length() == 0) {
            for (int i = 0; i < e.getHead().getOutEdges().size(); i++) {
                if (e.getHead().getOutEdges().elementAt(i).isExonic()
                        && !isNull(intersect(e.getHead().getOutEdges().elementAt(i).getTranscripts(), partition)))
                    return e.getHead().getOutEdges().elementAt(i);
            }
        } else if (Constants.verboseLevel >= Constants.VERBOSE_ERRORS)
            System.err.println("[ERROR] in graph recursion");

        return null;
    }

    private boolean isExon(AbstractEdge e) {
        if (e.getTail().getSite().isLeftFlank() && e.getHead().getSite().isRightFlank())
            return true;
        return false;
    }

    Path findPath(HashMap<AbstractEdge, Integer> endMap, Path p, int maxEnd, int lenMin, int lenMax, long[] validP) {

        // success
        //		if (p.isEmpty())
        //			System.currentTimeMillis();
        if (endMap.get(p.getSinkEdge()) != null) {    // p.getSinkEdge()
            int totLen = p.length() + endMap.get(p.getSinkEdge());
            if (totLen >= lenMin && totLen <= lenMax)
                return p;
            else
                return null;
        }

        Node startN = p.getSinkEdge().getHead();
        Iterator<SimpleEdge> outEdges = startN.outEdges.iterator();
        Path pp = null;
        while (outEdges.hasNext()) {
            SimpleEdge e = outEdges.next();
            if ((validP != null && !equalSet(intersect(e.transcripts, validP), validP))
                    || (isExon(e) && p.length() + length(e) > lenMax)
                    || e.getHead().getSite().getPos() > maxEnd)
                continue;    // die
            pp = p.clonePath();
            pp.addEdge(e, isExon(e) ? length(e) : 0);
            pp = findPath(endMap, pp, maxEnd, lenMin, lenMax, validP);
            if (pp != null)
                break;
        }
        return pp;
    }

    private int length(SimpleEdge e) {
        return (e.getHead().getSite().getPos() - e.getTail().getSite().getPos() + 1);
    }

    public static boolean isAcceptableIntrons() {
        return acceptableIntrons;
    }

    public static void setAcceptableIntrons(boolean acceptableIntrons) {
        SplicingGraph.acceptableIntrons = acceptableIntrons;
    }

    public static byte getIntronConfidenceLevel() {
        return intronConfidenceLevel;
    }

    public static void setIntronConfidenceLevel(byte intronConfidenceLevel) {
        SplicingGraph.intronConfidenceLevel = intronConfidenceLevel;
    }

    public static boolean isRetrieveASEvents() {
        return retrieveASEvents;
    }

    public static void setRetrieveASEvents(boolean retrieveASEvents) {
        SplicingGraph.retrieveASEvents = retrieveASEvents;
    }

    public static boolean isRetrieveDSEvents() {
        return retrieveDSEvents;
    }

    public static void setRetrieveDSEvents(boolean retrieveDSEvents) {
        SplicingGraph.retrieveDSEvents = retrieveDSEvents;
    }

    public static boolean isRetrieveVSEvents() {
        return retrieveVSEvents;
    }

    public static void setRetrieveVSEvents(boolean retrieveVSEvents) {
        SplicingGraph.retrieveVSEvents = retrieveVSEvents;
    }

    /**
     * adds junction edges to the graph
     *
     * @param readLen
     * @return
     */
    public int addEJ(int readLen) {
        AbstractEdge[] edges = getExonicEdgesInGenomicOrder();

        //System.out.println("edges EX "+edges.length);
        int nrSJ = 0;
        for (int i = 0; i < edges.length; i++) {
            //edges[i].setPossReadNr(Math.max(edges[i].length()- readLen+ 1,0));
            Vector<AbstractEdge> tailV = new Vector<AbstractEdge>(2);
            tailV.add(edges[i]);
            Iterator<SimpleEdge> iter = edges[i].getHead().getOutEdges().iterator();
            while (iter.hasNext()) {
                nrSJ = addEJ((Vector<SimpleEdge>) tailV.clone(), iter.next(), readLen, 0, edges[i].getTranscripts(), nrSJ);
            }
        }
        return nrSJ;
        //System.out.println("edges SJ "+g.edgeHash.size());
    }

    public long[] getSupport(AbstractEdge e, AbstractEdge f, int readLen, int[] insertMinMax, long[] part) {
        if (insertMinMax == null || isNull(part))
            return part;
        Transcript[] t = decodeTset(part);
        int cnt = 0;
        for (int i = 0; i < t.length; i++) {
            int[] ee = e.getFrac(t[i], readLen), ff = f.getFrac(t[i], readLen);
            int eMinDist = ff[0] - 1 - (ee[1] + readLen - 1); // -1 for real insert
            int eMaxDist = ff[1] - 1 - (ee[0] + readLen - 1);
            if (eMaxDist < insertMinMax[0] || eMinDist > insertMinMax[1])
                t[i] = null;
            else
                ++cnt;
        }
        Transcript[] tt = new Transcript[cnt];
        cnt = 0;
        for (int i = 0; i < t.length; i++)
            if (t[i] != null)
                tt[cnt++] = t[i];

        return encodeTset(tt);
    }

    int addEJ(Vector<SimpleEdge> tailV, SimpleEdge f, int readLen, int cumuLen, long[] supp, int nr) {

        // break conditions
        if (f.getHead().getSite().getPos() == Integer.MAX_VALUE)
            return nr;

        supp = SplicingGraph.intersect(supp, f.getTranscripts());    // eben net. s.u.
        //supp= Graph.intersect(supp, f.getTail().getTranscripts()); // dont understand the change
        //supp= Graph.intersect(supp, f.getHead().getTranscripts());	// not interesting for exon junction
        if (SplicingGraph.isNull(supp))
            return nr;

        if (f.isIntronic()) {    // pass thru introns
            Vector<SimpleEdge> outV = f.getHead().getOutEdges();
//			if (outV.size()!= 1)
//				System.currentTimeMillis();
            assert (outV.size() == 1);    // should hold, no sedges there, we iterate in genomic order
            supp = intersect(supp, f.getTranscripts());
            nr = addEJ(tailV, outV.elementAt(0), readLen, cumuLen, supp, nr);
            return nr;
        }

        // do
        SimpleEdge e = tailV.elementAt(0);    // if exon fragment, make exon junction to first exon
        int freePos = readLen - cumuLen - 1;
        int possible = freePos, elen = e.length(), flen = f.length();
        if (elen > 0)
            possible = Math.min(possible, elen);    // possible needed lateron for setPoss()
        if (flen > 0)
            possible = Math.min(possible, flen);    // but 0-length edges kill here the extension
        //assert(possible>0);	// there can be 0-length edges, see Edge.length()
        tailV.add(f);
        cumuLen += f.length();

        //if (tailV.elementAt(0).toString().startsWith("19418[19902^"))
        //System.currentTimeMillis();

        if (possible > 0 && cumuLen + e.length() >= 0) { // readLenMin, not too short, create super-edge
            SimpleEdge[] edges = new SimpleEdge[tailV.size()];
            for (int i = 0; i < edges.length; i++)
                edges[i] = tailV.elementAt(i);
            //Arrays.sort(edges, defaultEdgeCoordComparator);
            Transcript[] t = decodeTset(supp);

            int[] pp = SuperEdge.getFrac(t[0], readLen, edges);
            if (pp[0] <= pp[1]) {
                SuperEdge sedge = createSuperEdge(edges, supp, false);
                //sedge.setPossReadNr(Math.max(possible,0));
                ++nr;
            }
        }

        // recursive extension
        if (cumuLen < readLen - 1) {        // break condition: 1nt overlap on one side, 1nt on the other side
            Iterator<SimpleEdge> iter = f.getHead().getOutEdges().iterator();
            while (iter.hasNext()) {
                SimpleEdge g = iter.next();
                if (g.isIntronic() || g.isExonic()) {
                    Vector<SimpleEdge> v = (Vector<SimpleEdge>) tailV.clone();
                    nr = addEJ(v, g, readLen, cumuLen, supp, nr);
                }
            }
        }

        return nr;
    }

    /**
     * finds one/all pathes between start and end in the graph
     *
     * @param start
     * @param end
     * @param lenBounds
     * @param existsOne
     * @return
     */
    public int[] walker(DirectedRegion[] start, DirectedRegion[] end, int[] lenBounds, boolean existsOne) {

        AbstractEdge e = getEdge(start);
        assert (e != null);
        AbstractEdge f = getEdge(end);
        assert (f != null);
        if (e == f)
            return new int[]{end[0].get5PrimeEdge()
                    - start[start.length - 1].get3PrimeEdge() - 1};    // trivial

        long[] partition = intersect(e.getTranscripts(), f.getTranscripts());
        if (isNull(partition))
            return new int[0];    // also trivial

        // from here on, not so trivial
        int pfx = (e.getHead().getSite().getPos() - start[start.length - 1].get3PrimeEdge())
                + (end[0].get5PrimeEdge() - f.getTail().getSite().getPos());
        IntVector pathLengths = new IntVector();
        pathLengths = walkerRek(e, f, partition, pfx, lenBounds, existsOne, pathLengths);

        return pathLengths.toIntArray();
    }

    private IntVector walkerRek(AbstractEdge e, AbstractEdge f, long[] partition, int pfx,
                                int[] lenBounds, boolean existsOne, IntVector pathLengths) {

        if (lenBounds != null && pfx > lenBounds[1])
            return pathLengths;    // easy

        if (e == f) {
            if (lenBounds == null || (pfx >= lenBounds[0] && pfx <= lenBounds[1]))
                pathLengths.add(pfx);
            return pathLengths;    // also easy
        }

        AbstractEdge g = e;
        while ((g = nextEdge(g, partition, false)) != null) { // still easy
            int pfx2 = pfx;
            if (g != f && (!(g instanceof SuperEdge)))
                pfx2 += g.length();    // real length ok?
            pathLengths = walkerRek(g, f, intersect(g.transcripts, partition),
                    pfx2, lenBounds, existsOne, pathLengths);
            if (existsOne && pathLengths.length > 0)
                return pathLengths;
        }

        return pathLengths;

    }

    /**
     * iterates existing super-edges and returns the one with the
     * given edge set, creates new if necessary
     *
     * @param v
     * @param pe
     * @return
     */
    public SuperEdge getSuperEdge(Vector<AbstractEdge> v, boolean pe, long[] part) {

        SuperEdge se;
        AbstractEdge[] ee;
        Collections.sort(v, SimpleEdge.getDefaultPositionComparator());

        // search whether exists
        if (v.elementAt(0).getSuperEdges() != null) {
            Iterator<SuperEdge> iter = v.elementAt(0).getSuperEdges().iterator();
            while (iter.hasNext()) {
                se = iter.next();
                if (se.isPend() != pe)
                    continue;
                ee = se.getEdges();
                if (ee.length != v.size())
                    continue;
                int i = 0;
                for (; i < ee.length; i++) {
                    if (ee[i] != v.elementAt(i))
                        break;
                }
                if (i == ee.length)
                    return se;    // found
            }
        }

        // create new
        ee = new AbstractEdge[v.size()];
        long[] supp = part == null ? v.elementAt(0).getTranscripts().clone() : part;
        for (int i = 0; i < ee.length; i++) {
            ee[i] = v.elementAt(i);
            supp = intersect(supp, ee[i].getTranscripts());
        }
        if (isNull(supp))
            return null;    // paired-end without tx evidence
        se = createSuperEdge(ee, supp, pe);
        return se;
    }

    public SuperEdge getSuperEdge(SimpleEdge e1, SimpleEdge e2, boolean pe, long[] part) {

        SuperEdge se;
        AbstractEdge[] ee;
        if (SimpleEdge.getDefaultPositionComparator().compare(e1, e2) == 1) {
            SimpleEdge h = e1;
            e1 = e2;
            e2 = h;
        }

        // search whether exists
        if (e1.getSuperEdges() != null) {
            Iterator<SuperEdge> iter = e1.getSuperEdges().iterator();
            while (iter.hasNext()) {
                se = iter.next();
                if (se.isPend() != pe)
                    continue;
                ee = se.getEdges();
                if (ee.length != 2)
                    continue;
                if (ee.length == 2 && ee[0] == e1 && ee[1] == e2)
                    return se;    // found
            }
        }

        // create new
        long[] supp = part == null ? e1.getTranscripts().clone() : part;
        supp = intersect(supp, e1.getTranscripts());
        supp = intersect(supp, e2.getTranscripts());
        if (isNull(supp))
            return null;    // paired-end without tx evidence
        se = createSuperEdge(new SimpleEdge[]{e1, e2}, supp, pe);    // TODO: make pairs, no arrays
        return se;
    }


    protected void add(SuperEdge se, long[][] sig, Vector<Vector<AbstractEdge>> v) {
        for (int i = 0; i < sig.length; i++) {
            if (isNull(without(sig[i], se.getTranscripts()))) {    // 1st var completely included
                v.elementAt(i).add(se);
            }
        }
    }

    public static boolean checkEtype(byte eType, AbstractEdge e) {

        if (eType == ETYPE_AL)
            return true;
        if (eType == ETYPE_EX) {
            if (e.isExonic() && (!(e instanceof SuperEdge)))
                return true;
            return false;
        }
        if (eType == ETYPE_IN) {
            if (e.isIntronic() && (!(e instanceof SuperEdge)))
                return true;
            return false;
        }
        if (eType == ETYPE_SJ) {
            try {
                SuperEdge se = (SuperEdge) e;
                if (!se.isPend()) {
                    for (int i = 1; i < se.getEdges().length; i++) {
                        if (se.getEdges()[i - 1].getHead().getSite().isDonor()
                                && se.getEdges()[i].getTail().getSite().isAcceptor())
                            return true;    // at least one
                    }
                    return false;
                } else
                    return false;
            } catch (ClassCastException e2) {
                return false;
            }
        }
        if (eType == ETYPE_XJ) {
            try {
                SuperEdge se = (SuperEdge) e;
                if (!se.isPend()) {
                    for (int i = 1; i < se.getEdges().length; i++) {
                        if (se.getEdges()[i - 1].getHead().getSite() !=
                                se.getEdges()[i].getTail().getSite())
                            return false;    // at least one
                    }
                    return true;
                } else
                    return false;
            } catch (ClassCastException e2) {
                return false;
            }
        }
        if (eType == ETYPE_PE) {
            try {
                SuperEdge se = (SuperEdge) e;
                if (se.isPend())
                    return true;
                return false;
            } catch (ClassCastException e2) {
                return false;
            }

        }
        return false;
    }

    public AbstractEdge getEdge(DirectedRegion regs[]) {

        Vector<SimpleEdge> v = new Vector<SimpleEdge>();
        for (int i = 0; i < regs.length; i++) {

            // convert it to directionality of trpt, now with reads on antisense strand
            DirectedRegion rr = regs[i];
            if (regs[i].getStrand() != gene.getStrand()) {
                rr = new DirectedRegion(-regs[i].getStart(), -regs[i].getEnd(), gene.getStrand());
            }

            Node[] nodes = getNodesInGenomicOrder();
            Node dummy5prime = new Node(new SpliceSite(rr.get5PrimeEdge(), SpliceSite.TYPE_ACCEPTOR, this.gene),
                    SplicingGraph.createNullArray(taSize));

            int pStart = Arrays.binarySearch(nodes, dummy5prime, Node.getDefaultPositionTypeComparator());
            if (pStart < 0)
                pStart = (-(pStart + 1)) - 1; // +1, -1 for the node before ins point

            // collect all normal edges the regions align to
            // regs are exonic regions from 1 read, they are each contained in exactly one edge
            // NO !!! a read can span more than one exonic stretch !
            Node n = nodes[pStart];

            while (n.getSite().getPos() < rr.get3PrimeEdge() ||
                    (n.getSite().isLeftFlank() && n.getSite().getPos() <= rr.get3PrimeEdge())) {
                Vector<SimpleEdge> outEdges = n.getOutEdges();
                int sizeB4 = v.size();
                SimpleEdge iEdge = null;
                for (int j = 0; j < outEdges.size(); j++) {
                    if (outEdges.elementAt(j).isExonic()) {
                        v.add(outEdges.elementAt(j));    // there can only be one exonic outedge
                        n = outEdges.elementAt(j).getHead();
                        break;    // inner loop
                    } else
                        iEdge = outEdges.elementAt(j);
                }
                if (v.size() == sizeB4) {
                    // lemma: there can only be one
                    //						if (iEdge!= null) {
                    //							if (iEdge.getTail().getSite().getPos()< rr.get3PrimeEdge())	// < corrects for exon flank pos
                    //								v.add(iEdge);
                    //						}
                    break;
                }
            }

        }

        // with introns cannot longer be
        if (v.size() == 0)
            return null;    // can be now, due to antisense mapping !

        if (v.size() == 1) {
            int a = Math.abs(v.elementAt(0).getTail().getSite().getPos()), b = Math.abs(v.elementAt(0).getHead().getSite().getPos());
            int edgeLeft = Math.min(a, b), edgeRight = Math.max(a, b);
            a = Math.abs(regs[0].get5PrimeEdge());
            b = Math.abs(regs[regs.length - 1].get3PrimeEdge());
            int regLeft = Math.min(a, b), regRight = Math.max(a, b);
            if (edgeLeft > regLeft || edgeRight < regRight)
                return null; // exceeding gene
            return v.elementAt(0);    // finished, read spans only one edge
        }

        assert (v.size() != 0);
        SimpleEdge[] ve = new SimpleEdge[v.size()];
        for (int i = 0; i < ve.length; i++)
            ve[i] = v.elementAt(i);
        Arrays.sort(ve, SuperEdge.getDefaultEdgeByTailComparator());    // sort for comparison

        // else, check if there is already a superedge with this edge-set
        Vector<SuperEdge> seV = ve[0].getSuperEdges();
        SuperEdge se = null;
        for (int j = 0; seV != null && j < seV.size(); j++) {
            AbstractEdge[] e = seV.elementAt(j).getEdges();    // these are sorted
            if (e.length != ve.length)
                continue;
            int k;
            for (k = 0; k < e.length; k++) {
                if (e[k] != ve[k])    // neg strand not necesarily sorted..
                    break;
            }
            if (k == e.length) {
                se = seV.elementAt(j);
                break;    // se found
            }
        }
        if (se == null) {
            // TODO can be now due to antisense overlap
            return null;
        }
        int a = Math.abs(se.getTail().getSite().getPos()), b = Math.abs(se.getHead().getSite().getPos());
        int edgeLeft = Math.min(a, b), edgeRight = Math.max(a, b);
        a = Math.abs(regs[0].get5PrimeEdge());
        b = Math.abs(regs[regs.length - 1].get3PrimeEdge());
        int regLeft = Math.min(a, b), regRight = Math.max(a, b);
        if (edgeLeft > regLeft || edgeRight < regRight)
            return null;
        //		if (se== null) {
        //			Edge[] e= new Edge[v.size()];
        //			for (int j = 0; j < e.length; j++) 
        //				e[j]= v.elementAt(j);
        //			se= new SuperEdge(e);
        //		}
        return se;

    }

    // no primer/junction ovl <> minPovl= minPlen
    public void getVariations(int minAmplicon, int maxAmplicon,
                              int minPlen, int maxPlen, int minPovl, int seqLen, int minIntronSize) {

        if (nodeHash.size() == 0)
            return;
        Node[] nodes = getNodesInGenomicOrder();
        AbstractEdge[] edges = getExonicEdgesInGenomicOrder();
        //Vector<Node> openSrcVec= new Vector<Node>(nodes.length/ 2);

        HashMap<Transcript, Node> validSinceMap = new HashMap<Transcript, Node>(trpts.length, 1f);

        long[] inter, unity, without;
        HashMap<PartitionSet, Integer> map;
        //				if (trpts[0].getGene().getReferenceTranscript().getTranscriptID().equals("NM_001006658"))
        //					System.currentTimeMillis();

        /*
					 * 080821: parameter to skip edge bubbles, skip all bubbles with src/snk=root/leaf
					 * Theorem: the downproject exclusion cannot be fucked up by that, because
					 * as soon as an inner bubble contains root/leaf as an flank, there cannot be an
					 * outer bubble not including root/leaf.
					 */
        for (int i = 0; i < nodes.length; i++) {

            if (onlyInternal && (i == 0 || i == nodes.length - 1))
                continue;

            // get rev Edges
            Vector<SimpleEdge> revEdges = new Vector<SimpleEdge>(1, 1);
            int revIdxMin = getEdges(nodes[i].getOutEdges(), minPlen, maxPlen, minPovl, true, revEdges);

            // check bubble
            Vector<SimpleEdge> inEdges = nodes[i].getInEdges();
            if (inEdges.size() >= 2) {    // exon/intron or intron/intron


                Vector<Partition> partitions = new Vector<Partition>(inEdges.size());    // TODO hashmaps
                Vector<PartitionSet> partitionSets = new Vector<PartitionSet>(inEdges.size());
                for (int j = 0; j < inEdges.size(); j++) {
                    Partition p = new PartitionExt();
                    p.transcripts = inEdges.elementAt(j).getTranscripts();
                    PartitionSet s = new PartitionSet();
                    p.addParent(s);
                    partitions.add(p);
                    partitionSets.add(s);
                }

                //for (int j = openSrcVec.size()-1; j >= 0; --j) {	// openSrcVec.elementAt(j)
                for (int j = i - 1; j >= 0; --j) {

                    if (onlyInternal && j == 0)
                        continue;

                    if (!intersects(nodes[j].getTranscripts(), nodes[i].getTranscripts()))
                        continue;

                    Vector<SimpleEdge> fwEdges = new Vector<SimpleEdge>(1, 1);
                    int fwIdxMin = getEdges(nodes[j].getInEdges(), minPlen, maxPlen, minPovl, false, fwEdges);

                    long[] interSect = intersect(nodes[j].getTranscripts(), nodes[i].getTranscripts());


                    // check for path lengths/intron lengths
                    for (int m = 0; m < nodes[j].getOutEdges().size(); m++) {
                        SimpleEdge e = nodes[j].getOutEdges().elementAt(m);
                        int len = e.length();
                        for (int k = 0; k < partitions.size(); k++) {
                            inter = intersect(partitions.elementAt(k).transcripts, nodes[j].getOutEdges().elementAt(m).getTranscripts());
                            if (isNull(inter))
                                continue;
                            boolean removeP = false;
                            if (e.isExonic()) {
                                PartitionExt pe = (PartitionExt) partitions.elementAt(k);
                                pe.pathLen += len;
                                if (pe.pathLen > maxAmplicon)
                                    removeP = true;
                            } else if (minIntronSize > 0 && len < minIntronSize) {
                                removeP = true;
                            }

                            if (removeP) {
                                Iterator<PartitionSet> iter = partitions.elementAt(k).parents.keySet().iterator();
                                while (iter.hasNext()) {
                                    PartitionSet ps = iter.next();
                                    ps.partitions.remove(partitions.elementAt(k));
                                    if (ps.partitions.size() == 0)
                                        partitionSets.remove(ps);
                                }
                                partitions.remove(k--);
                            }
                        }
                        if (partitions.size() == 0)
                            break;
                    }
                    if (partitions.size() == 0)
                        break;

//								if (nodes[j].getOutEdges().size()<= 1)
//									continue;

                    // split partitions
                    Vector<Vector<Partition>> splitPathes = new Vector<Vector<Partition>>(partitions.size());
                    for (int k = 0; k < partitions.size(); k++) {
                        Vector<Partition> newPartitions = new Vector<Partition>();    // TODO size

                        for (int m = 0; m < nodes[j].getOutEdges().size(); m++) {

                            if (!nodes[j].getOutEdges().elementAt(m).valid)
                                continue;

                            // TODO check for equalset ?
                            inter = intersect(partitions.elementAt(k).transcripts, nodes[j].getOutEdges().elementAt(m).getTranscripts());
                            if (isNull(inter))
                                continue;

                            without = without(partitions.elementAt(k).transcripts, inter);
                            if (isNull(without)) {
                                Transcript[] tx = decodeTset(inter);    // DEBUG remove
                                newPartitions.add(partitions.remove(k--));    // just temporary remove, parent cannot disappear
                                break;
                            } else {
                                Partition newPartition = (Partition) partitions.elementAt(k).clonePartitionWithoutTx();
                                newPartition.transcripts = inter;
                                newPartitions.add(newPartition);    // new partition
                                partitions.elementAt(k).transcripts = without;
                            }
                        }
                        if (newPartitions.size() > 0)
                            splitPathes.add(newPartitions);
                    }

                    // now add new partitions
                    for (int k = 0; k < splitPathes.size(); k++) {
                        for (int h = 0; h < splitPathes.elementAt(k).size(); h++) {
                            partitions.add(splitPathes.elementAt(k).elementAt(h));
                        }
                    }

                    // combinations
                    if (splitPathes.size() >= 2) {    // DEBUG: later >=1
                        // flatten playground
                        int s = 0;
                        for (int y = 0; y < splitPathes.size(); y++) {
                            for (int z = 0; z < splitPathes.elementAt(y).size(); z++) {
                                PartitionExt pe = ((PartitionExt) splitPathes.elementAt(y).elementAt(z));
                                if (pe.pathLen >= minAmplicon && pe.pathLen <= maxAmplicon)
                                    ++s;
                            }
                        }
                        Partition[] playground = new Partition[s];
                        for (int y = 0; y < splitPathes.size(); y++) {
                            for (int z = 0; z < splitPathes.elementAt(y).size(); z++) {
                                PartitionExt pe = ((PartitionExt) splitPathes.elementAt(y).elementAt(z));
                                if (pe.pathLen >= minAmplicon && pe.pathLen <= maxAmplicon)
                                    playground[s++] = pe;
                            }
                        }
                        System.err.println();
                    }

                    // create a new partition set
                    if (splitPathes.size() > 1) {
                        PartitionSet newSet = new PartitionSet();
                        partitionSets.add(newSet);
                        map = new HashMap<PartitionSet, Integer>();
                        Iterator<PartitionSet> iter;
                        PartitionSet pset;
                        for (int k = 0; k < splitPathes.size(); k++) {
                            for (int h = 0; h < splitPathes.elementAt(k).size(); h++) {
                                iter = splitPathes.elementAt(k).elementAt(h).parents.keySet().iterator();
                                while (iter.hasNext()) {
                                    pset = iter.next();
                                    // 20101022: changed from .get()== null to !containsKey()
                                    if (!pset.partitions.containsKey(splitPathes.elementAt(k).elementAt(h)))
                                        continue;
                                    if (map.get(pset) == null)
                                        map.put(pset, new Integer(1));
                                    else
                                        map.put(pset, new Integer(map.get(pset).intValue() + 1));
                                }
                                splitPathes.elementAt(k).elementAt(h).addParent(newSet);
                            }
                        }

                        // remove un-needed partition-sets
                        iter = map.keySet().iterator();
                        while (iter.hasNext()) {
                            pset = iter.next();
                            if (pset.partitions.size() == map.get(pset).intValue()) {
                                Object[] o = pset.partitions.keySet().toArray();
                                for (int k = 0; k < o.length; k++)
                                    ((Partition) o[k]).parents.remove(pset);
                                pset.partitions = null;
                                partitionSets.remove(pset);
                            }
                        }

                    }

                    // stop
                    inter = intersect(nodes[j].getTranscripts(), nodes[i].getTranscripts());
                    if (equalSet(inter, nodes[i].getTranscripts())) {
                        break;
                    }
                }
            }


        } // for all nodes

    }

    private int getEdges(Vector<SimpleEdge> edgeV, int minPlen, int maxPlen,
                         int minPovl, boolean toRight, Vector<SimpleEdge> value) {
        Iterator<SimpleEdge> iter = edgeV.iterator();
        while (iter.hasNext()) {
            SimpleEdge e = iter.next();
            if (e.isExonic()) {
                assert (value.size() == 0);
                value.add(e);    // TODO break
            }
        }
        if (value.size() == 0)
            return 0;    // only introns adjacent, will be handled on next exonic edge

        // get min/max target edge set
        assert (value.size() == 1);
        int revLen = value.elementAt(0).length(),
                revLen2 = minPovl;    // independent right edge == donor/acceptor
        int revIdxMin = -1;
        if (minPovl > 0) {
            while (revLen2 < maxPlen) {
                if (revIdxMin < 0 && revLen >= minPlen)
                    revIdxMin = value.size() - 1;

                SimpleEdge e = null;
                while (iter.hasNext()) {
                    e = iter.next();
                    if (e.isExonic()) {
                        value.add(e);
                        break;
                    }
                }
                int len = e.length();
                revLen += len;
                revLen2 += len;
            }
        } else {    // minPovl== -1
            if (revLen < minPlen)
                revIdxMin = -1;
            else
                revIdxMin = 0;
        }
        return revIdxMin;
    }

    boolean checkValid(int[] idx, int[] combi, Vector<Vector<Partition>> buckets) {
        // ensure there is no partition containing all elements
        // => intersection of parents is empty
        HashMap<PartitionSet, PartitionSet> partitions = (HashMap<PartitionSet, PartitionSet>)
                buckets.elementAt(idx[0]).elementAt(combi[0]).parents.clone();
        for (int i = 1; i < idx.length; i++) {
            Object[] o = partitions.values().toArray();
            for (int j = 0; j < o.length; j++) {
                if (buckets.elementAt(idx[i]).elementAt(combi[i]).parents.get(o[j]) == null)
                    partitions.remove(o[j]);
                if (partitions.size() == 0)    // as soon as one is from a partition the others are not from
                    return true;
            }
        }
        return false;
    }

    boolean checkValid(int[] idx, int[] idx2, Partition[] playground) {
        // ensure there is no parent partition containing all elements
        // => intersection of parents is empty
        HashMap<PartitionSet, PartitionSet> partitions = (HashMap<PartitionSet, PartitionSet>)
                playground[idx[0]].parents.clone();
        for (int i = 1; i < idx.length; i++) {
            Object[] o = partitions.values().toArray();
            for (int j = 0; j < o.length; j++) {
                if (playground[idx[i]].parents.get(o[j]) == null)
                    partitions.remove(o[j]);
                if (partitions.size() == 0)
                    return true;
            }
        }

        for (int i = 0; i < idx2.length; i++) {
            Object[] o = partitions.values().toArray();
            for (int j = 0; j < o.length; j++) {
                if (playground[idx2[i]].parents.get(o[j]) == null)
                    partitions.remove(o[j]);
                if (partitions.size() == 0)
                    return true;
            }
        }

        return false;
    }

    public boolean contains(long[] transcripts, int txIdx) {

        assert (txIdx >= 0);
        int cnt = 0;
        while (txIdx >= 64) {
            txIdx -= 64;
            ++cnt;
        }

        long x = MyMath.pow(2, txIdx);
        return ((transcripts[cnt] & x) == x);
    }

    public void getEventsByPartitions(int n) {

        boolean outputCombinations = true;
        if (n <= 1) {
            outputCombinations = false;
            n = 2;
        }

        if (nodeHash.size() == 0)
            return;
        Node[] nodes = getNodesInGenomicOrder();
        long[] inter;

        /*
		 * 080821: parameter to skip edge bubbles, skip all bubbles with src/snk=root/leaf
		 * Theorem: the downproject exclusion cannot be fucked up by that, because
		 * as soon as an inner bubble contains root/leaf as an flank, there cannot be an
		 * outer bubble not including root/leaf.
		 */
        // iterate forward
        for (int i = 0; i < nodes.length; i++) {

            if (onlyInternal && (i == 0 || i == nodes.length - 1))
                continue;

            Vector<SimpleEdge> inEdges = nodes[i].getInEdges();
            if (inEdges.size() < n && !outputCDS)
                continue;    // speed-up, rechecked in initPartitions()

            // initialize partitions
            Vector<Partition> partitions = new Vector<Partition>(inEdges.size());    // TODO hashmaps
            Vector<PartitionSet> partitionSets = new Vector<PartitionSet>(inEdges.size());
            int nrPartitions = initPartitions(nodes[i], partitions, partitionSets);
            if ((!outputCDS) && (nrPartitions < n))    // for cds extend everything until 
                continue;

            // iterate backward
            for (int j = i - 1; j >= 0; --j) {

                if (onlyInternal && j == 0)
                    continue;

                if (!intersects(nodes[j].getTranscripts(), nodes[i].getTranscripts()))
                    continue;

                // check for invalid introns TODO outside the (indegree>= n) condition (?)
                if (acceptableIntrons) {
                    removeInvalidPartitions(nodes[j], partitions, partitionSets);
                }
                if (partitions.size() == 0)    // inside acceptableIntrons?
                    break;

                if (nodes[j].getOutEdges().size() <= 1)    // after removal!
                    continue;

                // split partitions wrt CDS? .. should not make any difference
                Vector<Vector<Partition>> splitPathes = splitPartitions(nodes[j], partitions);
                // now add new partitions -> inside splitPartitions?
                for (int k = 0; k < splitPathes.size(); k++) {
                    for (int h = 0; h < splitPathes.elementAt(k).size(); h++) {
                        partitions.add(splitPathes.elementAt(k).elementAt(h));
                    }
                }

                // generate events
                int nrValidSplits = splitPathes.size();    // splicing variants here
                Vector<PartitionCDS> predictedCDS = null;
                if (outputCDS) {

                    predictedCDS = predictCDS(nodes[j], nodes[i], splitPathes);

                    // count valid partitions
                    nrValidSplits = 0;
                    for (int ii = 0; ii < splitPathes.size(); ++ii)
                        for (int jj = 0; jj < splitPathes.elementAt(ii).size(); ++jj)
                            nrValidSplits += ((PartitionCDS) splitPathes.elementAt(ii).elementAt(jj)).cdsValid53 ? 1 : 0;
                }

//				if (nodes[j].getSite().getPos()== 261866
//						&& nodes[i].getSite().getPos()== 262192)
//					System.currentTimeMillis();

                if (nrValidSplits >= n) {
                    if (outputCombinations)
                        generateTuples(n, nodes[j], nodes[i], splitPathes);
                    else {
                        generateTuples(-1, nodes[j], nodes[i], splitPathes);
                    }
                }

                // create new pset, remove no longer needed parts
                if (!outputCDS)
                    updatePartitions(splitPathes, partitions, partitionSets);

                // stop
                inter = intersect(nodes[j].getTranscripts(), nodes[i].getTranscripts());
                if (equalSet(inter, nodes[i].getTranscripts())) {
                    boolean valid = true;
                    if (outputCDS) {
                        // this heuristic does not work
/*						int nrSplitPathes= 0;
						for (int k = 0; k < splitPathes.size(); k++) 
							nrSplitPathes+= splitPathes.elementAt(k).size();
						if (nrValidSplits>= nrSplitPathes)
							valid= true;
*/
                        // all with same 3frame must have same curr5frame		
                        HashMap<Integer, Integer> map3p5p = new HashMap<Integer, Integer>(nrValidSplits);
                        for (int k = 0; valid && k < splitPathes.size(); k++) {
                            for (int k2 = 0; valid && k2 < splitPathes.elementAt(k).size(); k2++) {
                                PartitionCDS p = (PartitionCDS) splitPathes.elementAt(k).elementAt(k2);
                                if (map3p5p.containsKey(p.frame3) && map3p5p.get(p.frame3) != p.currFrame5)
                                    valid = false;
                                else
                                    map3p5p.put(p.frame3, p.currFrame5);
                            }
                        }

                    }
                    if (valid)
                        break;    // exit inner loop
                }

                // re-set cds values for predicted ESTs
                if (outputCDS) {
                    // re-set predicted CDSs
                    for (int k = 0; k < predictedCDS.size(); k++) {
                        PartitionCDS p = predictedCDS.elementAt(k);
                        p.currFrame5 = p.frame3 = Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI];
                    }

                    // re-set vaild pairs
                    for (int k = 0; k < splitPathes.size(); k++)
                        for (int k2 = 0; k2 < splitPathes.elementAt(k).size(); k2++) {
                            PartitionCDS p = (PartitionCDS) splitPathes.elementAt(k).elementAt(k2);
                            p.cdsValid53 = false;     // re-set vaild pairs
                        }

                }

            }    // backward iteration nodes
        } // forward iteration nodes

    }

    private void updatePartitions(Vector<Vector<Partition>> splitPathes, Vector<Partition> partitions, Vector<PartitionSet> partitionSets) {
        // create a new partition set
        if (splitPathes.size() > 1) {
            PartitionSet newSet = new PartitionSet();
            partitionSets.add(newSet);
            HashMap<PartitionSet, Integer> map = new HashMap<PartitionSet, Integer>();
            Iterator<PartitionSet> iter;
            PartitionSet pset;
            for (int k = 0; k < splitPathes.size(); k++) {
                for (int h = 0; h < splitPathes.elementAt(k).size(); h++) {
                    iter = splitPathes.elementAt(k).elementAt(h).parents.keySet().iterator();
                    while (iter.hasNext()) {
                        pset = iter.next();
                        if (pset.partitions.get(splitPathes.elementAt(k).elementAt(h)) == null)
                            continue;
                        if (map.get(pset) == null)
                            map.put(pset, new Integer(1));
                        else
                            map.put(pset, new Integer(map.get(pset).intValue() + 1));
                    }
                    splitPathes.elementAt(k).elementAt(h).addParent(newSet);
                }
            }

            // remove un-needed partition-sets
            iter = map.keySet().iterator();
            while (iter.hasNext()) {
                pset = iter.next();
                if (pset.partitions.size() == map.get(pset).intValue()) {
                    Object[] o = pset.partitions.keySet().toArray();
                    for (int k = 0; k < o.length; k++)
                        ((Partition) o[k]).parents.remove(pset);
                    pset.partitions = null;
                    partitionSets.remove(pset);
                }
            }

        }

    }

    private Vector<Vector<Partition>> splitPartitions(Node nodeJ,
                                                      Vector<Partition> partitions) {

        long[] inter, without;
        SimpleEdge e;
        Partition p;
        Vector<SimpleEdge> outEdges = nodeJ.getOutEdges();
        Vector<Vector<Partition>> splitPathes = new Vector<Vector<Partition>>(partitions.size());
        for (int k = 0; k < partitions.size(); k++) {
            p = partitions.elementAt(k);
            Vector<Partition> newPartitions = new Vector<Partition>();    // TODO size
            for (int m = 0; m < outEdges.size(); m++) {

                e = outEdges.elementAt(m);
                if (!e.valid)
                    continue;

                // TODO check for equalset (?)
                inter = intersect(p.transcripts, e.getTranscripts());
                if (isNull(inter))
                    continue;

                without = without(p.transcripts, inter);
                if (isNull(without)) {
                    newPartitions.add(partitions.remove(k--));    // just temporary remove, parent cannot disappear
                    break;
                } else {
                    Partition newPartition = (Partition) p.clonePartitionWithoutTx();
                    newPartition.transcripts = inter;
                    newPartitions.add(newPartition);    // new partition
                    partitions.elementAt(k).transcripts = without;
                }
            }
            if (newPartitions.size() > 0)
                splitPathes.add(newPartitions);
        }

        return splitPathes;
    }

    private int removeInvalidPartitions(Node nodeJ, Vector<Partition> partitions, Vector<PartitionSet> partitionSets) {

        long[] inter, without;
        Partition p;
        SimpleEdge e;
        int nrPartRemoved = 0, nrPSetsRemoved = 0;
        Vector<SimpleEdge> outEdges = nodeJ.getOutEdges();
        for (int m = 0; m < outEdges.size(); m++) {
            e = outEdges.elementAt(m);
            if (!e.valid) {
                for (int k = 0; k < partitions.size(); k++) {
                    p = partitions.elementAt(k);
                    inter = intersect(p.transcripts, e.getTranscripts());
                    if (isNull(inter))
                        continue;
                    without = without(p.transcripts, inter);
                    if (isNull(without)) {
                        Iterator<PartitionSet> iter = p.parents.keySet().iterator();
                        while (iter.hasNext()) {
                            PartitionSet ps = iter.next();
                            ps.partitions.remove(p);
                            if (ps.partitions.size() == 0) {
                                partitionSets.remove(ps);
                                ++nrPSetsRemoved;
                            }
                        }
                        partitions.remove(k--);
                        ++nrPartRemoved;
                        if (partitions.size() == 0)
                            break;
                    } else
                        p.transcripts = without;

                }
                if (partitions.size() == 0)
                    break;
            }
            if (partitions.size() == 0)
                break;
        }

        return nrPartRemoved;
    }


    private Vector<PartitionCDS> predictCDS(Node nodeJ, Node nodeI, Vector<Vector<Partition>> splitPathes) {

        int posI = nodeI.getSite().getPos(), posJ = nodeJ.getSite().getPos();
        int ni = Translation.getCombinedFrame(Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI],
                Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI]),
                nc = Translation.FRAME_BYTEVAL[Translation.FRAME_BITNC];

        // iterate partitions: determine 5'frame, keep ESTs for later
        Vector<PartitionCDS> predV = new Vector<PartitionCDS>(2, 2);
        HashMap<Integer, Vector<PartitionCDS>> allCombinedFrames = new HashMap<Integer, Vector<PartitionCDS>>(1);
        for (int i = 0; i < splitPathes.size(); i++) {
            for (int j = 0; j < splitPathes.elementAt(i).size(); j++) {
                PartitionCDS p = (PartitionCDS) splitPathes.elementAt(i).elementAt(j);
                if (p.frame3 == ni) {
                    predV.add(p);
                    continue;
                }

                // else (not an EST)
                Transcript tx = trpts[getNextTxIdx(p.transcripts, -1)];    // enough to check 1
                if (tx.getTranslations() == null) {
                    p.setCurrFrame5(nc);
                } else
                    p.setCurrFrame5(tx.getTranslations()[0].getFrameOrRegion(posJ));

                int cds = Translation.getCombinedFrame(p.getCurrFrame5(), p.getFrame3());
                Vector<PartitionCDS> v = allCombinedFrames.get(cds);
                if (v == null) {
                    v = new Vector<PartitionCDS>(2, 2);
                    allCombinedFrames.put(cds, v);
                }
                v.add(p);
            } // for j
        } // for i


        // predict CDSs on ESTs
        Iterator<PartitionCDS> iter = predV.iterator();
        Integer[] a = new Integer[allCombinedFrames.size()];
        a = allCombinedFrames.keySet().toArray(a);
        while (a.length > 0 && iter.hasNext()) {
            PartitionCDS p = iter.next();
            int nrCDSfound = createCDSsingle(p, posJ, posI, a);
            int combi = Translation.getCombinedFrame(p.getCurrFrame5(), p.getFrame3());
            Vector<PartitionCDS> v = allCombinedFrames.get(combi);
            if (v == null) {
                v = new Vector<PartitionCDS>(2, 2);
                allCombinedFrames.put(combi, v);
            }
            v.add(p);
        }

        // mark all with sufficient support >1 as valid
        Iterator<Vector<PartitionCDS>> ii = allCombinedFrames.values().iterator();
        while (ii.hasNext()) {
            Vector<PartitionCDS> v = ii.next();
            for (int j = 0; j < v.size(); j++) {    // iterate ALL m*n
                PartitionCDS p1 = v.elementAt(j);
                Object[] o = p1.parents.values().toArray();
                int cnt = 0;
                for (int i = 0; i < v.size(); i++) {
                    if (i == j)
                        continue;
                    Partition p2 = v.elementAt(i);
                    for (int x = 0; x < o.length; ++x) {
                        if (p2.parents.get(o[x]) == null)
                            ++cnt;
                    }
                }
                if (cnt > 0)
                    v.elementAt(j).cdsValid53 = true;
            }
        }

        return predV;    // ESTs, where CDS has been predicted
    }


    private int createCDSsingle(PartitionCDS partEST, int posJ, int posI, Integer[] combinedCDS) {

        Transcript txEST = trpts[SplicingGraph.getNextTxIdx(partEST.getTranscripts(), -1)];
        int pos1 = txEST.getExonicPosition(Math.max(posJ, txEST.get5PrimeEdge())),
                pos2 = txEST.getExonicPosition(Math.min(posI, txEST.get3PrimeEdge()));    // exonic pos of event
        if (pos1 + 3 > pos2) {    // not possible
            partEST.setCurrFrame5(Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI]);
            partEST.setFrame3(Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI]);
            return 0;
        }
        String seq = txEST.getSplicedSequence();
        int[] cdsLen = new int[combinedCDS.length], cdsFrames = new int[combinedCDS.length];
        int maxLen = -1, cnt = 0;
        for (int i = 0; i < combinedCDS.length; i++) {
            cdsLen[i] = createCDSpartNew(seq, pos1, pos2, partEST, combinedCDS[i]);
            cdsFrames[i] = Translation.getCombinedFrame(partEST.currFrame5, partEST.frame3);
            if (maxLen < cdsLen[i]) {
                maxLen = cdsLen[i];
                cnt = 1;
            } else if (maxLen == cdsLen[i])
                ++cnt;
        }

        // decide for one frame
        if (cnt > 1) {    // try to repair, prefer the predicted with both flanks like ref to only 1 flank
            int cnt2 = 0;
            for (int i = 0; i < cdsFrames.length; ++i) {
                if (cdsLen[i] == maxLen && cdsFrames[i] == combinedCDS[i])
                    ++cnt2;
            }
            for (int i = 0; cnt2 > 0 && i < cdsFrames.length; ++i) {
                if (cdsLen[i] == maxLen && cdsFrames[i] != combinedCDS[i])
                    cdsLen[i] = -1;
            }
        }
        cnt = 0;
        for (int i = 0; i < cdsFrames.length; ++i) {
            if (cdsLen[i] != maxLen)
                continue;
            int frame5 = Translation.get5Frame(cdsFrames[i]),
                    frame3 = Translation.get3Frame(cdsFrames[i]);
            // first CDS
            if (cnt == 0) {
                partEST.frame3 = frame3;
                partEST.currFrame5 = frame5;
                ++cnt;
                // replcace NC with something else .. (UTR, ..)
            } else {
                if (Translation.getFrameVerbose(frame5, frame3).equals("NC"))
                    continue;    // keep the one b4
                else if (Translation.getFrameVerbose(partEST.getCurrFrame5(), partEST.getFrame3()).equals("NC")) {
                    partEST.setCurrFrame5(frame5);
                    partEST.setFrame3(frame3);
                } else
                    ++cnt;    // found another (possibly ident.) CDS for another reference combi
            }
        }

        if (cnt > 1)
            System.currentTimeMillis();

        return cnt;
    }

    private int createCDSpartNew(String seq, int pos1, int pos2, PartitionCDS partEST, int combined) {

        int delta = pos2 - pos1;
        int frame5 = Translation.get5Frame(combined),
                frame3 = Translation.get3Frame(combined);
        boolean coding5 = Translation.isCDS(frame5),
                coding3 = Translation.isCDS(frame3);
        if (!(coding5 || coding3)) {
            partEST.setCurrFrame5(frame5);
            partEST.setFrame3(frame3);
            return 0;    // non-coding
        }
        if (coding5) {
            if (frame5 <= pos1)
                pos1 -= frame5;
            else
                pos1 += 3 - frame5;
        }
        if (coding3) {
            if ((pos2 + 2 - frame3) < seq.length())
                pos2 += (3 - frame3);
            else
                pos2 -= frame3;
        }
        if (pos1 + 3 > pos2)
            return (-1);    // not enough sequence left

        assert (seq.length() % 3 == 0);
        seq = seq.substring(
                Math.min(Math.max(pos1, 0), seq.length()),
                Math.max(Math.min(pos2 + 1, seq.length()), 0)
        );

        // translate and create cds partition
        if (coding5) {
            int p = Translation.findStop(seq);
            partEST.setCurrFrame5(frame5);
            if (p == -1) {
                partEST.setFrame3(Translation.FRAME_BYTEVAL[(frame5 + delta) % 3]);
                return seq.length();
            } else {
                partEST.setFrame3(Translation.FRAME_BYTEVAL[Translation.FRAME_BIT3UTR]);
                return (p + 3);
            }

        } else {
            int p = Translation.findStart(seq);
            if (p == -2) {    // inframe stop encountered
                partEST.setFrame3(Translation.FRAME_BYTEVAL[Translation.FRAME_BITNC]);
                partEST.setCurrFrame5(Translation.FRAME_BYTEVAL[Translation.FRAME_BITNC]);
                return (-2);
            } else {
                partEST.setFrame3(frame3);
                if (p == -1) {    // no start found
                    partEST.setCurrFrame5(Translation.FRAME_BYTEVAL[(delta - frame3) % 3]);
                    return seq.length();
                } else {
                    partEST.setCurrFrame5(Translation.FRAME_BYTEVAL[Translation.FRAME_BIT5UTR]);
                    return (seq.length() - p);
                }
            }
        }
    }

    private boolean checkESTonly(long[] trpts2, int someTxIdx1) {
        byte t = trpts[someTxIdx1].getSourceType();
        boolean b = (t != Transcript.ID_SRC_REFSEQ && t != Transcript.ID_SRC_UCSC);
        while (b && (someTxIdx1 = getNextTxIdx(trpts2, someTxIdx1)) >= 0) {
            t = trpts[someTxIdx1].getSourceType();
            b &= (t != Transcript.ID_SRC_REFSEQ && t != Transcript.ID_SRC_UCSC); // trpts[someTxIdx1].getSourceType()== Transcript.ID_SRC_EST;
        }

        return b;
    }

    /**
     * Partitions tx by inedges and in case by frame of annotated CDS
     *
     * @param nodeI
     * @param partitions
     * @param partitionSets
     */
    private int initPartitions(Node nodeI, Vector<Partition> partitions, Vector<PartitionSet> partitionSets) {

        boolean debug = false;
        Vector<SimpleEdge> inEdges = nodeI.getInEdges();
        Partition p;
        PartitionSet s;
        int ni = Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI],
                nc = Translation.FRAME_BYTEVAL[Translation.FRAME_BITNC];

        // for each inEdge/splice structure variant
        HashMap<Integer, Integer> mapCDSall = null;
        if (outputCDS)
            mapCDSall = new HashMap<Integer, Integer>();
        int posI = nodeI.getSite().getPos(), cds;

        // iterate in-edges
        for (int j = 0; j < inEdges.size(); j++) {
            long[] tx = inEdges.elementAt(j).getTranscripts();
            // split same splice structures according to CDSs
            if (outputCDS) {
                HashMap<Integer, long[]> mapCDSpart = new HashMap<Integer, long[]>();

                // iterate over transcripts,
                // partition subsets with same 3'CDS
                for (int idx = -1; (idx = getNextTxIdx(tx, idx)) >= 0; ) {
                    if (trpts[idx].getTranslations() == null) {
                        if (trpts[idx].getSourceType() == Transcript.ID_SRC_MRNA ||
                                trpts[idx].getSourceType() == Transcript.ID_SRC_EST)
                            cds = ni;
                        else
                            cds = nc;
                    } else
                        cds = trpts[idx].getTranslations()[0].getFrameOrRegion(posI);

                    long[] tcds = mapCDSpart.get(cds);
                    if (tcds == null) {
                        tcds = encodeTx(idx);
                        mapCDSpart.put(cds, tcds);
                    } else
                        addTx(idx, tcds);
                }

                // reference CDS overrides EST (NI)
                if (mapCDSpart.containsKey(ni)) {
                    if (mapCDSpart.size() > 1)
                        mapCDSpart.remove(ni);    // no own partition for EST/mRNA subset
                }

                // create partitions, count how many of each kind
                Iterator<Integer> iter = mapCDSpart.keySet().iterator();
                while (iter.hasNext()) {
                    int b = iter.next();
                    if (mapCDSall.containsKey(b))
                        mapCDSall.put(b, mapCDSall.get(b) + 1);
                    else
                        mapCDSall.put(b, 1);
                    long[] part = mapCDSpart.get(b);
                    p = new PartitionCDS(part, b);
                    s = new PartitionSet();
                    p.addParent(s);
                    partitions.add(p);
                    partitionSets.add(s);
                }


            } else {    /* no cds consideration */
                p = new Partition(tx);
                s = new PartitionSet();
                p.addParent(s);
                partitions.add(p);
                partitionSets.add(s);
            }
        } // for all in-edges


        // count how many variants can generate an event
        // (~size of outer event)
        int nrEffPart = 0;
        if (outputCDS) {
            Iterator<Integer> iter = mapCDSall.keySet().iterator();
            while (iter.hasNext()) {
                int frame3 = iter.next();
                int nrPart = mapCDSall.get(frame3);
                // either NI (EST,mRNA) + something else, or, 2 of a kind to start event
                if (frame3 == ni || nrPart > 1)
                    nrEffPart += nrPart;
            }
        } else
            nrEffPart = partitions.size();    // = nrInEdges 

        return nrEffPart;

    }

    void appendCDS(ASEvent event, int[] idx, int[] idx2, Partition[] playground) {

        Transcript[][] tt = event.getTranscripts();
        HashMap<TxSet, IntVector> littleMap = new HashMap<TxSet, IntVector>(tt.length);
        for (int i = 0; i < tt.length; i++) {
            TxSet p = new TxSet(encodeTset(tt[i]));
            IntVector v = littleMap.get(p);    // can be more then once in different frames
            if (v == null) {
                v = new IntVector(1);
                littleMap.put(p, v);
            }
            v.add(i);
        }

        // collect frame5 x frame3
        PartitionCDS p;
        int evIdx;
        Vector<Partition> v;
        int dim = event.getDimension();
        HashMap<Integer, Vector<Partition>> mapCDSprojection = new HashMap<Integer, Vector<Partition>>(dim * dim);
        for (int i = 0; i < idx.length; i++) {
            p = (PartitionCDS) playground[idx[i]];
            TxSet tset = new TxSet(p.transcripts);
            evIdx = littleMap.get(new TxSet(p.transcripts)).remove(0);
            event.set53Valid(evIdx, p.cdsValid53);
            event.setFrame3(evIdx, p.getFrame3());
            event.setFrame5(evIdx, p.getCurrFrame5());
            if (!p.cdsValid53)
                continue;

            int combi = Translation.getCombinedFrame(p.getCurrFrame5(), p.getFrame3());
            v = mapCDSprojection.get(combi);
            if (v == null) {
                v = new Vector<Partition>(dim);
                mapCDSprojection.put(combi, v);
            }
/*			Iterator<Integer> iter= p.getMapConfirmedCDSflanks().keySet().iterator();
			while(iter.hasNext()) {
				//int combi= Translation.getCombinedFrame(p.getCurrFrame5(), p.getFrame3());
				int combi= iter.next();
				v= mapCDSprojection.get(combi);
				if (v== null) {
					v= new Vector<Partition>(dim);
					mapCDSprojection.put(combi, v);
				}
				v.add(p);
			}
*/
        }

        for (int i = 0; i < idx2.length; i++) {
            p = (PartitionCDS) playground[idx2[i]];
            TxSet tset = new TxSet(p.transcripts);
            evIdx = littleMap.get(tset).remove(0);
            event.setFrame3(evIdx, p.getFrame3());
            event.setFrame5(evIdx, p.getCurrFrame5());
            if (!p.cdsValid53)
                continue;


            int combi = Translation.getCombinedFrame(p.getCurrFrame5(), p.getFrame3());
            v = mapCDSprojection.get(combi);
            if (v == null) {
                v = new Vector<Partition>(dim);
                mapCDSprojection.put(combi, v);
            }
/*			Iterator<Integer> iter= p.getMapConfirmedCDSflanks().keySet().iterator();
			while(iter.hasNext()) {
				//int combi= Translation.getCombinedFrame(p.getCurrFrame5(), p.getFrame3());
				int combi= iter.next();
				v= mapCDSprojection.get(combi);
				if (v== null) {
					v= new Vector<Partition>(dim);
					mapCDSprojection.put(combi, v);
				}
				v.add(p);
			}
*/
        }

        // iterate different frames that are joined between i and j
        Iterator<Integer> iter = mapCDSprojection.keySet().iterator();
        //int[] cdsImpact= new int[mapCDSprojection.size()], temp= new int[mapCDSprojection.size()];
        //int[][] cdsVariants= new int[mapCDSprojection.size()][];
        StringBuilder sbImpact = new StringBuilder(), sbVariants = new StringBuilder();
        while (iter.hasNext()) {    // iterate frame53 combis

            int combi = iter.next();    // frame5-frame3
            v = mapCDSprojection.get(combi);    // vec(partition)
            assert (v.size() > 1);

            // make non-redundant set of psets
            HashMap<PartitionSet, MyLittleInt> mapPset = new HashMap<PartitionSet, MyLittleInt>();
            for (int i = 0; i < v.size(); i++) {
                Iterator<PartitionSet> iter2 = ((PartitionCDS) v.elementAt(i)).parents.keySet().iterator();
                while (iter2.hasNext()) {
                    mapPset.put(iter2.next(), null);
                }
            }

            // group by common parents
            int c = 0;
            Object[] oo = mapPset.keySet().toArray();
            HashMap<Partition, MyLittleInt> mapGroups = new HashMap<Partition, MyLittleInt>(v.size());
            for (int i = 0; i < oo.length; ++i, ++c) {    // iterate parents
                if (oo[i] == null)
                    continue;
                PartitionSet pset = (PartitionSet) oo[i];
                // iterate children, put new counter
                Iterator<Partition> iter2 = pset.partitions.keySet().iterator();
                while (iter2.hasNext()) {    // iterate partitions of parent, update group counter
                    Partition pp = iter2.next();
                    MyLittleInt val = mapGroups.get(pp);
                    if (val == null)
                        mapGroups.put(pp, new MyLittleInt(c));
                    else
                        val.val = c;
                }
            }

            // create map groupNr x vector
            HashMap<MyLittleInt, Vector<Partition>> mapFinal = new HashMap<SplicingGraph.MyLittleInt, Vector<Partition>>();
            Iterator<Partition> iter2 = mapGroups.keySet().iterator();
            while (iter2.hasNext()) {
                Partition p1 = iter2.next();
                MyLittleInt val = mapGroups.get(p1);
                Vector<Partition> v1 = mapFinal.get(val);
                if (v1 == null) {
                    v1 = new Vector<Partition>();
                    mapFinal.put(val, v1);
                }
                v1.add(p1);
            }

            Partition p1;
            Vector<Partition> v1;
            if (mapFinal.size() == 1)
                sbImpact.append("-");
            else {
                String s = Translation.getFrameVerbose(combi);
                if (s == null || s.contains("null")) {
                    int f5 = Translation.get5Frame(combi), f3 = Translation.get3Frame(combi);
                    int combi2 = Translation.getCombinedFrame(f5, f3);
                    System.currentTimeMillis();
                }
                sbImpact.append(s);
            }
            sbImpact.append(";");

/*			Iterator<MyLittleInt> iter3= mapFinal.keySet().iterator();
			while(iter3.hasNext()) {
				MyLittleInt key= iter3.next();
				v1= mapFinal.get(key);
				for (int i = 0; i < v1.size(); i++) {
					if (littleMap.get(new TxSet(v1.elementAt(i).transcripts))!= null)
						sbVariants.append(littleMap.get(new TxSet(v1.elementAt(i).transcripts))+ 1);
					sbVariants.append("/");
				}
				sbVariants.deleteCharAt(sbVariants.length()- 1);
				sbVariants.append(",");
			}
			if (sbVariants.length()> 0)
				sbVariants.deleteCharAt(sbVariants.length()- 1);
			sbVariants.append(";");
*/
        }
        if (sbImpact.length() > 0) {
            sbImpact.deleteCharAt(sbImpact.length() - 1);
//			sbVariants.deleteCharAt(sbVariants.length()- 1);
            event.cdsImpact = sbImpact.toString();
//			event.cdsVariants= sbVariants.toString();
        }
    }

    private void alternateCDS(PartitionCDS partEST, Transcript txEST, int posJ, int posI, Transcript txNEST, int frame5, int frame3) {
        boolean coding5 = Translation.isCDS(frame5),
                coding3 = Translation.isCDS(frame3);
        if (!(coding5 || coding3)) {
            partEST.setCurrFrame5(frame5);
            partEST.setFrame3(frame3);
            return;
        }
        int pos1 = txEST.getExonicPosition(Math.max(posJ, txEST.get5PrimeEdge())),
                pos2 = txEST.getExonicPosition(Math.min(posI, txEST.get3PrimeEdge()));
        int delta = pos2 - pos1;
        String seq = txEST.getSplicedSequence();
        if (coding5) {
            if (frame5 <= pos1)
                pos1 -= frame5;
            else
                pos1 += 3 - frame5;
        }
        if (coding3) {
            if ((pos2 + 2 - frame3) < seq.length())
                pos2 += (3 - frame3);
            else
                pos2 -= frame3;
        }
        if (pos1 + 3 > pos2)
            return;
        assert (seq.length() % 3 == 0);
        seq = seq.substring(
                Math.min(Math.max(pos1, 0), seq.length()),
                Math.max(Math.min(pos2 + 1, seq.length()), 0)
        );
        if (coding5) {
            int p = Translation.findStop(seq);
            partEST.setCurrFrame5(frame5);
            if (p == -1) {
                partEST.setFrame3(Translation.FRAME_BYTEVAL[(frame5 + delta) % 3]);
            } else
                partEST.setFrame3(Translation.FRAME_BYTEVAL[Translation.FRAME_BIT3UTR]);

        } else {
            int p = Translation.findStart(seq);
            if (p == -2) {
                partEST.setFrame3(Translation.FRAME_BYTEVAL[Translation.FRAME_BITNC]);
                partEST.setCurrFrame5(Translation.FRAME_BYTEVAL[Translation.FRAME_BITNC]);
            } else {
                partEST.setFrame3(frame3);
                if (p == -1) {
                    partEST.setCurrFrame5(Translation.FRAME_BYTEVAL[(delta - frame3) % 3]);
                } else
                    partEST.setCurrFrame5(Translation.FRAME_BYTEVAL[Translation.FRAME_BIT5UTR]);
            }
        }
    }

    private Vector<Vector<Partition>> createCDSpartitions(Node nodeJ, Node nodeI, Vector<Vector<Partition>> splitPathes) {

        int posI = nodeI.getSite().getPos(), posJ = nodeJ.getSite().getPos();
        if (posJ == 267227 && posI == 268280) {
            Transcript[][][] tt = new Transcript[splitPathes.size()][][];
            for (int i = 0; i < tt.length; i++) {
                tt[i] = new Transcript[splitPathes.elementAt(i).size()][];
                for (int j = 0; j < tt[i].length; j++) {
                    tt[i][j] = decodeTset(splitPathes.elementAt(i).elementAt(j).transcripts);
                }
            }
            System.currentTimeMillis();
        }

        int nini = Translation.getCombinedFrame(Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI],
                Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI]),
                ncnc = Translation.getCombinedFrame(Translation.FRAME_BYTEVAL[Translation.FRAME_BITNC],
                        Translation.FRAME_BYTEVAL[Translation.FRAME_BITNC]);
        Vector<Vector<Partition>> cdsPathes = new Vector<Vector<Partition>>(splitPathes.size());

        // iterate partitions
        HashMap<Integer, long[]> mapCDSpart = new HashMap<Integer, long[]>(1); // local map when splitting cds'
        HashMap<long[], Integer> mapPart2trans = new HashMap<long[], Integer>(1); // partitions to be translated
        HashMap<Integer, Vector<PartitionCDS>> allCombinedFrames = new HashMap<Integer, Vector<PartitionCDS>>(1);
        for (int i = 0; i < splitPathes.size(); i++) {
            cdsPathes.add(new Vector<Partition>(splitPathes.elementAt(i).size()));
            for (int j = 0; j < splitPathes.elementAt(i).size(); j++) {
                Partition p = splitPathes.elementAt(i).elementAt(j);

                // iterate all transcripts in this partition
                int cds;
                mapCDSpart.clear();
                for (int idx = -1; (idx = getNextTxIdx(p.transcripts, idx)) >= 0; ) {
                    if (trpts[idx].getTranslations() == null) {
                        if (trpts[idx].getSourceType() == Transcript.ID_SRC_MRNA ||
                                trpts[idx].getSourceType() == Transcript.ID_SRC_EST)
                            cds = nini;
                        else
                            cds = ncnc;
                    } else
                        cds = Translation.getCombinedFrame(
                                trpts[idx].getTranslations()[0].getFrameOrRegion(posJ),
                                trpts[idx].getTranslations()[0].getFrameOrRegion(posI));
                    long[] tcds = mapCDSpart.get(cds);
                    if (tcds == null) {
                        tcds = encodeTx(idx);
                        mapCDSpart.put(cds, tcds);
                    } else
                        addTx(idx, tcds);
                }

                // check for non-inited CDSs
                if (mapCDSpart.containsKey(nini)) {
                    if (mapCDSpart.size() == 1)
                        mapPart2trans.put(mapCDSpart.get(nini), i);
                    mapCDSpart.remove(nini);
                }

                // create all partitions that do not have to be predictied
                Iterator<Integer> iter = mapCDSpart.keySet().iterator();
                while (iter.hasNext()) {
                    int combined = iter.next();
                    assert (combined != nini);
                    long[] tx = mapCDSpart.get(combined);

                    PartitionCDS nuCDS = new PartitionCDS();
                    nuCDS.transcripts = tx;
                    nuCDS.frame3 = Translation.get3Frame(combined);
                    nuCDS.currFrame5 = Translation.get5Frame(combined);

                    cdsPathes.elementAt(i).add(nuCDS);
                    Vector<PartitionCDS> v = allCombinedFrames.get(combined);
                    if (v == null) {
                        v = new Vector<PartitionCDS>(1);
                        allCombinedFrames.put(combined, v);
                    }
                    v.add(nuCDS);
                }

            } // for j
        } // for i


        Iterator<Vector<PartitionCDS>> ii = allCombinedFrames.values().iterator();
        while (ii.hasNext()) {    // mark all with sufficient support >1 as valid
            Vector<PartitionCDS> v = ii.next();
            for (int j = 0; v.size() > 1 && j < v.size(); j++)
                v.elementAt(j).cdsValid53 = true;
        }

        if (isRoot(nodeJ) &&
                nodeI.getSite().getPos() == 27764470)
            System.currentTimeMillis();

        // create new cds
        Iterator<long[]> iter = mapPart2trans.keySet().iterator();
        Integer[] a = new Integer[allCombinedFrames.size()];
        a = allCombinedFrames.keySet().toArray(a);
        for (int j = 0; j < a.length; j++) {
            int frame5 = Translation.get5Frame(a[j]),
                    frame3 = Translation.get3Frame(a[j]);
            System.currentTimeMillis();
        }
        while (iter.hasNext()) {
            long[] tx = iter.next();
            int x = mapPart2trans.get(tx);

            Vector<Partition> baseV = cdsPathes.elementAt(x);
            if (a.length == 0) {    // no reference CDSs to extend
                PartitionCDS p = new PartitionCDS();
                p.transcripts = tx;
                int combi = Translation.getCombinedFrame(p.currFrame5, p.frame3);
                Vector<PartitionCDS> v = allCombinedFrames.get(combi);
                if (v == null) {
                    v = new Vector<PartitionCDS>(1);
                    allCombinedFrames.put(combi, v);
                } else {
                    p.cdsValid53 = true;
                    for (int j = 0; j < v.size(); j++)
                        v.elementAt(j).cdsValid53 = true;
                }
                v.add(p);
                baseV.add(p);

            } else {    // there are references


                Vector<PartitionCDS> vcds = createCDSpartNew(tx, posJ, posI, a);
                for (int c = 0; c < vcds.size(); c++) {
                    PartitionCDS p = vcds.elementAt(c);
                    int combi = Translation.getCombinedFrame(p.currFrame5, p.frame3);
                    Vector<PartitionCDS> v = allCombinedFrames.get(combi);
                    if (v == null) {
                        v = new Vector<PartitionCDS>(1);
                        allCombinedFrames.put(combi, v);
                    } else {
                        p.cdsValid53 = true;
                        for (int j = 0; j < v.size(); j++)
                            v.elementAt(j).cdsValid53 = true;
                    }
                    v.add(p);
                    baseV.add(p);
                }
            } // end else references
        } // end iter partitions

        return cdsPathes;
    }

    private Vector<PartitionCDS> createCDSpartNew(long[] tx, int posJ, int posI, Integer[] combinedCDS) {

        Transcript txEST = trpts[SplicingGraph.getNextTxIdx(tx, -1)];
        int pos1 = txEST.getExonicPosition(Math.max(posJ, txEST.get5PrimeEdge())),
                pos2 = txEST.getExonicPosition(Math.min(posI, txEST.get3PrimeEdge()));    // exonic pos of event
        String seq = txEST.getSplicedSequence();
        int[] cdsLen = new int[combinedCDS.length], cdsFrames = new int[combinedCDS.length];
        PartitionCDS partEST = new PartitionCDS();
        partEST.transcripts = tx;
        int maxLen = -1;
        for (int i = 0; i < combinedCDS.length; i++) {
            cdsLen[i] = createCDSpartNew(seq, pos1, pos2, partEST, combinedCDS[i]);
            cdsFrames[i] = Translation.getCombinedFrame(partEST.currFrame5, partEST.frame3);
            if (cdsFrames[i] == 131328) {
                int frame5 = Translation.get5Frame(cdsFrames[i]),
                        frame3 = Translation.get3Frame(cdsFrames[i]);
                System.currentTimeMillis();
            }
            if (maxLen < cdsLen[i]) {
                maxLen = cdsLen[i];
            }
        }

        // decide for one frame
        Vector<PartitionCDS> vv = new Vector<PartitionCDS>(1, 1);
        for (int i = 0; i < cdsFrames.length; i++) {
            if (cdsLen[i] < maxLen)
                continue;
            int frame5 = Translation.get5Frame(cdsFrames[i]),
                    frame3 = Translation.get3Frame(cdsFrames[i]);
            if (frame5 == 0 || frame3 == 0)
                System.currentTimeMillis();
            if (vv.size() == 0) {
                partEST.frame3 = frame3;
                partEST.currFrame5 = frame5;
                vv.add(partEST);
            } else {
                if (Translation.getFrameVerbose(frame5, frame3).equals("NC"))
                    continue;
                else {
                    for (int j = 0; j < vv.size(); j++) {    // replace nc
                        PartitionCDS p2 = vv.elementAt(j);
                        if (Translation.getFrameVerbose(p2.currFrame5, p2.frame3).equals("NC")) {
                            p2.currFrame5 = frame5;
                            p2.frame3 = frame3;
                            break;
                        }
                    }
                }
            }
        }

        // DEBUG
        if (vv.size() > 1) {
            for (int i = 0; i < vv.size(); i++) {
                int frame5 = Translation.get5Frame(vv.elementAt(i).currFrame5),
                        frame3 = Translation.get3Frame(vv.elementAt(i).frame3);
                System.currentTimeMillis();
            }
        }

        return vv;
    }


}
