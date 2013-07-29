package barna.astalavista;

import barna.model.*;
import barna.model.commons.IntVector;
import barna.model.commons.MyTime;
import barna.model.splicegraph.*;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/16/12
 * Time: 5:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class EventExtractor extends SplicingGraph implements Runnable {


    public static int n = 2;
    boolean output = false;

    protected Vector<ASEvent> eventV = null;

    boolean outputCDS = false;

    static long cumulGC = 0l, cumulGF = 0l, cumulGT = 0l, cumulEV = 0l;


    /**
     * Number of AS events found.
     */
    public static int counter= 0;   // TODO not static


    public static WriterThread writerThread = null;    // let it null here, see addEvent()


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

    AStalavistaSettings settings;

    public EventExtractor(Gene g, AStaSettings settings) {
        super(g);
        this.settings= settings;

        // EVENTS_FILE OPTIONS
        if (!settings.get(AStaSettings.EVENTS).isEmpty()) {

            this.onlyInternal= !(settings.get(AStaSettings.EVENTS).contains(AStaSettings.EventTypes.ASE)
                    || settings.get(AStaSettings.EVENTS).contains(AStaSettings.EventTypes.DSP)
                    || settings.get(AStaSettings.EVENTS).contains(AStaSettings.EventTypes.VST));
            this.retrieveASEvents= (settings.get(AStaSettings.EVENTS).contains(AStaSettings.EventTypes.ASE)
                    || settings.get(AStaSettings.EVENTS).contains(AStaSettings.EventTypes.ASI));
            this.retrieveDSEvents= settings.get(AStaSettings.EVENTS).contains(AStaSettings.EventTypes.DSP);
            this.retrieveVSEvents= settings.get(AStaSettings.EVENTS).contains(AStaSettings.EventTypes.VST);

            // Dimension
            int v= (Integer) settings.get(AStaSettings.EVENTS_DIMENSION);
            if (v< 2)
                n= (-1); // complete events
            else
                n= 2;

            acceptableIntrons= settings.get(AStaSettings.EVENTS_ATR).contains(AStaSettings.EventOptions.IOK);
            // consider canonical sites only
            SplicingGraph.canonicalSS= settings.get(AStaSettings.EVENTS_ATR).contains(AStaSettings.EventOptions.CSS);

            // intron confidence
            SplicingGraph.intronConfidenceLevel= (byte) settings.get(AStaSettings.INTRON_CONFIDENCE).intValue();

            // edge confidence
            Transcript.setEdgeConfidenceLevel((byte) settings.get(AStaSettings.EDGE_CONFIDENCE).intValue());

            // Events: predict 3'complete transcripts
            if (settings.get(AStaSettings.EVENTS_ATR).contains(AStaSettings.EventOptions.CP3))
                ASEvent.check3Pcomplete= true;
            // Events: predict NMD
            if (settings.get(AStaSettings.EVENTS_ATR).contains(AStaSettings.EventOptions.NMD))
                ASEvent.checkNMD= true;
            // Events: output flank type (constitutive/alternative)
            if (settings.get(AStaSettings.EVENTS_ATR).contains(AStaSettings.EventOptions.FLT))
                ASEvent.setOutputFlankMode(true);
            if (settings.get(AStaSettings.EVENTS_ATR).contains(AStaSettings.EventOptions.SEQ))
                ; // TODO SplicingGraph.outputSeq= true;

        }


    }

    @Override
    public void run() {

            String tmp = null;
            try {
                long t1 = System.currentTimeMillis();

                // DEBUG
//							Transcript[] trpts= g[i].getTranscripts();
//							Arrays.sort(trpts, defaultTranscriptByNameComparator);
//							if (!trpts[0].getTranscriptID().equals("AA002101"))
//								continue;

                //SplicingGraph gr = new SplicingGraph(g[i]);
                if (output) {
                    System.err.println(trpts[0].getTranscriptID() + "  transcripts " + gene.getTranscriptCount());
                }
                //gr.init(n);
                constructGraph();
                Node[] no = getNodesInGenomicOrder();
                long t2 = System.currentTimeMillis();
                long dA = (t2 - t1);
                cumulGC += dA;
                if (output) {
                    tmp = trpts[0].getTranscriptID() + " construct " + (dA / 1000) + " n " + nodeHash.size() + " e " + edgeHash.size();
                    System.err.println(tmp);
                    System.err.flush();
                }

                // in again
                // DEBUG: take out fuzzy flanks
                collapseFuzzyFlanks(); // formerly called in init()
                //gr.cleanESTborderEdges();	// for Claudia, kill events linked to EST edge variations

                long t22 = System.currentTimeMillis();
                long dF = (t22 - t2);
                cumulGF += dF;
                if (output) {
                    tmp = trpts[0].getTranscriptID() + " flanks " + (dF / 1000) + " n " + nodeHash.size() + " e " + edgeHash.size();
                    System.err.println(tmp);
                    System.err.flush();
                }

                //gr.contractGraph(2);	// not n !!!
                no = getNodesInGenomicOrder();
                long t3 = System.currentTimeMillis();
                long dT = (t3 - t22);
                cumulGT += dT;
                if (output) {
                    tmp = trpts[0].getTranscriptID() + " contract " + (dT / 1000) + " n " + nodeHash.size() + " e " + edgeHash.size();
                    System.err.println(tmp);
                    System.err.flush();
                }
                //int oldSize= evMap.size();

                getEventsByPartitions(n);

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
                System.err.println(gene.getTranscripts()[0].getTranscriptID() + " " + gene.getTranscriptCount() + " " + gene.getSpliceSites().length);
                if (tmp != null)
                    System.err.println(tmp);
            } finally {
                writerThread.interrupt();
            }

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

    public Vector<ASEvent> getEventV() {
        if (eventV == null) {
            eventV = new Vector<ASEvent>();

        }

        return eventV;
    }


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

}
