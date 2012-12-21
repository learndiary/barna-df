package barna.astalavista;

import barna.model.Transcript;
import barna.model.splicegraph.*;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 12/21/12
 * Time: 3:49 PM
 */
public class PrimerDesignerOld extends SplicingGraph {
    
    public PrimerDesignerOld() {
        super(null);
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
                                pe.setPathLen(pe.getPathLen()+ len);
                                if (pe.getPathLen() > maxAmplicon)
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

                            if (!nodes[j].getOutEdges().elementAt(m).isValid())
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
                                if (pe.getPathLen() >= minAmplicon && pe.getPathLen() <= maxAmplicon)
                                    ++s;
                            }
                        }
                        Partition[] playground = new Partition[s];
                        for (int y = 0; y < splitPathes.size(); y++) {
                            for (int z = 0; z < splitPathes.elementAt(y).size(); z++) {
                                PartitionExt pe = ((PartitionExt) splitPathes.elementAt(y).elementAt(z));
                                if (pe.getPathLen() >= minAmplicon && pe.getPathLen() <= maxAmplicon)
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

}
