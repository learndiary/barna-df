package fbi.genome.model.splicegraph;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class GraphTest {

//    static void _240808_test_multithread(String[] args) {
//        SpliceGraphIO io = new SpliceGraphIO();
//        MyFile inputFile= new MyFile(io.parseArguments(args).getAbsolutePath());
//        //
//        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf
//        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716.gtf
//        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716.gtf
//        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr11.gtf
//        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr6.gtf
//        //
//        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807.gtf
//        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807_mRNAs_fromUCSC070919.gtf
//        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807_mRNAs_fromUCSC070919_splicedESTs_fromUCSC070919.gtf
//        boolean output= false, output2= true;
////		if (rusc)
////			outputFname= "delme.asta";
//
//        Graph.WriterThread writerThread= new Graph.WriterThread();
//        writerThread.start();
//
//            // init and start threads
//        long t0= System.currentTimeMillis();
//        if (output2) {
//            io.outputStats(new OutputStreamWriter(System.err));
//            //Date ti= new Date(t0);
//            //System.out.println("["+ti+"]  started, k= "+EventExtractorThread.n+" species "+EventExtractorThread.species+", input file "+inputFile.getAbsolutePath()+", output file= "+outputFname);
//        }
//        //GTFChrReader reader= new GTFChrReader(file.getAbsolutePath());
//        //ChromosomeReaderThread readerThread= new ChromosomeReaderThread(reader);
//        GFFReader reader= new GFFReader(inputFile.getAbsolutePath());
//        if (Graph.readAheadLimit> 0)
//            reader.setReadAheadLimit(Graph.readAheadLimit);
//        reader.setNoIDs(null);
//        //reader.sweepToChromosome("chr17");
//        GeneAheadReaderThread readerThread= new GeneAheadReaderThread(reader);
//        readerThread.setOutput(output);
//        readerThread.setOutput2(output2);
//        readerThread.start();
//        try {
//            readerThread.join();
//            readerThread.getDownstreamThread().join();
//        } catch (InterruptedException e1) {
//            ;	// :)
//        }
//
//        System.err.println("took "+((System.currentTimeMillis()- t0)/1000)+" sec.");
//        try {
//            writerThread.setKill(true);
//            writerThread.interrupt();
//            writerThread.join();
//        } catch (InterruptedException e) {
//            // TODO Auto-generated catch block
//        }
//        System.err.println("found "+Graph.counter+" events.");
//        if (Graph.acceptableIntrons) {
//            System.err.println("discarded "+Graph.invalidIntrons+" introns, " +
//                    "found "+(Graph.totalIntrons- Graph.invalidIntrons)+" valid ones when checking splice sites: " +
//                            "ratio (invalid/total) = "+ MyFormatter.fprint(((double) Graph.invalidIntrons) / Graph.totalIntrons, 2));
//        }
//    }

}
