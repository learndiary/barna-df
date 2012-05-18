package barna.astalavista;

import barna.io.GeneAheadReaderThread;
import barna.io.gtf.GTFwrapper;
import barna.model.ASEvent;
import barna.model.Graph;
import barna.model.Species;
import barna.model.Transcript;
import barna.model.commons.MyFile;
import barna.model.constants.Constants;
import barna.model.splicegraph.SplicingGraph;

import java.io.BufferedWriter;
import java.io.File;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.Date;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 5/16/12
 * Time: 10:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class Asta {


    private static MyFile inputFile;
    public static int counter= 0;
    private static int readAheadLimit= -1;
    static boolean onlyInternal= true;
    public static boolean DEBUG= false;
    static long invalidIntrons= 0, totalIntrons= 0;
    static boolean acceptableIntrons= false; // schmu-buh, immer true sobald intronConfidence gesetzt


    public static void main(String[] args) {
        //_070808_test();
        //gphase.Constants.DATA_DIR= "/home/msammeth";
        //System.out.println(gphase.Constants.DATA_DIR);


        _240808_test_multithread(args);

        /*IntronModel iModel= new IntronModel();
          iModel.read(new File("test.model"));
          genome.model.Graph.overrideSequenceDirPath= "c:\\genomes\\H.sapiens\\golden_path_200603\\chromFa";
          extractSpliceJunctions(30, 30, iModel,
                  new File("M:\\annotations\\hg18_all-mrna_UCSC090507.gtf"),
                  new File("C:\\testJunctions.fasta"));
          */
    }

    static void _240808_test_multithread(String[] args) {

        SplicingGraph.writerThread= new SplicingGraph.WriterThread();
        inputFile= new MyFile(parseArguments(SplicingGraph.writerThread, args).getAbsolutePath());
        //
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr11.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr6.gtf
        //
        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807.gtf
        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807_mRNAs_fromUCSC070919.gtf
        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807_mRNAs_fromUCSC070919_splicedESTs_fromUCSC070919.gtf
        boolean output= false, output2= true;
//		if (rusc)
//			outputFname= "delme.asta";

        SplicingGraph.writerThread.start();

        // init and start threads
        long t0= System.currentTimeMillis();
        if (output2) {
            // writerThread
            outputStats(SplicingGraph.writerThread, new OutputStreamWriter(System.err));
            //Date ti= new Date(t0);
            //System.out.println("["+ti+"]  started, k= "+EventExtractorThread.n+" species "+EventExtractorThread.species+", input file "+inputFile.getAbsolutePath()+", output file= "+outputFname);
        }
        //GTFChrReader reader= new GTFChrReader(file.getAbsolutePath());
        //ChromosomeReaderThread readerThread= new ChromosomeReaderThread(reader);
        GTFwrapper reader= new GTFwrapper(inputFile.getAbsolutePath());
        if (readAheadLimit> 0)
            reader.setReadAheadLimit(readAheadLimit);
        reader.setNoIDs(null);
        //reader.sweepToChromosome("chr17");
        GeneAheadReaderThread readerThread= new GeneAheadReaderThread(reader);
        readerThread.setOutput(output);
        readerThread.setOutput2(output2);
        readerThread.start();
        try {
            readerThread.join();
            readerThread.getDownstreamThread().join();
        } catch (InterruptedException e1) {
            ;	// :)
        }

        System.err.println("took "+((System.currentTimeMillis()- t0)/1000)+" sec.");
        try {
            SplicingGraph.writerThread.setKill(true);
            SplicingGraph.writerThread.interrupt();
            SplicingGraph.writerThread.join();
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
        }
        System.err.println("found "+counter+" events.");
        if (acceptableIntrons) {
            DecimalFormat df = new DecimalFormat("#.##");
            System.err.println("discarded " + invalidIntrons + " introns, " +
                    "found " + (totalIntrons - invalidIntrons) + " valid ones when checking splice sites: " +
                    "ratio (invalid/total) = " + df.format(((double) invalidIntrons) / totalIntrons));
        }
    }

    /**
     * @deprecated recalled method
     * @param args
     * @return
     */
    static File parseArguments(SplicingGraph.WriterThread writerThread, String[] args) {

        System.err.println("\nThis is ASta"
                +", graph-based AS event retriever of the AStalavista package.");

        boolean helpRequested= false;
        for (int i = 0; args!= null&& i < args.length; i++) {
            if (args[i].equals("--help"))
                helpRequested= true;
        }

        if (helpRequested|| args== null|| args.length< 2) {	// -i input file
            System.err.println("Here is a list of the options I understand:\n");
            System.err.println("-i, --input <input file>");
            System.err.println("This is a bit important, I cannot work without an input annotation. I want a GTF file with " +
                    "transcript annotations (exon features, with a mandatory optional attribute named \'transcript_id\') " +
                    "IN THE SAME COLUMN (i.e., if the transcript identifier of the 1st line is in column #10, it has to be in " +
                    "all lines of the file in column #10. The rest of the file should comply with the standard as specified " +
                    "at http://mblab.wustl.edu/GTF2.html.\n"+
                    "There may also be CDS features, but they become only interesting when checking for additional things " +
                    "as NMD probability etc.."+
                    "\n");
            System.err.println("-o, --output <output file|\'stdout\'>");
            System.err.println("Optional, the name of the output file (fully qualified path) OR the keyword \'stdout\' for " +
                    "writing everything to the standard output stream. " +
                    "If nothing is specified, the output will be written to a file \'<input file>_astalavista.gtf.gz\'. " +
                    "\n");
            System.err.println("-g, --genome <path to directory>");
            System.err.println("Path to the directory containing sequence files corresponding to the <seqname> field " +
                    "in the input GTF. A genome directory is required if a intron confidence value is specified." +
                    "\n");
            System.err.println("-k, --dimension <int value>");
            System.err.println("Dimension >1 of the events to be extracted. Default is 2 (i.e., \'pairwise events\'). " +
                    "\n");
            System.err.println("-tmp");
            System.err.println("Set temporary directory" +
                    "\n");
            System.err.println("-ext, +ext");
            System.err.println("(De-)activate external events, i.e. events that include the transcript start or the " +
                    "poly-adenylation site" +
                    "\n");
            System.err.println("-ic, --intronConfidence [int value]");
            System.err.println("Level of intron confidence. The default is to trust all introns. Introns are assigned " +
                    "a confidency class:\n" +
                    "\t 0 if 'RefSeq' appears in the source field of the annotation\n" +
                    "\t 1 if 'mRNA' appears in the source field of the annotation\n" +
                    "\t 2 if 'EST' appears in the source field of the annotation\n" +
                    "\t 3 if if none of the above applies\n" +
                    "all introns of confidency level > intronConfidence will be checked for proper splice sites when extracting events." +
                    "\n");
            System.err.println("-as, +as");
            System.err.println("Deactivate (\'-as\') or activate (\'+as\') the retrieval of Alternative Splicing events. See documentation " +
                    "for the definition of events that suffice an alternative splicing event." +
                    "\n");
            System.err.println("-ds, +ds");
            System.err.println("Deactivate (\'-ds\') or activate (\'+ds\') the retrieval of aDditional splicing events. See documentation " +
                    "for the definition of events that suffice an additional splicing event." +
                    "\n");
            System.err.println("-s, --seqsite");
            System.err.println("Output splice site sequences with events. Requires a reference genome."+
                    "\n");
            System.err.println("--flankType");
            System.err.println("Output the type of the event flanks, i.e., \'constitutive\' or \'alternative\'."+
                    "\n");

            // reactivated on 20100112
            System.err.println("-ec, --edgeConfidence [int value]");
            System.err.println("Level of confidence for edges (i.e., annotated transcription starts/poly-adenylation sites). " +
                    "The default is to trust no annotated edge and to extend overlapping first/last exons of a transcript to " +
                    "their most extreme position. :\n" +
                    "\t 0 if 'RefSeq' appears in the source field of the annotation\n" +
                    "\t 1 if 'mRNA' appears in the source field of the annotation\n" +
                    "\t 2 if 'EST' appears in the source field of the annotation\n" +
                    "\t 3 if if none of the above applies\n" +
                    "all transcript edges of confidency level > edgeConfidence will be extended in case the annotation shows " +
                    "another exon with the same adjacent splice site and an earlier/later start/end." +
                    "\n");
            System.err.println("AStalavista.");
            System.exit(0);
        }

        java.io.File file= null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-i")|| args[i].equals("--input")) {
                if (i+1>= args.length) {
                    System.err.println("Hey, you forgot to specify the input file!");
                    System.exit(-1);
                }
                file= new java.io.File(args[++i]);
                continue;
            }
            if (args[i].equals("-s")|| args[i].equals("--seqsite")) {
                SplicingGraph.outputSeq= true;
                continue;
            }
            if (args[i].equals("-o")|| args[i].equals("--output")) {
                String s= args[++i];
                if (s.equalsIgnoreCase("stdout"))
                    SplicingGraph.writeStdOut= true;
                else
                    writerThread.outputFname= s+".gz";
                continue;
            }
            if (args[i].equals("-g")|| args[i].equals("--genome")) {
                if (i+1>= args.length) {
                    System.err.println("You forgot to specify the genome!");
                    System.exit(-1);
                }
                MyFile checkFile= new MyFile(args[i+1]);
                if (checkFile.exists()&& checkFile.isDirectory())
                    Graph.overrideSequenceDirPath= checkFile.getAbsolutePath();
                else {
                    String[] s= args[++i].split("_");
                    if (s.length!= 2) {
                        System.err.println("Invalid genome directory: "+args[i+1]);
                        System.exit(-1);
                    }
                    Species spe= new Species(s[0]);
                    spe.setGenomeVersion(s[1]);
                    SplicingGraph.EventExtractorThread.setSpecies(spe);
                }
                acceptableIntrons= true;
                continue;
            }

            if (args[i].equals("-k")|| args[i].equals("--dimension")) {
                try {
                    int x= Integer.parseInt(args[++i]);
                    if (x< -1|| ((x> -1)&& (x< 2)))
                        System.err.println(args[i]+" is not a valid dimension, ignored");
                    SplicingGraph.EventExtractorThread.n= x;
                } catch (NumberFormatException e) {
                    System.err.println(args[i]+" is not a valid dimension, ignored"); // :)
                }
                continue;
            }
            // -c is now RESERVED for "command" multicaster
//				if (args[i].equals("-c")|| args[i].equals("--canonical")) {
//					canonicalSS= true;
//					continue;
//				}
            //			if (args[i].equals("-c")|| args[i].equals("--cluster")) {
            //				readAheadLimit= Integer.parseInt(args[++i]);
            //				continue;
            //			}

            if (args[i].equalsIgnoreCase("-3pc")) {
                ASEvent.check3Pcomplete= true;
                continue;
            }
            if (args[i].equals("-tmp")) {
                System.setProperty(Constants.PROPERTY_TMPDIR, args[++i]);
                continue;
            }
            if (args[i].equals("-ds")) {
                SplicingGraph.retrieveDSEvents= false;
                continue;
            }
            if (args[i].equals("+ds")) {
                SplicingGraph.retrieveDSEvents= true;
                continue;
            }
            if (args[i].equals("-as")) {
                SplicingGraph.retrieveASEvents= false;
                continue;
            }
            if (args[i].equals("+as")) {
                SplicingGraph.retrieveASEvents= true;
                continue;
            }
            if (args[i].equalsIgnoreCase("-nmd")) {
                ASEvent.checkNMD= true;
                continue;
            }

            if (args[i].equals("+vs")) {
                SplicingGraph.retrieveVSEvents= true;
                continue;
            }
            if (args[i].equals("-vs")) {
                SplicingGraph.retrieveVSEvents= false;
                continue;
            }
            if (args[i].equals("-ext")) {
                onlyInternal= true;	// net false
                continue;
            }
            if (args[i].equals("+ext")) {
                onlyInternal= false;
                continue;
            }
            if (args[i].equalsIgnoreCase("-ic")|| args[i].equalsIgnoreCase("--intronConfidence")) {
                if (i+1== args.length)
                    System.err.println("You did not provide an intron confidence.");
                try {
                    SplicingGraph.intronConfidenceLevel= Byte.parseByte(args[i+1]);
                    ++i;	// ignore if missing
                    acceptableIntrons= true;
                } catch (NumberFormatException e) {
                    System.err.println("Intron confidence must be an integer value, you gave me "+args[i+1]);
                }
                continue;
            }

            // reactivated 20100112
            if (args[i].equalsIgnoreCase("-ec")|| args[i].equalsIgnoreCase("--edgeConfidence")) {
                if (i+1== args.length)
                    System.err.println("You did not provide an edge confidence.");
                try {
                    Transcript.setEdgeConfidenceLevel(Byte.parseByte(args[i + 1]));
                    ++i;	// ignore if missing
                } catch (NumberFormatException e) {
                    System.err.println("Exon confidence must be an integer value, you gave me "+args[i+1]);
                }
                continue;
            }

            if (args[i].equals("--flankType")) {
                ASEvent.setOutputFlankMode(true);
                continue;
            }


        }

        if (ASEvent.isOutputFlankMode()&& Graph.overrideSequenceDirPath== null) {
            System.err.println("[OOOPS] Parameters require genomic sequence. Use the option -g.");
            System.exit(-1);
        }

        // check for necessary components
        if (file== null|| !file.exists()) {
            System.err.print("Hey, you forgot to specify a valid input file! ");
            if (file!= null)
                System.err.println("Cannot find: "+file.getAbsolutePath());
            else
                System.err.println();
            System.exit(-1);
        }

        if ((SplicingGraph.canonicalSS|| acceptableIntrons)&& Graph.overrideSequenceDirPath== null) {
            System.err.println("You want me to check introns for valid/canonical splice sites, but you did not provide a valid sequence directory");
            System.exit(-1);
        }

        if (writerThread.outputFname== null&& (!SplicingGraph.writeStdOut)) {
            writerThread.outputFname= file.getAbsolutePath()+"_astalavista.gtf.gz";
        }
        if (writerThread.outputFname!= null&& new MyFile(writerThread.outputFname).exists()) {
            // Confirm o..+"\n by typing \'yes\':"
            System.err.println("Overwriting output file "+writerThread.outputFname+".");
//				try {
//					StringBuffer sb= new StringBuffer(3);
//					for (int i = System.in.read(); i != '\n';i= System.in.read()) {
//						sb.append((char) i);
//					}
//					sb.deleteCharAt(sb.length()-1);
//					if (sb.toString().trim().equalsIgnoreCase("yes")) {
//						System.err.println("AStalavista.");
//						System.exit(0);
//					}
//				} catch (Exception e) {
//					e.printStackTrace();
//					System.err.println("Output file exists, I give up.\nAStalavista.");
//					System.exit(0);
//				}
            while (new MyFile(writerThread.outputFname).exists())	// TODO strange access probs under win32
                try {
                    new MyFile(writerThread.outputFname).delete();
                    Thread.sleep(100);
                } catch (InterruptedException e) {
                    ; // e.printStackTrace();
                }
            System.err.println(writerThread.outputFname+ " deleted.");
            //System.out.println("File "+outputFname+" exists, check - I give up now.");
            //System.exit(0);
        }


        barna.io.gtf.GTFwrapper checkReader= new GTFwrapper(file.getAbsolutePath());
        //System.err.println("DEBUG -- temporarily deactivated file check");
        if (!checkReader.isApplicable()) {
            System.err.println("sorting input file, temporary directory "+System.getProperty(Constants.PROPERTY_TMPDIR));
            file= checkReader.sort();  // WAS: GFFreader.createSortedFile();
            System.err.println("Here is a sorted version of your file: "+file.getAbsolutePath());
        }


        return file;
    }

    static void outputStats(SplicingGraph.WriterThread writerThread, Writer writer) {
        BufferedWriter buffy= new BufferedWriter(writer);
        try {
            buffy.write("# started\t"+new Date(System.currentTimeMillis())+"\n");
            buffy.write("# input\t"+inputFile.getAbsolutePath()+"\n");
            buffy.write("# output");
            if (!SplicingGraph.writeStdOut)
                buffy.write("\t"+writerThread.outputFname+"\n");
            else
                buffy.write("\tstdout\n");
            if (Graph.overrideSequenceDirPath== null) {
                if (DEBUG)
                    buffy.write("# genome\t"+SplicingGraph.EventExtractorThread.species+"\n");
            } else
                buffy.write("# genome\t"+Graph.overrideSequenceDirPath+"\n");
            buffy.write("# dimension\t"+SplicingGraph.EventExtractorThread.n+"\n");
            buffy.write("# internalOnly\t"+ SplicingGraph.onlyInternal+ "\n");
            //buffy.write("# canonicalSS "+canonicalSS+"\n");
            //buffy.write("# acceptableIntrons "+acceptableIntrons+"\n");
            if (acceptableIntrons)
                buffy.write("# intronConfidenceLevel "+ SplicingGraph.intronConfidenceLevel+"\n");
            if (!onlyInternal)
                buffy.write("# edgeConfidenceLevel "+Transcript.getEdgeConfidenceLevel()+"\n");
            buffy.write("# as_events\t");
            if (SplicingGraph.retrieveASEvents)
                buffy.write("true");
            else
                buffy.write("false");
            buffy.write("\n");
            if (SplicingGraph.retrieveDSEvents) {
                buffy.write("# ds_events\t");
                if (SplicingGraph.retrieveDSEvents)
                    buffy.write("true");
                else
                    buffy.write("false");
            }
            buffy.write("\n");
            buffy.write("# vs_events\t");
            if (SplicingGraph.retrieveVSEvents)
                buffy.write("true");
            else
                buffy.write("false");
            buffy.write("\n");

            buffy.write("\n");
            buffy.flush();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
