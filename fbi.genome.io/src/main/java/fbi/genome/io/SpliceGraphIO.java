/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package fbi.genome.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.Writer;
import java.util.Date;

import fbi.genome.io.gtf.GTFwrapper;
import fbi.genome.model.ASEvent;
import fbi.genome.model.Gene;
import fbi.genome.model.Graph;
import fbi.genome.model.IntronModel;
import fbi.genome.model.Species;
import fbi.genome.model.Transcript;
import fbi.genome.model.commons.MyFile;
import fbi.genome.model.constants.Constants;
import fbi.genome.model.splicegraph.SpliceGraph;

/**
 * IO methods to create a spliece graph
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class SpliceGraphIO {
    public static boolean DEBUG= false;
    private String outputFname= null;
    private MyFile inputFile= null;
    private Graph graph;


    public static void extractSpliceJunctions(int eFlankDon, int eFlankAcc, IntronModel iModel, File inF, File outF) {
        if (outF!= null&& outF.exists())
            outF.delete();
        GTFwrapper reader= new GTFwrapper(inF.getAbsolutePath());
        Gene[] g= null;
        try {
            for (reader.read(); (g= reader.getGenes())!= null; reader.read()) {
                for (int i = 0; i < g.length; i++) {
                	SpliceGraph gr= new SpliceGraph(g[i]);
                    gr.constructGraph();
                    gr.writeSpliceJunctionSeqs(eFlankDon, eFlankAcc, iModel, outF);
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public SpliceGraphIO(){

    }


    public void outputStats(Writer writer) {
        BufferedWriter buffy= new BufferedWriter(writer);
        try {
            buffy.write("# started\t"+new Date(System.currentTimeMillis())+"\n");
            buffy.write("# input\t"+inputFile.getAbsolutePath()+"\n");
            buffy.write("# output");
            if (!SpliceGraph.writeStdOut)
                buffy.write("\t"+outputFname+"\n");
            else
                buffy.write("\tstdout\n");
            if (fbi.genome.model.Graph.overrideSequenceDirPath== null) {
                if (DEBUG)
                    buffy.write("# genome\t"+ SpliceGraph.EventExtractorThread.species+"\n");
            } else
                buffy.write("# genome\t"+ fbi.genome.model.Graph.overrideSequenceDirPath+"\n");
            buffy.write("# dimension\t"+ SpliceGraph.EventExtractorThread.n+"\n");
            buffy.write("# internalOnly\t"+ SpliceGraph.onlyInternal +"\n");
            //buffy.write("# canonicalSS "+canonicalSS+"\n");
            //buffy.write("# acceptableIntrons "+acceptableIntrons+"\n");
            if (SpliceGraph.acceptableIntrons)
                buffy.write("# intronConfidenceLevel "+SpliceGraph.intronConfidenceLevel+"\n");
            if (!SpliceGraph.onlyInternal)
                buffy.write("# edgeConfidenceLevel "+ Transcript.getEdgeConfidenceLevel()+"\n");
            buffy.write("# as_events\t");
            if (SpliceGraph.retrieveASEvents)
                buffy.write("true");
            else
                buffy.write("false");
            buffy.write("\n");
            if (SpliceGraph.retrieveDSEvents) {
                buffy.write("# ds_events\t");
                if (SpliceGraph.retrieveDSEvents)
                    buffy.write("true");
                else
                    buffy.write("false");
            }
            buffy.write("\n");
            buffy.write("# vs_events\t");
            if (SpliceGraph.retrieveVSEvents)
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


    public java.io.File parseArguments(String[] args) {

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
                    SpliceGraph.outputSeq= true;
                    continue;
                }
                if (args[i].equals("-o")|| args[i].equals("--output")) {
                    String s= args[++i];
                    if (s.equalsIgnoreCase("stdout"))
                        SpliceGraph.writeStdOut= true;
                    else
                        outputFname= s+".gz";
                    continue;
                }
                if (args[i].equals("-g")|| args[i].equals("--genome")) {
                    if (i+1>= args.length) {
                        System.err.println("You forgot to specify the genome!");
                        System.exit(-1);
                    }
                    MyFile checkFile= new MyFile(args[i+1]);
                    if (checkFile.exists()&& checkFile.isDirectory())
                        fbi.genome.model.Graph.overrideSequenceDirPath= checkFile.getAbsolutePath();
                    else {
                        String[] s= args[++i].split("_");
                        if (s.length!= 2) {
                            System.err.println("Invalid genome directory: "+args[i+1]);
                            System.exit(-1);
                        }
                        Species spe= new Species(s[0]);
                        spe.setGenomeVersion(s[1]);
                        SpliceGraph.EventExtractorThread.setSpecies(spe);
                    }
                    SpliceGraph.acceptableIntrons= true;
                    continue;
                }

                if (args[i].equals("-k")|| args[i].equals("--dimension")) {
                    try {
                        int x= Integer.parseInt(args[++i]);
                        if (x< -1|| ((x> -1)&& (x< 2)))
                            System.err.println(args[i]+" is not a valid dimension, ignored");
                        SpliceGraph.EventExtractorThread.n= x;
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
                    SpliceGraph.retrieveDSEvents= false;
                    continue;
                }
                if (args[i].equals("+ds")) {
                    SpliceGraph.retrieveDSEvents= true;
                    continue;
                }
                if (args[i].equals("-as")) {
                    SpliceGraph.retrieveASEvents= false;
                    continue;
                }
                if (args[i].equals("+as")) {
                    SpliceGraph.retrieveASEvents= true;
                    continue;
                }
                if (args[i].equalsIgnoreCase("-nmd")) {
                    ASEvent.checkNMD= true;
                    continue;
                }

                if (args[i].equals("+vs")) {
                    SpliceGraph.retrieveVSEvents= true;
                    continue;
                }
                if (args[i].equals("-vs")) {
                    SpliceGraph.retrieveVSEvents= false;
                    continue;
                }
                if (args[i].equals("-ext")) {
                    SpliceGraph.onlyInternal= true;	// net false
                    continue;
                }
                if (args[i].equals("+ext")) {
                    SpliceGraph.onlyInternal= false;
                    continue;
                }
                if (args[i].equalsIgnoreCase("-ic")|| args[i].equalsIgnoreCase("--intronConfidence")) {
                    if (i+1== args.length)
                        System.err.println("You did not provide an intron confidence.");
                    try {
                        SpliceGraph.intronConfidenceLevel= Byte.parseByte(args[i+1]);
                        ++i;	// ignore if missing
                        SpliceGraph.acceptableIntrons= true;
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
                        Transcript.setEdgeConfidenceLevel(Byte.parseByte(args[i+1]));
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

            if (ASEvent.isOutputFlankMode()&& fbi.genome.model.Graph.overrideSequenceDirPath== null) {
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

            if ((SpliceGraph.canonicalSS|| SpliceGraph.acceptableIntrons)&& fbi.genome.model.Graph.overrideSequenceDirPath== null) {
                System.err.println("You want me to check introns for valid/canonical splice sites, but you did not provide a valid sequence directory");
                System.exit(-1);
            }

            if (outputFname== null&& (!SpliceGraph.writeStdOut)) {
                outputFname= file.getAbsolutePath()+"_astalavista.gtf.gz";
            }
            if (outputFname!= null&& new MyFile(outputFname).exists()) {
                // Confirm o..+"\n by typing \'yes\':"
                System.err.println("Overwriting output file "+outputFname+".");
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
                while (new MyFile(outputFname).exists())	// TODO strange access probs under win32
                    try {
                        new MyFile(outputFname).delete();
                        Thread.sleep(100);
                    } catch (InterruptedException e) {
                        ; // e.printStackTrace();
                    }
                System.err.println(outputFname+ " deleted.");
                //System.out.println("File "+outputFname+" exists, check - I give up now.");
                //System.exit(0);
            }


            GTFwrapper checkReader= new GTFwrapper(file.getAbsolutePath());
            //System.err.println("DEBUG -- temporarily deactivated file check");
            if (!checkReader.isApplicable()) {
                System.err.println("sorting input file, temporary directory "+System.getProperty(Constants.PROPERTY_TMPDIR));
                file= checkReader.sort();
                System.err.println("Here is a sorted version of your file: "+file.getAbsolutePath());
            }


            return file;
        }


}
