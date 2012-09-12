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

package barna.io;

import barna.io.gtf.GTFwrapper;
import barna.model.*;
import barna.model.commons.MyFile;
import barna.model.constants.Constants;
import barna.model.splicegraph.SplicingGraph;

import java.io.BufferedWriter;
import java.io.File;
import java.io.Writer;
import java.util.Date;

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
                	SplicingGraph gr= new SplicingGraph(g[i]);
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
            buffy.write("# started\t"+new Date(System.currentTimeMillis())+barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# input\t"+inputFile.getAbsolutePath()+barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# output");
            if (!SplicingGraph.writeStdOut)
                buffy.write("\t"+outputFname+barna.commons.system.OSChecker.NEW_LINE);
            else
                buffy.write("\tstdout\n");
            if (barna.model.Graph.overrideSequenceDirPath== null) {
                if (DEBUG)
                    buffy.write("# genome\t"+ SplicingGraph.EventExtractorThread.species+barna.commons.system.OSChecker.NEW_LINE);
            } else
                buffy.write("# genome\t"+ barna.model.Graph.overrideSequenceDirPath+barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# dimension\t"+ SplicingGraph.EventExtractorThread.n+barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# internalOnly\t"+ SplicingGraph.onlyInternal +barna.commons.system.OSChecker.NEW_LINE);
            //buffy.write("# canonicalSS "+canonicalSS+barna.commons.system.OSChecker.NEW_LINE);
            //buffy.write("# acceptableIntrons "+acceptableIntrons+barna.commons.system.OSChecker.NEW_LINE);
            if (SplicingGraph.acceptableIntrons)
                buffy.write("# intronConfidenceLevel "+SplicingGraph.intronConfidenceLevel+barna.commons.system.OSChecker.NEW_LINE);
            if (!SplicingGraph.onlyInternal)
                buffy.write("# edgeConfidenceLevel "+ Transcript.getEdgeConfidenceLevel()+barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# as_events\t");
            if (SplicingGraph.retrieveASEvents)
                buffy.write("true");
            else
                buffy.write("false");
            buffy.write(barna.commons.system.OSChecker.NEW_LINE);
            if (SplicingGraph.retrieveDSEvents) {
                buffy.write("# ds_events\t");
                if (SplicingGraph.retrieveDSEvents)
                    buffy.write("true");
                else
                    buffy.write("false");
            }
            buffy.write(barna.commons.system.OSChecker.NEW_LINE);
            buffy.write("# vs_events\t");
            if (SplicingGraph.retrieveVSEvents)
                buffy.write("true");
            else
                buffy.write("false");
            buffy.write(barna.commons.system.OSChecker.NEW_LINE);

            buffy.write(barna.commons.system.OSChecker.NEW_LINE);
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
                        barna.commons.system.OSChecker.NEW_LINE);
                System.err.println("-o, --output <output file|\'stdout\'>");
                System.err.println("Optional, the name of the output file (fully qualified path) OR the keyword \'stdout\' for " +
                        "writing everything to the standard output stream. " +
                        "If nothing is specified, the output will be written to a file \'<input file>_astalavista.gtf.gz\'. " +
                        barna.commons.system.OSChecker.NEW_LINE);
                System.err.println("-g, --genome <path to directory>");
                System.err.println("Path to the directory containing sequence files corresponding to the <seqname> field " +
                        "in the input GTF. A genome directory is required if a intron confidence value is specified." +
                        barna.commons.system.OSChecker.NEW_LINE);
                System.err.println("-k, --dimension <int value>");
                System.err.println("Dimension >1 of the events to be extracted. Default is 2 (i.e., \'pairwise events\'). " +
                        barna.commons.system.OSChecker.NEW_LINE);
                System.err.println("-tmp");
                System.err.println("Set temporary directory" +
                        barna.commons.system.OSChecker.NEW_LINE);
                System.err.println("-ext, +ext");
                System.err.println("(De-)activate external events, i.e. events that include the transcript start or the " +
                        "poly-adenylation site" +
                        barna.commons.system.OSChecker.NEW_LINE);
                System.err.println("-ic, --intronConfidence [int value]");
                System.err.println("Level of intron confidence. The default is to trust all introns. Introns are assigned " +
                        "a confidency class:\n" +
                        "\t 0 if 'RefSeq' appears in the source field of the annotation\n" +
                        "\t 1 if 'mRNA' appears in the source field of the annotation\n" +
                        "\t 2 if 'EST' appears in the source field of the annotation\n" +
                        "\t 3 if if none of the above applies\n" +
                        "all introns of confidency level > intronConfidence will be checked for proper splice sites when extracting events." +
                        barna.commons.system.OSChecker.NEW_LINE);
                System.err.println("-as, +as");
                System.err.println("Deactivate (\'-as\') or activate (\'+as\') the retrieval of Alternative Splicing events. See documentation " +
                        "for the definition of events that suffice an alternative splicing event." +
                        barna.commons.system.OSChecker.NEW_LINE);
                System.err.println("-ds, +ds");
                System.err.println("Deactivate (\'-ds\') or activate (\'+ds\') the retrieval of aDditional splicing events. See documentation " +
                        "for the definition of events that suffice an additional splicing event." +
                        barna.commons.system.OSChecker.NEW_LINE);
                System.err.println("-s, --seqsite");
                System.err.println("Output splice site sequences with events. Requires a reference genome."+
                        barna.commons.system.OSChecker.NEW_LINE);
                System.err.println("--flankType");
                System.err.println("Output the type of the event flanks, i.e., \'constitutive\' or \'alternative\'."+
                        barna.commons.system.OSChecker.NEW_LINE);

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
                        barna.commons.system.OSChecker.NEW_LINE);
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
                        barna.model.Graph.overrideSequenceDirPath= checkFile.getAbsolutePath();
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
                    SplicingGraph.acceptableIntrons= true;
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
                    SplicingGraph.onlyInternal= true;	// net false
                    continue;
                }
                if (args[i].equals("+ext")) {
                    SplicingGraph.onlyInternal= false;
                    continue;
                }
                if (args[i].equalsIgnoreCase("-ic")|| args[i].equalsIgnoreCase("--intronConfidence")) {
                    if (i+1== args.length)
                        System.err.println("You did not provide an intron confidence.");
                    try {
                        SplicingGraph.intronConfidenceLevel= Byte.parseByte(args[i+1]);
                        ++i;	// ignore if missing
                        SplicingGraph.acceptableIntrons= true;
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

            if (ASEvent.isOutputFlankMode()&& barna.model.Graph.overrideSequenceDirPath== null) {
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

            if ((SplicingGraph.canonicalSS|| SplicingGraph.acceptableIntrons)&& barna.model.Graph.overrideSequenceDirPath== null) {
                System.err.println("You want me to check introns for valid/canonical splice sites, but you did not provide a valid sequence directory");
                System.exit(-1);
            }

            if (outputFname== null&& (!SplicingGraph.writeStdOut)) {
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
