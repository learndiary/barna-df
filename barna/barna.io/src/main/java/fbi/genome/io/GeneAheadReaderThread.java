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

import java.util.Date;

import fbi.genome.io.gtf.GTFwrapper;
import fbi.genome.model.Gene;
import fbi.genome.model.splicegraph.SplicingGraph;

/**
* @author Thasso Griebel (Thasso.Griebel@googlemail.com)
*/
public class GeneAheadReaderThread extends Thread {

        GTFwrapper reader;
        Gene[] g;
        SplicingGraph.EventExtractorThread downstreamThread;
        boolean output= false, output2= true, checkIntrons= true;

        public GeneAheadReaderThread(GTFwrapper newReader) {
            super();
            this.reader= newReader;
            reader.setSilent(true);
            setName("gene_ahead_reader_thread");
        }

        @Override
        public void run() {

            while (true) {
                long t0= System.currentTimeMillis();
                try {
                    reader.read();
                } catch (Exception e1) {
                    e1.printStackTrace();
                }
                if (reader.getUnclusteredGeneNb()== 0)
                    break;
                if (output2) {
                    Date ti= new Date(System.currentTimeMillis());
                    System.err.println("["+ti+"] read: "+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
                }
                t0= System.currentTimeMillis();
                g= reader.getGenes();
                if (g== null) {
                    System.err.println(" => no genes, stopping.");	// just in case
                }

//					if (output2) {
//						Date ti= new Date(System.currentTimeMillis());
//						System.err.println("["+ti+"] clustered: "+g[0].getChromosome()+" "+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
//					}

                    // quickly init splice sites, no vale la pena
//				t0= System.currentTimeMillis();
//				for (int i = 0; i < g.length; i++) {
//					if (g[i].getTranscriptCount()== 1)
//						continue;
//					for (int j = 0; j < g[i].getTranscripts().length; j++) {
//						SpliceSite[] ss= g[i].getTranscripts()[j].getSpliceSitesAll();
//						for (int x = 0; ss!= null&& x < ss.length; x++) {
//							ss[x].getTranscripts();
//						}
//					}
//				}
//				if (output2) {
//					Date ti= new Date(System.currentTimeMillis());
//					System.err.println("["+ti+"] ss init: "+g[0].getChromosome()+" "+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
//				}

            if (downstreamThread!= null&& downstreamThread.isAlive())
                    try {
                        downstreamThread.join();
                    } catch (InterruptedException e1) {
                        ; //:)
                    }

//					if (!reader.getGenes()[0].getChromosome().equals("chr1"))
//						break; 	// DEBUG !!!

                downstreamThread= new SplicingGraph.EventExtractorThread(this);
                downstreamThread.setOutput(output);
                downstreamThread.setOutput2(output2);
                downstreamThread.setG(g);
                downstreamThread.start();
                //downstreamThread.run();
                g= null;
                //break;
            }

        }

        public SplicingGraph.EventExtractorThread getDownstreamThread() {
            return downstreamThread;
        }

        public void setDownstreamThread(SplicingGraph.EventExtractorThread downstreamThread) {
            this.downstreamThread = downstreamThread;
        }

        public Gene[] getG() {
            return g;
        }

        public void setG(Gene[] g) {
            this.g = g;
        }

        public GTFwrapper getReader() {
            return reader;
        }

        public void setReader(GTFwrapper reader) {
            this.reader = reader;
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
    }
