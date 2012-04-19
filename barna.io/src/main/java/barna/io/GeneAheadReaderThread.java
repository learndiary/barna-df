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
import barna.model.Gene;
import barna.model.splicegraph.SplicingGraph;

import java.util.Date;

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
