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

package barna.flux.simulator.tools;

import barna.commons.Execute;
import barna.commons.log.Log;
import barna.io.BEDMappingIterator;
import barna.io.FileHelper;
import barna.io.MSIterator;
import barna.io.bed.BEDFileReader;
import barna.io.bed.BEDFileReader;
import barna.io.gtf.GTFwrapper;
import barna.io.state.MappingReaderState;
import barna.model.Gene;
import barna.model.Graph;
import barna.model.Mapping;
import barna.model.Transcript;
import barna.model.bed.BEDMapping;
import barna.model.bed.BEDobject2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Collection;
import java.util.Collections;

/**
 * Count and print breakpoint distribution
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
//@Cli(name = "extractPWM", description = "Extract PWM from GTF annotation and BED file")
public class PWMExtractor {  //implements FluxTool {

	public static void main(String[] args) {
		
		PWMExtractor ex= new PWMExtractor();
		Graph.overrideSequenceDirPath= "/Users/micha/genomes/hg19";
		ex.setBedFile(new File("/Users/micha/tmp/simulator/FRT_STD/STD.bed"));
		ex.setGtfFile(new File("/Users/micha/annotation/hg19/hg19_RefSeq_fromUCSC100615_sorted.clean.gtf"));
		ex.setOutFile(new File("/Users/micha/tmp/simulator/FRT_STD/STD.pwm"));
		
		Execute.initialize(2);
		//Future f= Execute.getExecutor().submit(ex);
		try {
			ex.call();
		} catch (Exception e) {
			e.printStackTrace();
		}
		Execute.shutdown();
	}
	
    private File gtfFile;
    private File bedFile;
    private File outFile;


    public File getGtfFile() {
        return gtfFile;
    }

//    @Option(name = "g", longName = "gtf", description = "The GTF file", required = true)
    public void setGtfFile(final File gtfFile) {
        this.gtfFile = gtfFile;
    }

    public File getBedFile() {
        return bedFile;
    }

//    @Option(name = "b", longName = "bed", description = "The .bed file", required = true)
    public void setBedFile(final File bedFile) {
        this.bedFile = bedFile;
    }

    public File getOutFile() {
        return outFile;
    }

//    @Option(name = "o", longName = "out", description = "The output file", required = true)
    public void setOutFile(final File outFile) {
        this.outFile = outFile;
    }

//    //@Override
//    public boolean validateParameters(final HelpPrinter printer, final ArgumentProcessor toolArguments) {
//        return true;
//    }

    //@Override
    public Object call() throws Exception {
        int flank5 = 10, flank3 = 20;
        GTFwrapper anoReader = new GTFwrapper(gtfFile.getAbsolutePath());
        if (!anoReader.isApplicable()) {
            Log.message("\tsorting GTF file");
            File f = anoReader.sort();
            Log.message("\tsorted file in " + f.getAbsolutePath());
            gtfFile = f;
            anoReader = new GTFwrapper(f.getAbsolutePath());
        }
        Log.message("");

        File ff = new File(bedFile.getAbsolutePath() + "_sorted");
        BEDFileReader bedReader = null;
        if (ff.exists()) {
            Log.message("\tusing sorted file " + ff.getName());
            bedFile = ff;
            bedReader = new BEDFileReader(bedFile.getAbsolutePath());
        } else {
            bedReader = new BEDFileReader(bedFile.getAbsolutePath());
            if (!bedReader.isApplicable()) {
                Log.message("\tsorting BED file");

                File f = FileHelper.createTempFile("bed-sort", "bed");
                bedReader.sort(f);
                if (FileHelper.move(f, ff, null)) {
                    bedFile = ff;
                } else {
                    bedFile = f;
                }
                Log.message("\tsorted file in " + bedFile.getAbsolutePath());
                bedReader = new BEDFileReader(bedFile.getAbsolutePath());
            }
        }

        Log.message("\tprocessing ");

        int cntTrpt = 0, cntReads = 0;
        int[][] sense = new int[flank5 + flank3][], asense = new int[flank5 + flank3][];
        for (int i = 0; i < asense.length; i++) {
            sense[i] = new int[4];
            asense[i] = new int[4];
            for (int j = 0; j < asense[i].length; j++) {
                sense[i][j] = 0;
                asense[i][j] = 0;
            }
        }
        Gene[] genes = null;
        for (anoReader.read(); (genes = anoReader.getGenes()) != null; anoReader.read()) {
            for (int i = 0; i < genes.length; i++) {
                if (genes[i].getTranscriptCount() > 1)    // non-AS genes
                {
                    continue;
                }
                ++cntTrpt;
                Transcript t = genes[i].getTranscripts()[0];
                
                /*MappingReaderState state= bedReader.read(t.getChromosome(), t.getStart(), t.getEnd());
                BEDobject2[] beds = (BEDobject2[]) state.result;
                if (beds == null) {
                    continue;
                } */
                String s = t.getSplicedSequence().toUpperCase();
                for (MSIterator<Mapping> mappingIterator = bedReader.read(t.getChromosome(), t.getStart(), t.getEnd());mappingIterator.hasNext();) {
                //for (int j = 0; j < beds.length; j++) {
                    // get t-coordinates
                    Mapping mapping = mappingIterator.next();
                    int tstart = t.getExonicPosition(mapping.getStart() + 1),
                            tend = t.getExonicPosition(mapping.getEnd());    // t-coordinates, 0-based
                    if (tstart < 0 || tstart >= s.length()) {
                        continue;
                    }

                    // count on subsequence
                    ++cntReads;
                    boolean sens = mapping.getStrand() == t.getStrand();
                    int[][] a = sens ? sense : asense;
                    if (sens) {
                        for (int k = 0; k < a.length; ++k) {
                            int p = tstart - flank5 + k;
                            if (p < 0) {
                                continue;
                            }
                            if (p >= s.length()) {
                                break;
                            }
                            if (s.charAt(p) == 'A') {
                                ++a[k][0];
                            } else if (s.charAt(p) == 'C') {
                                ++a[k][1];
                            } else if (s.charAt(p) == 'G') {
                                ++a[k][2];
                            } else if (s.charAt(p) == 'T') {
                                ++a[k][3];
                            }
                        }

                    } else {    // read asense to s
                        for (int k = 0; k < a.length; ++k) {
                            int p = tend + flank5 - k;
                            if (p >= s.length()) {
                                continue;
                            }
                            if (p < 0) {
                                break;
                            }
                            // reverse complement?
/*								if (s.charAt(p)== 'A')
                                ++a[k][3];	// A -> count T
                            else if (s.charAt(p)== 'C')
                                ++a[k][2];
                            else if (s.charAt(p)== 'G')
                                ++a[k][1];
                            else if (s.charAt(p)== 'T')
                                ++a[k][0];
*/
                            if (s.charAt(p) == 'A') {
                                ++a[k][0];
                            } else if (s.charAt(p) == 'C') {
                                ++a[k][1];
                            } else if (s.charAt(p) == 'G') {
                                ++a[k][2];
                            } else if (s.charAt(p) == 'T') {
                                ++a[k][3];
                            }
                        }
                    }

                }
            }
        }
        Log.message("\tOK");
        Log.message("\tFound " + cntTrpt + " non-AS tx with " + cntReads + " reads.");
        Log.message("");

        // output
        BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
        for (int i = 0; i < sense.length; i++) {
            int pos = (i >= flank5 ? i - flank5 + 1 : i - flank5);
            writer.write(pos + "\t" + sense[i][0] + "\t" + sense[i][1] + "\t" + sense[i][2] + "\t" + sense[i][3] + "\n");
        }
        writer.flush();
        writer.close();


//        writer= new BufferedWriter(new FileWriter(fileBed+"_asense.pwm"));
//        for (int i = 0; i < asense.length; i++) {
//            int pos= (i>= flank5? i- flank5+ 1: i- flank5);
//            writer.write(pos+ "\t"+ asense[i][0]+ "\t"+ asense[i][1]+ "\t"+ asense[i][2]+ "\t"+ asense[i][3]+ "\n");
//        }
//        writer.flush();
//        writer.close();
//        System.err.println("wrote "+ fileBed+"_sense.pwm, and "+ fileBed+"_asense.pwm.");
        return null;
    }
}
