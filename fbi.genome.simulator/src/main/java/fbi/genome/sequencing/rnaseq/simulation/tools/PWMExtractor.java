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

package fbi.genome.sequencing.rnaseq.simulation.tools;

import fbi.commons.Log;
import fbi.commons.file.FileHelper;
import fbi.commons.flux.FluxTool;
import fbi.commons.options.HelpPrinter;
import fbi.genome.io.bed.BEDwrapper;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.Gene;
import fbi.genome.model.Transcript;
import fbi.genome.model.bed.BEDobject2;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

/**
 * Count and print breakpoint distribution
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
//@Cli(name = "extractPWM", description = "Extract PWM from GTF annotation and BED file")
public class PWMExtractor {  //implements FluxTool {

    private File gtfFile;
    private File bedFile;
    private File outFile;


    public File getGtfFile() {
        return gtfFile;
    }

    @Option(name = "g", longName = "gtf", description = "The GTF file", required = true)
    public void setGtfFile(final File gtfFile) {
        this.gtfFile = gtfFile;
    }

    public File getBedFile() {
        return bedFile;
    }

    @Option(name = "b", longName = "bed", description = "The .bed file", required = true)
    public void setBedFile(final File bedFile) {
        this.bedFile = bedFile;
    }

    public File getOutFile() {
        return outFile;
    }

    @Option(name = "o", longName = "out", description = "The output file", required = true)
    public void setOutFile(final File outFile) {
        this.outFile = outFile;
    }

    //@Override
    public boolean validateParameters(final HelpPrinter printer, final ArgumentProcessor toolArguments) {
        return true;
    }

    //@Override
    public Object call() throws Exception {
        int flank5 = 10, flank3 = 20;
        GFFReader anoReader = new GFFReader(gtfFile.getAbsolutePath());
        if (!anoReader.isApplicable()) {
            Log.message("\tsorting GTF file");
            File f = anoReader.createSortedFile();
            Log.message("\tsorted file in " + f.getAbsolutePath());
            gtfFile = f;
            anoReader = new GFFReader(f.getAbsolutePath());
        }
        Log.message("");

        File ff = new File(bedFile.getAbsolutePath() + "_sorted");
        BEDwrapper bedReader = null;
        if (ff.exists()) {
            Log.message("\tusing sorted file " + ff.getName());
            bedFile = ff;
            bedReader = new BEDwrapper(bedFile.getAbsolutePath());
        } else {
            bedReader = new BEDwrapper(bedFile.getAbsolutePath());
            if (!bedReader.isApplicable()) {
                Log.message("\tsorting BED file");
                File f = bedReader.sortBED(bedFile);
                if (FileHelper.move(f, ff, null)) {
                    bedFile = ff;
                } else {
                    bedFile = f;
                }
                Log.message("\tsorted file in " + bedFile.getAbsolutePath());
                bedReader = new BEDwrapper(bedFile.getAbsolutePath());
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
                BEDobject2[] beds = bedReader.read(t.getChromosome(), t.getStart(), t.getEnd());
                if (beds == null) {
                    continue;
                }
                String s = t.getSplicedSequence().toUpperCase();
                for (int j = 0; j < beds.length; j++) {
                    // get t-coordinates
                    int tstart = t.getExonicPosition(beds[j].getStart() + 1),
                            tend = t.getExonicPosition(beds[j].getEnd());    // t-coordinates, 0-based
                    if (tstart < 0 || tstart >= s.length()) {
                        continue;
                    }

                    // count on subsequence
                    ++cntReads;
                    boolean sens = beds[j].getStrand() == t.getStrand();
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
