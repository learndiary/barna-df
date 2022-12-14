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

package barna.flux.capacitor.reconstruction;

import barna.commons.Execute;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.CommandLine;
import barna.commons.launcher.Tool;
import barna.commons.log.Log;
import barna.commons.parameters.ParameterException;
import barna.commons.system.OSChecker;
import barna.commons.utils.StringUtils;
import barna.flux.capacitor.graph.AnnotationMapper;
import barna.flux.capacitor.graph.MappingsInterface;
import barna.flux.capacitor.profile.BiasProfiler;
import barna.flux.capacitor.profile.MappingStats;
import barna.flux.capacitor.profile.Profile;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping;
import barna.genome.lpsolver.LPSolverLoader;
import barna.io.*;
import barna.io.bed.BEDReader;
import barna.io.gtf.GTFwrapper;
import barna.model.rna.UniversalReadDescriptor;
import barna.io.sam.SAMReader;
import barna.model.*;
import barna.model.commons.MyFile;
import barna.model.constants.Constants;
import barna.model.gff.GFFObject;
import barna.model.splicegraph.AbstractEdge;
import barna.model.splicegraph.SimpleEdge;
import barna.model.splicegraph.SplicingGraph;
import barna.model.splicegraph.SuperEdge;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import lpsolve.LpSolve;
import lpsolve.VersionInfo;
import net.sf.samtools.SAMFileReader;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;


/**
 * The Flux Capacitor class performes a deconvolution of reads falling into common areas of transcripts.
 *
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class FluxCapacitor implements Tool<MappingStats>, ReadStatCalculator {

    // global debug switch
    public static boolean DEBUG= false;

    /**
     * Enumerates possible tasks for the FluxCapacitor
     */
    private enum Task {
        COUNT_INTRONS, COUNT_SJ, PROFILE, DECOMPOSE
    };

    /**
     * Task to be executed in the current run - initialized empty
     */
    private EnumSet<Task> currentTasks = EnumSet.noneOf(Task.class);

    /**
     * Possible features for the output
     */
    private enum OutputFlag {
        BALANCED, OBSERVATION, PREDICTION, GENE, TRANSCRIPT, EXON, SJ, EVENT, UNKNOWN, LP
    };

    private EnumSet<OutputFlag> output = EnumSet.of(OutputFlag.BALANCED, OutputFlag.TRANSCRIPT);

    /**
     * Stores a reference to the parsed command line arguments
     * to update settings
     */
    private JSAPResult commandLineArgs;

    /**
     * Variable to store the stats about input data
     */
    private MappingStats stats;

    /**
     * Thread to parallelize annotation reading.
     */
    class GtfReaderThread extends Thread {

        public GtfReaderThread() {
            super("GTF_reader");
        }

        @Override
        public void run() {
            try {
                getWrapperGTF().read();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }


    /**
     * A class that encapsulates all information necessary to carry out the deconvolution
     * of the reads in a locus.
     */
    static class LocusSolver implements Callable<MappingStats> {

        /**
         * The locus that is to be solved.
         */
        private Gene gene = null;

        /**
         * AStalavista events found in the locus.
         */
        private ASEvent[] events = null;

        /**
         * Iterator over mappings.
         */
		private MSIterator mappings = null;

        /**
         * EnumSet indicating which task(s) has(have) to be performed in the current run.
         */
        private EnumSet<Task> tasks;

        /**
         * EnumSet indicating which feature(s) has(ve) to be put on the ouput.
         */
        private EnumSet<OutputFlag> output;

        /**
         * Whether reads are paired
         */
        private boolean pairedEnd;

        /**
         * Whether reads are stranded
         */
        private boolean stranded;

        /**
         * Used for chaining threads to be executed sequentially.
         */
        private Thread threadBefore = null;

        /**
         * The number of annotation-mapped mappings or reads.
         */
        private int nrMappingsReadsOrPairs;

        /**
         * Variable to store invariant of observed split frequency.
         */
        private float invariantTestObsSplitFreq = 0;

        /**
         * Variable to store invariant of predicted split frequency.
         */
        private float invariantTestPredSplitFreq = 0;

        /**
         * Variable to exchange stats
         */
        private MappingStats stats = null;

        /**
         * Settings for the current FluxCapacitor run
         */
        private FluxCapacitorSettings settings;

        /**
         * Reads profile
         */
        private Profile profile;

        /**
         * Cost model for deconvolution.
         */
        private byte costModel = GraphLPsolver.COSTS_LINEAR;

        /**
         * Granularity of deconvolution cost function.
         */
        private byte costSplit = 1;

        /**
         * Lower and upper bound of how much of the original observation can be substracted respectively added.
         */
        private float[] costBounds = new float[]{0.95f, Float.NaN};

        /**
         * Constructor providing reads and mappings for deconvolution.
         * The mode of the run can be switched between profiling and deconvolution.
         *
         * @param newGene   the locus model
         * @param newMappings the mappings that fall in the locus
         * @param tasks     tasks to be preformed
         */
		public LocusSolver(Gene newGene, MSIterator newMappings, EnumSet tasks, EnumSet<OutputFlag> output, boolean pairedEnd, boolean stranded, FluxCapacitorSettings settings, Profile profile) {

            this.gene = newGene;
			this.mappings = newMappings;
            this.tasks = tasks;
            this.output = output;
            this.pairedEnd = pairedEnd;
            this.stranded = stranded;
            this.settings = settings;
            this.profile = profile;
            this.stats = new MappingStats();
            this.stats.add(profile.getMappingStats());
            this.stats.reset();

            nrMappingsReadsOrPairs = 0;
        }


        /**
         * Launches the profiling/deconvolution on the locus
         */
        /*@Override
        public void run() {

            AnnotationMapper mapper = null;

            if(tasks.isEmpty())
                return;

                mapper = new AnnotationMapper(this.gene);
                mapper.map(this.mappings, settings);

                nrReadsLoci += mapper.nrMappingsLocus;
                nrReadsMapped += mapper.getNrMappingsMapped();
                nrMappingsReadsOrPairs += mapper.getNrMappingsMapped() / 2;
                nrPairsNoTxEvidence += mapper.getNrMappingsNotMappedAsPair();
                nrPairsWrongOrientation += mapper.getNrMappingsWrongPairOrientation();
            }


            //TODO move this on call()
            //Execute tasks
            for (Task t : this.tasks) {
                switch (t) {
                    case COUNT_INTRONS:
                        outputIntronsGFF(this.gene, mapper);
                        break;
                    case COUNT_SJ:
                        outputSJGFF(this.gene, mapper);
                        break;
                    case DECOMPOSE:
                        GraphLPsolver mySolver = null;
                        if (mapper.nrMappingsMapped > 0 && this.gene.getTranscriptCount() > 1) {    // OPTIMIZE
                            mySolver = getSolver(mapper, (int) (mapper.nrMappingsMapped * 2)); // not: getMappedReadcount()
                            mySolver.run();
                        }
                        outputGFF(mapper, events, mySolver);
                        break;
                }
            }

//            if (decompose) {
//
//                // BUG 110301: do not use mapTrivial, problems with split-maps
////				if (this.gene.getTranscriptCount()== 1) {
////					mapTrivial(gene.getTranscripts()[0], beds);
////					outputGFF(null, null, null);
////				} else {
//
//                GraphLPsolver mySolver = null;
//                if (mapper.nrMappingsMapped > 0 && this.gene.getTranscriptCount() > 1) {    // OPTIMIZE
//                    mySolver = getSolver(mapper, (int) (mapper.nrMappingsMapped * 2)); // not: getMappedReadcount()
//                    mySolver.run();
//                }
//                outputGFF(mapper, events, mySolver);
////					}
//
//            } else {
//                // map all reads
//                if (this.gene.getTranscriptCount() == 1) {
//                    ++nrSingleTranscriptLearn;
//                    learn(this.gene.getTranscripts()[0], beds);
//                }
//            }
            if (mappings != null) {
                mappings.clear();
                mappings = null;
            }
            gene = null;
            // makes it terribly slow
            //System.gc();

//			synchronized(FluxCapacitor.this.threadPool) {
            //FluxCapacitor.this.threadPool.remove(this);
//			}
        }*/

        /**
         * Count reads to splice junction within the current locus and output them in GTF format.
         *
         * @param gene current locus
         * @param mapper mapping graph for the current locus
         */
        private void outputSJGFF(Gene gene, AnnotationMapper mapper) {
            Map<String, Integer> m = mapper.getSJReads(settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).isPaired() ? true : false);
            StringBuilder sb = new StringBuilder();
            for (String s : m.keySet()) {
                String[] junction = s.split("\\^");
                sb.append(gene.getChromosome());
                sb.append("\t");
                sb.append("flux");
                sb.append("\t");
                sb.append(FluxCapacitorConstants.GFF_FEATURE_JUNCTION);
                sb.append("\t");
                sb.append(junction[1].contains("-") ? junction[2].replace("-", "") : junction[1]);
                sb.append("\t");
                sb.append(junction[2].contains("-")?junction[1].replace("-",""):junction[2]);
                sb.append("\t");
                sb.append(".");
                sb.append("\t");
                sb.append(gene.getStrand() > 0 ? "+" : "-");
                sb.append("\t");
                sb.append(".");
                sb.append("\t");
                sb.append("gene_id \""+junction[0]+"\";");
                sb.append(" ");
                sb.append("locus_id \""+gene.getLocusID()+"\";");
                sb.append(" ");
                sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_READS+" "+String.format("%1$f", (float) m.get(s)));// +";");
                sb.append(barna.commons.system.OSChecker.NEW_LINE);
            }
            Log.print(sb.toString());
        }

        /**
         * Count reads to all-intronic regions within the current locus and output them in GTF format.
         *
         * @param gene current locus
         * @param mapper mapping graph for the current locus
         */
        private void outputIntronsGFF(Gene gene, AnnotationMapper mapper) {
            Map<String, Float[]> m = mapper.getAllIntronicReads(settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).isPaired() ? true : false);
            StringBuilder sb = new StringBuilder();
            for (String s : m.keySet()) {
                String[] intron = s.split("\\^");
                sb.append(gene.getChromosome());
                sb.append("\t");
                sb.append("flux");
                sb.append("\t");
                sb.append(FluxCapacitorConstants.GFF_FEATURE_INTRON);
                sb.append("\t");
                sb.append(intron[1].contains("-")?intron[2].replace("-",""):intron[1]);
                sb.append("\t");
                sb.append(intron[2].contains("-")?intron[1].replace("-",""):intron[2]);
                sb.append("\t");
                sb.append(".");
                sb.append("\t");
                sb.append(gene.getStrand() > 0?"+":"-");
                sb.append("\t");
                sb.append(".");
                sb.append("\t");
                sb.append("gene_id \"" + intron[0] + "\";");
                sb.append(" ");
                sb.append("locus_id \"" + gene.getLocusID() + "\";");
                sb.append(" ");
                sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_READS+" "+String.format("%1$f", (float) m.get(s)[0].intValue()) +";");
                sb.append(" ");
                sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_FRAC_COVERED+" "+m.get(s)[1]);//+";");
                sb.append(barna.commons.system.OSChecker.NEW_LINE);
            }
            Log.print(sb.toString());
        }


        /**
         * Writes the GTF output file
         *
         * @param g      the splcing graph with annotation-mapped reads
         * @param events AStalavista events
         * @param solver the corresponding <code>GraphLPSolver</code> instance
         */
        private void outputGFF(AnnotationMapper g, ASEvent[] events, GraphLPsolver solver) {

            //TODO decompose in multiple methods

            // check locus
            stats.incrLoci(1);
            if (solver != null || nrMappingsReadsOrPairs > 0)
                stats.incrLociExp(1);
            double valOF = solver == null ? 0 : solver.getValObjFunc();
            if (valOF > FluxCapacitorConstants.BIG) {
                stats.incrLociUnsolved(1);
                Log.warn("Unsolved system: " + gene.getLocusID());
            }

            // pre-build rpkm hash
            HashMap<String, Double> rpkmMap = null;
            long nrReads = stats.getReadsTotal();
            long base = (nrReads <= 0 ? 1 : nrReads);
            Transcript[] tt = gene.getTranscripts();
            if (output.contains(OutputFlag.BALANCED)) {
                rpkmMap = new HashMap<String, Double>(tt.length, 1f);
                for (int i = 0; i < tt.length; i++) {
                    Transcript tx = tt[i];
                    String tid = tt[i].getTranscriptID();

                    double val = 0d;
                    if (solver == null)
                        val = nrMappingsReadsOrPairs * 2;
                    else {
                        val = solver.getTrptExprHash().get(tid);
                        if (val < 1 - costBounds[0]) // 1- 0.95
                            val = 0;
                    }

                    if (val > 0 && !(output.contains(OutputFlag.OBSERVATION) || output.contains(OutputFlag.PREDICTION)))
                        stats.incrTxsExp(1);

                    double rpkm = FluxCapacitor.calcRPKM(val, tx.getExonicLength(), base);
                    if (Double.isNaN(rpkm))
                        Log.warn("NaN RPKM produced: " + val + " / " + base + " = " + rpkm);

                    rpkmMap.put(tid, rpkm);
                }
            }


            // reproduce original
            boolean foundExons = true, foundTranscripts = false;
            /*if (((GTFwrapper) getWrapperGTF()).isKeepOriginalLines() && origLines != null) {
                foundTranscripts = outputGFForiginalLines(rpkmMap);
            }*/

            StringBuilder sb = new StringBuilder();
            // LOCUS TODO genes
            if (output.contains(OutputFlag.GENE)) {
                if (output.contains(OutputFlag.OBSERVATION) || output.contains(OutputFlag.PREDICTION)) {
                    //getGTF(sb, g.trpts[0].getGene(), g, solver, perM, pv);
                    try {
                        assert (testInvariant(invariantTestObsSplitFreq,
                                pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs, 0.05));
                    }    // min: 0.01
                    catch (AssertionError e) {
                        Log.warn(getClass().getName() + ".outputGFF():\n\tinvariantTestObsSplitFreq= "
                                + invariantTestObsSplitFreq + ", nrMappingsReadsOrPairs= " + (pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs)
                                + "\n\tlocus: " + g.trpts[0].getTranscriptID());
                    }
                    ;
                    try {
                        assert (testInvariant(invariantTestPredSplitFreq,
                                pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs, 0.1));
                    } catch (AssertionError e) {
                        Log.warn(getClass().getName() + ".outputGFF():\n\tinvariantTestPredSplitFreq= "
                                + invariantTestPredSplitFreq + ", nrMappingsReadsOrPairs= " + (pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs)
                                + "\n\tlocus: " + g.trpts[0].getTranscriptID());
                    }
                    ;

                } else if (output.contains(OutputFlag.BALANCED)) {
                }
            }


            // TRANSCRIPTS
            if (output.contains(OutputFlag.TRANSCRIPT) || output.contains(OutputFlag.EXON) || output.contains(OutputFlag.SJ)) {
                float invariantObsAllTx = 0, invariantPredAllTx = 0,
                        invariantObsAllEx = 0, invariantPredAllEx = 0;
                for (int i = 0; i < tt.length; i++) {
                    stats.incrTxs(1);
//					float invariantObsTx= invariantTestObsSplitFreq,
//					invariantPredTx= invariantTestPredSplitFreq;
                    String tid = tt[i].getTranscriptID();
                    float invariantObsTx = 0, invariantPredTx = 0;
                    if (output.contains(OutputFlag.TRANSCRIPT) && !foundTranscripts) {
                        if (output.contains(OutputFlag.OBSERVATION) || output.contains(OutputFlag.PREDICTION)) {
                            //getGTF(sb, g.trpts[i], solver, g, perM, null, false);	// writer.write
                            invariantObsAllTx += invariantTestObsSplitFreq; //invariantObsTx;
                            invariantPredAllTx += invariantTestPredSplitFreq; // invariantPredTx;
                            invariantObsTx = invariantTestObsSplitFreq;
                            invariantPredTx = invariantTestPredSplitFreq;
                            if (invariantPredTx > 0)
                                stats.incrTxsExp(1);

                        } else if (output.contains(OutputFlag.BALANCED)) {

                            GFFObject obj = GFFObject.createGFFObject(tt[i]);
                            sb.append(obj.toString());
                            int x = sb.length();
                            while (Character.isWhitespace(sb.charAt(--x)))
                                sb.delete(x, x + 1);
                            if (sb.charAt(x) != ';')
                                sb.append("; ");
                            else
                                sb.append(" ");

                            // deconvoluted reads
                            sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_READS);
                            sb.append(" ");
                            sb.append(String.format("%1$f",
                                    (float) (rpkmMap.get(tid) * tt[i].getExonicLength() * ((double)base / 1000000000l))));
                            sb.append("; ");

                            // spliced length
                            sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_LENGTH);
                            sb.append(" ");
                            sb.append(Integer.toString(tt[i].getExonicLength()));
                            sb.append("; ");

                            // rpkm
                            sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
                            sb.append(" ");
                            //sb.append(rpkmMap.get(g.trpts[i].getTranscriptID()));
                            // avoid scientific notation
                            sb.append(String.format("%1$f", rpkmMap.get(tid).floatValue()));
                            sb.append(barna.commons.system.OSChecker.NEW_LINE);
                        }
                    }
                    // EXONS
                    float invariantObsEx = 0, invariantPredEx = 0;
                    if (output.contains(OutputFlag.EXON) && !foundExons) {
                        Exon[] exons = tt[i].getExons();
                        for (int j = 0; j < exons.length; j++) {
                            //getGTF(sb, exons[j], tt[i], g, solver, unsolvedSystem, perM, null, false);
                            invariantObsEx += invariantTestObsSplitFreq;
                            invariantPredEx += invariantTestPredSplitFreq;
                        }
                    }

                    // SJ
                    if (output.contains(OutputFlag.SJ)) {
                        Vector<Vector<AbstractEdge>> eeV = new Vector<Vector<AbstractEdge>>(5, 5);
                        eeV.add(new Vector<AbstractEdge>());
                        g.getRPK(tt[i], pairedEnd, SplicingGraph.ETYPE_SJ, eeV);
                        //long[][] sig= new long[][]{g.encodeTset(tt[i])};
                        for (int j = 0; j < eeV.elementAt(0).size(); j++) {
                            //getGTF(sb, eeV.elementAt(0).elementAt(j), sig, g, solver, perM);
                            invariantObsEx += invariantTestObsSplitFreq;
                            invariantPredEx += invariantTestPredSplitFreq;
                        }
                    }
                    invariantObsAllEx += invariantObsEx;
                    invariantPredAllEx += invariantPredEx;

                    if (output.contains(OutputFlag.EXON) && output.contains(OutputFlag.SJ) && output.contains(OutputFlag.TRANSCRIPT)) {
                        try {
                            assert (testInvariant(invariantObsEx, invariantObsTx, 0.05));
                        }    // min: 0.02
                        catch (AssertionError e) {
                            Log.warn(getClass().getName() + ".outputGFF():\n\tinvariantObsEx= "
                                    + invariantObsEx + ", invariantObsTx= " + invariantObsTx
                                    + "\n\tlocus: " + tt[0].getTranscriptID());
                        }
                        ;
                        try {
                            assert (testInvariant(invariantPredEx, invariantPredTx, 0.1));
                        } catch (AssertionError e) {
                            Log.warn(getClass().getName() + ".outputGFF():\n\tinvariantPredEx= "
                                    + invariantPredEx + ", invariantPredTx= " + invariantPredTx
                                    + "\n\tlocus: " + tt[0].getTranscriptID());
                        }
                        ;
                    }
                }
                if (output.contains(OutputFlag.TRANSCRIPT)) {
                    try {
                        assert (testInvariant(invariantObsAllTx,
                                pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs, 0.05));
                    }    // min: 0.01
                    catch (AssertionError e) {
                        Log.warn(getClass().getName() + ".outputGFF():\n\tinvariantObsAllTx= "
                                + invariantObsAllTx + ", nrMappingsReadsOrPairs= " + (pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs)
                                + "\n\tlocus: " + tt[0].getTranscriptID());
                    }
                    try {
                        assert (testInvariant(invariantPredAllTx,
                                pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs, 0.1));
                    } catch (AssertionError e) {
                        Log.warn(getClass().getName() + ".outputGFF():\n\tinvariantPredAllTx= "
                                + invariantPredAllTx + ", nrMappingsReadsOrPairs= " + (pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs)
                                + "\n\tlocus: " + tt[0].getTranscriptID());
                    }
                }
                if (output.contains(OutputFlag.EXON) && output.contains(OutputFlag.SJ)) {
                    try {
                        assert (testInvariant(invariantObsAllEx,
                                pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs, 0.05));
                    }    // min: 0.02
                    catch (AssertionError e) {
                        Log.warn(getClass().getName() + ".outputGFF():\n\tinvariantObsAllEx= "
                                + invariantObsAllEx + ", nrMappingsReadsOrPairs= " + (pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs)
                                + "\n\tlocus: " + tt[0].getTranscriptID());
                    }
                    try {
                        assert (testInvariant(invariantPredAllEx,
                                pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs, 0.1));
                    } catch (AssertionError e) {
                        Log.warn(getClass().getName() + ".outputGFF():\n\tinvariantPredAllEx= "
                                + invariantPredAllEx + ", nrMappingsReadsOrPairs= " + (pairedEnd ? nrMappingsReadsOrPairs * 2 : nrMappingsReadsOrPairs)
                                + "\n\tlocus: " + tt[0].getTranscriptID());
                    }
                }
            }

            // EVENTS
            if (output.contains(OutputFlag.EVENT)) {
                HashMap<Object, Double> tExpMap = null;
                if (solver != null) {
                    tExpMap = solver.getTrptExprHash();
                    Object[] keys = tExpMap.keySet().toArray();
                    for (int i = 0; i < keys.length; i++) {
                        if (!(keys[i] instanceof String))
                            continue;
                        if (tExpMap.get(keys[i]) < 0)
                            tExpMap.put((String) keys[i], 0d);    // TODO ugly
                    }
                }
                for (int i = 0; events != null && i < events.length; i++) {
                    if (output.contains(OutputFlag.OBSERVATION) || output.contains(OutputFlag.PREDICTION))
                        ; //getGTF(sb, events[i], g, solver, unsolvedSystem, perM, pv, tExpMap);
                    else
                        stats.incrEvents(1);
                    if (output.contains(OutputFlag.BALANCED)) {
                        sb.append(events[i].toStringGTF());
                        sb.append(" ");
                        sb.append("\"");
                        boolean allPos = true;
                        for (int j = 0; j < events[i].getTranscripts().length; j++) {
                            float sum = 0;
                            for (int k = 0; k < events[i].getTranscripts()[j].length; k++)
                                sum += rpkmMap.get(events[i].getTranscripts()[j][k].getTranscriptID());
                            sb.append(sum);
                            sb.append(",");
                            allPos &= (sum > 0);
                        }
                        if (allPos && !(output.contains(OutputFlag.OBSERVATION) || output.contains(OutputFlag.PREDICTION)))
                            stats.incrEventsExp(1);
                        sb.replace(sb.length() - 1, sb.length(), "\";"+OSChecker.NEW_LINE);
                    }

                }
            }

            // FRAGMENTS and XJUNCTIONS
            if (false && solver != null) {
                ArrayList<AbstractEdge> cc = new ArrayList<AbstractEdge>();
                if (solver != null) {
                    Iterator<Object> iter = solver.getConstraintHash().keySet().iterator();
                    while (iter.hasNext()) {
                        Object o = iter.next();
                        if (o instanceof SimpleEdge)
                            cc.add((SimpleEdge) o);
                    }
                }
                Collections.sort(cc, SimpleEdge.getDefaultPositionComparator());

                Iterator<AbstractEdge> iter = cc.iterator();
                while (iter.hasNext()) {
                    AbstractEdge e = iter.next();
                    // no INTRONS
                    if ((!(e instanceof SuperEdge)) && (!e.isExonic()))
                        continue;
                    //getGTF(sb, e, new long[][]{e.getTranscripts()}, g, solver, perM);
                }
            }

            Log.print(sb.toString());

        }

        /**
         * Re-produces the original lines read in the GTF input file
         *
         * @param rpkmMap hash to map transcript ID to an deconvoluted expression value
         * @return <code>true</code> if <code>transcript</code> features were found in
         *         the input, <code>false</code> otherwise
         */
        /*private boolean outputGFForiginalLines(HashMap<String, Double> rpkmMap) {

            Transcript[] tt = gene.getTranscripts();
            boolean foundTranscripts = false;
            for (int i = 0; i < origLines.size(); i++) {
                String s = origLines.elementAt(i);
                String feat = GFFObject.getField(3, s);
                String tid = GFFObject.getTranscriptID(s);
                int tx = 0;
                if ((feat.equals(feat.equals(Transcript.GFF_FEATURE_TRANSCRIPT))
                        || feat.equals(Exon.GFF_FEATURE_EXON))
                        && (output.contains(OutputFlag.OBSERVATION) || output.contains(OutputFlag.PREDICTION)))
                    for (tx = 0; tx < tt.length; tx++)
                        if (tt[tx].getTranscriptID().equals(tid))
                            break;
                if (tx >= tt.length) {
                    System.err.println("\nTranscript " + tid + " not found in: ");
                    for (int j = 0; j < tt.length; j++)
                        System.err.println("\t" + tt[j].getTranscriptID());
                    System.err.println();
                }

                if (feat.equals(Transcript.GFF_FEATURE_TRANSCRIPT) && output.contains(OutputFlag.TRANSCRIPT)) {
                    foundTranscripts = true;
                    StringBuilder sb = new StringBuilder(s);
                    int x = sb.length();    // trim
                    while (Character.isWhitespace(sb.charAt(--x)))
                        sb.delete(x, x + 1);
                    if (sb.charAt(x) != ';')
                        sb.append("; ");
                    else
                        sb.append(Constants.SPACE);

                    if ((output.contains(OutputFlag.OBSERVATION) || output.contains(OutputFlag.PREDICTION)) && tx < tt.length)
                        ; //getGTF(sb, tt[tx], solver, g, perM, pv, true);
                    else if (output.contains(OutputFlag.BALANCED)) {
                        sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
                        sb.append(Constants.SPACE);
                        if (rpkmMap.containsKey(tid))
                            sb.append(String.format("%1$f", rpkmMap.get(tid).floatValue()));    // rgasp parser does not like scientific notation
                        else
                            sb.append(Constants.NULL);
                        sb.append(";"+ OSChecker.NEW_LINE);
                    }

                    Log.print(sb.toString());

                } else if (feat.equals(Exon.GFF_FEATURE_EXON) && output.contains(OutputFlag.EXON)) {

                    StringBuilder sb = new StringBuilder(s);
                    int x = sb.length();
                    while (Character.isWhitespace(sb.charAt(--x)))
                        sb.delete(x, x + 1);
                    if (sb.charAt(x) != ';')
                        sb.append("; ");
                    else
                        sb.append(' ');


                    if ((output.contains(OutputFlag.OBSERVATION) || output.contains(OutputFlag.PREDICTION)) && tx < tt.length) {
                        int start = Integer.parseInt(GFFObject.getField(4, s));
                        int end = Integer.parseInt(GFFObject.getField(5, s));
                        int j = 0;
                        for (; j < tt[x].getExons().length; j++) {
                            int begin = Math.abs(tt[x].getExons()[j].getStart()),
                                    ende = Math.abs(tt[x].getExons()[j].getEnd());
                            if (begin > ende) {
                                int h = begin;
                                begin = ende;
                                ende = h;
                            }
                            if (begin == start && ende == end)
                                break;
                        }
                        //getGTF(sb, tt[x].getExons()[j], tt[i], g, solver, unsolvedSystem, perM, pv, true);
                    } else if (output.contains(OutputFlag.BALANCED)) {
                        sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
                        sb.append(Constants.SPACE);
                        if (rpkmMap.containsKey(tid))
                            sb.append(String.format("%1$f", rpkmMap.get(tid).floatValue()));    // rgasp parser does not like scientific notation
                        else
                            sb.append(Constants.NULL);
                        sb.append(";"+OSChecker.NEW_LINE);
                    }


                    Log.print(sb.toString());
                } else if (output.contains(OutputFlag.UNKNOWN)) {
                    Log.print(s + System.getProperty("line.separator"));
                }
            }

            return foundTranscripts;
        }*/


        /**
         * Compares an invariant to the reference value and shouts if the difference
         * is larger than the specified tolerance.
         *
         * @param invariant  the tested value
         * @param reference  the reference value
         * @param stringency the relative tolerance level
         * @return <code>true</code> if the invariant passes the test, <code>false</code>
         *         otherwise
         */
        private boolean testInvariant(double invariant, double reference, double stringency) {
            double delta = Math.abs(reference == 0 ? invariant : (invariant - reference) / reference);
            if (delta > stringency) {
                if (invariant <= 0)
                    return true;     // catch 0-predictions
                return false;
            }

            return true;
        }

        /**
         * Sets the dependency on a thread that is to finish before execution can start.
         *
         * @param threadBefore the thread that is waited for
         */
        void setThreadBefore(Thread threadBefore) {
            this.threadBefore = threadBefore;
        }

        /**
         * Creates a <code>GraphLPsolver</code> instance.
         *
         * @param mapper      the annotation mapper
         * @param mappedReads the number of mapped reads
         * @return a <code>GraphLPsolver</code> instance
         */
        private GraphLPsolver getSolver(AnnotationMapper mapper, int mappedReads) {

            GraphLPsolver solver = new GraphLPsolver(
                    settings,
                    mapper,
                    profile.getMappingStats(),
                    //pairedEnd ? insertMinMax : null,
                    mappedReads,
                    stranded,
                    pairedEnd);
            /*if (output.contains(OutputFlag.LP))
                solver.setFileLPdir(getFileLP());*/
            solver.costModel = costModel;    // COSTS_LINEAR
            solver.setCostSplit(costSplit);
            solver.setProfile(profile);
            solver.setMappingStats(profile.getMappingStats());
            solver.costBounds = costBounds;

            return solver;
        }


        /*
         * Learns systematic biases along a transcript
         *
         * @param tx   the Transcript
         * @param mappings the mappings
         *
         */
        /*@Deprecated
        private void learn(Transcript tx, MSIterator<Mapping> mappings) {

            if (mappings== null)
                return;
            if (!mappings.hasNext())
                return;

            Mapping bed1, bed2;
            UniversalReadDescriptor.Attributes
                    attributes = settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).createAttributes(),
                    attributes2 = settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).createAttributes();
            int elen = tx.getExonicLength();    // this is the "effective" length, modify by extensions
//				if (elen< readLenMin)
//					return;	// discards reads

            UniversalMatrix m = profile.getMatrix(elen);
            if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
                if (coverage == null)
                    coverage = new Coverage(elen);
                else
                    coverage.reset(elen);
            }

            while (mappings.hasNext()) {
                stats.incrReadsSingleTxLoci(1);
                bed1= mappings.next();
                CharSequence tag = bed1.getName();
                attributes = settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).getAttributes(tag, attributes);
                if (pairedEnd) {
                    if (attributes.flag < 1)
                        Log.warn("Read ignored, error in readID: " + tag);
                    if (attributes.flag == 2)    // don't iterate second read
                        continue;
                }

                if (stranded) {
                    if ((tx.getStrand() == bed1.getStrand() && attributes.strand == 2)
                            || (tx.getStrand() != bed1.getStrand() && attributes.strand == 1)) {
                        stats.incrMappingsWrongStrand(1);
                        continue;
                    }
                }

                int bpoint1 = getBpoint(tx, bed1);
                if (bpoint1 < 0 || bpoint1 >= elen) {    // outside tx area, or intron (Int.MIN_VALUE)
                    stats.incrMappingsSingleTxLociNoAnn(1);
                    continue;
                }

                stats.incrMappingsSingleTxLoci(1);  // the (first) read maps

                if (pairedEnd) {

//                    mappings.mark();
                    Iterator<Mapping> mates = mappings.getMates(bed1,settings.get(FluxCapacitorSettings.READ_DESCRIPTOR));
                    while(mates.hasNext()) {
                        bed2= mates.next();
//                        attributes2 = settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).getAttributes(bed2.getName(), attributes2);
//                        if (attributes2 == null)
//                            continue;
//                        if (!attributes.id.equals(attributes2.id))
//                            break;
//                        if (attributes2.flag == 1)    // not before break, inefficient
//                            continue;

                        int bpoint2 = getBpoint(tx, bed2);
                        if (bpoint2 < 0 || bpoint2 >= elen) {
                            stats.incrMappingsSingleTxLociNoAnn(1);
                            continue;
                        }

                        // check again strand in case one strand-info had been lost
                        if (stranded) {
                            if ((tx.getStrand() == bed2.getStrand() && attributes2.strand == 2)
                                    || (tx.getStrand() != bed2.getStrand() && attributes2.strand == 1)) {
                                stats.incrMappingsWrongStrand(1);
                                continue;
                            }
                        }

                        // check directionality (sequencing-by-synthesis)
                        if ((bed1.getStrand() == bed2.getStrand())
                                || ((bed1.getStart() < bed2.getStart()) && (bed1.getStrand() != DirectedRegion.STRAND_POS))
                                || ((bed2.getStart() < bed1.getStart()) && (bed2.getStrand() != DirectedRegion.STRAND_POS))) {
                            stats.incrPairsWrongOrientation(2);
                            continue;
                        }

                        m.add(bpoint1, bpoint2, -1, -1, elen);    // 5TODO rlen currently not used
                        // update coverage
                        if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
                            if (bpoint1 < bpoint2) {
                                for (int i = bpoint1; i < bpoint1 + bed1.getLength(); i++)
                                    coverage.increment(i);
                                for (int i = bpoint2 - bed2.getLength() + 1; i <= bpoint2; i++)
                                    coverage.increment(i);
                            } else {
                                for (int i = bpoint2; i < bpoint2 + bed2.getLength(); i++)
                                    coverage.increment(i);
                                for (int i = bpoint1 - bed1.getLength() + 1; i <= bpoint1; i++)
                                    coverage.increment(i);
                            }
                        }
                        //addInsertSize(Math.abs(bpoint2- bpoint1)+ 1);	// TODO write out insert size distribution

                        stats.incrMappingPairsSingleTxLoci(2);

                    }
//                    mappings.reset();

                } else {    // single reads
                    m.add(bpoint1, -1, elen,
                            bed1.getStrand() == tx.getStrand() ? Constants.DIR_FORWARD : Constants.DIR_BACKWARD);
                    // update coverage
                    if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
                        if (bed1.getStrand() == tx.getStrand()) {
                            for (int i = bpoint1; i < bpoint1 + bed1.getLength(); i++)
                                coverage.increment(i);
                        } else {
                            for (int i = bpoint1 - bed1.getLength() + 1; i <= bpoint1; i++)
                                coverage.increment(i);
                        }
                    }
                }

            } // iterate bed objects


            // output coverage stats
            if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
                writeCoverageStats(
                        tx.getGene().getLocusID(),
                        tx.getTranscriptID(),
                        tx.isCoding(),
                        tx.getExonicLength(),
                        pairedEnd ? stats.getMappingPairsSingleTxLoci() : stats.getMappingsSingleTxLoci(),
                        coverage.getFractionCovered(),
                        coverage.getChiSquare(true),
                        coverage.getCV(true));
            }
        }*/


        /*
         * Returns the breakpoint indicated by a mapping within a transcript.
         *
         * @param tx  transcript to which a read maps
         * @param bed genomic mappping
         * @return transcript coordinate of the breakpoint indicated by the mapping
         */
        /*private int getBpoint(Transcript tx, Mapping bed) {

            // TODO add check whether complete read is contained in transcript

            // just depends on genomic position, not on sense/antisense!
            int gpos = bed.getStrand() >= 0 ? bed.getStart() + 1 : bed.getEnd();
            int epos = tx.getExonicPosition(gpos);

            // security check, get distance between both exonic coordinates
            int epos2 = tx.getExonicPosition(bed.getStrand() >= 0 ? bed.getEnd() : bed.getStart() + 1);
            int len = bed.getLength();
            if (readLenMin < 0 || len < readLenMin)
                readLenMin = len;
            if (len > readLenMax)
                readLenMax = len;

            if (len != Math.abs(epos - epos2) + 1)
                return Integer.MIN_VALUE;
            return epos;
        }*/ /* End of LocusSolver*/

        @Override
        public MappingStats call() throws Exception {
            AnnotationMapper mapper = null;

            if(tasks.isEmpty())
                return null;

            AnnotationMapping am = settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING);

            mapper = new AnnotationMapper(this.gene, am.isPaired(), am.isStranded(), !settings.get(FluxCapacitorSettings.DISABLE_MULTIMAP_WEIGHTING), settings.get(FluxCapacitorSettings.READ_STRAND));
            mapper.setStats(stats);
            mapper.map(this.mappings, settings.get(FluxCapacitorSettings.INSERT_FILE));

            /*stats.incrReadsLoci(mapper.nrMappingsLocus);
            stats.incrMappingsMapped(mapper.getNrMappingsMapped());
            nrMappingsReadsOrPairs += mapper.getNrMappingsMapped() / 2;
            stats.incrMappingPairsNoTx(mapper.getNrMappingsNotMappedAsPair());
            stats.incrPairsWrongOrientation(mapper.getNrMappingsWrongPairOrientation());*/
            stats.setReadsLoci(mapper.nrMappingsLocus);
            stats.setMappingsMapped(mapper.getNrMappingsMapped());
            nrMappingsReadsOrPairs += mapper.getNrMappingsMapped() / 2;
            stats.setMappingPairsNoTx(mapper.getNrMappingsNotMappedAsPair());
            stats.setPairsWrongOrientation(mapper.getNrMappingsWrongPairOrientation());
            stats.setCtrHits(mapper.ctrHits);
            stats.setCtrHitsMultiGenome(mapper.ctrHitsMultiGenome);
            stats.setCtrHitsMultiLocus(mapper.ctrHitsMultiLocus);
            stats.setCtrHitsMultiLocusAndGenome(mapper.ctrHitsMultiLocusAndGenome);
            stats.setCtrHitsNone(mapper.ctrHitsNone);
            stats.setCtrHitsNoneMultiGenome(mapper.ctrHitsNoneMultiGenome);
            // complete profile, if necessary
            if (profile.getMappingStats().getReadLenMin()< 2)
                profile.getMappingStats().setReadLenMin(stats.getReadLenMin());
            if (profile.getMappingStats().getReadLenMax()< 2)
                profile.getMappingStats().setReadLenMax(stats.getReadLenMax());

            //Execute tasks
            for (Task t : this.tasks) {
                switch (t) {
                    case COUNT_INTRONS:
                        outputIntronsGFF(this.gene, mapper);
                        break;
                    case COUNT_SJ:
                        outputSJGFF(this.gene, mapper);
                        break;
                    case DECOMPOSE:
                        GraphLPsolver mySolver = null;
                        if (mapper.nrMappingsMapped > 0 && this.gene.getTranscriptCount() > 1) {    // OPTIMIZE
                            mySolver = getSolver(mapper, (int) (mapper.nrMappingsMapped * 2)); // not: getMappedReadcount()
                            mySolver.run();
                        }
                        outputGFF(mapper, events, mySolver);
                        break;
                }
            }

//            if (decompose) {
//
//                // BUG 110301: do not use mapTrivial, problems with split-maps
////				if (this.gene.getTranscriptCount()== 1) {
////					mapTrivial(gene.getTranscripts()[0], beds);
////					outputGFF(null, null, null);
////				} else {
//
//                GraphLPsolver mySolver = null;
//                if (mapper.nrMappingsMapped > 0 && this.gene.getTranscriptCount() > 1) {    // OPTIMIZE
//                    mySolver = getSolver(mapper, (int) (mapper.nrMappingsMapped * 2)); // not: getMappedReadcount()
//                    mySolver.run();
//                }
//                outputGFF(mapper, events, mySolver);
////					}
//
//            } else {
//                // map all reads
//                if (this.gene.getTranscriptCount() == 1) {
//                    ++nrSingleTranscriptLearn;
//                    learn(this.gene.getTranscripts()[0], beds);
//                }
//            }
            if (mappings != null) {
                mappings.clear();
                mappings = null;
            }
            gene = null;
            // makes it terribly slow
            //System.gc();

//			synchronized(FluxCapacitor.this.threadPool) {
            //FluxCapacitor.this.threadPool.remove(this);
//			}
            return stats;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    /**
     * Workaround for wrapper script, provides first argument as parameter file.
     *
     * @param args the program's arguments
     */
    public static void main(String[] args) {

        Execute.initialize(2);

        try {

            int ok = loadLibraries();
            if (ok < 0)
                exit(-1);


            final FluxCapacitor myCapacitor = new FluxCapacitor();
            myCapacitor.setFile(new File(args[0]));

            // run
            try {
                myCapacitor.call();
            } catch (Throwable XXX) {
                XXX.printStackTrace();
                System.exit(0);
            }

        } catch (Throwable t) {
            if (t instanceof Exception)
                ((Exception) t).printStackTrace();
            else if (t instanceof Error)
                ((Error) t).printStackTrace();
            else
                System.err.println(t.getMessage());
            if (cheatDoNotExit) {
                int in = 0;
                while (in != '\n')
                    try {
                        in = System.in.read();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
            }

        } finally {
            Execute.shutdown();
        }
    }


    /**
     * File name extensions of supported file formats.
     *
     * @deprecated marked for removal
     */
    public static enum SupportedFormatExtensions {
        GTF, GFF, BED, BAM
    }

    /**
     * The parameter file.
     */
    protected File file = null;

    /**
     * Mode of strand consideration for LP constraints:
     * consider XXX (STRAND_NONE), both strands (STRAND_ENABLED), or XXX (STRAND_SPECIFIC).
     */
    byte strand = FluxCapacitorConstants.STRAND_ENABLED;

    /**
     * Delta-region around gene (up-/downstream) within which
     * mapped reads are still considered to solve the locus.
     */
    int tolerance = 1000;

    /**
     * File to which the learned bias profiles are written.
     */
    public File fileProfile = null;

    /**
     * Directory to which the linear programs are written.
     */
    public File fileLPdir = null;

    /**
     * The minimum read length.
     */
    int readLenMin = 75;

    /**
     * The maximum read length.
     */
    int readLenMax = -1;

    /**
     * @deprecated marked for removal
     */
    public static boolean cheatDoNotExit = false;
    /**
     * @deprecated marked for removal
     */
    public static boolean cheatDisableFCheck = false;

    /**
     * The number of non-redundant read IDs in the mappings.
     */
//    int nrBEDreads = -1;

    /**
     * The total number of mappings.
     */
//    int nrBEDmappings = -1;

    /**
     * The number of mappings read in a second run, invariant check.
     */
    private int checkBEDscanMappings = 0;

    /**
     * Vector of threads that are carried out in parallel.
     */
    private Vector<Thread> threadPool = new Vector<Thread>();

    /**
     * Maximum number of parallel threads.
     */
    int maxThreads = 1;

    /**
     * Vector of Strings representing the original lines of the annotation read annotation file.
     *
     * @deprecated marked for removal
     */
    Vector<String> origLines = null;  // TODO check that this vector is no longer initiallized somewhere
    int checkGTFscanExons = 0;

    /**
     * Flag whether to ouput the observation (i.e., values before deconvolution).
     */
//    boolean outputObs = false;

    /**
     * Flag whether to ouput the prediction (i.e., values after deconvolution).
     */
//    boolean outputPred = false;

    /**
     * Flag whether to output balanced values, i.e., values that should be observed as cast back from deconvoluted values.
     */
//    boolean outputBalanced = true;

    /**
     * Flag whether to output exon features.
     */
//    boolean outputExon = false;

    /**
     * Flag whether to output unknown features, e.g., reproduced ones from the input annotation.
     */
//    boolean outputUnknown = false;

    /**
     * Flag whether to output splice-junction features.
     */
//    boolean outputSJunction = false;

    /**
     * Flag whether to output gene/loci features.
     */
//    boolean outputGene = false;

    /**
     * Flag whether to output transcript features.
     */
//    boolean outputTranscript = true;

    /**
     * Flag whether to output AStalavista event features.
     */
//    boolean outputEvent = false;

    /**
     * Flag whether to output linear programs.
     */
//    boolean outputLP = false;

    /**
     * Default compression for output files.
     */
    byte compressionOut = FileHelper.COMPRESSION_NONE;

    /**
     * Dimension of AStalavista events retrieved from the gene structures.
     */
    int eventDim = 2;

    /**
     * Flag indicating whether read biases should be ignored (<code>true</code>),
     * i.e. uniform read distribution is assumed.
     */
    boolean uniform = false;

    /*
     * Volatile:
     * The value of this variable will never be cached thread-locally:
     * all reads and writes will go straight to "main memory";
     * Access to the variable acts as though it is enclosed in a
     * synchronized block, synchronized on itself.
     */

    /**
     * Number of loci in the annotation.
     */
//    volatile int nrLoci = 0;

    /**
     * Number of loci expressed, i.e., with positive quantification values.
     */
//    volatile int nrLociExp = 0;

    /**
     * Number of transcripts in the annotation.
     */
//    volatile int nrTx = 0;

    /**
     * Number of transcripts expressed, i.e., with positive quantification values.
     */
//    volatile int nrTxExp = 0;

    /**
     * Number of AStalavista events in the annotation.
     */
//    volatile int nrEvents = 0;

    /**
     * Number of AStalavista events expressed, i.e., with positive quantification values.
     */
//    volatile int nrEventsExp = 0;

    /**
     * Number of mappings that could be mapped to the annotation.
     */
//    volatile int nrReadsMapped = 0;

    /**
     * Number of mappings that were found in annotated genes/loci.
     */
//    volatile int nrReadsLoci = 0;

    /**
     * Number of not alternatively processed loci in the annotation.
     */
//    volatile int nrSingleTranscriptLoci = 0;    // TODO counted in learn AND decompose redundantly

    /**
     * Number of not alternatively processed loci that have been used for learning biases.
     */
    volatile int nrSingleTranscriptLearn = 0;  // TODO what is the difference between this and nrSingleTranscriptLoci

    /**
     * Number of mappings within the boundaries of not alternatively processed loci.
     */
//    volatile int nrReadsSingleLoci = 0;

    /**
     * Number of mappings that mapped to not alternatively processed loci.
     */
//    volatile int nrReadsSingleLociMapped = 0;

    /**
     * Number of mapping pairs that mapped to not alternatively processed loci.
     */
//    volatile int nrReadsSingleLociPairsMapped = 0;

    /**
     * Number of mappings in not alternatively processed loci that did not match
     * the expected gene structure provided in the annotation.
     */
//    volatile int nrReadsSingleLociNoAnnotation = 0;

    /**
     * Number of loci that could not be solved.
     */
//    volatile int nrUnsolved = 0;

    /**
     * Number of mappings that have an unexpected length.
     *
     * @deprecated should no longer be used
     */
//    volatile int nrReadsWrongLength = 0; // TODO (how) can this still occur?

    /**
     * Number of mappings that do not match the expected strand.
     */
//    volatile int nrMappingsWrongStrand = 0;

    /**
     * Number of mapping pairs that do not match the annotation.
     */
//    volatile int nrPairsNoTxEvidence = 0;

    /**
     * Number of mapping pairs with the wrong orientation of both mates.
     */
//    volatile int nrPairsWrongOrientation = 0;

    /**
     * The minimum and the maximum insert size found.
     */
    int[] insertMinMax = null;   // TODO check if correctly used

    /**
     * Reader to read mappings from a file.
     */
    private MappingReader mappingReader;

    /**
     * Wrapper to read the annotation from a GTF file format.
     *
     * @deprecated marked for removal
     */
    private GTFwrapper gtfReader;   // TODO pull up to AnnotationWrapper

    /**
     * A profile instance representing bias matrices.
     */
    Profile profile;

    /**
     * Object describing the settings of the run.
     */
    FluxCapacitorSettings settings = null;

    /**
     * Flag indicating whether paramters and their description are to be output.
     */
    protected boolean printParameters;


    /**
     * Dummy constructor.
     */
    public FluxCapacitor() {
    }

    /**
     * Intialize the capacitor with settings
     */
    public FluxCapacitor(FluxCapacitorSettings settings) {
        this.settings = settings;
    }

    /**
     * Check operating system and load the native libraries. Exceptions are catched and logged here.
     * Use the return value to check whether loading was successfull.
     *
     * @return 0 if native libraries could be loaded successfully, (-1) otherwise
     */
    public static int loadLibraries() {
        Log.info("PRE-CHECK", "I am checking availability of the required lpsolve JNI libs.");
        VersionInfo lpVer = null;
        try {
            LPSolverLoader.load();
            lpVer = LpSolve.lpSolveVersion();
            if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP) {
                Log.info("PRE-CHECK", "\t* successfully loaded lpsolve JNI (version " + lpVer.getMajorversion() + "." + lpVer.getMinorversion()
                        + ",release " + lpVer.getRelease() + ",build " + lpVer.getBuild() + ")\n");
            }
            return 0;
        } catch (Exception e) {
            Log.error("Error while loading native libraries: " + e.getMessage());
            Log.error("You can try to set the environment variables " + LPSolverLoader.ENV_JNI + " for \n" +
                    "the JNI library and " + LPSolverLoader.ENV_LIB + " for the shared object library to support \n" +
                    "your operating system. The sources for the LPSolver can be found @ http://sourceforge.net/projects/lpsolve");
        }
        return -1;
    }

    // don't load class if LP libs are not present
    static {
        if (loadLibraries() < 0)
            System.exit(-1);
    }

    /**
     * Says goodbye and exits providing the given return code to the environment.
     *
     * @param code the return code
     */
    static void exit(int code) {
        String pfx = "[ASTALAVISTA] ";
        if (code < 0)
            pfx = "[CIAO] ";
        if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP)
            System.err.println(pfx + "I'm exiting.");
        System.exit(code);
    }

    /**
     * Writes the initialized settings of a run to a logging stream.
     */
    private void printSettings() {
        Log.info("HEHO", "We are set, so let's go!");
        // TODO
        // settings.write(Log.logStream);
        StringBuilder sb;

        // INPUT
        Log.info(FluxCapacitorSettings.ANNOTATION_FILE.getName(),
                settings.get(FluxCapacitorSettings.ANNOTATION_FILE).getAbsolutePath());
        Log.info(FluxCapacitorSettings.MAPPING_FILE.getName(),
                settings.get(FluxCapacitorSettings.MAPPING_FILE).getAbsolutePath());
        Log.info(FluxCapacitorSettings.READ_DESCRIPTOR.getName(),
                settings.getReadDescriptor().toString());
        String ext = FileHelper.getExtension(settings.get(FluxCapacitorSettings.MAPPING_FILE)).toUpperCase();
        if(ext.equals("BED") || ext.equals("GZ")) {
            Log.info(FluxCapacitorSettings.READ_DESCRIPTOR.getName(),
                    settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).toString());
        }
        Log.info("\tminimum intron length "+ Transcript.maxLengthIntronIsGap);
        if (settings.get(FluxCapacitorSettings.PROFILE_FILE) != null) {
            Log.info(FluxCapacitorSettings.PROFILE_FILE.getName(),
                    settings.get(FluxCapacitorSettings.PROFILE_FILE).toString());
        }

        Log.info(settings.SORT_IN_RAM.getName(),
                Boolean.toString(settings.get(FluxCapacitorSettings.SORT_IN_RAM)));

        // OUTPUT
        Log.info(FluxCapacitorSettings.TMP_DIR.getName(),
                settings.get(FluxCapacitorSettings.TMP_DIR).getAbsolutePath());
        sb = new StringBuilder();
        if (settings.get(FluxCapacitorSettings.STDOUT_FILE) == null)
            sb.append("stdout");
        else {
            sb.append(settings.get(FluxCapacitorSettings.STDOUT_FILE).getAbsolutePath());
            if (compressionOut != FileHelper.COMPRESSION_NONE)
                sb.append("\t" + FluxCapacitorConstants.CLI_LONG_COMPRESSION + "\t" + FileHelper.COMPRESSION_KEYWORDS[compressionOut]);
        }
        Log.info(settings.STDOUT_FILE.getName(), sb.toString());

        if (settings.isStranded())
            Log.info("\tstrand information considered.");
        if (settings.isPaired())
            Log.info("\tmate pairing information considered");

        if (settings.get(FluxCapacitorSettings.INSERT_FILE) != null) {
            Log.info("\twriting insert sizes to " + settings.get(FluxCapacitorSettings.INSERT_FILE).getAbsolutePath());
        }

    }

    /**
     * Finishes all pending file I/O operations and closes handles.
     */
    void fileFinish() {

        // TODO close input should occur by reader or interface method
		mappingReader.close();
        gtfReader.close();

        // TODO close files for non-/mapped reads, insert sizes, LPs, profiles

        // profiles
        /*if (settings.get(FluxCapacitorSettings.PROFILE_FILE) != null)
            writeProfiles();*/

        // close output
        if (Log.outputStream != System.out && Log.outputStream != System.err)
            Log.outputStream.close();

    }

    /**
     * Retrieves statistics about the provided annotation.
     *
     * @param wrapper annotation reader
     */
    private void fileStats(AnnotationWrapper wrapper) {

        // (3) scan
        ((AbstractFileIOWrapper) wrapper).scanFile();
        if (((AbstractFileIOWrapper) wrapper).getNrInvalidLines() > 0)
            Log.warn("Skipped " + ((AbstractFileIOWrapper) wrapper).getNrInvalidLines() + " lines.");

        Log.info(Constants.TAB + wrapper.getNrGenes() + " loci, "
                + wrapper.getNrTranscripts() + " transcripts, "
                + wrapper.getNrExons() + " exons.");
    }


    /**
     * Executes a FC run.
     *
     * @return statistics about the run
     * @throws Exception if something went wrong
     */
    @Override
    public MappingStats call() throws Exception {

        // parameters
        if (settings== null) // if not pre-inited by API..
            settings= createSettings(file, commandLineArgs);
        FileHelper.tempDirectory = settings.get(FluxCapacitorSettings.TMP_DIR).getAbsoluteFile();
        currentTasks= getTasks(settings);

        // pre-processing
        gtfReader= createAnnotationReader(settings.get(FluxCapacitorSettings.ANNOTATION_FILE), settings);
        mappingReader= createMappingReader(settings.get(FluxCapacitorSettings.MAPPING_FILE), settings);
        if (!settings.get(FluxCapacitorSettings.DISABLE_FILE_CHECK)) {
            if (stats == null)
                stats = new MappingStats(); //Initialize stats
            fileStats(mappingReader);   // dont deactivate, rpkm will be 0
            Log.info("Annotation and mapping input checked");
        } else {
            stats= new MappingStats();
            stats.setReadsTotal(100007838); // TODO
        }

        //Gene[] oGenes= gtfReader.getGenes();
        Gene[] genes= null; //PreProcessor.collapse(oGenes);
        //Log.info("Collapsed "+ oGenes.length+ " genes into "+ genes.length+ " loci.");
        //oGenes= null;
        System.gc();

        // process tasks
        long t0 = System.currentTimeMillis();
        printSettings();         // settings only complete after pre-processing
        // 1 profiling to file
        if (currentTasks.contains(Task.PROFILE)) {
            BiasProfiler profiler = new BiasProfiler(this, strand, settings.isPaired(),
                    !settings.get(FluxCapacitorSettings.DISABLE_MULTIMAP_WEIGHTING), gtfReader,mappingReader);
            stats=profiler.call().getMappingStats();
            printProfile((System.currentTimeMillis() - t0) / 1000);
        }

        if (!currentTasks.isEmpty()) {

            fileProfile = settings.get(FluxCapacitorSettings.PROFILE_FILE);
            profile = getProfile();
            if (profile == null) {
                throw new RuntimeException("Cannot evaluate profile");
            }
            profile.getMappingStats().add(stats);
            stats=profile.getMappingStats();
            printProfile((System.currentTimeMillis() - t0) / 1000);

            explore(genes);

            if (settings.get(FluxCapacitorSettings.STATS_FILE)!=null)
                stats.writeStats(settings.get(FluxCapacitorSettings.STATS_FILE));
        }

        fileFinish();
        Log.info("\n[TICTAC] I finished flux in "
                + ((System.currentTimeMillis() - t0) / 1000) + " sec.\nCheers!");

        return stats;
    }

    /**
     * Creates or returns the tasks for this run.
     * @param settings parameters for setting up tasks
     * @return task list
     */
    public EnumSet<Task> getTasks(FluxCapacitorSettings settings) {

        if (currentTasks== null || currentTasks.isEmpty()) {
            currentTasks= createTasks(settings);
        }

        return currentTasks;
    }

    /**
     * Composes a list of tasks from the provided settings.
     * @param settings parameters for setting up tasks
     * @return task list
     */
    public static EnumSet<Task> createTasks(FluxCapacitorSettings settings) {

        EnumSet<Task> currentTasks= EnumSet.noneOf(Task.class);

        //Get from settings the tasks to be executed in the current run
        if (!settings.get(FluxCapacitorSettings.COUNT_ELEMENTS).isEmpty()) {
            for (FluxCapacitorSettings.CountElements e : settings.get(FluxCapacitorSettings.COUNT_ELEMENTS)) {
                switch (e) {
                    case SPLICE_JUNCTIONS:
                        currentTasks.add(Task.COUNT_SJ);
                        break;
                    case INTRONS:
                        currentTasks.add(Task.COUNT_INTRONS);
                        break;
                }
            }
        }
        if (!settings.get(FluxCapacitorSettings.DISABLE_DECONVOLUTION)) {
            currentTasks.add(Task.DECOMPOSE);
        }

        return currentTasks;
    }

    public static MappingReader createMappingReader(File mappingFile, FluxCapacitorSettings settings) {
        // reads
        Log.progressStart("Scanning mapping file");
        MappingReader mappingReader =
                (MappingReader) fileInit(settings.get(FluxCapacitorSettings.MAPPING_FILE), settings);
        Log.progressFinish("OK", true);

        return mappingReader;
    }


    public static GTFwrapper createAnnotationReader(File annotationFile, FluxCapacitorSettings settings) {
        // annotation
        Log.progressStart("Scanning annotation file");
        GTFwrapper gtfReader = (GTFwrapper) fileInit(annotationFile, settings);
        gtfReader.setBasic(true);
        //gtfReader.loadAllGenes();
        //fileStats((gtfReader));
        Log.progressFinish("OK", true);
        return gtfReader;
    }

    private void printProfile(long secs) {
        System.err.println("\tfirst round finished .. took " + secs + " sec.\n\n\t"
                + stats.getSingleTxLoci() + " single transcript loci\n\t"
                + mappingReader.getCountMappings()+ " mappings in file\n\t"  //    removed getNrLines()
                + stats.getReadsSingleTxLoci() + " mappings fall in single transcript loci\n\t"    // these loci(+/-"+tolerance+"nt)\n\t"
                // counter un-reliable, /2 read is skipped in paired-end mode
                + ((strand == FluxCapacitorConstants.STRAND_SPECIFIC) ? stats.getMappingsWrongStrand() + " mappings map to annotation in antisense direction,\n\t" : "")
                //+ (pairedEnd?(nrReadsSingleLociPotentialPairs+ " mappings form potential pairs,\n\t"):"")
                + (settings.isPaired() ? (stats.getMappingPairsSingleTxLoci()) + " mappings in annotation-mapped pairs\n\t"
                : stats.getMappingsSingleTxLoci() + " mappings map to annotation\n\t")
                //+ nrReadsSingleLociNoAnnotation+ " mappings do NOT match annotation,\n\t"
                //+ (uniform?"":func.profiles.size()+" profiles collected\n\t")
                + (stats.getReadLenMin() == Integer.MAX_VALUE ? "NA" : stats.getReadLenMin()) + "," + (stats.getReadLenMax() == -1 ? "NA" : stats.getReadLenMax()) + " min/max read length\n\t"
                + (settings.isPaired() && insertMinMax != null ? insertMinMax[0] + "," + insertMinMax[1] + " min/max insert size\n\t" : ""));
        //nrUniqueReads= getBedReader().getNrUniqueLinesRead();
        //System.err.println("\ttotal lines in file "+nrUniqueReads);
//        System.err.println();
    }

    /**
     * Performs the RPKM normalization, given a number of reads mapping to a transcript of a certain length.
     *
     * @param reads the number of mappings to the transcript
     * @param len   the length of the transcript
     * @return the RPKM value
     */
    public static double calcRPKM(double reads, long len, long totalReads) {
        return ((reads / (double) len) * (1000000000l / (double) (totalReads <= 0 ? 1 : totalReads)));
    }

    /**
     * Creates a composite name between both input files.
     *
     * @return a string representation the composite name
     */
    public String getCompositeFName() {
        File f = settings.get(FluxCapacitorSettings.ANNOTATION_FILE),
                g = settings.get(FluxCapacitorSettings.MAPPING_FILE);
        return MyFile.stripExtension(f.getName()) + "__"
                + MyFile.stripExtension(g.getName());
    }

    /**
     * Creates a temporary file in the location provided, iff write access is
     * available there. Otherwise the file is created in the custom or system
     * temporary directory.
     *
     * @param location     a file in the target directory or the directory itself,
     *                     may be <code>null</code>
     * @param name         prefix of the file to be created, class name is appended
     *                     at the beginning
     * @param extension    (optional) suffix of the temporary file that is created
     * @param deleteOnExit flag for calling the <code>deleteOnExit()</code>
     *                     method for the file
     * @return a temporary file according to the specifications
     */
    protected static File createTempFile(File location, File tmpDir, String name, String extension, boolean deleteOnExit) {

        // get location
        if (location == null)
            location = tmpDir;
        else {
            if (!location.isDirectory())
                location = location.getParentFile();
            if (!location.canWrite())
                location = tmpDir;
        }

        // get name
        if (name == null)
            name = FluxCapacitor.class.getSimpleName();
        else
            name = FluxCapacitor.class.getSimpleName() + "_" + name;

        File f = null;
        try {
            f = FileHelper.createTempFile(name, extension, location);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return createFile(f, deleteOnExit);
    }

    /**
     * Control gateway for file creation from the main class,
     * adds a hook for delete on exit in case.
     *
     * @param f            the file that has been created
     * @param deleteOnExit flag to mark for deletion on exit
     * @return
     */
    protected static File createFile(File f, boolean deleteOnExit) {
        if (deleteOnExit)
            f.deleteOnExit();

        return f;
    }

    /**
     * Returns the file with linear programs, respectively creates
     * a file from the name from the names of the input files.
     *
     * @return the name of the file with systematic biases
     */
    public File getFileLP() {
        if (fileLPdir == null) {
            fileLPdir = createTempFile(null,
                    settings.get(FluxCapacitorSettings.TMP_DIR),
                    getCompositeFName() + FluxCapacitorConstants.SFX_LP,
                    null,
                    false);
            if (fileLPdir.exists())
                fileLPdir.delete();
            fileLPdir.mkdir();
        }

        return fileLPdir;
    }

    /**
     * Set the parameter file
     *
     * @param file parameter file
     */
    public void setFile(File file) {
        this.file = file;
    }


    @Override
    public String getName() {
        return "capacitor";
    }

    @Override
    public String getDescription() {
        return "The Flux Capacitor";
    }

    @Override
    public String getLongDescription() {
        return null;
    }

    @Override
    public List<Parameter> getParameter() {
        List<Parameter> parameters = JSAPParameters.getJSAPParameter(new FluxCapacitorSettings());
        parameters.add(JSAPParameters.switchParameter("profile").type(String.class).help("do the profiling").get());
        parameters.add(JSAPParameters.flaggedParameter("parameter", 'p').type(File.class).help("specify parameter file (PAR file)").valueName("file").get());
        parameters.add(JSAPParameters.switchParameter("printParameters").help("Print default parameters").get());
        return parameters;
    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        commandLineArgs = args;
        setPrintParameters(args.userSpecified("printParameters"));
        setFile(args.getFile("parameter"));

        Iterator errors = args.getErrorMessageIterator();
        boolean haserror = false;
        while(errors.hasNext()){
            Log.error("Error parsing command line argument : " + errors.next());
            haserror = true;
        }
        if(haserror){
            return false;
        }

        if (isPrintParameters()) {
            FluxCapacitorSettings settings = new FluxCapacitorSettings();
            settings.write(System.out);
            return false;
        }

        if (args.userSpecified("profile")) {
            currentTasks.add(Task.PROFILE);
        }

        if (getFile() != null && !getFile().canRead()) {
            Log.error("");
            Log.error("Parameter file " + getFile().getAbsolutePath() + " does not exist or I can not read it!");
            Log.error(barna.commons.system.OSChecker.NEW_LINE);
            return false;
        }

        return true;
    }

    /**
     * Returns the flag wheter parameter names and explanations are to be output.
     *
     * @return flag wheter parameter names and explanations are to be output.
     */
    private boolean isPrintParameters() {
        return printParameters;
    }

    /**
     * Processes one locus encapsulated as a thread.
     *
     * @param gene      the locus
     * @param mappings the mappings in the locus
     * @param tasks     the tasks to be performed
     */
    /*private void solve(Gene gene, MSIterator mappings, EnumSet tasks) {

        // create LP and solve
        LocusSolver lsolver = new LocusSolver(gene, mappings, tasks, output, pairedEnd, stranded, settings, profile);
        if (maxThreads > 1) {
            //Thread outThread= new Thread(lsolver);
            Thread lastThread = getLastThread();
            //		int retry= 0;
            while (threadPool.size() >= maxThreads)
                //			|| (retry< maxThreads&& Runtime.getRuntime().freeMemory()< (0.10* Runtime.getRuntime().maxMemory())))
                try {
                    //System.err.println(Runtime.getRuntime().freeMemory()+"<"+ (0.25* Runtime.getRuntime().maxMemory()));
                    //				++retry;
                    //System.gc();

                    //Thread.currentThread().sleep(10); // polling bad
                    lastThread.join();

                    //				if (threadPool.size()< maxThreads&& retry> maxThreads)
                    //					break;
                } catch (InterruptedException e) {
                    ; // :)
                }

            lastThread = getLastThread();
            lsolver.setThreadBefore(lastThread);
            //		synchronized(FluxCapacitor.this.threadPool) {
            threadPool.add(lsolver);
            //		}
            lsolver.start(); //;

        } else
            lsolver.run();
    }*/

    /**
     * Returns the thread last added to the thread pool,
     * or <code>null</code> if the pool is empty.
     *
     * @return
     */
    private Thread getLastThread() {
        synchronized (FluxCapacitor.this.threadPool) {
            if (this.threadPool.size() > 0)
                return this.threadPool.get(this.threadPool.size() - 1);
            else
                return null;
        }
    }

    /**
     * Access the capacitor settings
     *
     * @return settings the capacitor settings
     */
    public FluxCapacitorSettings getSettings() {
        return settings;
    }

    /**
     * Constructs default settings, or settings based on the provided parameter file and command line arguments.
     * @param file handle for the parameter file, might be <code>null</code>
     * @param commandLineArgs parsed command line arguments
     * @return new settings instance
     * @throws ParameterException if the parameters provided do not correspond to the specification
     */
    public static FluxCapacitorSettings createSettings(File file, JSAPResult commandLineArgs) throws ParameterException {

        FluxCapacitorSettings settings= null;

        // load parameters
        if (file == null) {
            if (settings == null) {
                // create default settings
                settings = new FluxCapacitorSettings();
                FluxCapacitorSettings.relativePathParser.setParentDir(new File(""));
            }
        } else {
            if (!file.exists())
                throw new RuntimeException("I have no parameter file and I want to scream!");
            try {
                settings = FluxCapacitorSettings.createSettings(file);
            } catch (Exception e) {
                throw new RuntimeException("Unable to load settings from " + file + "\n\n " + e.getMessage(), e);
            }
        }

        // add command line parameters
        if (commandLineArgs != null) {
            for (barna.commons.parameters.Parameter p : settings.getParameters().values()) {
                if (p.getLongOption()!=null && commandLineArgs.userSpecified(p.getLongOption())) {
                    String value = commandLineArgs.getObject(p.getLongOption()).toString();
                    Object parsed = p.parse(value);
                    settings.set(p, parsed);
                }
            }
        }

        // validate the settings
        settings.validate();

        // prepare output file
        File f = settings.get(FluxCapacitorSettings.STDOUT_FILE);
        if (f!= null&& f.exists() && !CommandLine.confirm(
                "[CAUTION] I overwrite the output file " +
                        settings.get(FluxCapacitorSettings.STDOUT_FILE).getName() +
                        ", please confirm:\n\t(Yes,No,Don't know)")) {
            exit(-1);
        }
        if (settings.get(FluxCapacitorSettings.STDOUT_FILE) != null) {
            try {
                Log.outputStream = new PrintStream(new FileOutputStream(f));
            } catch (FileNotFoundException e) {
                Log.warn("Cannot write log file to " + f.getAbsolutePath());    // let it on stderr?!
            }
        }

        return settings;
    }

    /**
     * Merges and smoothens bias profiles.
     *
     */
    Profile getProfile() {
        Log.info("PROFILE","Loading profile");
        if (uniform) {
            profile = new Profile();
            profile.fill();
        } else {
            if (fileProfile != null && fileProfile.exists()) {
                profile = BiasProfiler.readProfile(fileProfile, true);
                /*if (profile != null) {
                    System.err.println("\tsmoothing..");
                    for (int i = 0; i < profile.masters.length; i++) {
                        int w = profile.masters[i].sense.length / 5;
                        profile.masters[i].sums =
                                Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV,
                                        w, profile.masters[i].sense);
                        profile.masters[i].suma =
                                Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV,
                                        w, profile.masters[i].asense);
                    }
                }*/
                }
            if (profile == null) {
                profile = new Profile();
                try {
                    BiasProfiler profiler = new BiasProfiler(this, strand, settings.isPaired(), !settings.get(FluxCapacitorSettings.DISABLE_MULTIMAP_WEIGHTING),gtfReader, mappingReader);
                    profile = profiler.call();
                } catch (Throwable e) {
                    e.printStackTrace();
                    if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP)
                        System.err.println("[FATAL] Error occured during scanning\n\t" + e.getMessage());
                }
                //writeProfiles(); // TODO
            }
        }

        if (profile == null)
            return null;
        // check
        for (int i = 0; i < profile.getMasters().length; i++) {
            if (profile.getMasters()[i].hasEmptyPositions())
                profile.getMasters()[i].fill();
        }


        return profile;
    }


    /**
     * Return the number of mappings (sense and/or anti-sense) falling in exonic edges of
     * a certain transcript signature, before or after deconvolution.
     *
     * @param v          edge set of a locus
     * @param dir        directionality of the mappings to be retrieved, sense or anti-sense
     * @param sig        signature of the transcript(s) with which exons the mappings have to
     *                   coincide
     * @param normalized flag to indicate whether observed mappings (<code>false</code>) or
     *                   predicted mappings after deconvolution (<code>true</code>) are counted
     * @return the number of mappings according to the specified paramters
     * @deprecated marked for removal
     */
    @Override
    public double getReads(Vector<AbstractEdge> v, byte dir, long[] sig, boolean normalized) {
        int sum = 0;
        for (int i = 0; i < v.size(); i++) {
            AbstractEdge e = v.elementAt(i);
            long[] inter = SplicingGraph.intersect(e.getTranscripts(), sig);
            if (SplicingGraph.isNull(inter) || !e.isExonic())
                continue;

            if (settings.isPaired()) {
                for (int j = 0; e.getSuperEdges() != null && j < v.elementAt(i).getSuperEdges().size(); j++) {
                    SuperEdge se = e.getSuperEdges().elementAt(j);
                    if (!se.isPend())
                        continue;
                    int cnt = 0;
                    for (int k = 0; k < se.getEdges().length; k++)
                        if (se.getEdges()[k] == v.elementAt(i))
                            ++cnt;
                    if (dir >= 0)
                        sum += cnt * ((MappingsInterface) se).getMappings().getReadNr();
                    if (dir <= 0)
                        sum += cnt * ((MappingsInterface) se).getMappings().getRevReadNr();
                }
            } else {
                if (dir >= 0)
                    sum += ((MappingsInterface) e).getMappings().getReadNr();
                if (dir <= 0)
                    sum += ((MappingsInterface) e).getMappings().getRevReadNr();
            }
        }
        return sum;
    }

    /**
     * Return the number of mappings (sense and/or anti-sense) falling in exonic edges of
     * a certain transcript signature, before or after deconvolution.
     *
     * @param v          edge set of a locus
     * @param dir        directionality of the mappings to be retrieved, sense or anti-sense
     * @param g          splicing graph underlying the locus
     * @param sig        signature of the transcript(s) with which exons the mappings have to
     *                   coincide
     * @param excl       flag to indicate whether exclusively edges with exactly the provided transcript
     *                   signature are considered (<code>true</code>)
     * @param normalized flag to indicate whether observed mappings (<code>false</code>) or
     *                   predicted mappings after deconvolution (<code>true</code>) are counted
     * @return the number of mappings according to the specified paramters
     * @deprecated marked for removal
     */
    public double getReadsAvg(Vector<AbstractEdge> v, byte dir, SplicingGraph g, long[] sig, boolean excl, boolean normalized) {
        double sum = 0;
        for (int i = 0; i < v.size(); i++) {
            AbstractEdge e = v.elementAt(i);
            long[] trpts = v.elementAt(i).getTranscripts();
            long[] inter = SplicingGraph.intersect(trpts, sig);
            if (SplicingGraph.isNull(inter) || (excl && !SplicingGraph.equalSet(sig, trpts)) || !e.isExonic())
                continue;
            double sf = (double) g.decodeCount(v.elementAt(i).getTranscripts());
            int mult = g.decodeCount(inter);

            if (settings.isPaired()) {
                for (int j = 0; e.getSuperEdges() != null &&
                        j < e.getSuperEdges().size(); j++) {
                    SuperEdge se = e.getSuperEdges().elementAt(j);
                    if (!se.isPend())
                        continue;
                    int cnt = 0;
                    for (int k = 0; k < se.getEdges().length; k++)
                        if (se.getEdges()[k] == e)
                            ++cnt;
                    if (dir >= 0)
                        sum += (((MappingsInterface) se).getMappings().getReadNr() * mult * cnt) / sf;
                    if (dir <= 0)
                        sum += (((MappingsInterface) se).getMappings().getRevReadNr() * mult * cnt) / sf;
                }
            } else {
                if (dir >= 0)
                    sum += (((MappingsInterface) e).getMappings().getReadNr() * mult) / sf;
                if (dir <= 0)
                    sum += (((MappingsInterface) e).getMappings().getRevReadNr() * mult) / sf;
            }

            System.currentTimeMillis();
        }

        return sum;
    }


    /**
     * If instantiated, the method returns the GTF file wrapper,
     * otherwise <code>null</code>.
     *
     * @return a wrapper instance for GTF files, or <code>null</code>
     */
    private AbstractFileIOWrapper getWrapperGTF() {
        if (gtfReader == null) {

            return  getWrapperGTF(settings.get(FluxCapacitorSettings.ANNOTATION_FILE));
        }
        return gtfReader;
    }

    /**
     * If instantiated, the method returns the GTF file wrapper,
     * otherwise a new instance is created using the provided file
     * argument.
     *
     * @param inputFile the GTF file's handle
     * @return a wrapper instance for the GTF file
     */
    public static AbstractFileIOWrapper getWrapperGTF(File inputFile) {

        GTFwrapper gtfReader = new GTFwrapper(inputFile.getAbsolutePath());
        gtfReader.setNoIDs(null);
        gtfReader.setReadGene(true);
        gtfReader.setReadFeatures(new String[]{"exon", "CDS"});
        gtfReader.setReadAheadTranscripts(1000000);    // only one locus a time
//		gtfReader.setReadAheadTranscripts(-1);
//		gtfReader.setReadAll(true);
        gtfReader.setGeneWise(true);    // when setting false, clusterLoci sorts genes in another ordering
        gtfReader.setPrintStatistics(false);
        gtfReader.setReuse(true);

        //gtfReader.setSourceExclude();
        Transcript.removeGaps = false;

        return gtfReader;
    }

    /**
     * Creates a wrapper to read mapping as BED objects from a file.
     *
     * @param inputFile file with mappings in BED format
     * @return a wrapper instance providing read access to the specified file
     * @deprecated marked for removal
     */
    @Deprecated
    private AbstractFileIOWrapper getWrapperBED(File inputFile) {
		/*mappingReader= new BEDReader(inputFile.getAbsolutePath());
		return mappingReader;*/  // TODO pull up to MappingReader
        mappingReader = new BEDReader(inputFile, settings.get(FluxCapacitorSettings.SORT_IN_RAM),settings.getReadDescriptor(),settings.get(FluxCapacitorSettings.TMP_DIR));
        return null;// removed mappingReader;
    }


    /**
     * Finalizes an iteration over all mappings: closes readers, outputs stats, etc.
     *
     * @param secs  time measured for the iteration, in [sec]
     * @param stats object capturing statistics of the run
     */
    void exploreFinish(long secs, MappingStats stats) {

        try {

            while (threadPool.size() > 0 && threadPool.elementAt(0).isAlive())
                try {
                    threadPool.elementAt(0).join();
                } catch (Exception e) {
                    ; //:)
                }

            //					if (fileMappings!= null)
            //						getSammy().close();

            //assert(nrUniqueReads==getBedReader().getNrUniqueLinesRead());	// take out for cheat
            if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP) {
                System.err.println();
                System.err.println("\treconstruction finished .. took " + secs + " sec.\n\n\t"
                        + mappingReader.getCountMappings()+" mappings read from file\n\t"
                        // no info, reads in redundantly many reads
                        //+ nrReadsLoci+" mappings in annotated loci regions\n\t"
                        + stats.getMappingsMapped() + " mapping" + (settings.isPaired() ? " pairs" : "s") + " map to annotation\n"
                        + (settings.isPaired() ?
                        "\t" + stats.getMappingPairsNoTx() + " mapping"+ (settings.isPaired() ? " pairs" : "s") + " without transcript evidence\n"
                                + "\t" + stats.getPairsWrongOrientation() + " mapping"+ (settings.isPaired() ? " pairs" : "s") + " in wrong orientation\n"
                        //+ "\t"+ nrMappingsForced+ " single mappings forced\n"
                        : "")
                        + ((strand == FluxCapacitorConstants.STRAND_SPECIFIC) ? stats.getMappingsWrongStrand() + " mappings map to annotation in antisense direction\n\t" : "")
                        //+ nrMultiMaps+" mapped multiply.\n\n\t"
                        + (output.contains(OutputFlag.GENE) ? "\n\t" + stats.getLoci() + " loci, " + stats.getLociExp() + " detected" : "")
                        + (output.contains(OutputFlag.TRANSCRIPT) ? "\n\t" + stats.getTxs() + " transcripts, " + stats.getTxsExp() + " detected" : "")
                        + (output.contains(OutputFlag.EVENT) ? "\n\t" + stats.getEvents() + " ASevents of dimension " + eventDim + ", " + stats.getEventsExp() + " detected" : "")
                        + barna.commons.system.OSChecker.NEW_LINE
                        //+ nrUnsolved+" unsolved systems."
                );
            }

            // output stats
            /*if (stats != null) {
                stats.setMappingsTotal(mappingReader.getCountMappings());
                stats.setMappingsMapped(nrReadsMapped);
                if (pairedEnd) {
                    stats.setMappingPairsNoTx(nrPairsNoTxEvidence);
                    stats.setPairsWrongOrientation(nrPairsWrongOrientation);
                }
                if (strand == FluxCapacitorConstants.STRAND_SPECIFIC) {
                    stats.setMappingsWrongStrand(nrMappingsWrongStrand);
                }
                if (outputGene) {
                    stats.setLociExp(nrLociExp);
                }
                if (outputTranscript) {
                    stats.setTxsExp(nrTxExp);
                }
                if (outputEvent) {

                    stats.setEventsExp(nrEvents);
                }
            }*/

        } catch (Exception e) {
            throw new RuntimeException(e);
        }


    }


    /**
     * A complete iteration cycle over all mappings, either for profiling or for deconvolution.
     *
     */
    public boolean explore(Gene[] genePreclustered) {

        /*if (!settings.get(FluxCapacitorSettings.STATS_FILE_APPEND))
            stats.reset();*/

        long t0 = System.currentTimeMillis();
        try {

            Transcript.setEdgeConfidenceLevel(Transcript.ID_SRC_MOST_INCONFIDENT);
            //this.gtfReader= null;
            //GFFReader gtfReader= getGTFreader();
            gtfReader.reset();
            gtfReader.setSourceInclude(null);
            gtfReader.setSourceExclude(null);
			mappingReader.reset();

            if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP) {
                    if (currentTasks.contains(Task.COUNT_INTRONS)||currentTasks.contains(Task.COUNT_SJ)) {
                        StringBuilder message = new StringBuilder();
                        message.append("Counting reads to the following elements: ");
                        if (currentTasks.contains(Task.COUNT_SJ)) {
                            message.append("splice junctions");
                            if (currentTasks.contains(Task.COUNT_INTRONS))
                                message.append(", ");
                        }
                        if (currentTasks.contains(Task.COUNT_INTRONS))
                            message.append("all-intronic regions");
                        message.append(".");
                        Log.info("","");
                        Log.info("COUNT", message.toString());
                    }
                    if (currentTasks.contains(Task.DECOMPOSE)) {
                        Log.info("","");
                        Log.info("SOLVE", "Deconvolving reads of overlapping transcripts.");
                    }
                }

            Log.progressStart("deconvolving");


            // TODO BARNA-112 disable keeping original lines

            gtfReader.read();
            Gene[] genes = gtfReader.getGenes();

            long tlast = System.currentTimeMillis();
            boolean output = false;

            String lastChr = null;
            byte lastStr = 0;
            int lastEnd = -1;

            if (genes != null) {
                lastChr = genes[0].getChromosome();
                lastStr = genes[0].getStrand();
            }

            Thread readerThread = null;
            int readObjects = 0;
            for (;lastChr != null&& genes!= null;genes = gtfReader.getGenes()) {    // MAIN LOOP

                /*if (gtfReader.getVLines() != null)
                    origLines = (Vector<String>) gtfReader.getVLines().clone();    // TODO make array, trim..*/

                // http://forums.sun.com/thread.jspa?threadID=5171135&tstart=1095
                if (readerThread == null)
                    readerThread = new GtfReaderThread();
                //readerThread.start();
                readerThread.run();

                for (int i = 0; (genes != null) && i < genes.length; i++) {


                    // flop strand
                    if (lastChr.equals(genes[i].getChromosome())) {
                        if (lastStr != genes[i].getStrand()) {
                            //System.err.println(lastChr+" "+lastStr+ " "+ readObjects+ " wrote "+ dbgCntWriteMap +" not "+ dbgCntWriteNonmap);
                            readObjects = 0;
                            // jump back
								mappingReader.reset(genes[i].getChromosome());
                            lastStr = genes[i].getStrand();
                            lastEnd = -1;
                        }
                    } else {                        // flop chr
                        //System.err.println(lastChr+" "+lastStr+ " "+ readObjects+ " wrote "+ dbgCntWriteMap +" not "+ dbgCntWriteNonmap);
                        readObjects = 0;
                        lastChr = genes[i].getChromosome();
                        lastStr = genes[i].getStrand();
                        lastEnd = -1;
                    }

                    // do it..
                    quantify(genes[i], mappingReader, currentTasks);

                    if (output) {
                        System.out.println(genes[i].getChromosome() + " " +
                                genes[i].getStrand() +
                                " cluster " + genes[i].getLocusID());
                        // TODO beds.size() no longer available
                        //", "+beds.size()+" reads.");
                        if ((lastStr != genes[i].getStrand()
                                || !(lastChr.equals(genes[i].getChromosome())))) {
                            long t = System.currentTimeMillis();
                            if (lastStr != 0 && (!lastChr.equals("")))
                                System.out.println(lastChr + " " + lastStr +
                                        " " + ((t - tlast) / 1000) + " sec.");
                            tlast = t;
                            lastStr = genes[i].getStrand();
                            lastChr = genes[i].getChromosome();
                        }
                    }

                }
                //getWriter().flush();


            }    // end iterate GTF

				//mappingReader.finish(); //TODO check

            while (threadPool.size() > 0 && threadPool.elementAt(0).isAlive())
                try {
                    threadPool.elementAt(0).join();
                } catch (Exception e) {
                    ; //:)
                }
            Log.progressFinish(StringUtils.OK, true);
            if (checkGTFscanExons > 0 && checkGTFscanExons != gtfReader.getNrExons())
                System.err.println("[ERROR] consistency check failed in GTF reader: " + checkGTFscanExons + "<>" + gtfReader.getNrExons());
            checkGTFscanExons = gtfReader.getNrExons();
				if (checkBEDscanMappings> 0&& checkBEDscanMappings!= mappingReader.getCountMappings())
					System.err.println("[ERROR] consistency check failed in BED reader "+ checkBEDscanMappings+ "<>"+ mappingReader.getCountMappings());
            //checkBEDscanMappings= getBedReader().getNrLines();

            // close readers, output
            exploreFinish(((System.currentTimeMillis() - t0) / 1000), stats);

        } catch (Exception e1) {
            Log.error("Error while iterating loci:", e1);
            throw new RuntimeException(e1);
        }

        return true;
    }

    protected void quantify(Gene gene, MappingReader mappingReader, EnumSet<Task> currentTasks) {

        // boundaries
        int start = gene.getStart();
        int end = gene.getEnd();

        if (gene.getStrand() < 0) {
            start = -start;
            end = -end;
        }
        int tol = this.tolerance; // 1000
        tol = 0;
        start = Math.max(1, start - tol);
        end = end + tol;

        MSIterator<Mapping> mappings= mappingReader.read(gene.getChromosome(), start, end);
        LocusSolver lsolver = new LocusSolver(gene, mappings, currentTasks, this.output, settings.isPaired(),
                settings.isStranded(), settings, profile);

        try {
            stats.addLocus(lsolver.call());
        } catch (Exception e) {
            Log.error("Error during deconvolution: "+ e.getMessage());
            throw new RuntimeException(e);
        }

        if (mappings != null) {
            mappings.clear();
        }
    }

    /**
     * Checks whether an input file is to be uncompressed and/or sorted before reading.
     *
     * @param inputFile file which is to be read
     * @return a wrapper instance that accesses the file, after potential decompression and/or sorting
     */
    public static AbstractFileIOWrapper fileInit(File inputFile, FluxCapacitorSettings settings) {

        // (1) unpack, if compressed
        byte cb = FileHelper.COMPRESSION_NONE;
        if (!FileHelper.getExtension(inputFile).toUpperCase().equals("BAM")) {
            cb = FileHelper.getCompression(inputFile);
            if (cb != FileHelper.COMPRESSION_NONE) {
                File f = new File(FileHelper.stripExtension(inputFile.getAbsolutePath()));
                if (f.exists()) {
                    Log.println("Assuming file " + f.getName() + " is a decompressed version of " + inputFile.getName());
                } else {
                    f = createTempFile(null, settings.get(FluxCapacitorSettings.TMP_DIR),
                            FileHelper.getFileNameWithoutExtension(f), FileHelper.getExtension(f), true);
                    try {
                        FileHelper.inflate(inputFile, f, cb);
                    } catch (Exception e) {
                        throw new RuntimeException(e);
                    }
                }
                inputFile = f;
            }
        }
        // (2) sort, if needed
        AbstractFileIOWrapper wrapper = getWrapper(inputFile, settings);
        if (!wrapper.isApplicable()) {
            File sortedDir = settings.get(FluxCapacitorSettings.KEEP_SORTED);
            File f;
            if (sortedDir!=null)
                f = FileHelper.getSortedFile(new File(sortedDir, inputFile.getName()));
            else
                f = FileHelper.getSortedFile(inputFile);
            File lock = FileHelper.getLockFile(f);

            if (f.exists() && !lock.exists()) {

                Log.warn("Assuming file " + f.getName() + " is a sorted version of " + inputFile.getName());

            } else {    // we have to sort

                boolean lockCreated = false;
                if (sortedDir!=null) {//settings.get(FluxCapacitorSettings.KEEP_SORTED)) {    // try to store in original

                    if (lock.exists()) {    // switch to sorting to temp
                        Log.warn("Seems that another process is just sorting file " + inputFile +
                                "\nremove lock file " + lock.getName() + " if dead leftover." +
                                "\nContinuing with sorting to temporary file " +
                                (f = createTempFile(f,     // access to non-Temp
                                    settings.get(FluxCapacitorSettings.TMP_DIR),
                                        FileHelper.getFileNameWithoutExtension(f),
                                        FileHelper.getExtension(f),
                                        false)).getAbsolutePath());

                    } else if (!f.getParentFile().canWrite()) {    // sort to temp, but do not delete (parameter)
                        Log.warn("Cannot write sorted file to " + f.getAbsolutePath() +
                                "\nContinuing with sorting to temporary file " +
                                (f = createTempFile(f, // access to non-Temp
                                    settings.get(FluxCapacitorSettings.TMP_DIR),
                                        FileHelper.getFileNameWithoutExtension(f),
                                        FileHelper.getExtension(f),
                                        false)).getAbsolutePath());

                    } else {    // sort to default sorted file
                        try {
                            lock.createNewFile();
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }
                        lockCreated = true;
                    }

                } else {    // do not keep sorted files, sort to temp and delete on exit
                    f = createTempFile(null,
                        settings.get(FluxCapacitorSettings.TMP_DIR),
                            FileHelper.getFileNameWithoutExtension(f),
                            FileHelper.getExtension(f),
                            true);
                }

                // doit
                if (wrapper.getInputFile() != null)
                    Log.info("Sorting " + wrapper.getInputFile().getAbsolutePath());
                wrapper.sort(f);

                // if locked
                if (lockCreated)
                    lock.delete();

            }
            inputFile = f;
            wrapper = getWrapper(inputFile, settings);
        }

        return wrapper;
    }

    /**
     * Obtains global statistics from the mapping file, e.g., number of total mappings etc.
     *
     * @param reader mapping file reader
     */
	private void fileStats(MappingReader reader) {

        // (3) scan
        ((AbstractFileIOWrapper) reader).scanFile();
        if (((AbstractFileIOWrapper) reader).getNrInvalidLines() > 0)
            Log.warn("Skipped " + ((AbstractFileIOWrapper) reader).getNrInvalidLines() + " lines.");

        // ensure sync between paired and stranded annotation mapping:
        // mappingReader only knows now whether there are paired reads
        if (settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).equals(AnnotationMapping.AUTO)) {
            settings.setAnnotationMappingAuto(mappingReader.isPaired());
        }
        if (settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).isPaired() && !mappingReader.isPaired())
            throw new RuntimeException("Annotation mapping " + settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING) +" requires paired reads");

        checkBEDscanMappings = reader.getCountMappings();
        stats.setReadsTotal(reader.getCountReads());
        stats.setMappingsTotal(reader.getCountMappings());

		Log.info("\t"+ stats.getReadsTotal() + " mapped reads"
                + (stats.getMappingsTotal() > 0 ? ", " + stats.getMappingsTotal() + " mappings: R-factor " + (reader.getCountMappings() / (float) reader.getCountReads()) : ""));
        if (stats.getMappingsTotal() > 0)
            Log.info("\t" + reader.getCountContinuousMappings() + " entire, " + reader.getCountSplitMappings()
                    + " split mappings (" + (reader.getCountSplitMappings() * 100f / reader.getCountMappings()) + "%)");

        // SAM/BAM wrapper prohibits to specify read descriptors
        //if (mappingReader instanceof SAMReader) {
        //    if (mappingReader.isPaired())
        //        settings.set(FluxCapacitorSettings.READ_DESCRIPTOR, new UniversalReadDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED));
        //    else
        //        settings.set(FluxCapacitorSettings.READ_DESCRIPTOR, new UniversalReadDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMPLE));
        //}


        // (4) check if read descriptor is applicable
        String ext = FileHelper.getExtension(settings.get(FluxCapacitorSettings.MAPPING_FILE)).toUpperCase();
        if(ext.equals("BED") || ext.equals("GZ")) {
            if (reader.isApplicable(settings.get(FluxCapacitorSettings.READ_DESCRIPTOR)))
                Log.info("\tRead descriptor seems OK");
            else {
                String msg = "Read Descriptor " + settings.get(FluxCapacitorSettings.READ_DESCRIPTOR)
                        + " incompatible with read IDs";
                Log.error(msg);
                throw new RuntimeException(msg);
            }
        }
        if (reader.isApplicable(settings.getReadDescriptor()))
            Log.info("\tRead descriptor seems OK");
        else {
            String msg = "Read Descriptor " + settings.getReadDescriptor()
                    + " incompatible with read IDs";
            Log.error(msg);
            throw new RuntimeException(msg);
        }

    }

    /**
     * Initializes a wrapper for an input file, depending on the file extension.
     *
     * @param inputFile the file from which is read
     * @return a wrapper instance for reading from the specified file
     */
    public static AbstractFileIOWrapper getWrapper(File inputFile, FluxCapacitorSettings settings) {

        String ext = FileHelper.getExtension(inputFile).toUpperCase();
        SupportedFormatExtensions sup = null;
        try {
            sup = SupportedFormatExtensions.valueOf(ext);
        } catch (Exception e) {
            throw new RuntimeException("Unsupported file format " + ext);
        }


        switch (sup) {
            case GTF:
            case GFF:
                return getWrapperGTF(inputFile);
            case BED:
                return new BEDReader(inputFile, settings.get(FluxCapacitorSettings.SORT_IN_RAM),settings.get(FluxCapacitorSettings.READ_DESCRIPTOR),settings.get(FluxCapacitorSettings.TMP_DIR), settings.get(FluxCapacitorSettings.MIN_SCORE));
            case BAM:
                SAMReader r = new SAMReader(inputFile, true, settings.get(FluxCapacitorSettings.SORT_IN_RAM), settings.get(FluxCapacitorSettings.MIN_SCORE), !settings.get(FluxCapacitorSettings.IGNORE_SAM_FLAGS), settings.get(FluxCapacitorSettings.SAM_PRIMARY_ONLY), !settings.get(FluxCapacitorSettings.IGNORE_SAM_PAIRING_INFORMATION), settings.get(FluxCapacitorSettings.SAM_UNIQUE_ONLY));
                if (!settings.get(FluxCapacitorSettings.SAM_VALIDATION_STRINGENCY).equals(SAMFileReader.ValidationStringency.DEFAULT_STRINGENCY)) {
                    Log.info("SAM","Setting validation stringency to " + settings.get(FluxCapacitorSettings.SAM_VALIDATION_STRINGENCY));
                    r.setValidationStringency(settings.get(FluxCapacitorSettings.SAM_VALIDATION_STRINGENCY));



                }
                return r;
            default:
                return null;
        }

        //return null;    // make compiler happy
    }

    /**
     * Returns the parameter file.
     *
     * @return the parameter file
     */
    public File getFile() {
        return file;
    }

    /**
     * Enable or disable default parameter printing.
     *
     * @param printParameters if <code>true</code>, parameter printing is enables,
     *                        otherwise it is disabled
     */
    public void setPrintParameters(final boolean printParameters) {
        this.printParameters = printParameters;
    }



}
 