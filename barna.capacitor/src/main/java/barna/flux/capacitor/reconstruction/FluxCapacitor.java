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
import barna.commons.launcher.FluxTool;
import barna.commons.log.Log;
import barna.commons.utils.StringUtils;
import barna.flux.capacitor.graph.AnnotationMapper;
import barna.flux.capacitor.graph.MappingsInterface;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping;
import barna.genome.lpsolver.LPSolverLoader;
import barna.io.*;
import barna.io.bed.BEDDescriptorComparator;
import barna.io.bed.BEDwrapper;
import barna.io.gtf.GTFwrapper;
import barna.io.rna.UniversalReadDescriptor;
import barna.io.state.MappingWrapperState;
import barna.model.*;
import barna.model.bed.BEDMapping;
import barna.model.commons.Coverage;
import barna.model.commons.MyFile;
import barna.model.constants.Constants;
import barna.model.gff.GFFObject;
import barna.model.splicegraph.AbstractEdge;
import barna.model.splicegraph.SimpleEdge;
import barna.model.splicegraph.SplicingGraph;
import barna.model.splicegraph.SuperEdge;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import lpsolve.LpSolve;
import lpsolve.VersionInfo;

import java.io.*;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;


/**
 * The Flux Capacitor class performes a deconvolution of reads falling into common areas of transcripts.
 *
 * @author Micha Sammeth (gmicha@gmail.com)
 */
public class
        FluxCapacitor implements FluxTool<FluxCapacitorStats>, ReadStatCalculator {

    /**
     * Enumerates possible tasks for the FluxCapacitor
     */
    private enum Task {
        COUNT_INTRONS, COUNT_SJ, LEARN, DECOMPOSE
    };

    /**
     * Task to be executed in the current run - initialized empty
     */
    private EnumSet<Task> currentTasks = EnumSet.noneOf(Task.class);

    /**
     * Store a reference to the parsed command line arguments
     * to update settings
     */
    private JSAPResult commandLineArgs;

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
     * Comparator for comparing read identifiers according to the provided descriptor.
     */
    MappingComparator comp = null;

    /**
     * Returns an instance for comparing read identifiers according to the provided descriptor.
     *
     * @return instance for comparing read identifiers according to the provided descriptor
     */
    private Comparator<? super Mapping> getDescriptorComparator() {
        if (comp == null) {
            comp = new MappingComparator(
                    settings.get(FluxCapacitorSettings.READ_DESCRIPTOR));
        }

        return comp;
    }

    /**
     * A class that encapsulates all information necessary to carry out the deconvolution
     * of the reads in a locus.
     */
    class LocusSolver extends Thread { //TODO implement Callable interface and make top level class (or static member class)

        /**
         * The locus that is to be solved.
         */
        Gene gene = null;

        /**
         * AStalavista events found in the locus.
         */
        ASEvent[] events = null;

        /**
         * Iterator over mappings.
         */
        BufferedIterator beds = null;

        /**
         * EnumSet indicating which task(s) has(have) to be performed in the current run.
         */
        EnumSet<Task> tasks;

        /**
         * The coverage profile of systematic biases.
         */
        Coverage coverage = null;

        /**
         * Used for chaining threads to be executed sequentially.
         */
        Thread threadBefore = null;

        /**
         * The number of annotation-mapped mappings or reads.
         */
        int nrMappingsReadsOrPairs;

        /**
         * Variable to store invariant of observed split frequency.
         */
        private float invariantTestObsSplitFreq = 0;

        /**
         * Variable to store invariant of predicted split frequency.
         */
        private float invariantTestPredSplitFreq = 0;


        /**
         * Constructor providing reads and mappings for deconvolution.
         * The mode of the run can be switched between profiling and deconvolution.
         *
         * @param newGene   the locus model
         * @param newBeds   the mappings that fall in the locus
         * @param tasks     tasks to be preformed
         */
        public LocusSolver(Gene newGene, BufferedIterator newBeds, EnumSet tasks) {

            this.gene = newGene;
            this.beds = newBeds;
            this.tasks = tasks;

            nrMappingsReadsOrPairs = 0;
        }


        /**
         * Launches the profiling/deconvolution on the locus
         */
        @Override
        public void run() {

            AnnotationMapper mapper = null;

            if (!tasks.contains(Task.LEARN)) {
                mapper = new AnnotationMapper(this.gene);
                mapper.map(this.beds, settings);

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
                    case LEARN:
                        if (this.gene.getTranscriptCount() == 1) {
                            ++nrSingleTranscriptLearn;
                            learn(this.gene.getTranscripts()[0], beds);
                        }
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

            beds = null;
            gene = null;
            // makes it terribly slow
            //System.gc();

//			synchronized(FluxCapacitor.this.threadPool) {
            FluxCapacitor.this.threadPool.remove(this);
//			}
        }

        /**
         * Count reads to splice junction within the current locus and output them in GTF format.
         *
         * @param gene current locus
         * @param mapper mapping graph for the current locus
         */
        private void outputSJGFF(Gene gene, AnnotationMapper mapper) {
            Map<String, Integer> m = mapper.getSJReads(settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).equals(AnnotationMapping.PAIRED) ? true : false);
            StringBuilder sb = new StringBuilder();
            for (String s : m.keySet()) {
                String[] junction = s.split("\\^");
                sb.append(gene.getChromosome());
                sb.append("\t");
                sb.append("flux");
                sb.append("\t");
                sb.append(FluxCapacitorConstants.GFF_FEATURE_JUNCTION);
                sb.append("\t");
                sb.append(junction[0].contains("-") ? junction[1].replace("-", "") : junction[0]);
                sb.append("\t");
                sb.append(junction[1].contains("-")?junction[0].replace("-",""):junction[1]);
                sb.append("\t");
                sb.append(".");
                sb.append("\t");
                sb.append(gene.getStrand() > 0 ? "+" : "-");
                sb.append("\t");
                sb.append(".");
                sb.append("\t");
                sb.append("gene_id \""+gene.getGeneID()+"\";");
                sb.append(" ");
                sb.append("locus_id \""+gene.getLocusID()+"\";");
                sb.append(" ");
                sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_READS+" "+String.format("%1$f", (float) m.get(s)));// +";");
                sb.append("\n");
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
            Map<String, Float[]> m = mapper.getAllIntronicReads(settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).equals(AnnotationMapping.PAIRED) ? true : false);
            StringBuilder sb = new StringBuilder();
            for (String s : m.keySet()) {
                String[] intron = s.split("\\^");
                sb.append(gene.getChromosome());
                sb.append("\t");
                sb.append("flux");
                sb.append("\t");
                sb.append(FluxCapacitorConstants.GFF_FEATURE_INTRON);
                sb.append("\t");
                sb.append(intron[0].contains("-")?intron[1].replace("-",""):intron[0]);
                sb.append("\t");
                sb.append(intron[1].contains("-")?intron[0].replace("-",""):intron[1]);
                sb.append("\t");
                sb.append(".");
                sb.append("\t");
                sb.append(gene.getStrand() > 0?"+":"-");
                sb.append("\t");
                sb.append(".");
                sb.append("\t");
                sb.append("gene_id \""+gene.getGeneID()+"\";");
                sb.append(" ");
                sb.append("locus_id \""+gene.getLocusID()+"\";");
                sb.append(" ");
                sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_READS+" "+String.format("%1$f", (float) m.get(s)[0].intValue()) +";");
                sb.append(" ");
                sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_FRAC_COVERED+" "+m.get(s)[1]);//+";");
                sb.append("\n");
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
            ++nrLoci;
            if (solver != null || nrMappingsReadsOrPairs > 0)
                ++nrLociExp;
            double valOF = solver == null ? 0 : solver.getValObjFunc();
            if (valOF > FluxCapacitorConstants.BIG) {
                ++nrUnsolved;
                Log.warn("Unsolved system: " + gene.getLocusID());
            }

            // pre-build rpkm hash
            HashMap<String, Double> rpkmMap = null;
            double base = (nrBEDreads < 0 ? 1 : nrBEDreads);
            Transcript[] tt = gene.getTranscripts();
            if (outputBalanced) {
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

                    if (val > 0 && !(outputObs || outputPred))
                        ++nrTxExp;

                    double rpkm = (float) ((val / (double) tx.getExonicLength()) * (1000000000l / base));
                    if (Double.isNaN(rpkm))
                        Log.warn("NaN RPKM produced: " + val + " / " + base + " = " + rpkm);

                    rpkmMap.put(tid, rpkm);
                }
            }


            // reproduce original
            boolean foundExons = true, foundTranscripts = false;
            if (((GTFwrapper) getWrapperGTF()).isKeepOriginalLines() && origLines != null) {
                foundTranscripts = outputGFForiginalLines(rpkmMap);
            }

            StringBuilder sb = new StringBuilder();
            // LOCUS TODO genes
            if (outputGene) {
                if (outputObs || outputPred) {
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

                } else if (outputBalanced) {
                }
            }


            // TRANSCRIPTS
            if (outputTranscript || outputExon || outputSJunction) {
                float invariantObsAllTx = 0, invariantPredAllTx = 0,
                        invariantObsAllEx = 0, invariantPredAllEx = 0;
                for (int i = 0; i < tt.length; i++) {
                    ++nrTx;
//					float invariantObsTx= invariantTestObsSplitFreq,
//					invariantPredTx= invariantTestPredSplitFreq;
                    String tid = tt[i].getTranscriptID();
                    float invariantObsTx = 0, invariantPredTx = 0;
                    if (outputTranscript && !foundTranscripts) {
                        if (outputObs || outputPred) {
                            //getGTF(sb, g.trpts[i], solver, g, perM, null, false);	// writer.write
                            invariantObsAllTx += invariantTestObsSplitFreq; //invariantObsTx;
                            invariantPredAllTx += invariantTestPredSplitFreq; // invariantPredTx;
                            invariantObsTx = invariantTestObsSplitFreq;
                            invariantPredTx = invariantTestPredSplitFreq;
                            if (invariantPredTx > 0)
                                ++nrTxExp;

                        } else if (outputBalanced) {

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
                                    (float) (rpkmMap.get(tid) * tt[i].getExonicLength() * (base / 1000000000l))));
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
                            sb.append("\n");
                        }
                    }
                    // EXONS
                    float invariantObsEx = 0, invariantPredEx = 0;
                    if (outputExon && !foundExons) {
                        Exon[] exons = tt[i].getExons();
                        for (int j = 0; j < exons.length; j++) {
                            //getGTF(sb, exons[j], tt[i], g, solver, unsolvedSystem, perM, null, false);
                            invariantObsEx += invariantTestObsSplitFreq;
                            invariantPredEx += invariantTestPredSplitFreq;
                        }
                    }

                    // SJ
                    if (outputSJunction) {
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

                    if (outputExon && outputSJunction && outputTranscript) {
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
                if (outputTranscript) {
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
                if (outputExon && outputSJunction) {
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
            if (outputEvent) {
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
                    if (outputObs || outputPred)
                        ; //getGTF(sb, events[i], g, solver, unsolvedSystem, perM, pv, tExpMap);
                    else
                        ++nrEvents;
                    if (outputBalanced) {
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
                        if (allPos && !(outputObs || outputPred))
                            ++nrEventsExp;
                        sb.replace(sb.length() - 1, sb.length(), "\";\n");
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
        private boolean outputGFForiginalLines(HashMap<String, Double> rpkmMap) {

            Transcript[] tt = gene.getTranscripts();
            boolean foundTranscripts = false;
            for (int i = 0; i < origLines.size(); i++) {
                String s = origLines.elementAt(i);
                String feat = GFFObject.getField(3, s);
                String tid = GFFObject.getTranscriptID(s);
                int tx = 0;
                if ((feat.equals(feat.equals(Transcript.GFF_FEATURE_TRANSCRIPT))
                        || feat.equals(Exon.GFF_FEATURE_EXON))
                        && (outputObs || outputPred))
                    for (tx = 0; tx < tt.length; tx++)
                        if (tt[tx].getTranscriptID().equals(tid))
                            break;
                if (tx >= tt.length) {
                    System.err.println("\nTranscript " + tid + " not found in: ");
                    for (int j = 0; j < tt.length; j++)
                        System.err.println("\t" + tt[j].getTranscriptID());
                    System.err.println();
                }

                if (feat.equals(Transcript.GFF_FEATURE_TRANSCRIPT) && outputTranscript) {
                    foundTranscripts = true;
                    StringBuilder sb = new StringBuilder(s);
                    int x = sb.length();    // trim
                    while (Character.isWhitespace(sb.charAt(--x)))
                        sb.delete(x, x + 1);
                    if (sb.charAt(x) != ';')
                        sb.append("; ");
                    else
                        sb.append(Constants.SPACE);

                    if ((outputObs || outputPred) && tx < tt.length)
                        ; //getGTF(sb, tt[tx], solver, g, perM, pv, true);
                    else if (outputBalanced) {
                        sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
                        sb.append(Constants.SPACE);
                        if (rpkmMap.containsKey(tid))
                            sb.append(String.format("%1$f", rpkmMap.get(tid).floatValue()));    // rgasp parser does not like scientific notation
                        else
                            sb.append(Constants.NULL);
                        sb.append(";\n");
                    }

                    Log.print(sb.toString());

                } else if (feat.equals(Exon.GFF_FEATURE_EXON) && outputExon) {

                    StringBuilder sb = new StringBuilder(s);
                    int x = sb.length();
                    while (Character.isWhitespace(sb.charAt(--x)))
                        sb.delete(x, x + 1);
                    if (sb.charAt(x) != ';')
                        sb.append("; ");
                    else
                        sb.append(' ');


                    if ((outputObs || outputPred) && tx < tt.length) {
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
                    } else if (outputBalanced) {
                        sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
                        sb.append(Constants.SPACE);
                        if (rpkmMap.containsKey(tid))
                            sb.append(String.format("%1$f", rpkmMap.get(tid).floatValue()));    // rgasp parser does not like scientific notation
                        else
                            sb.append(Constants.NULL);
                        sb.append(";\n");
                    }


                    Log.print(sb.toString());
                } else if (outputUnknown) {
                    Log.print(s + System.getProperty("line.separator"));
                }
            }

            return foundTranscripts;
        }


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

            GraphLPsolver solver = new GraphLPsolver(mapper, readLenMin,
                    pairedEnd ? insertMinMax : null, mappedReads,
                    strand == FluxCapacitorConstants.STRAND_ENABLED,
                    pairedEnd);
            if (outputLP)
                solver.setFileLPdir(getFileLP());
            solver.costModel = costModel;    // COSTS_LINEAR
            solver.setCostSplit(costSplit);
            solver.setProfile(profile);
            solver.setReadLen(readLenMin);
            solver.costBounds = costBounds;

            return solver;
        }


        /**
         * Learns systematic biases along a transcript
         *
         * @param tx   the Transcript
         * @param beds the mappings
         */
        private void learn(Transcript tx, BufferedIterator beds) {

            if (beds == null)
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

            while (beds.hasNext()) {

                ++nrReadsSingleLoci;
                bed1 = new BEDMapping(beds.next());
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
                        ++nrMappingsWrongStrand;
                        continue;
                    }
                }

                int bpoint1 = getBpoint(tx, bed1);
                if (bpoint1 < 0 || bpoint1 >= elen) {    // outside tx area, or intron (Int.MIN_VALUE)
                    ++nrReadsSingleLociNoAnnotation;
                    continue;
                }

                ++nrReadsSingleLociMapped;    // the (first) read maps

                if (pairedEnd) {

                    beds.mark();
                    while (beds.hasNext()) {
                        bed2 = new BEDMapping(beds.next());
                        attributes2 = settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).getAttributes(bed2.getName(), attributes2);
                        if (attributes2 == null)
                            continue;
                        if (!attributes.id.equals(attributes2.id))
                            break;
                        if (attributes2.flag == 1)    // not before break, inefficient
                            continue;

                        int bpoint2 = getBpoint(tx, bed2);
                        if (bpoint2 < 0 || bpoint2 >= elen) {
                            ++nrReadsSingleLociNoAnnotation;
                            continue;
                        }

                        // check again strand in case one strand-info had been lost
                        if (stranded) {
                            if ((tx.getStrand() == bed2.getStrand() && attributes2.strand == 2)
                                    || (tx.getStrand() != bed2.getStrand() && attributes2.strand == 1)) {
                                ++nrMappingsWrongStrand;
                                continue;
                            }
                        }

                        // check directionality (sequencing-by-synthesis)
                        if ((bed1.getStrand() == bed2.getStrand())
                                || ((bed1.getStart() < bed2.getStart()) && (bed1.getStrand() != DirectedRegion.STRAND_POS))
                                || ((bed2.getStart() < bed1.getStart()) && (bed2.getStrand() != DirectedRegion.STRAND_POS))) {
                            nrPairsWrongOrientation += 2;
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

                        nrReadsSingleLociPairsMapped += 2;

                    }
                    beds.reset();

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
                        pairedEnd ? nrReadsSingleLociPairsMapped : nrReadsSingleLociMapped,
                        coverage.getFractionCovered(),
                        coverage.getChiSquare(true),
                        coverage.getCV(true));
            }
        }


    } /* End of LocusSolver*/


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
        GTF, GFF, BED,
    }

    /**
     * The parameter file.
     */
    protected File file = null;

    /**
     * Flag that indicates whether annotation mapping enforces read paring (or not).
     */
    public boolean pairedEnd = false;

    /**
     * Flag that indicates whether annotation mapping enforces correct strand (or not).
     */
    public boolean stranded = false;

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
    int nrBEDreads = -1;

    /**
     * The total number of mappings.
     */
    int nrBEDmappings = -1;

    /**
     * The number of mappings read in a second run, invariant check.
     */
    private int checkBEDscanMappings = 0;

    /**
     * Temporary file for coverage statistics of the 5' to 3' read distribution.
     */
    File fileTmpCovStats = null;

    /**
     * Writer of the coverage statistics of the 5' to 3' read distribution.
     */
    private BufferedWriter writerTmpCovStats = null;

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
    boolean outputObs = false;

    /**
     * Flag whether to ouput the prediction (i.e., values after deconvolution).
     */
    boolean outputPred = false;

    /**
     * Flag whether to output balanced values, i.e., values that should be observed as cast back from deconvoluted values.
     */
    boolean outputBalanced = true;

    /**
     * Flag whether to output exon features.
     */
    boolean outputExon = false;

    /**
     * Flag whether to output unknown features, e.g., reproduced ones from the input annotation.
     */
    boolean outputUnknown = false;

    /**
     * Flag whether to output splice-junction features.
     */
    boolean outputSJunction = false;

    /**
     * Flag whether to output gene/loci features.
     */
    boolean outputGene = false;

    /**
     * Flag whether to output transcript features.
     */
    boolean outputTranscript = true;

    /**
     * Flag whether to output AStalavista event features.
     */
    boolean outputEvent = false;

    /**
     * Flag whether to output linear programs.
     */
    boolean outputLP = false;

    /**
     * Default compression for output files.
     */
    byte compressionOut = FileHelper.COMPRESSION_NONE;

    /**
     * Dimension of AStalavista events retrieved from the gene structures.
     */
    int eventDim = 2;

    /**
     * Cost model for deconvolution.
     */
    byte costModel = GraphLPsolver.COSTS_LINEAR;

    /**
     * Granularity of deconvolution cost function.
     */
    byte costSplit = 1;

    /**
     * Flag indicating whether read biases should be ignored (<code>true</code>),
     * i.e. uniform read distribution is assumed.
     */
    boolean uniform = false;

    /**
     * Lower and upper bound of how much of the original observation can be substracted respectively added.
     */
    float[] costBounds = new float[]{0.95f, Float.NaN};


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
    volatile int nrLoci = 0;

    /**
     * Number of loci expressed, i.e., with positive quantification values.
     */
    volatile int nrLociExp = 0;

    /**
     * Number of transcripts in the annotation.
     */
    volatile int nrTx = 0;

    /**
     * Number of transcripts expressed, i.e., with positive quantification values.
     */
    volatile int nrTxExp = 0;

    /**
     * Number of AStalavista events in the annotation.
     */
    volatile int nrEvents = 0;

    /**
     * Number of AStalavista events expressed, i.e., with positive quantification values.
     */
    volatile int nrEventsExp = 0;

    /**
     * Number of mappings that could be mapped to the annotation.
     */
    volatile int nrReadsMapped = 0;

    /**
     * Number of mappings that were found in annotated genes/loci.
     */
    volatile int nrReadsLoci = 0;

    /**
     * Number of not alternatively processed loci in the annotation.
     */
    volatile int nrSingleTranscriptLoci = 0;    // TODO counted in learn AND decompose redundantly

    /**
     * Number of not alternatively processed loci that have been used for learning biases.
     */
    volatile int nrSingleTranscriptLearn = 0;  // TODO what is the difference between this and nrSingleTranscriptLoci

    /**
     * Number of mappings within the boundaries of not alternatively processed loci.
     */
    volatile int nrReadsSingleLoci = 0;

    /**
     * Number of mappings that mapped to not alternatively processed loci.
     */
    volatile int nrReadsSingleLociMapped = 0;

    /**
     * Number of mapping pairs that mapped to not alternatively processed loci.
     */
    volatile int nrReadsSingleLociPairsMapped = 0;

    /**
     * Number of mappings in not alternatively processed loci that did not match
     * the expected gene structure provided in the annotation.
     */
    volatile int nrReadsSingleLociNoAnnotation = 0;

    /**
     * Number of loci that could not be solved.
     */
    volatile int nrUnsolved = 0;

    /**
     * Number of mappings that have an unexpected length.
     *
     * @deprecated should no longer be used
     */
    volatile int nrReadsWrongLength = 0; // TODO (how) can this still occur?

    /**
     * Number of mappings that do not match the expected strand.
     */
    volatile int nrMappingsWrongStrand = 0;

    /**
     * Number of mapping pairs that do not match the annotation.
     */
    volatile int nrPairsNoTxEvidence = 0;

    /**
     * Number of mapping pairs with the wrong orientation of both mates.
     */
    volatile int nrPairsWrongOrientation = 0;

    /**
     * The minimum and the maximum insert size found.
     */
    int[] insertMinMax = null;   // TODO check if correctly used

    /**
     * Wrapper to read mappings from a BED file format.
     *
     * @deprecated marked for removal
     */
    private BEDwrapper bedWrapper;  // TODO pull up to MappingWrapper

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
    private void printStats() {
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
                settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).toString());
        // TODO
        //p.println("\t"+CLI_LONG_VERBOSE+"\t"+Constants.VERBOSE_KEYWORDS[Constants.verboseLevel]);
        //if (copyLocal)
        //	p.println("\t"+FluxCapacitorParameters.PAR_COPY_INPUT);
        Log.info(settings.SORT_IN_RAM.getName(),
                Boolean.toString(settings.get(FluxCapacitorSettings.SORT_IN_RAM)));

        // OUTPUT
        Log.info(FluxCapacitorSettings.TMP_DIR.getName(),
                settings.get(FluxCapacitorSettings.TMP_DIR).getAbsolutePath());
        // TODO
        //if (Constants.globalPfx!= null)
        //p.println("\t"+FluxCapacitorConstants.CLI_LONG_TPX+"\t"+ Constants.globalPfx);
        sb = new StringBuilder();
        if (settings.get(FluxCapacitorSettings.STDOUT_FILE) == null)
            sb.append("stdout");
        else {
            sb.append(settings.get(FluxCapacitorSettings.STDOUT_FILE).getAbsolutePath());
            if (compressionOut != FileHelper.COMPRESSION_NONE)
                sb.append("\t" + FluxCapacitorConstants.CLI_LONG_COMPRESSION + "\t" + FileHelper.COMPRESSION_KEYWORDS[compressionOut]);
        }
        Log.info(settings.STDOUT_FILE.getName(), sb.toString());
/*			p.print("\tfeatures:\t");
			if (outputExon)
				p.print("Exons ");
			if (outputSJunction)
				p.print("splice-Junctions ");
			if (outputTranscript)
				p.print("Transcripts ");
			if (outputGene)
				p.print("Genes ");
			if (outputEvent)
				p.print("eVents ");
			p.println();
*/
/*			p.print("\tbases:\t");
			if (outputObs)
				p.print("Observed ");
			if (outputPred)
				p.print("preDicted ");
			if (outputBalanced)
				p.print("Balanced ");
			p.println();
*/
/*			p.print("\tscopes:\t");
			if (outputAll)
				p.print("All ");
			if (outputSplit)
				p.print("Split ");
			if (outputUnique)
				p.print("Unique ");
			p.println();
*/
/*			p.print("\tmeasure:\t");
			if (outputFreq)
				p.print("Freq ");
			if (outputSplit)
				p.print("Split ");
			if (outputUnique)
				p.print("Unique ");
			p.println();
*/
/*			if (outputMapped|| outputNotmapped|| outputProfiles|| outputLP) {
				p.print("\tsave:\t");
				if (outputMapped)
					p.print("Mapped-alignments ");
				if (outputNotmapped)
					p.print("Notmapped-alignments ");
				if (outputProfiles)
					p.print("Profiles ");
				if (outputLP)
					p.print("Linear-programs ");
				p.println();
			}
*/
        // ALGORITHM
        //p.println("\t"+ CLI_LONG_THREAD+" "+ maxThreads);
//		sb= new StringBuilder("Read Distribution\t");
//		if (uniform)
//			sb.append("uniform");
//		else if (fileProfile!= null&& fileProfile.exists())
//			sb.append("from profiles in "+ fileProfile.getAbsolutePath());
//		else {
//			sb.append("profiling is carried out");
//			if (fileProfile!= null)
//				sb.append(" and stored in "+ fileProfile.getAbsolutePath());
//		}
//		Log.info(sb.toString());

        if (stranded)
            Log.info("\tstrand information considered.");
        if (pairedEnd)
            Log.info("\tmate pairing information considered");


/*			p.print("\t"+ CLI_LONG_COST_MODEL+" "+GraphLPsolver.COSTS_NAMES[costModel]);
			if (!Double.isNaN(costModelPar))
				p.print(" "+Double.toString(costModelPar));
			p.println();
			p.println("\t"+ CLI_LONG_COST_SPLIT+ " "+costSplit);
			p.print("\t"+ CLI_LONG_COST_BOUNDS+ " ");
			if (costBounds== null)
				p.println(" none.");
			else
				p.println(costBounds[0]+","+costBounds[1]);
			//p.println("\tread length " + readLen);
			p.println("\t"+ CLI_LONG_STRAND+" "+ strandSpecific);
*/
        if (settings.get(FluxCapacitorSettings.INSERT_FILE) != null) {
            Log.info("\twriting insert sizes to " + settings.get(FluxCapacitorSettings.INSERT_FILE).getAbsolutePath());
        }
//		if (pairedEnd)
//			p.println("\t"+CLI_LONG_PAIR+"\t"+insertMinMax[0]+","+insertMinMax[1]);
        //System.err.println("\t"+CLI_LONG_NOISE+"\t"+Float.toString(1- GraphLPsolver.min_read_rest_frac));

    }

    /**
     * Finishes all pending file I/O operations and closes handles.
     */
    void fileFinish() {

        // TODO close input should occur by reader or interface method
        bedWrapper.close();
        gtfReader.close();


        if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
            if (FileHelper.move(
                    fileTmpCovStats,
                    settings.get(FluxCapacitorSettings.COVERAGE_FILE))) {
                fileTmpCovStats.delete();
                Log.info("Coverage statistics in " + settings.get(FluxCapacitorSettings.COVERAGE_FILE).getAbsolutePath());
            } else
                Log.warn("Failed to move coverage statistics to " +
                        settings.get(FluxCapacitorSettings.COVERAGE_FILE).getAbsolutePath() + "\n"
                        + "\tinformation in " + fileTmpCovStats.getAbsolutePath());
        }

        // TODO close files for non-/mapped reads, insert sizes, LPs, profiles

        // profiles
        if (settings.get(FluxCapacitorSettings.PROFILE_FILE) != null)
            writeProfiles();

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
    public FluxCapacitorStats call() throws Exception {

        // TODO not here
        if (loadLibraries() < 0)
            System.exit(-1);

        // load parameters
        if (file != null && !file.exists()) {
            throw new RuntimeException("I have no parameter file and I want to scream!");
        }

        if (file != null) {
            try {
                settings = FluxCapacitorSettings.createSettings(file);
            } catch (Exception e) {
                throw new RuntimeException("Unable to load settings from " + file + "\n\n " + e.getMessage(), e);
            }
        } else {
            // create default settings
            settings = new FluxCapacitorSettings();
            FluxCapacitorSettings.relativePathParser.setParentDir(new File(""));
        }

        // add command line parameter
        if (commandLineArgs != null) {
            if (commandLineArgs.userSpecified("annotation")) {
                settings.set(FluxCapacitorSettings.ANNOTATION_FILE, commandLineArgs.getFile("annotation"));
            }
            if (commandLineArgs.userSpecified("input")) {
                settings.set(FluxCapacitorSettings.MAPPING_FILE, commandLineArgs.getFile("input"));
            }
            if (commandLineArgs.userSpecified("output")) {
                settings.set(FluxCapacitorSettings.STDOUT_FILE, commandLineArgs.getFile("output").getAbsoluteFile());
            }
            if (commandLineArgs.userSpecified("annotation-mapping")) {
                try {
                    settings.set(FluxCapacitorSettings.ANNOTATION_MAPPING, AnnotationMapping.valueOf(commandLineArgs.getString("annotation-mapping")));
                } catch (Exception e) {
                    throw new RuntimeException("Invalid Annotation Mapping : " + commandLineArgs.getString("annotation-mapping"));
                }
            }
            if (commandLineArgs.userSpecified("read-descriptor")) {
                settings.set(FluxCapacitorSettings.READ_DESCRIPTOR, new UniversalReadDescriptor());
                settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).init(commandLineArgs.getString("read-descriptor"));
            }
            if (commandLineArgs.userSpecified("sort-in-ram")) {
                settings.set(FluxCapacitorSettings.SORT_IN_RAM, true);
            }
        }

        // validate the settings
        settings.validate();


        FileHelper.tempDirectory = settings.get(FluxCapacitorSettings.TMP_DIR);

        // prepare output files
        if (settings.get(FluxCapacitorSettings.STDOUT_FILE) != null) {
            File f = settings.get(FluxCapacitorSettings.STDOUT_FILE);
            if (f.exists() && !CommandLine.confirm(
                    "[CAUTION] I overwrite the output file " +
                            settings.get(FluxCapacitorSettings.STDOUT_FILE).getName() +
                            ", please confirm:\n\t(Yes,No,Don't know)")) {
                exit(-1);
            }

            try {
                Log.outputStream = new PrintStream(new FileOutputStream(f));
            } catch (FileNotFoundException e) {
                Log.warn("Cannot write log file to " + f.getAbsolutePath());    // let it on stderr?!
            }
        }

        // prepare input files
        AbstractFileIOWrapper wrapperAnnotation;
        AbstractFileIOWrapper wrapperMappings;
        if (cheatDisableFCheck) {
            Log.warn("Development run, file check disabled !!!");
            wrapperAnnotation = getWrapper(settings.get(FluxCapacitorSettings.ANNOTATION_FILE));
            wrapperMappings = getWrapper(settings.get(FluxCapacitorSettings.MAPPING_FILE));
        } else {

            Log.progressStart("Scanning annotation file");
            wrapperAnnotation =
                    fileInit(settings.get(FluxCapacitorSettings.ANNOTATION_FILE));
            fileStats((AnnotationWrapper) wrapperAnnotation);
            Log.progressFinish("OK", true);

            Log.progressStart("Scanning mapping file");
            wrapperMappings =
                    fileInit(settings.get(FluxCapacitorSettings.MAPPING_FILE));
            fileStats((MappingWrapper) wrapperMappings);
            Log.progressFinish("OK", true);
            Log.info("Annotation and mapping input checked");
        }


        // TODO parameters
        pairedEnd = settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).equals(AnnotationMapping.PAIRED)
                || settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).equals(AnnotationMapping.COMBINED);
        stranded = settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).equals(AnnotationMapping.STRANDED)
                || settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).equals(AnnotationMapping.COMBINED);


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
        //Get from settings the tasks to be executed in the current run
        if (!settings.get(FluxCapacitorSettings.NO_DECOMPOSE)) {
            currentTasks.add(Task.DECOMPOSE);
        }

        //print current run stats
        printStats();

        // run
        long t0 = System.currentTimeMillis();

        FluxCapacitorStats stats = new FluxCapacitorStats();
        if (currentTasks.contains(Task.DECOMPOSE)) {
            profile = getProfile(stats);
            if (profile == null) {
                exit(-1);
            }
        }

        explore(FluxCapacitorConstants.MODE_RECONSTRUCT, stats);

        // BARNA-103 : write stats to file
        File statsFile = settings.get(FluxCapacitorSettings.STATS_FILE);
        if (statsFile != null) {
            FluxCapacitorStats statsToWrite = stats;
            BufferedWriter writer = null;
            BufferedReader reader = null;
            Boolean append = settings.get(FluxCapacitorSettings.STATS_FILE_APPEND);

            File lockFile = new File(statsFile.getAbsolutePath() + ".lock");
//            if(!lockFile.exists()) lockFile.createNewFile();
            FileChannel channel = new RandomAccessFile(lockFile, "rw").getChannel();
            FileLock lock = channel.lock();

            try {
                Gson gson = new GsonBuilder().setPrettyPrinting().create();
                if (statsFile.exists() && append) {
                    // read stats file and append
                    reader = new BufferedReader(new FileReader(statsFile));
                    FluxCapacitorStats existingStats = gson.fromJson(reader, FluxCapacitorStats.class);
                    reader.close();
                    existingStats.add(stats);
                    statsToWrite = existingStats;
                }
                Log.info((append ? "Appending stats to " : "Writing stats to ") + statsFile.getAbsolutePath());
                writer = new BufferedWriter(new FileWriter(statsFile));
                gson.toJson(statsToWrite, writer);
                writer.close();
            } catch (Exception e) {
                Log.error("Unable to " + (append ? "append stats to " : "write stats to ") + statsFile.getAbsolutePath() + " : " + e.getMessage(), e);
            } finally {
                if (reader != null) reader.close();
                if (writer != null) writer.close();
                // release the lock
                try {
                    lock.release();
                } catch (IOException e) {
                    Log.error("Unable to release lock");
                }
                channel.close();
            }
        }

        fileFinish();

        Log.info("\n[TICTAC] I finished flux in "
                + ((System.currentTimeMillis() - t0) / 1000) + " sec.\nCheers!");

        //System.err.println("over "+ GraphLPsolver.nrOverPredicted+", under "+GraphLPsolver.nrUnderPredicted);

        return stats;
    }


    /**
     * Reads bias profiles from the provided source file.
     *
     * @return an instance describing the bias distribution
     */
    private Profile readProfiles() {

        try {
            profile = new Profile(this);

            ZipFile zf = new ZipFile(fileProfile);
            Enumeration entries = zf.entries();
            String line;
            Vector<Integer> v = new Vector<Integer>();
            Vector<UniversalMatrix> w = new Vector<UniversalMatrix>();
            System.err.println("[LOAD] getting profiles");
            while (entries.hasMoreElements()) {
                ZipEntry ze = (ZipEntry) entries.nextElement();
                BufferedReader buffy = new BufferedReader(
                        new InputStreamReader(zf.getInputStream(ze)));
                int lcount = 0;
                while ((line = buffy.readLine()) != null)
                    ++lcount;
                buffy.close();
                v.add(lcount);
                UniversalMatrix m = new UniversalMatrix(lcount);
                buffy = new BufferedReader(
                        new InputStreamReader(zf.getInputStream(ze)));
                lcount = 0;
                while ((line = buffy.readLine()) != null) {
                    String[] ss = line.split("\t");
                    assert (ss.length == 2);
                    m.sense[lcount] = Integer.parseInt(ss[0]);
                    m.sums += m.sense[lcount];
                    m.asense[lcount] = Integer.parseInt(ss[1]);
                    m.suma += m.asense[lcount];
                    ++lcount;
                }
                buffy.close();
                assert (lcount == m.sense.length);
                w.add(m);
            }
            zf.close();

            int[] len = new int[v.size()];
            for (int i = 0; i < len.length; i++)
                len[i] = v.elementAt(i);
            Arrays.sort(len);
            profile.masters = new UniversalMatrix[w.size()];
            for (int i = 0; i < len.length; i++) {
                for (int j = 0; j < len.length; j++) {
                    if (len[i] == v.elementAt(j)) {
                        profile.masters[i] = w.elementAt(j);
                        // check
                        for (int n = 0; n < profile.masters[i].getLength(); n++) {
                            if (profile.masters[i].asense[n] == 0 || profile.masters[i].sense[n] == 0) {
                                if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP)
                                    System.err.println("\tprofile with 0-count positions");
                                return null;
                            }
                        }
                    }
                }
            }
            System.err.println("\tfound " + profile.masters.length + " profiles.");

            return profile;

        } catch (Exception e) {
            e.printStackTrace();
        }

        return null;
    }


    /**
     * Performs the RPKM normalization, given a number of reads mapping to a transcript of a certain length.
     *
     * @param reads the number of mappings to the transcript
     * @param len   the length of the transcript
     * @return the RPKM value
     */
    public float calcRPKM(float reads, int len) {
        float rpkm = (float) ((reads / (double) len) * (1000000000l / (double) (nrBEDreads < 0 ? 1 : nrBEDreads)));
        return rpkm;
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
    protected File createTempFile(File location, String name, String extension, boolean deleteOnExit) {

        // get location
        if (location == null)
            location = settings.get(FluxCapacitorSettings.TMP_DIR);
        else {
            if (!location.isDirectory())
                location = location.getParentFile();
            if (!location.canWrite())
                location = settings.get(FluxCapacitorSettings.TMP_DIR);
        }

        // get name
        if (name == null)
            name = getClass().getSimpleName();
        else
            name = getClass().getSimpleName() + "_" + name;

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
    protected File createFile(File f, boolean deleteOnExit) {
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
     * Writes bias profiles to disk.
     */
    private void writeProfiles() {
        try {
            final String MSG_WRITING_PROFILES = "writing profiles";

            Log.progressStart(MSG_WRITING_PROFILES);

            BufferedWriter buffy = new BufferedWriter(new FileWriter(settings.get(FluxCapacitorSettings.PROFILE_FILE)));

            UniversalMatrix[] mm = profile.getMasters();
            for (int i = 0; i < mm.length; i++) {
                String lenString = Integer.toString(mm[i].getLength());
                buffy.write(lenString);
                buffy.write("\n");
                buffy.write(mm[i].toStringBuilder(20).toString());
            }
            buffy.close();
            Log.progressFinish(StringUtils.OK, true);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Writes coverage statistics of a transcript to disk.
     *
     * @param geneID       locus identifier
     * @param transcriptID transcript identifier
     * @param cds          flag to indicate whether transcript has an annotated ORF
     * @param length       (processed) length of the transcript
     * @param nrReads      number of Mappings
     * @param fracCov      the fraction covered
     * @param chiSquare    chi-square value for the coverage profile
     * @param cv           coefficient of variation for the coverage profile
     */
    private void writeCoverageStats(String geneID, String transcriptID,
                                    boolean cds, int length, int nrReads,
                                    float fracCov, long chiSquare, double cv) {

        try {
            if (fileTmpCovStats == null)
                fileTmpCovStats = FileHelper.createTempFile("tmpCovStats", ".pro");
            if (writerTmpCovStats == null)
                writerTmpCovStats = new BufferedWriter(new FileWriter(fileTmpCovStats));
            writerTmpCovStats.write(
                    geneID + "\t" + transcriptID + "\t" + (cds ? "CDS" : "NC") + "\t"
                            + Integer.toString(length) + "\t" + Integer.toString(nrReads) + "\t"
                            + Float.toString(fracCov) + "\t" + Long.toString(chiSquare) + "\t"
                            + Float.toString((float) cv) + "\n"
            );

        } catch (Exception e) {
            throw new RuntimeException(e);
        }
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
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        parameters.add(JSAPParameters.flaggedParameter("parameter", 'p').type(File.class).help("specify parameter file (PAR file)").valueName("file").get());
        parameters.add(JSAPParameters.flaggedParameter("annotation", 'a').type(File.class).help("Path to the annotation file").valueName("gtf").get());
        parameters.add(JSAPParameters.flaggedParameter("input", 'i').type(File.class).help("Path to the mapping file").valueName("bed").get());
        parameters.add(JSAPParameters.flaggedParameter("output", 'o').type(File.class).help("Path to the output file").valueName("gtf").get());
        parameters.add(JSAPParameters.flaggedParameter("annotation-mapping", 'm').type(String.class).help("Annotation Mapping (default PAIRED)").valueName("mapping").defaultValue("PAIRED").get());
        parameters.add(JSAPParameters.flaggedParameter("read-descriptor", 'd').type(String.class).help("Read Descriptor (default PAIRED)").valueName("descriptor").defaultValue("PAIRED").get());
        parameters.add(JSAPParameters.switchParameter("sort-in-ram", 'r').help("Sort in RAM").get());

        parameters.add(JSAPParameters.switchParameter("printParameters").help("Print default parameters").get());
        return parameters;
    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        commandLineArgs = args;
        setPrintParameters(args.userSpecified("printParameters"));
        setFile(args.getFile("parameter"));

        if (isPrintParameters()) {
            FluxCapacitorSettings settings = new FluxCapacitorSettings();
            settings.write(System.out);
            return false;
        }

        if (getFile() != null && !getFile().canRead()) {
            Log.error("");
            Log.error("Parameter file " + getFile().getAbsolutePath() + " does not exist or I can not read it!");
            Log.error("\n");
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
     * @param beds      the mappings in the locus
     * @param tasks     the tasks to be performed
     */
    private void solve(Gene gene, BufferedIterator beds, EnumSet tasks) {

        // create LP and solve
        LocusSolver lsolver = new LocusSolver(gene, beds, tasks);
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
    }

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
     * Merges and smoothens bias profiles.
     *
     * @param stats statistics object for storing attributes
     * @return a model for the bias profile
     */
    Profile getProfile(FluxCapacitorStats stats) {
        if (uniform) {
            profile = new Profile(this);
            profile.fill();
        } else {
            if (fileProfile != null && fileProfile.exists()) {
                profile = readProfiles();
                if (profile != null) {
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
                }
            }
            if (profile == null) {
                profile = new Profile(this);
                try {
                    explore(FluxCapacitorConstants.MODE_LEARN, stats);
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

            if (pairedEnd) {
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

            if (pairedEnd) {
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
     * Returns the breakpoint indicated by a mapping within a transcript.
     *
     * @param tx  transcript to which a read maps
     * @param bed genomic mappping
     * @return transcript coordinate of the breakpoint indicated by the mapping
     */
    private int getBpoint(Transcript tx, Mapping bed) {

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
    }


    /**
     * If instantiated, the method returns the GTF file wrapper,
     * otherwise <code>null</code>.
     *
     * @return a wrapper instance for GTF files, or <code>null</code>
     */
    private AbstractFileIOWrapper getWrapperGTF() {
        if (gtfReader == null) {

            return getWrapperGTF(settings.get(FluxCapacitorSettings.ANNOTATION_FILE));
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
    private AbstractFileIOWrapper getWrapperGTF(File inputFile) {

        gtfReader = new GTFwrapper(inputFile.getAbsolutePath());
        gtfReader.setNoIDs(null);
        gtfReader.setReadGene(true);
        gtfReader.setReadFeatures(new String[]{"exon", "CDS"});
        gtfReader.setReadAheadTranscripts(1);    // only one locus a time
//		gtfReader.setReadAheadTranscripts(-1);
//		gtfReader.setReadAll(true);
        gtfReader.setGeneWise(true);
        gtfReader.setPrintStatistics(false);
        gtfReader.setReuse(true);
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
    private AbstractFileIOWrapper getWrapperBED(File inputFile) {
        bedWrapper = new BEDwrapper(inputFile.getAbsolutePath());
        return bedWrapper;  // TODO pull up to MappingWrapper
    }

    /**
     * Retrieves all mappings in a certain region from a BED input file.
     *
     * @param gene the locus for which reads are to be read
     * @param from start coordinate on chromosome
     * @param to   end coordinate on chromosome
     * @return an iterator instance that enumerates all mappings in the specified region
     */
    private BufferedIterator readBedFile(Gene gene, int from, int to) {
        return readBedFile(gene, from, to, 0, 1);
    }

    /**
     * Out-of-memory-proof method that retrieves all mappings in a certain region from a BED input file.
     * The strategy is try-and-see, first it is attempted to try to load all requested reads into memory (RAM);
     * if latter attempt fails, the iterator is initialized on disk. Method retries if disk/filesystem blocks.
     * in the latter case.
     *
     * @param gene          the locus for which reads are to be read
     * @param from          start coordinate on chromosome
     * @param to            end coordinate on chromosome
     * @param retryCount    number of retries that are attempted in the case of disk/filesystem temporarily unreachable
     * @param timeInSeconds time between retries
     * @return an iterator instance that enumerates all mappings in the specified region
     */
    private BufferedIterator readBedFile(Gene gene, int from, int to, int retryCount, long timeInSeconds) {
        if (settings.get(FluxCapacitorSettings.SORT_IN_RAM)) {
            try {
                return readBedFileRAM(gene, from, to);
            } catch (OutOfMemoryError memoryError) {
                System.gc();
                Thread.yield();
                Log.warn("Not enough memory to sort BED entries in RAM. Switching to disk sorting. This run is NOT failed!\n " +
                        "You can increase the amount of memory used " +
                        "by the capacitor using the FLUX_MEM environment variable. For example: export FLUX_MEM=\"6G\"; flux-capacitor ... to use" +
                        "6 GB of memory.");
                return readBedFileDisk(gene, from, to, retryCount, timeInSeconds);
            }
        } else {
            return readBedFileDisk(gene, from, to, retryCount, timeInSeconds);
        }
    }


    /**
     * Loads all mappings in the respective region into RAM.
     *
     * @param gene the locus for which reads are to be read
     * @param from start coordinate on chromosome
     * @param to   end coordinate on chromosome
     * @return an iterator instance that enumerates elements of an array stored in RAM
     */
    private BufferedIterator readBedFileRAM(Gene gene, int from, int to) {

        if (from > to || from < 0 || to < 0)
            throw new RuntimeException("BED reading range error: " + from + " -> " + to);
        // init iterator
        BufferedIterator iter = null;
        // memory
        MappingWrapperState state = bedWrapper.read(gene.getChromosome(), from, to);
        if (state.result == null)
            return null;
        BEDMapping[] beds = (BEDMapping[]) state.result;//TODO move to Mapping
        Arrays.sort(beds, getDescriptorComparator());
        iter = new BufferedIteratorRAM(beds);

        return iter;

    }

    /**
     * Writes all mappings in the respective region to disk, retries if disk/filesystem blocks.
     *
     * @param gene the locus for which reads are to be read
     * @param from start coordinate on chromosome
     * @param to   end coordinate on chromosome
     * @return an iterator instance that enumerates elements of an array stored in RAM
     */
    private BufferedIterator readBedFileDisk(Gene gene, int from, int to, int retryCount, long timeInSeconds) {

        if (from > to || from < 0 || to < 0)
            throw new RuntimeException("BED reading range error: " + from + " -> " + to);

        // init iterator
        BufferedIterator iter = null;

        try {
            // read, maintain main thread
            PipedInputStream pin = new PipedInputStream();
            PipedOutputStream pout = new PipedOutputStream(pin);
            Comparator<CharSequence> c = new BEDDescriptorComparator(settings.get(FluxCapacitorSettings.READ_DESCRIPTOR));
            File tmpFile = createTempFile(null, gene.getChromosome() + ":" + from + "-" + to + ".", "bed", true);
            BufferedIteratorDisk biter = new BufferedIteratorDisk(pin, tmpFile, c);
            biter.init();
            iter = biter;
            MappingWrapperState state = bedWrapper.read(pout, gene.getChromosome(), from, to);
            pout.flush();
            pout.close();
            if (state.count == 0)
                return null;
        } catch (IOException e) {
            /*
             * "Resource temporarily unavailable"
             * Catch this exception and try again after sleeping for a while
             */
            if (e.getMessage().contains("Resource temporarily unavailable")) {
                if (retryCount < 6) {
                    Log.warn("Filesystem reports : 'Resource temporarily unavailable', I am retrying (" + (retryCount + 1) + ")");
                    try {
                        Thread.sleep(1000 * (timeInSeconds));
                    } catch (InterruptedException e1) {
                    }
                    return readBedFileDisk(gene, from, to, retryCount + 1, timeInSeconds * 6);
                }
            }
            throw new RuntimeException(
                    "Could not get reads for locus " + gene.getChromosome() + ":" + from + "-" + to + ", retried " + retryCount + " times", e);
        }

        return iter;

    }


    /**
     * Finalizes an iteration over all mappings: closes readers, outputs stats, etc.
     *
     * @param mode  task carried out during the iteration, profiling or deconvolution
     * @param secs  time measured for the iteration, in [sec]
     * @param stats object capturing statistics of the run
     */
    void exploreFinish(byte mode, long secs, FluxCapacitorStats stats) {

        try {


            if (mode == FluxCapacitorConstants.MODE_LEARN) {

                // TODO calculate insert size distribution parameters

                // close coverage writer
                if (settings.get(FluxCapacitorSettings.COVERAGE_STATS))
                    writerTmpCovStats.close();
                writerTmpCovStats = null;


                if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP) {

                    //System.err.println(" OK.");
                    System.err.println("\tfirst round finished .. took " + secs + " sec.\n\n\t"
                            + nrSingleTranscriptLoci + " single transcript loci\n\t"
                            + bedWrapper.getNrLines() + " mappings in file\n\t"
                            + nrReadsSingleLoci + " mappings fall in single transcript loci\n\t"    // these loci(+/-"+tolerance+"nt)\n\t"
                            // counter un-reliable, /2 read is skipped in paired-end mode
                            + ((strand == FluxCapacitorConstants.STRAND_SPECIFIC) ? nrMappingsWrongStrand + " mappings map to annotation in antisense direction,\n\t" : "")
                            //+ (pairedEnd?(nrReadsSingleLociPotentialPairs+ " mappings form potential pairs,\n\t"):"")
                            + (pairedEnd ? (nrReadsSingleLociPairsMapped) + " mappings in annotation-mapped pairs\n\t"
                            : nrReadsSingleLociMapped + " mappings map to annotation\n\t")
                            //+ nrReadsSingleLociNoAnnotation+ " mappings do NOT match annotation,\n\t"
                            //+ (uniform?"":func.profiles.size()+" profiles collected\n\t")
                            + readLenMin + "," + readLenMax + " min/max read length\n\t"
                            + (pairedEnd && insertMinMax != null ? insertMinMax[0] + "," + insertMinMax[1] + " min/max insert size\n\t" : ""));
                    //nrUniqueReads= getBedReader().getNrUniqueLinesRead();
                    //System.err.println("\ttotal lines in file "+nrUniqueReads);
                    System.err.println();
                }

                // output stats
                if (stats != null) {
                    stats.setLociSingle(nrSingleTranscriptLoci);
                    stats.setMappingsTotal(bedWrapper.getNrLines());
                    stats.setMappingsSingle(nrReadsSingleLoci);
                    if (strand == FluxCapacitorConstants.STRAND_SPECIFIC) {
                        stats.setMappingsNotSens(nrMappingsWrongStrand);
                    }
                    if (pairedEnd) {
                        stats.setMappingsSinglePairs(nrReadsSingleLociPairsMapped);
                        stats.setMappingsSinglePairsMapped(nrReadsSingleLociMapped);
                    }
                }

            } else if (mode == FluxCapacitorConstants.MODE_RECONSTRUCT) {
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
                            + bedWrapper.getNrLines() + " mappings read from file\n\t"
                            // no info, reads in redundantly many reads
                            //+ nrReadsLoci+" mappings in annotated loci regions\n\t"
                            + nrReadsMapped + " mappings" + (pairedEnd ? " in pairs" : "s") + " map to annotation\n"
                            + (pairedEnd ?
                            "\t" + nrPairsNoTxEvidence + " mappings without tx evidence\n"
                                    + "\t" + nrPairsWrongOrientation + " mappings with wrong orientation\n"
                            //+ "\t"+ nrMappingsForced+ " single mappings forced\n"
                            : "")
                            + ((strand == FluxCapacitorConstants.STRAND_SPECIFIC) ? nrMappingsWrongStrand + " mappings map to annotation in antisense direction\n\t" : "")
                            //+ nrMultiMaps+" mapped multiply.\n\n\t"
                            + (outputGene ? "\n\t" + nrLoci + " loci, " + nrLociExp + " detected" : "")
                            + (outputTranscript ? "\n\t" + nrTx + " transcripts, " + nrTxExp + " detected" : "")
                            + (outputEvent ? "\n\t" + nrEvents + " ASevents of dimension " + eventDim + ", " + nrEventsExp + " detected" : "")
                            + "\n"
                            //+ nrUnsolved+" unsolved systems."
                    );
                }

                // output stats
                if (stats != null) {
                    stats.setMappingsTotal(bedWrapper.getNrLines());
                    stats.setMappingsMapped(nrReadsMapped);
                    if (pairedEnd) {
                        stats.setMappingsPairsNa(nrPairsNoTxEvidence);
                        stats.setMappingsPairsWo(nrPairsWrongOrientation);
                    }
                    if (strand == FluxCapacitorConstants.STRAND_SPECIFIC) {
                        stats.setMappingsNotSens(nrMappingsWrongStrand);
                    }
                    if (outputGene) {
                        stats.setLociExp(nrLociExp);
                    }
                    if (outputTranscript) {
                        stats.setTxExp(nrTxExp);
                    }
                    if (outputEvent) {

                        stats.setEventsExp(nrEvents);
                    }
                }
            }

        } catch (Exception e) {
            throw new RuntimeException(e);
        }


    }


    /**
     * A complete iteration cycle over all mappings, either for profiling or for deconvolution.
     *
     * @param mode  task carried out during the iteration, profiling or deconvolution
     * @param stats object capturing statistics of the run
     */
    public boolean explore(byte mode, FluxCapacitorStats stats) {

        nrSingleTranscriptLoci = 0;
        nrReadsLoci = 0;
        nrReadsMapped = 0;
        nrReadsWrongLength = 0;
        nrMappingsWrongStrand = 0;

        if (mode == FluxCapacitorConstants.MODE_LEARN) {
            nrReadsSingleLoci = 0;
            nrReadsSingleLociMapped = 0;
        }

        //System.out.println(System.getProperty("java.library.path"));
        long t0 = System.currentTimeMillis();
        try {

            Transcript.setEdgeConfidenceLevel(Transcript.ID_SRC_MOST_INCONFIDENT);
            //this.gtfReader= null;
            //GFFReader gtfReader= getGTFreader();
            gtfReader.reset();
            bedWrapper.reset();

            if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP) {
                if (mode == FluxCapacitorConstants.MODE_LEARN) {
                    Log.info("","");
                    Log.info("LEARN", "Scanning the input and getting the attributes.");
                }
                else if (mode == FluxCapacitorConstants.MODE_RECONSTRUCT) {
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
            }
            final String profiling = "profiling ", decomposing = "decomposing ";

            if (mode == FluxCapacitorConstants.MODE_LEARN)
                Log.progressStart(profiling);
            else if (mode == FluxCapacitorConstants.MODE_RECONSTRUCT)
                Log.progressStart(decomposing);


            // BARNA-112 disable keeping original lines
//				if (mode== FluxCapacitorConstants.MODE_LEARN)
//					gtfReader.setKeepOriginalLines(false);
//				else if (mode== FluxCapacitorConstants.MODE_RECONSTRUCT)
//					gtfReader.setKeepOriginalLines(true);

            gtfReader.read();
            Gene[] gene = null, geneNext = gtfReader.getGenes();

            long tlast = System.currentTimeMillis();
            boolean output = false;

            String lastChr = null;
            byte lastStr = 0;
            int lastEnd = -1;
            int tol = this.tolerance; // 1000

            if (geneNext != null) {
                lastChr = geneNext[0].getChromosome();
                lastStr = geneNext[0].getStrand();
            }

            Thread readerThread = null;
            int readObjects = 0;
            while (lastChr != null) {    // MAIN LOOP


                if ((gene = geneNext) == null)
                    break;
                if (mode == FluxCapacitorConstants.MODE_RECONSTRUCT && gtfReader.getVLines() != null)
                    origLines = (Vector<String>) gtfReader.getVLines().clone();    // TODO make array, trim..

                // http://forums.sun.com/thread.jspa?threadID=5171135&tstart=1095
                if (readerThread == null)
                    readerThread = new GtfReaderThread();
                //readerThread.start();
                readerThread.run();
                geneNext = gtfReader.getGenes();

                for (int i = 0; (gene != null) && i < gene.length; i++) {


                    // flop strand
                    if (lastChr.equals(gene[i].getChromosome())) {
                        if (lastStr != gene[i].getStrand()) {
                            //System.err.println(lastChr+" "+lastStr+ " "+ readObjects+ " wrote "+ dbgCntWriteMap +" not "+ dbgCntWriteNonmap);
                            readObjects = 0;
                            // jump back
                            bedWrapper.reset(gene[i].getChromosome());
                            lastStr = gene[i].getStrand();
                            lastEnd = -1;
                        }
                    } else {                        // flop chr
                        //System.err.println(lastChr+" "+lastStr+ " "+ readObjects+ " wrote "+ dbgCntWriteMap +" not "+ dbgCntWriteNonmap);
                        readObjects = 0;
                        lastChr = gene[i].getChromosome();
                        lastStr = gene[i].getStrand();
                        lastEnd = -1;
                    }

                    if (gene[i].getTranscriptCount() == 1)
                        ++nrSingleTranscriptLoci;
                    else if (mode == FluxCapacitorConstants.MODE_LEARN)
                        continue;    // performance for not reading beds

                    BufferedIterator beds = null;

                    /*					File f= File.createTempFile("fluxpfx", ".bed");
                             FileOutputStream fos= new FileOutputStream(f);
                             handler.addStream(fos);
                             fileBED= f;
                             bedWrapper= null;
         */
                    // boundaries
                    int start = gene[i].getStart();
                    int end = gene[i].getEnd();
                    assert (geneNext == null || geneNext.length == 1);

                    if (gene[i].getStrand() < 0) {
                        start = -start;
                        end = -end;
                    }
                    tol = 0;
                    start = Math.max(1, start - tol);
                    end = end + tol;

                    beds = readBedFile(gene[i], start, end);

                    if (mode == FluxCapacitorConstants.MODE_LEARN && beds != null) {
                        solve(gene[i], beds, EnumSet.of(Task.LEARN));
                    } else if (mode == FluxCapacitorConstants.MODE_RECONSTRUCT) {
                        solve(gene[i], beds, currentTasks);
                    }

                    if (beds != null)
                        beds.clear();

                    if (output) {
                        System.out.println(gene[i].getChromosome() + " " +
                                gene[i].getStrand() +
                                " cluster " + gene[i].getLocusID());
                        // TODO beds.size() no longer available
                        //", "+beds.size()+" reads.");
                        if ((lastStr != gene[i].getStrand()
                                || !(lastChr.equals(gene[i].getChromosome())))) {
                            long t = System.currentTimeMillis();
                            if (lastStr != 0 && (!lastChr.equals("")))
                                System.out.println(lastChr + " " + lastStr +
                                        " " + ((t - tlast) / 1000) + " sec.");
                            tlast = t;
                            lastStr = gene[i].getStrand();
                            lastChr = gene[i].getChromosome();
                        }
                    }
                }
                //getWriter().flush();


            }    // end iterate GTF

            bedWrapper.finish();

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
            if (checkBEDscanMappings > 0 && checkBEDscanMappings != bedWrapper.getNrLines())
                System.err.println("[ERROR] consistency check failed in BED reader " + checkBEDscanMappings + "<>" + bedWrapper.getNrLines());
            //checkBEDscanMappings= getBedReader().getNrLines();

            // close readers, output
            exploreFinish(mode, ((System.currentTimeMillis() - t0) / 1000), stats);

        } catch (Exception e1) {
            Log.error("Error while iterating loci:", e1);
            throw new RuntimeException(e1);
        }

        return true;
    }

    /**
     * Checks whether an input file is to be uncompressed and/or sorted before reading.
     *
     * @param inputFile file which is to be read
     * @return a wrapper instance that accesses the file, after potential decompression and/or sorting
     */
    public AbstractFileIOWrapper fileInit(File inputFile) {

        // (1) unpack, if compressed
        byte cb = FileHelper.getCompression(inputFile);
        if (cb != FileHelper.COMPRESSION_NONE) {
            File f = new File(FileHelper.stripExtension(inputFile.getAbsolutePath()));
            if (f.exists()) {
                Log.println("Assuming file " + f.getName() + " is a decompressed version of " + inputFile.getName());
            } else {
                f = createTempFile(null, FileHelper.getFileNameWithoutExtension(f), FileHelper.getExtension(f), true);
                try {
                    FileHelper.inflate(inputFile, f, cb);
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            }
            inputFile = f;
        }

        // (2) sort, if needed
        AbstractFileIOWrapper wrapper = getWrapper(inputFile);
        if (!wrapper.isApplicable()) {
            File f = FileHelper.getSortedFile(inputFile);
            File lock = FileHelper.getLockFile(f);

            if (f.exists() && !lock.exists()) {

                Log.warn("Assuming file " + f.getName() + " is a sorted version of " + inputFile.getName());

            } else {    // we have to sort

                boolean lockCreated = false;
                if (settings.get(FluxCapacitorSettings.KEEP_SORTED_FILES)) {    // try to store in original

                    if (lock.exists()) {    // switch to sorting to temp
                        Log.warn("Seems that another process is just sorting file " + inputFile +
                                "\nremove lock file " + lock.getName() + " if dead leftover." +
                                "\nContinuing with sorting to temporary file " +
                                (f = createTempFile(f,     // access to non-Temp
                                        FileHelper.getFileNameWithoutExtension(f),
                                        FileHelper.getExtension(f),
                                        false)).getAbsolutePath());

                    } else if (!f.getParentFile().canWrite()) {    // sort to temp, but do not delete (parameter)
                        Log.warn("Cannot write sorted file to " + f.getAbsolutePath() +
                                "\nContinuing with sorting to temporary file " +
                                (f = createTempFile(f, // access to non-Temp
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
            wrapper = getWrapper(inputFile);
        }

        return wrapper;
    }

    /**
     * Obtains global statistics from the mapping file, e.g., number of total mappings etc.
     *
     * @param wrapper mapping file reader
     */
    private void fileStats(MappingWrapper wrapper) {

        if (settings.get(FluxCapacitorSettings.NR_READS_MAPPED) <= 0) {
            // (3) scan
            ((AbstractFileIOWrapper) wrapper).scanFile();
            if (((AbstractFileIOWrapper) wrapper).getNrInvalidLines() > 0)
                Log.warn("Skipped " + ((AbstractFileIOWrapper) wrapper).getNrInvalidLines() + " lines.");

            checkBEDscanMappings = wrapper.getCountMappings();
            nrBEDreads = wrapper.getCountReads();
            nrBEDmappings = wrapper.getCountMappings();
        } else {
            checkBEDscanMappings = -1;
            nrBEDreads = settings.get(FluxCapacitorSettings.NR_READS_MAPPED);
            nrBEDmappings = -1;
        }

        Log.info("\t" + nrBEDreads + " reads, "
                + (nrBEDmappings > 0 ? nrBEDmappings + " mappings: R-factor " + (wrapper.getCountMappings() / (float) wrapper.getCountReads()) : ""));
        if (nrBEDmappings > 0)
            Log.info("\t" + wrapper.getCountContinuousMappings() + " entire, " + wrapper.getCountSplitMappings()
                    + " split mappings (" + (wrapper.getCountSplitMappings() * 10f / wrapper.getCountMappings()) + "%)");

        // (4) check if read descriptor is applicable
        if (wrapper.isApplicable(settings.get(FluxCapacitorSettings.READ_DESCRIPTOR)))
            Log.info("\tRead descriptor seems OK");
        else {
            String msg = "Read Descriptor " + settings.get(FluxCapacitorSettings.READ_DESCRIPTOR)
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
    private AbstractFileIOWrapper getWrapper(File inputFile) {

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
                return getWrapperBED(inputFile);

        }

        return null;    // make compiler happy
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
 