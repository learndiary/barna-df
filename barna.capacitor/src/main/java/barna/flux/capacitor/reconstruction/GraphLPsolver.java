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

import barna.commons.log.Log;
import barna.commons.system.OSChecker;
import barna.flux.capacitor.graph.AnnotationMapper;
import barna.flux.capacitor.graph.MappingsInterface;
import barna.flux.capacitor.matrix.UniversalMatrix;
import barna.flux.capacitor.profile.MappingStats;
import barna.flux.capacitor.profile.Profile;
import barna.io.FileHelper;
import barna.model.SpliceSite;
import barna.model.Transcript;
import barna.model.commons.DoubleVector;
import barna.model.commons.IntVector;
import barna.model.constants.Constants;
import barna.model.splicegraph.AbstractEdge;
import barna.model.splicegraph.SimpleEdge;
import barna.model.splicegraph.SplicingGraph;
import barna.model.splicegraph.SuperEdge;
import lpsolve.LpSolve;
import lpsolve.LpSolveException;

import java.io.*;
import java.util.*;

import static lpsolve.LpSolve.LE;

/**
 * A class that takes a splicing graph with annotation mapped read counts
 * and transforms it into a system of linear equations that subsequently is
 * solved to yield a deconvolution of reads in common transcript areas.
 * @author Micha Sammeth (micha@sammeth.net)
 *
 */
public class GraphLPsolver {

    /**
     * Set of in-equation symbols.
     */
    static final String[] COND_SYMBOLS= new String[] {"", " <= ", " >= ", " = "};

    /**
     * Constant for linear cost function.
     */
    public final static byte COSTS_LINEAR= 0;

    /**
     * Constant for logarithmic cost function.
     */
    public final static byte COSTS_LOG= 1;

    /**
     * File extension for linear program files written to disk.
     */
    public final static String SFX_LPOUT= ".lp";

    /**
     * Flag for debug purposes.
     * @deprecated marked for removal
     */
    static boolean debug= true;    // TODO replace by JUnit tests

    /**
     * Splicing graph with mapped reads.
     */
    protected AnnotationMapper aMapper= null;

    /**
     * Linear Program solver.
     */
    protected LpSolve lpSolve= null;

    /**
     * Hash storing edges and transcripts that are mapped to integer constraint numbers.
     */
    Hashtable<Object,int[]> constraintHash= null;	// for synchronizing different accesses to same constraints

    /**
     * Variable to count up constraint numbers as they successively get added.
     */
    public int constraintCtr= 0;

    /**
     * Variable to count up restriciton numbers as they successively get added.
     */
    int restrNr= 0;

    /**
     * The (minimal) length of mapped reads.
     */
    int readLen= 0;

    /**
     * Statistics from annotation mapping.
     */
    MappingStats mappingStats= null;

    /**
     * Value of the objective function after solving the linear system (i.e., deconvolution).
     */
    double valObjFunc= -1d;

    /**
     * Sum of all deviations from observations in positive direction,
     * i.e., all reads added before deconvolution.
     * @deprecated not active
     */
    double valDeltaPos= -1d;

    /**
     * Sum of all deviations from observations in negative direction,
     * i.e., all reads removed before deconvolution.
     * @deprecated not active
     */
    double valDeltaNeg= -1d;

    /**
     * Array with the primal solution.
     */
    double[] result= null;

    double[] resultR= null;
    double[] resultC= null;

    /**
     * Hash that maps transcripts to their
     * predicted expression level after deconvolution.
     */
    HashMap<Object,Double> trptExprHash= null;

    /**
     * Profile of mapping distribution biases along transcript, from 5' to 3'.
     */
    Profile profile= null;

    /**
     * The number of mappings initially observed, for normalization purposes after deconvolution.
     */
    int nrMappingsObs= 0;

    /**
     * Normalization factor for deconvoluted expression levels.
     */
    double nFactor= Double.NaN;

    /**
     * Flag whether cost function is to be split into Watson and Crick mappings separately.
     */
    boolean costSplitWC= false;

    /**
     * Flag whether multimaps should be dis-/considered in the cost function.
     */
    boolean costUseMultimap= false;

    /**
     * Flag for paired-end mode.
     */
    boolean pairedEnd= false;

    /**
     * Model for cost function.
     */
    public byte costModel= COSTS_LOG;

    /**
     * Granularity of cost function.
     */
    public byte costSplit= 1;

    /**
     * Boundaries of cost function.
     */
    public float[] costBounds= null;

    /**
     * Flag whether LP programs should be documented in files, or not.
     */
    boolean writeFile= false;

    /**
     * Handle describing the directory where LP files are documented.
     */
    File fileLPdir= null;   // TODO move to settings

    /**
     *
     */
    File fileLPinput= null;
    IntVector costIdx;
    DoubleVector costVal;
    int[] insertMinMax;




    /**
     * Basic constructor, providing a graph with mappings, the minimum read length, and the total number of reads
     * observed in the locus.
     * @param aMapper splicing graph with annotation mapped reads
     * @param mapStats mapping statistics
     * @param realReads number of reads observed after annotation mapping
     */
    private GraphLPsolver(AnnotationMapper aMapper, MappingStats mapStats, int realReads) {
        this.aMapper= aMapper;
        this.mappingStats= mapStats;
        this.nrMappingsObs= realReads;
    }

    /**
     * Extended constructor, with read and insert size attributes.
     * @param aMapper splicing graph with annotation mapped reads
     * @param mapStats mapping statistics
     * @param realReads number of reads observed after annotation mapping
     * @param considerBothStrands flag for consideration of reads mapping in anti-/sense
     * @param pairedEnd flag for paired-end reads
     */
    public GraphLPsolver(AnnotationMapper aMapper, MappingStats mapStats, int realReads,
                         boolean considerBothStrands, boolean pairedEnd) {
        this(aMapper, mapStats, realReads);
        this.costSplitWC= considerBothStrands;
        this.pairedEnd= pairedEnd;
        // TODO this.insertMinMax= insertMinMax;
    }

    /**
     * Setter method for profile with mapping distribution biases.
     * @param profile profile of the biases
     */
    public void setProfile(Profile profile) {
        this.profile= profile;
    }

    /**
     * Retrieves or creates the linear program (LP) solver model.
     * @return the linear program (LP) solver model
     */
    LpSolve getLPsolve() {
        if (lpSolve == null) {
            if (writeFile)
                try {
                    getLPWriter().flush();
                    getLPWriter().close();
                    lpSolve= LpSolve.readLp(fileLPinput.getAbsolutePath(), LpSolve.IMPORTANT, null);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            else
                try {
//				Iterator iter= getConstraintHash().keySet().iterator();
//				int edgeCtr= 0, varCtr= 0;
//				while (iter.hasNext()) {
//					Object o= iter.next();
//					if (o instanceof Edge)
//						++edgeCtr;
//					else if (o instanceof Variation)
//						++varCtr;
//				}
//				
                    lpSolve = LpSolve.makeLp(0, constraintCtr);	// no 0-column here

                } catch (LpSolveException e) {
                    e.printStackTrace();
                }

        }

        return lpSolve;
    }

    /**
     * Outputs the constraint map for debugging purposes.
     * @param p stream to which the output is requested
     * @deprecated currently not in use
     */
    strictfp void printConstraintMap(PrintStream p) {

        // edges
        Object[] o= getConstraintHash().keySet().toArray();

        // sort a bit
        String[] ids= new String[o.length];
        HashMap<String,Object> map= new HashMap<String,Object>(o.length,1f);
        for (int i = 0; i < o.length; i++) {
            String s= o[i].toString();
            ids[i]= s;
            map.put(s,o[i]);
        }
        Arrays.sort(ids);
        for (String id : ids) {
            Object oo = map.get(id);
            int[] c = getConstraintHash().get(oo);
            p.print(id + "\tC" + c[0]);
            for (int j = 1; j < c.length; j++)
                p.print("\tC" + c[j]);
            if (oo instanceof SimpleEdge) {
                SimpleEdge e = (SimpleEdge) oo;
                int a = ((MappingsInterface) e).getMappings().getReadNr();
                //int b= e.getPossReadNr();
                p.print("\t" + Integer.toString(a));
                //p.print("\t"+Integer.toString(b));
                //p.print("\t"+Double.toString(((double) a)/((double) b)));
                Transcript[] t = aMapper.decodeTset(e.getTranscripts());
                for (Transcript aT : t)
                    p.print("\t" + aT.getTranscriptID());
            } else if (oo instanceof Transcript) {
                Transcript t = (Transcript) oo;
                p.print("\t" + t.getExonicLength());
            }
            p.print("\n");
        }
    }

    /**
     * Calculates and returns the normalization factor according to the read deviation
     * before/after deconvolution.
     * @return normalization factor to be applied to raw predicitions after deconvolutions
     */
    public double getNFactor() {
        //if (Double.isNaN(nFactor)|| true) {
        double fictReads= 0;
        for (Object o : getTrptExprHash().keySet()) {
            if (!(o instanceof String))
                continue;
            fictReads += getTrptExprHash().get(o);
        }

        if (fictReads< nrMappingsObs)
            ++nrUnderPredicted;
        else
            ++nrOverPredicted;
        // can happen, when reads are where none expected
//            if (fictReads== 0^ nrMappingsObs== 0)
//				System.currentTimeMillis();
        // avoid large scaling; was 0, avoid NaN; 0.5 too large
        if (fictReads> 0.000001) {
            nFactor= nrMappingsObs/ fictReads;
        } else
            nFactor= 1d;
        //}
        return nFactor;
    }


    /**
     * Creates a consecutive array where the specified indices
     * are set to the specified values.
     *
     * @param constIdx indices of constraints
     * @param constVal values of constraints
     * @return consecutive array implementing the specified values
     * at the specified indices
     */
    private double[] createArray(int[] constIdx, double[] constVal) {

        double[] a= new double[constraintCtr+ 1];
        for (int i = 0; i < a.length; i++)
            a[i]= 0d;
        for (int i = 0; i < constIdx.length; i++)
            a[constIdx[i]]= constVal[i];
        return a;
    }



    /**
     * Assigns the integer constraint numbers associated with the
     * restriction that is implied by the specified graph edge.
     *
     * @return array of constraint indices assigned to the specified
     * edge
     */
    int[] getConstraintIDs() {

        int size= 2; // plus minus
        if (costUseMultimap)	//
            ++size;		// only substract at 0-cost

        // init
        int[] a= new int[size];
        for (int i = 0; i < a.length; i++)
            a[i]= ++constraintCtr;

        return a;

    }


    private int[] allTrptIdx= null;

    int[] getAllTrptIdx() {
        if (allTrptIdx == null) {
            allTrptIdx = new int[aMapper.trpts.length];
            for (int i = 0; i < aMapper.trpts.length; i++)
                allTrptIdx[i]= getConstraintHash().get(aMapper.trpts[i])[0];
            Arrays.sort(allTrptIdx);
        }

        return allTrptIdx;
    }

    /**
     * Array to be reused for index values.
     */
    private int[] reuseIdx;

    /**
     * Array to be reused for values of the indices.
     */
    private double[] reuseVal;

    /**
     * Adds the constraints for a certain edge of the graph
     * to the LP system. The model size must already have been
     * determined.<br>
     * constraints[] organization:<br>
     * +splitW1 ... +splitWN -splitW1 ... -splitW2 (+multiW)<br>
     * +splitC1 ... +splitCN -splitC1 ... -splitC2 (+multiC)<br>
     * <br>
     * offset(+split.X,-split.X)= costSplit<br>
     *
     * @param e the edge for which the constraints are to be set
     * @deprecated currently not in use
     */
    void setConstraints(SimpleEdge e) {

        int[] a= getConstraintHash().get(e);
        assert(a!= null);
        int[] t= getAllTrptIdx();
        int restrEdgeConstr= costSplitWC? a.length/ 2: a.length;
        if (reuseIdx== null) {
            reuseIdx= new int[restrEdgeConstr+ t.length];
            System.arraycopy(t, 0, reuseIdx, reuseIdx.length - t.length, t.length);
            reuseVal= new double[reuseIdx.length];
            reuseVal[0]= 1d;	// add
            reuseVal[1]= -1d;	// substract
            if (costUseMultimap)
                reuseVal[2]= 1d;	// 0-cost substraction
        } else {
            assert(reuseIdx.length== a.length+ t.length);
            for (int i = 0; i < t.length; i++)
                reuseVal[i+ restrEdgeConstr]= 0d;	// reset trpt frac
        }

        int nrRestr= costSplitWC? 2: 1, offset= 0;
        Transcript[] tt= aMapper.decodeTset(e.getTranscripts());
        for (int i = 0; i < nrRestr; i++, offset+= restrEdgeConstr) {

            byte dir= i== 0? Constants.DIR_FORWARD: Constants.DIR_BACKWARD;

            // set individual edge constraints
            System.arraycopy(a, offset, reuseIdx, 0, restrEdgeConstr);

            // fill in transcript constraints
            double totVal= 0d;
            for (Transcript aTt : tt) {
                int tlen = aTt.getExonicLength();
                UniversalMatrix m = getMatrixMap().get(aTt.getTranscriptID());
                int[] area = e.getFrac(aTt, mappingStats.getReadLenMin(), dir);
                long reads = (long) m.get(area[0], area[1], tlen, dir);
                long sum = m.getSum(dir);
                double val = reads / (double) sum;
                totVal += val;
                if (val < 0 || Double.isNaN(val) || Double.isInfinite(val)) {
                    System.err.println("invalid val " + Double.toString(val) + " " + aMapper.trpts[0] + " " + e);
                }
                int tOffset = getConstraintHash().get(aTt)[0] - t[0];    // gotta be consecutive
                reuseVal[restrEdgeConstr + tOffset] = val;
            }
            if (totVal== 0) {
                System.err.println("edge with 0 expectation! "+aMapper.trpts[0]+" "+e);
                continue; // TODO check whether edge really should be disregarded
            }


            // set costs, bounds and cross-links for constraints
            if (costUseMultimap) {
                costIdx.add(reuseIdx[reuseIdx.length- 1]);
                costVal.add(0d);	// multi
            }
            setConstraintsCostsConstant(a, restrEdgeConstr, e);

            // create restrictions
            if (writeFile)
                writeRowLP(reuseIdx, reuseVal, LpSolve.EQ, i==1?
                        ((MappingsInterface) e).getMappings().getRevReadNr(): ((MappingsInterface) e).getMappings().getReadNr());
            else
                addConstraintToLp(reuseIdx, reuseVal, LpSolve.EQ, i==1?
                        ((MappingsInterface) e).getMappings().getRevReadNr(): ((MappingsInterface) e).getMappings().getReadNr());
            ++restrNr;

        }	// restriction iterator: i

    }

    /**
     * Sets the costs given the constraint indices of a certain
     * edge and the number of bins for increasing cost functions.<br>
     * Organization of the indices:<br>
     * +splitW1 ... +splitWN -splitW1 ... -splitW2 (+multiW)<br>
     * +splitC1 ... +splitCN -splitC1 ... -splitC2 (+multiC)<br>
     * @param a indices of the constraints for the edge
     * @param baseSize size of W/C cluster
     * @param e the edge for which
     */
    private void setConstraintsCostsConstant(int[] a, int baseSize, SimpleEdge e) {

        // costs
        for (int i = 0; i < (a.length/ baseSize); i++) { // w/o multimap, commonly handled elsewhere
            for (int j = 0; j < costSplit; j++) {
                costIdx.add(a[i*baseSize+j]);
                costVal.add(1d);	// everything costs 1
                costIdx.add(a[i*baseSize+j+costSplit]);
                costVal.add(1d);	// everything costs 1

                // adjust bounds
                if (costBounds!= null&& j== costSplit- 1) 	// last split, use stabi
                    try {

                        long obs= i== 0? ((MappingsInterface) e).getMappings().getReadNr():
                                ((MappingsInterface) e).getMappings().getRevReadNr();
                        double obs1= obs== 0? 1: obs;
                        double div= obs1/ costSplit;
                        double maxSub;
                        if(!Double.isNaN(costBounds[0])) {
                            // maxSub= div;
                            // maxSub= Math.max(0, (obs1/ costBounds[0])- ((costSplit-1)*div));
                            maxSub= Math.max(0, (obs1* costBounds[0])- ((costSplit-1)*div));
                            getLPsolve().setUpbo(a[i*baseSize+j], maxSub);	// ub of substracting values, +1 in matrix
                        }

//                      double maxAdd;
//						if (!Double.isNaN(costBounds[1])) {
//                          maxAdd= div;
//							maxAdd= Math.max(0, (obs1*costBounds[1])- ((costSplit-1)*div));
//							getLPsolve().setUpbo(a[i+ j+ costSplit], maxAdd);		// ub of adding reads, -1 in matrix
//						}

                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }

            }
        }

    }

    /**
     * Sets cost values for constraint bins that are logarithmically
     * increasing.
     *
     *
     * @param a	return value, array size is elementar block
     *          (either plus or minus)
     * @param e the edge for which
     * @deprecated orphaned
     */
    private void setConstraintsCostsLogarithmical(int[] a, SimpleEdge e) {

        long obs= ((MappingsInterface) e).getMappings().getReadNr();
        double obs1= (obs==0)?1:obs;
        double incr= obs1/ costSplit;

        double lastX= 0, lastY= 0d;
        double m;
        int costOffset= costSplit* 2- 1;
        for (int j = 0; j < costSplit; ++j) {	// m0,t1,m1,t2,m2..

            // calculate costs
            double newX= lastX+ incr;
            double newY= (-1d* Math.log((obs1- newX)/obs1));				// + weights <-> substr reads

            if (Double.isInfinite(newY)) {
                //m*= 2;	// good approximation? what bout Integer.MAXVAL?
                m= 9999;	// m still infinite
            } else {
                assert(!Double.isNaN(newY));
                assert(!Double.isInfinite(newY));
                assert(newY> 0);
                m= (newY- lastY)/ incr;
            }

            // add costs

            int mIdx= j* 2;
            costIdx.add(a[mIdx]);
            costVal.add(m);
            costIdx.add(a[mIdx+costOffset]);
            costVal.add(m);	// linear rise of segmental costs
            double max= 9999;
            if (j< costSplit- 1) {
                max= incr;
                try {
                    getLPsolve().setUpbo(a[mIdx], max);
                    getLPsolve().setUpbo(a[mIdx+costOffset], max);
                } catch (LpSolveException exx) {
                    exx.printStackTrace();
                }
            }
            if (j> 0) {
                costIdx.add(a[mIdx-1]);
                costVal.add(lastY);
                costIdx.add(a[mIdx-1+costOffset]);
                costVal.add(lastY);	// constant offset for segmental costs
                // connect bools
                try {
                    getLPsolve().setInt(a[mIdx-1], true);
                    getLPsolve().setUpbo(a[mIdx-1], 1);
                    addConstraintToLp(
                            new int[]{a[mIdx-1],a[mIdx]},
                            new double[] {-max, 1d}, LE, 0d);	// connect bool
                    ++restrNr;
                    getLPsolve().setInt(a[mIdx+costOffset-1], true);
                    getLPsolve().setUpbo(a[mIdx+costOffset-1], 1);
                    addConstraintToLp(
                            new int[]{a[mIdx+costOffset-1],a[mIdx+costOffset]},
                            new double[] {-max, 1d}, LE, 0d);	// connect bool
                    ++restrNr;
                } catch (LpSolveException exx) {
                    exx.printStackTrace();
                }
            }


            lastX= newX;
            lastY= newY;

        }	// iterate basesize (wo multi)


    }

    /**
     * Returns the mapping statistics.
     */
    public MappingStats getMappingStats() {
        return mappingStats;
    }

    /**
     * Sets statistics collected during annotation mapping.
     */
    public void setMappingStats(MappingStats mappingStats) {
        this.mappingStats= mappingStats;
    }

    /**
     * Retruns a map of transcript objects mapped to
     * their respective expression levels after
     * deconvolution and normalization.
     *
     * @return predicted expression levels after
     * deconvolution and normalization
     */
    public HashMap<Object, Double> getTrptExprHash() {
        return trptExprHash;
    }

    /**
     * Returns the value of the objective function
     * after deconvolution has been performed,
     * or (-1) if the linear system has not yet been
     * solved.
     *
     * @return value -1 before deconvolution,
     * otherwise the value of the objective function
     */
    public double getValObjFunc() {
        return valObjFunc;
    }

    /**
     * Estimates the error for the segment of a transcript
     * based on the total edge error.
     *
     * @param e edge representing a segment in the transcript
     * @param t transcript supporting the edge
     * @return error along the segment of
     * @deprecated orphaned
     */
    private double calcEdgeError(SimpleEdge e, Transcript t) {
		/*double tcov= getTrptExprHash().get(t.getTranscriptID()).doubleValue();
		int possReads= e.getPossReadNr();
		int[] c= getConstraintHash().get(e);
		double mod= 0d;
		if (c[0]!= 0d)
			mod= result[restrNr+c[0]];
		else
			mod= -result[restrNr+c[1]];
		
		double ovSig= 0d;	// sum of other trpts there
		Transcript[] tt= g.decodeTset(e.getTranscripts());				
		for (int i = 0; i < tt.length; i++) {
			if (tt[i]!= t)
				ovSig+= getTrptExprHash().get(tt[i].getTranscriptID()).doubleValue();
		}
		
		double edgeError= mod* possReads* ((((double) tcov)/(tcov+ ovSig)));
		return edgeError; */
        return 0d;
    }

    /**
     * Hash that maps transcript IDs to bias matrix.
     */
    HashMap<String, UniversalMatrix> matrixMap;

    /**
     * Constructs a map that maps the transcripts of the locus
     * by their ID to the appropriate bias matrix.
     *
     * @return hash that maps transcript IDs to bias matrix
     * @see #setConstraints(barna.model.splicegraph.SimpleEdge)
     * @deprecated orphaned
     */
    HashMap<String, UniversalMatrix> getMatrixMap() {
        if (matrixMap == null) {
            matrixMap = new HashMap<String, UniversalMatrix>(aMapper.trpts.length, 1f);
            int lenSum= 0;
            for (int i = 0; i < aMapper.trpts.length; i++)
                lenSum+= aMapper.trpts[i].getExonicLength();
            float avgLen= lenSum/ aMapper.trpts.length;
            double rpk= aMapper.getNrMappingsMapped()* 1000f/ avgLen;
            for (int i = 0; i < aMapper.trpts.length; i++) {
                int tlen= aMapper.trpts[i].getExonicLength();
                UniversalMatrix m= profile.getMatrix(tlen, rpk);
                matrixMap.put(aMapper.trpts[i].getTranscriptID(), m);
            }
        }

        return matrixMap;
    }

    /**
     * Adds a row to the LP system.
     *
     * @param idx array with indices of constraints
     * @param val array with factor values for constraints
     * @param eq equation type
     * @param cap right hand side of the equation
     */
    private void addConstraintToLp(int[] idx, double[] val, int eq, double cap) {

//		if (writeFile)
//			writeRowLP(idx, val, eq, cap);
//		else if (columnWise) {
//			;
//		} else 
        //double[] constraint= createArray(idx, val);
        try {
            getLPsolve().addConstraintex(idx.length, val, idx, eq, cap);
            //getLPsolve().addConstraint(constraint, eq, cap);
        } catch (LpSolveException e) {
            e.printStackTrace();
        }

    }

    static int nrUnderPredicted= 0, nrOverPredicted= 0;


    /**
     * Writes a line of the LP system to disk.
     *
     * @param idx array with indices of constraints
     * @param vals array with factor values for constraints
     * @param cond condition, equation type
     * @param rhs right hand side of the equation
     * @see #getLPWriter()
     */
    void writeRowLP(int[] idx, double[] vals, int cond, double rhs) {
        try {
            BufferedWriter buffy= getLPWriter();
            assert(idx.length== vals.length);
            for (int i = 0; i < vals.length; i++) {
                if (vals[i]% 1== 0)
                    buffy.write(Integer.toString((int) vals[i]));
                else
                    buffy.write(Double.toString(vals[i]));
                buffy.write(" x");
                buffy.write(Integer.toString(idx[i]));
                if (i< vals.length- 1&& vals[i+1]>= 0)
                    buffy.write(" +");
            }
            buffy.write(COND_SYMBOLS[cond]+Double.toString(rhs)+";"+ OSChecker.NEW_LINE);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Path to a file to output the LP of the current locus to.
     */
    String lpOutFName= null;

    /**
     * Returns the name and path for a file to output the LP system for the current locus to in the folder for
     * outputting LP report files. This is NOT to be confused with the linear programs written by
     * <code>getLPWriter()</code>!
     * @return path to output the LP of the current locus to, or <code>null</code> if no folder for LP files
     * has been specified
     * @see #getLPWriter()
     */
    String getLPoutFileName() {

        fileLPdir= null; //settings.get(FluxCapacitorSettings.TMP_DIR).getAbsoluteFile();
        if (lpOutFName== null&& fileLPdir!= null) {
            try {
                lpOutFName = FileHelper.createTempFile(aMapper.trpts[0].getGene().getLocusID().replace(":", "_"), SFX_LPOUT, fileLPdir).getAbsolutePath();
            } catch (IOException e) {
                throw new RuntimeException(e.getMessage());
            }
        }

        return lpOutFName;
    }

    /**
     * Writer for linear programs.
     */
    BufferedWriter lpWriter;

    /**
     * Returns the LP writer, eventually initializes it.
     * @return writer for linear programs
     * @see #writeFile
     * @deprecated init missing
     */
    BufferedWriter getLPWriter() {
        if (lpWriter == null) {
            try {
                // TODO link to settings
//				File dir= new File("I:\\solexa\\simulation\\lp");
//				fileLPinput= File.createTempFile(g.trpts[0].getTranscriptID(), ".lp", dir);
//				lpWriter = new BufferedWriter(new FileWriter(fileLPinput));
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        return lpWriter;
    }
    /**
     * Wrapper to solve the system of linear equations and return the status of the solver. Possible values are:
     *
     * <table>
     * <th><td>Message</td><td>Value</td><td>Explanation</td></th>
     * <tr><td>  NOMEMORY</td>    <td>(-2)</td>   <td>Out of memory</td></tr>
     * <tr><td>   OPTIMAL</td>    <td> (0)</td>   <td>An optimal solution was obtained</td></tr>
     * <tr><td>SUBOPTIMAL</td>    <td> (1)</td>   <td><p>The model is sub-optimal.
     * Only happens if there are integer variables and there is already an integer solution found.
     * The solution is not guaranteed the most optimal one.</p><ul>
     * <li>A timeout occured (set via set_timeout or with the -timeout option in lp_solve)</li>
     * <li>set_break_at_first was called so that the first found integer solution is found (-f option in lp_solve)</li>
     * <li>set_break_at_value was called so that when integer solution is found that is better than the specified value
     * that it stops (-o option in lp_solve)</li>
     * <li>set_mip_gap was called (-g/-ga/-gr options in lp_solve) to specify a MIP gap</li>
     * <li>An abort function is installed (put_abortfunc) and this function returned TRUE</li>
     * <li>At some point not enough memory could not be allocated</li></ul></td></tr>
     * <tr><td>INFEASIBLE</td>    <td> (2)</td>   <td>The model is infeasible</td></tr>
     * <tr><td>UNBOUNDED</td>     <td> (3)</td>   <td>The model is unbounded</td></tr>
     * <tr><td>DEGENERATE</td>    <td> (4)</td>   <td>The model is degenerative</td></tr>
     * <tr><td>NUMFAILURE</td>    <td> (5)</td>   <td>Numerical failure encountered</td></tr>
     * <tr><td>USERABORT</td>     <td> (6)</td>   <td>The abort routine returned TRUE. See put_abortfunc</td></tr>
     * <tr><td>TIMEOUT</td>       <td> (7)</td>   <td>A timeout occurred. A timeout was set via set_timeout</td></tr>
     * <tr><td>PRESOLVED</td>     <td> (9)</td>   <td>The model could be solved by presolve.
     * This can only happen if presolve is active via set_presolve</td></tr>
     * <tr><td>PROCFAIL</td>      <td>(10)</td>   <td>The B&B routine failed</td></tr>
     * <tr><td>PROCBREAK</td>     <td>(11)</td>   <td>The B&B was stopped because of a break-at-first (see
     * set_break_at_first) or a break-at-value (see set_break_at_value)</td></tr>
     * <tr><td>FEASFOUND</td>     <td>(12)</td>   <td>A feasible B&B solution was found</td></tr>
     * <tr><td>NOFEASFOUND</td>   <td>(13)</td>   <td>No feasible B&B solution found</td></tr>
     *
     * @param outFName name and absolute path of the file to which linear programs are written,
     *                 or <code>null</code> if no output to disk is to be performed
     * @return value specifying the status of the solver (cf. table in class documentation)
     */
    int solve(String outFName) {

        // TODO measure cumulative time for JNI calls
        //long t0= System.currentTimeMillis();
        //getLPsolve().printLp();
        //getLPsolve().setScaling(LpSolve.SCALE_CURTISREID);

        // call JNI
        int res= -1;
        try {
            if (outFName== null)
                // shut up! only IMPORTANT, SEVERE, CRITICAL
                getLPsolve().setVerbose(LpSolve.CRITICAL);
            else
                getLPsolve().setOutputfile(outFName);
            //getLPsolve().setPresolve(LpSolve.PRESOLVE_ROWS, 10);
            getLPsolve().setScaling(LpSolve.SCALE_DYNUPDATE);
            res= getLPsolve().solve();
        } catch (LpSolveException e) {
            throw new RuntimeException(e);
        }
        return res;
    }

    /**
     * @deprecated remove ASAP and replace by something more christian
     */
    public static boolean DEBUG= false;

    private void writeDEBUG(String s) {
        File f= new File("/home/micha/DEBUG_TEST_flux.out");
        try {
            BufferedWriter buffy= new BufferedWriter(new FileWriter(f, true));
            buffy.write(s+ "\n");
            buffy.close();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Algorithm to set up system of linear equations
     * and call the LP solver.
     */
    public strictfp void run() {

        // TODO init time management
        //long t0= System.currentTimeMillis();

        // initialize LP program
        setConstraints((byte) 0, null);
        getLPsolve();
        HashMap<String, Integer> tMap= setConstraints((byte) 1, null);

        // solve
        int ret= solve(getLPoutFileName());

        // append additional debug info
        if (DEBUG&& ret!= 0) {
            try {
                String fname= getLPoutFileName();

                getLPsolve().writeLp(fname+ "_wlp");
                getLPsolve().writeMps(fname+ "_mps");

                getLPsolve().setOutputfile(fname+ "_lp");
                getLPsolve().printLp();

                getLPsolve().setOutputfile(fname+ "_of");
                getLPsolve().printObjective();

                getLPsolve().setOutputfile(fname+ "_solv");
                getLPsolve().printSolution(1);

                // additional stream only afterwards
                PrintStream p= new PrintStream(new FileOutputStream(getLPoutFileName()+"_const", true));
                setConstraints((byte) 0, p);
                Log.warn("There was an issue with the linear problem. The linear system has been written to " + fname);
            } catch (Exception e) {
                Log.error("Failed to set lp output to:\n\t"+ getLPoutFileName(), e);
            }
        }

        // get transcription expression levels
        trptExprHash= getResult(tMap);
        // output
        if (DEBUG)
            setConstraints((byte) 2, null);

        //normalizeBack2LocusExpr(trptExprHash);
        getLPsolve().deleteLp();	// closes file outFName

        // output debug info
        if ( ret!=0) {
            Iterator<Object> idIter= trptExprHash.keySet().iterator();
            while(idIter.hasNext()) {
                Object o= idIter.next();
                if (!(o instanceof String))
                    continue;
                String id= (String) o;
                Log.debug(id + " " + trptExprHash.get(id));
            }
            Log.debug("\n");
            Log.debug("Settings:");
            Log.debug("paired-end\t" + pairedEnd);
            if (costBounds!= null) {
                if (!Double.isNaN(costBounds[0]))
                    Log.debug("cost boundaries:\tlower /" + costBounds[0]);
                if (!Double.isNaN(costBounds[1]))
                    Log.debug(" upper *" + costBounds[1]);
            }
            if (costModel== COSTS_LINEAR)
                Log.debug("costfunc\tlinear");
            else if(costModel== COSTS_LOG)
                Log.debug("costfunc\tlog");
            Log.debug(toStringConstraints());

            Log.debug("\n");
            Log.debug("Transcripts:");
            for (int i = 0; i < aMapper.trpts.length; i++) {
                StringBuilder sb= new StringBuilder(aMapper.trpts[i].getTranscriptID()+"\t");
                SpliceSite[] ss= aMapper.trpts[i].getSpliceSitesAll();
                for (SpliceSite s : ss) sb.append(s.toString());
                Log.debug(sb.toString());
            }
        }
    }

    /**
     * Generates a string representation of the constraint hash.
     *
     * @return string representing the hash of constraints.
     */
    public String toStringConstraints() {
        Iterator iter= getConstraintHash().keySet().iterator();
        StringBuilder sb= new StringBuilder();
        while (iter.hasNext()) {
            Object o= iter.next();
            if (!(o instanceof SimpleEdge))
                continue;
            sb.append(o);
            sb.append(":\t\t");
            int[] c= getConstraintHash().get(o);
            for (int aC : c) sb.append("C").append(Integer.toString(aC)).append(" ");
            sb.append("\n");
        }

        iter= getConstraintHash().keySet().iterator();
        while (iter.hasNext()) {
            Object o= iter.next();
            if (!(o instanceof Transcript))
                continue;
            sb.append(o);
            sb.append(":\t\t");
            int[] c= getConstraintHash().get(o);
            for (int aC : c) sb.append("C").append(Integer.toString(aC)).append(" ");
            sb.append("\n");
        }

        return sb.toString();
    }

    /**
     * Returns a hash with transcripts and their predicted
     * expression levels as extracted from the results of the
     * deconvolution process. Additional normalization steps
     * are performed to eliminate artifacts of the
     * deconvolution.
     *
     * @param tMap map of transcript IDs to constraint number
     * @return normalized transcript expression levels
     */
    protected HashMap<Object, Double> getResult(HashMap<String, Integer> tMap) {

        // get info about solution
        result= new double[1+ restrNr+ constraintCtr];

        try {
            getLPsolve().getPrimalSolution(result);
        } catch (LpSolveException e1) {
            e1.printStackTrace();
        }
        valObjFunc= result[0];
        resultR= new double[restrNr];
        System.arraycopy(result, 1, resultR, 0, restrNr);
        resultC= new double[constraintCtr];
        System.arraycopy(result, restrNr+ 1, resultC, 0, constraintCtr);
        trptExprHash= new HashMap<Object,Double>();

        // transcripts
        Transcript[] trpts= aMapper.trpts;
        double sum= 0;
        for (Transcript trpt : trpts) {
            int c = tMap.get(trpt.getTranscriptID());
            if (Double.isNaN(c))
                System.currentTimeMillis();
            double x = result[restrNr + c];
            double tot = mapCCheck.get(trpt.getTranscriptID());
            tot /= 2d;
            x /= tot;
            sum += x;
            assert (!Double.isNaN(x));
            trptExprHash.put(trpt.getTranscriptID(), x);
        }

        // normalizaton factor
        double nfac= nrMappingsObs/ sum;
        for (Transcript trpt : trpts) {
            double x = trptExprHash.get(trpt.getTranscriptID());
            x *= nfac;
            assert (!Double.isNaN(x));
            trptExprHash.put(trpt.getTranscriptID(), x);
        }

        // locus normalization
//        for (Transcript trpt : trpts) {
//            int tlen = trpt.getExonicLength();
//            UniversalMatrix m = profile.getMatrix(tlen);
//            double f = m.getNfactor(0.2d);
//            double x = trptExprHash.get(trpt.getTranscriptID());
//            x *= f;
//            assert (!Double.isNaN(x));
//            trptExprHash.put(trpt.getTranscriptID(), x);
//        }

        // apppend edge solutions
//        Object[] keys= constraintHash.keySet().toArray();
//        for (int i = 0; i < keys.length; i++) {
//            if (!(keys[i] instanceof SimpleEdge))
//                continue;
//            SimpleEdge e= (SimpleEdge) keys[i];
//            int[] c= constraintHash.get(e);
//            trptExprHash.put(e, ((MappingsInterface) e).getMappings().getReadNr()
//                    + ((MappingsInterface) e).getMappings().getRevReadNr()+result[restrNr+ c[0]]- result[restrNr+ c[1]]);
//        }

        return trptExprHash;
    }

    /**
     * Sets the resolution for increasing cost functions, i.e., how many cost intervals are distinguished between
     * 0 and maximum deviation.
     * @param costSplit number of intervals with different cost values
     */
    public void setCostSplit(byte costSplit) {
        if (costSplit<= 0)
            throw new IllegalArgumentException("Cost intervals have to be > 0!");
        this.costSplit = costSplit;
    }

    /**
     * Sets up a hash with that maps <code>Edge</code> respectively <code>Transcript</code> instances to an array
     * with the corresponding indices for constraints in the linear system.
     * @return a hash mapping <code>Edge</code> and <code>Transcript</code> instances to <code>int[]</code> instances
     * storing the constraint indices
     */
    public Hashtable<Object,int[]> getConstraintHash() {

        if (constraintHash == null) {

            constraintHash= new Hashtable<Object,int[]>();

            // edges
            AbstractEdge[] edges= aMapper.getExonicEdgesInGenomicOrder();
            for (AbstractEdge edge : edges) {
                if (!edge.isExonic())
                    continue;
                // add the edge itself
                if ((!pairedEnd) && edge.length() >= mappingStats.getReadLenMin()) {

                    constraintHash.put(edge, getConstraintIDs());
                }
                for (int j = 0; edge.getSuperEdges() != null && j < edge.getSuperEdges().size(); j++) {
                    SuperEdge se = edge.getSuperEdges().elementAt(j);
                    if (se.getEdges()[0] != edge)
                        continue;
                    assert (!constraintHash.containsKey(se)); // paired ends are iterated twice, ej maybe more
                    if (se.isPend()) {
                        if (!pairedEnd)
                            continue;    // should not incur, PEs never added
                        // else add
                    } else {
                        if (pairedEnd) {    // paired-end of EJ
                            for (int k = 0; se.getSuperEdges() != null && k < se.getSuperEdges().size(); k++) {
                                SuperEdge se2 = se.getSuperEdges().elementAt(k);
                                assert (se2.isPend());
                                if (se2.getEdges()[0] != se)
                                    continue;
                                assert (!constraintHash.containsKey(se2));
                                constraintHash.put(se2, getConstraintIDs());
                            }
                            continue;
                        } // else add
                    }

                    // add EJs eventually PEs
                    constraintHash.put(se, getConstraintIDs());

                }
            }

            // transcripts
            //System.err.println("hash "+constraintHash.size());
            Transcript[] trpts= aMapper.trpts;
            for (Transcript trpt : trpts) constraintHash.put(trpt, new int[]{++constraintCtr});

        }

        return constraintHash;
    }

    /**
     * Specifies the folder to store LP debug information files in, one per locus.
     * @param dir the folder to store LP files
     * @see #getLPoutFileName()
     */
    public void setFileLPdir(File dir) {
        this.fileLPdir = dir;
    }


    /**
     * Subroutine to retrieve constraints for a specific edge (and its super-edges) and a certain transcript supporting
     * that edge.
     *
     * @param e the base edge
     * @param sig signature of a transcript supporting that edge
     * @param v vector to cross-link transcript fraction of edge with super-edges
     * @param mapE hash that maps edges to vectors storing their constraint indices; the vectors for the base edge and
     *             super-edges are extended within the method
     * @param sense flag to distinguish between anti-/sense deconvolution along that edge
     * @param count flag to indicate whether only counting of constraint indices is performed
     */
    private void getConstraints(AbstractEdge e, long[] sig, IntVector v,
                                HashMap<AbstractEdge, IntVector> mapE, boolean sense, byte count) {

        // for the edge itself
        IntVector w= mapE.get(e);
        if (count== 0)
            ++constraintCtr;
        else if (count== 1|| count== 2) {
            v.add(++constraintCtr);	// for transcript fraction
            if (w== null)
                w= new IntVector();
            w.add(constraintCtr);	// edge consistency
        }
        mapE.put(e, w);

        // iterate super-edges
        for (int j = 0; e.getSuperEdges()!= null&& j < e.getSuperEdges().size(); j++) {
            SuperEdge se= e.getSuperEdges().elementAt(j);
            // skip intronic super-edges (wrong deconvolution)
            if (se.isPend()&& se.isIntronic())
                continue;
            if (SplicingGraph.isNull(SplicingGraph.intersect(se.getTranscripts(), sig)))
                continue;
            // sense/anti-sense.. e must be first/last in super-edge
            if ((sense&& se.getEdges()[0]!= e)|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
                continue;
            if (count== 0) {
                ++constraintCtr;
            } else {
                v.add(++constraintCtr);	// for transcript fraction
            }
            w= mapE.get(se);
            if (count== 1|| count== 2) {
                if (w== null)
                    w= new IntVector();
                w.add(constraintCtr); // for edge consistency
            }
            mapE.put(se, w);

            if (se.isPend())
                continue;	// no super-edges

            for (int k = 0; se.getSuperEdges()!= null&& k < se.getSuperEdges().size(); k++) {
                SuperEdge se2= se.getSuperEdges().elementAt(k);
                assert(se2.isPend());
                if (se2.isPend()&& se2.isIntronic())
                    continue;
                if (SplicingGraph.isNull(SplicingGraph.intersect(se2.getTranscripts(), sig)))
                    continue;
                // sense/anti-sense.. e must be first/last in super-edge
                if ((sense&& se2.getEdges()[0]!= se)|| ((!sense)&& se2.getEdges()[se2.getEdges().length- 1]!= se))
                    continue;
                if (count== 0) {
                    ++constraintCtr;
                } else {
                    v.add(++constraintCtr);	// tx
                }
                w= mapE.get(se2);
                if (count== 1|| count== 2) {
                    if (w== null)
                        w= new IntVector();
                    w.add(constraintCtr); // for edge consistency
                }
                mapE.put(se2, w);
            }
        }


    }

    /**
     * <code>true</code> to optimize flux,
     * <code>false</code> to optimize flow
     */
    boolean flux= true;

    HashMap<String, Double> mapCCheck= null;

    /**
     * Iterates the constraints for all edges and, counts them (<code>count</code> is <code>true</code>) or adds them to
     * the system of linear equations (<code>count</code> is <code>false</code>). Also the cost weights in the
     * objective function are set. Finally, a hash that maps transcript IDs to their corresponding expression value
     * constraint indices is provided.
     *
     * @param count flag, if <code>count</code> is <code>0</code> only counting of indices is performed.
     *              For value <code>1</code>, restrictions on the respective contraints are added to the linear program.
     *              For value <code>2</code>, output is performed.
     * @param p stream to output additional debug info for LP reports
     * @return hash that maps transcript IDs to the contraint indices corresponding to their expression levels
     */
    public HashMap<String, Integer> setConstraints(byte count, PrintStream p) {

        // DEBUG stuff
        int segmentCounter= 0;
        // DEBUG: maps edges to their segment number
        HashMap<SimpleEdge, Integer> segmentHash= null;
        HashMap<Integer, Integer> hashCxTx= null;
        HashMap<Transcript, Integer> hashTxNr= null;
        HashMap<Integer, double[]> txError= null;

        constraintCtr= 0;
        restrNr= 0;

        // transcript constraint variables
        Transcript[] trpts= aMapper.trpts;
        // w cost vector
        IntVector w= null;
        HashMap<String, Integer> tMap= null;
        if (count== 0)
            constraintCtr+= trpts.length;
        else {
            try {
                getLPsolve().setLpName(trpts[0].getTranscriptID());
            } catch (LpSolveException e1) {
                e1.printStackTrace();
            }
            tMap= new HashMap<String, Integer>(trpts.length* 2);
            for (Transcript trpt : trpts) {
                tMap.put(trpt.getTranscriptID(), ++constraintCtr);
                if (p != null)
                    p.println(trpt.getTranscriptID() + "\t" + constraintCtr);
            }
            w= new IntVector();	// indices for cost function
        }

        // iterate edges
        AbstractEdge[] edges= aMapper.getExonicEdgesInGenomicOrder();
        HashMap<Transcript, StringBuilder> txSegments= null, txEdges= null;
        if (count== 1) {
            mapCCheck= new HashMap<String, Double>();
            for (Transcript trpt : trpts) mapCCheck.put(trpt.getTranscriptID(), 0d);
        } else if (count== 2) {
            hashCxTx= new HashMap<Integer, Integer>();
            hashTxNr= new HashMap<Transcript, Integer>();
            txSegments= new HashMap<Transcript, StringBuilder>();
            txEdges= new HashMap<Transcript, StringBuilder>();
            segmentHash= new HashMap<SimpleEdge, Integer>();
            txError= new HashMap<Integer, double[]>();
            int i= 0;
            for(Transcript tx: trpts) {
                txSegments.put(tx, new StringBuilder());
                txEdges.put(tx,new StringBuilder());
                hashTxNr.put(tx, i);
                txError.put(i, new double[6]); // (sense,asense) x (OBS,SUB,ADD)
                ++i;
            }
            for (AbstractEdge e : edges) {
                if ((!e.isExonic())|| (!(e instanceof  SimpleEdge)))
                    continue;
                segmentHash.put((SimpleEdge) e, ++segmentCounter);   // DEBUG output
            }
        }

        // iterates only exonic segments
        for (AbstractEdge e : edges) {
            if (!e.isExonic())
                continue;

            // the base edge
            Transcript[] tt = aMapper.decodeTset(e.getTranscripts());
            // sense/anti
            for (int sa = 0; sa < 2; ++sa) {

                HashMap<AbstractEdge, IntVector> mapE = new HashMap<AbstractEdge, IntVector>();

                /* === Transcript Contributions === */

                for (Transcript aTt : tt) {
                    setConstraints(e, aTt, sa, count, tMap, mapE, txSegments, segmentHash, hashTxNr, hashCxTx);
                } // iterate transcripts


                /* === Edge Deconvolution === */
                double[] sumSEG= null;   // call-by-return: [0] sumObs, [1] sumPlus, [2] sumMinus
                if (count== 2) {
                    sumSEG= new double[3];
                    Arrays.fill(sumSEG, 0d);
                }
                AbstractEdge[] ee = new AbstractEdge[mapE.size()];
                mapE.keySet().toArray(ee);
                for (AbstractEdge f : ee) {
                    setConstraints(e, f, sa, count, sumSEG, tt, w,
                            mapE, txSegments, txEdges, segmentHash, hashTxNr, hashCxTx, txError);
                } // all edges in segment

                //close segment
                if (count== 2) {
                    Iterator<Transcript> i2= txSegments.keySet().iterator();
                    while (i2.hasNext()) {
                        Transcript tx= i2.next();
                        if (!SplicingGraph.intersects(aMapper.encodeTset(tx), e.getTranscripts()))
                            continue;
                        StringBuilder sb= txSegments.get(tx);
                        StringBuilder sb2= new StringBuilder(sb.length());
                        String[] ss= sb.toString().split("\t");
                        for (int i = 0; i < ss.length- 1; i++) {
                            sb2.append(ss[i]+ "\t");
                        }
                        StringTokenizer st= new StringTokenizer(ss[ss.length- 1], ";");
                        sb2.append(st.nextToken());
                        for (int i = 1; i < 4; i++)
                            sb2.append(";"+ st.nextToken());
                        // TODO can use as a check
                        //sb2.append(";DEC2="+ Math.round(sumTx[hashTxNr.get(tx)]* 100.00)/ 100.00);
                        sb2.append(";OBS="+ (Math.round(sumSEG[0]* 100.00)/ 100.00)
                                + ";SUB="+ (Math.round(sumSEG[1]* 100.00)/ 100.00)
                                + ";ADD="+ Math.round(sumSEG[2]* 100.00)/ 100.00);
                        while (st.hasMoreElements())
                            sb2.append(";"+ st.nextToken());
                        txSegments.put(tx, sb2);
                    }

                }

            }   // anti/sense
        } // end all edges

        if (count== 0)
            return null;

        if (count== 2) {
            double[] locusErr= new double[6];
            Arrays.fill(locusErr, 0d);
            Iterator<double[]> ii= txError.values().iterator();
            while(ii.hasNext()) {
                double[] tmp= ii.next();
                for (int i = 0; i < tmp.length; i++)
                    locusErr[i]+= tmp[i];
            }
            double locusTotSense= locusErr[0]- locusErr[1]+ locusErr[2],
                    locusTotAsense= locusErr[3]- locusErr[4]+ locusErr[5],
                    locusTot= locusTotSense+ locusTotAsense;
            double locusTotNor= 0d;
            for (int i = 0; i < trpts.length; i++)
                locusTotNor+= getTrptExprHash().get(trpts[i].getTranscriptID());
            String locusPfx= trpts[0].getGene().getLocusID()+ "\t"
                    + "LNOR="+ locusTotNor
                    + ";LEXP="+ Math.round(locusTot* 100.00)/ 100.00
                    + ";OBS="+ Math.round((locusErr[0]+ locusErr[3])* 100.00)/ 100.00
                    + ";SUB="+ Math.round((locusErr[1]+ locusErr[4])* 100.00)/ 100.00
                    + ";ADD="+ Math.round((locusErr[2]+ locusErr[5])* 100.00)/ 100.00
                    + ";OSN="+ Math.round(locusErr[0]* 100.00)/ 100.00
                    + ";SSN="+ Math.round(locusErr[1]* 100.00)/ 100.00
                    + ";ASN="+ Math.round(locusErr[2]* 100.00)/ 100.00
                    + ";OAS="+ Math.round(locusErr[3]* 100.00)/ 100.00
                    + ";SAS="+ Math.round(locusErr[4]* 100.00)/ 100.00
                    + ";AAS="+ Math.round(locusErr[5]* 100.00)/ 100.00
                    + "\t";

            Iterator<Transcript> i2= txSegments.keySet().iterator();
            while (i2.hasNext()) {
                Transcript tx= i2.next();
                StringBuilder sb= txSegments.get(tx);
                sb.append(txEdges.get(tx));

                int c = tMap.get(tx.getTranscriptID());
                int i = hashTxNr.get(tx);
                double[] err= txError.get(i);
                sb.insert(0, locusPfx+ tx.getTranscriptID()+ "\t"
                        + "TX="+ (i+ 1)
                        + ";TNOR="+ Math.round(getTrptExprHash().get(tx.getTranscriptID())* 100.00)/ 100.00
                        + ";TEXP="+ Math.round(resultC[c- 1]* 200.00)/ 100.00   // *2 to get to read-base
                        + ";OBS="+ Math.round((err[0]+ err[3])* 100.00)/ 100.00
                        + ";SUB="+ Math.round((err[1]+ err[4])* 100.00)/ 100.00
                        + ";ADD="+ Math.round((err[2]+ err[5])* 100.00)/ 100.00
                        + ";OSN="+ Math.round(err[0]* 100.00)/ 100.00
                        + ";SSN="+ Math.round(err[1]* 100.00)/ 100.00
                        + ";ASN="+ Math.round(err[2]* 100.00)/ 100.00
                        + ";OAS="+ Math.round(err[3]* 100.00)/ 100.00
                        + ";SAS="+ Math.round(err[4]* 100.00)/ 100.00
                        + ";AAS="+ Math.round(err[5]* 100.00)/ 100.00
                );

                try {
                    BufferedWriter buffy= new BufferedWriter(
                            new FileWriter("/Volumes/Raptor/scratch/simulations/hg19_gencode_paired_sorted_tx_dbug.txt", true));
                    buffy.write(sb.toString()+ "\n");
                    buffy.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            return null;
        }

        // set objective function/costs
        //double[] a= createArray(w.toIntArray());	// linear costs
        //double min= Math.exp(-nrMappingsObs);
        double[] costs= new double[w.size()];
        for (int i = 0; i < costs.length; i++) {
            //double x= w.get(i);
            //costs[i]= 1+ Math.exp(-Math.log(1+ x));	// neutralizes
            //costs[i]= 1d+ Math.pow(x+1d, -1d/2d);
            //costs[i]= 0;	//
            //costs[i]= 1d+ Math.log(x+ 1);
            //costs[i]= 1d+ Math.exp(-x);
            //costs[i]= 100d/ (x+ 1);	// last
            costs[i]= 1d;

            //costs[i]= 1d+ Math.log(x+ 1);
            //costs[i]= 1d+ Math.sqrt(x+ 1);	// best
            //costs[i]= 1d+ Math.pow(x+ 1, 1d/ 3d);

            //costs[i]= (Math.log(x+ 1)/ (x+ 1));	// logdiv
            //costs[i]= 1d/ (1d+ Math.log(x+ 1d));	// divlog
        }
        double[] a= createArray(w.toIntArray(), costs);

        try {
            getLPsolve().setObjFn(a);
            getLPsolve().setMinim();
        } catch (LpSolveException e) {
            e.printStackTrace();
        }

        // consistency check
        Object[] oo= mapCCheck.keySet().toArray();
        for (Object anOo : oo) {
            double val = mapCCheck.get(anOo);
            if (Math.abs(2d - val) > 0.2d)
                Log.warn("Fraction inconsistency "+ anOo+ "\t"+ val);
        }

        return tMap;
    }

    protected void setConstraints(AbstractEdge e, Transcript aTt, int sa, byte count,
                                  HashMap<String, Integer> tMap,
                                  HashMap<AbstractEdge, IntVector> mapE,
                                  HashMap<Transcript, StringBuilder> txSegments,
                                  HashMap<SimpleEdge, Integer> segmentHash,
                                  HashMap<Transcript, Integer> hashTxNr,
                                  HashMap<Integer, Integer> hashCxTx) {

        IntVector v= new IntVector(); // re-use to avoid GC?
        StringBuilder txSegmentBuilder= null;
        if (count== 2) {
            txSegmentBuilder= txSegments.get(aTt);
            txSegmentBuilder.append("\tSEG=" + segmentHash.get(e) + ";DIR=" + (sa == 0 ? "sense" : "anti"));
        }


        long[] sig = aMapper.encodeTset(aTt);
        int saveCCtr= constraintCtr;
        getConstraints(e, sig, v, mapE, sa == 0, count);
        if (count== 2) {
            for (int i = saveCCtr+ 1; i <= constraintCtr; i++)
                hashCxTx.put(i, hashTxNr.get(aTt));
        }


        if (count== 1|| count== 2) {
            int[] idx = new int[v.length + 1]; // obs parts+ tx frac
            System.arraycopy(v.vector, 0, idx, 0, v.length);
            idx[idx.length - 1] = tMap.get(aTt.getTranscriptID());
            double[] val = new double[idx.length];
            Arrays.fill(val, 1d);
            int tlen = aTt.getExonicLength();
            UniversalMatrix m = profile.getMatrix(tlen);

            int e1= aTt.getExonicPosition(e.getDelimitingPos(true));
            int e2= aTt.getExonicPosition(e.getDelimitingPos(false));
            byte dir= sa == 0 ? Constants.DIR_FORWARD : Constants.DIR_BACKWARD;
            double f = m.getFrac(e1, e2, tlen, dir);
            if (count== 2) {
                txSegmentBuilder.append(";TEX=" + Math.round(f * 10000.0) / 100.0);
                double sum= 0d;
                for (int i = 0; i < idx.length- 1; i++)
                    sum+= resultC[idx[i]- 1];
                txSegmentBuilder.append(";DEC=" + Math.round(sum * 100.00) / 100.0);
            }
            assert (!(Double.isInfinite(f) || Double.isNaN(f)));
            mapCCheck.put(aTt.getTranscriptID(),
                    mapCCheck.get(aTt.getTranscriptID()) + f);
            val[val.length - 1] = -f;
            if (debug && count== 1) {
                StringBuilder sb = new StringBuilder(e.toString());
                sb.append(": ");
                for (int k = 0; k < idx.length; k++) {
                    sb.append(val[k] > 0 ? "+" : "");
                    sb.append(val[k] % 1 == 0 ? ((int) val[k]) : val[k]);
                    sb.append("C");
                    sb.append(idx[k]).append(" ");
                }
                sb.append("= 0");
                Log.debug(sb.toString());
            }
            if (count== 1)
                addConstraintToLp(idx, val, LpSolve.EQ, 0);
            ++restrNr;

        }

    }

    protected void setConstraints(AbstractEdge e, AbstractEdge f, int sa, byte count,
                                  double[] sumSEG,
                                  Transcript[] tt,
                                  IntVector w,
                                  HashMap<AbstractEdge, IntVector> mapE,
                                  HashMap<Transcript, StringBuilder> txSegments,
                                  HashMap<Transcript, StringBuilder> txEdges,
                                  HashMap<SimpleEdge, Integer> segmentHash,
                                  HashMap<Transcript, Integer> hashTxNr,
                                  HashMap<Integer, Integer> hashCxTx,
                                  HashMap<Integer, double[]> txError) {

        boolean paird = (f instanceof SuperEdge) && ((SuperEdge) f).isPend();

        int nr = ((paird || sa == 0) ? ((MappingsInterface) f).getMappings().getReadNr()
                : ((MappingsInterface) f).getMappings().getRevReadNr());

        IntVector v = mapE.remove(f);
        if (count== 0)
            constraintCtr += 2;
        else if (count== 1|| count== 2) {

            int effLen= f.getEffLength((sa==0? Constants.DIR_FORWARD: Constants.DIR_BACKWARD), mappingStats.getReadLenMax());
            effLen= (effLen== 0? 1: effLen);    // prevent from div-by-0
            double obs= (flux? (nr/ (double) effLen): nr);
            assert((!Double.isInfinite(obs))&& (!Double.isNaN(obs)));

            int[] idx = new int[v.length + 2];    // +/-
            System.arraycopy(v.vector, 0, idx, 0, v.length);
            int c = ++constraintCtr;

            // plus, substracts from obs
            if (paird || !pairedEnd) {
                w.add(c);
            }
            idx[idx.length - 2] = c;

            // prevent from substracting complete observation
            double lim = (paird || !pairedEnd) ? Math.max(nr - (1d), 0) : nr;
            //assert(effLen> 0|| nr== 0); // might occur for clipped mappings
            if (flux)
                lim/= effLen;
            assert(lim>= 0&& (!Double.isInfinite(lim))&& (!Double.isNaN(lim)));
            // TODO
            if (count== 1)
                try {
                    getLPsolve().setUpbo(constraintCtr, lim);
                } catch (LpSolveException e1) {
                    e1.printStackTrace();
                }


            // minus, adds reads
            // do not limit adding, might cause unsolvable systems
            c = ++constraintCtr;
            w.add(c);
            idx[idx.length - 1] = c;

            // output
            if (count== 2) {
                for (Transcript tx: tt) {
                    double  sub= resultC[constraintCtr- 2],
                            add= resultC[constraintCtr- 1];  // resultC is 0-based

                    StringBuilder sb= txEdges.get(tx);
                    String eid= getEID(f, segmentHash);
                    sb.append("\tEID=E"+ eid);
                    sb.append(";SEG="+ segmentHash.get(e)+ ";DIR="+ (sa == 0 ? "sense" : "anti"));

                    sb.append(";OBS="+ Math.round(obs* 100.00)/ 100.00);
                    sb.append(";SUB="+ Math.round(sub* 100.00)/ 100.00);
                    sb.append(";ADD="+ Math.round(add* 100.00)/ 100.00); // +1 -1

                    double tot= obs+ add- sub;
                    for (int i = 0; i < idx.length- 2; i++) {
                        int tc= hashCxTx.get(idx[i]);

                        double dec= resultC[idx[i]- 1],
                                frac= (tot==0? 0d: (dec/ tot));
                        sb.append(";TX"+ (tc+ 1)+ "="+ Math.round(dec* 100.00)/ 100.00);

                        if (tc!= hashTxNr.get(tx))
                            continue;   // only this transcript

                        double[] err= txError.get(tc);
                        err[sa* 3+ 0]+= obs* frac;
                        err[sa* 3+ 1]+= sub* frac;
                        err[sa* 3+ 2]+= add* frac;

                        //else
                        StringBuilder segb= txSegments.get(tx);
                        segb.append(";E"+ eid+ "="+ Math.round(resultC[idx[i]- 1]* 100.00)/ 100.00);
                    }

                }
                sumSEG[0]+= obs;
                sumSEG[1]+= resultC[constraintCtr- 2];
                sumSEG[2]+= resultC[constraintCtr- 1];  // resultC is 0-based
            }

            // contribution weights
            double[] val = new double[idx.length];
            double x= (flux? 1d/ effLen: 1d);
            assert(x> 0&& (!Double.isInfinite(x))&& (!Double.isNaN(x)));

            Arrays.fill(val, x);
            val[val.length - 2] = 1d;
            val[val.length - 1] = -1d;
            if (count== 1) {
                addConstraintToLp(idx, val, LpSolve.EQ, obs);
            }
            ++restrNr;
        }   // count > 0


    }

    private String getEID(AbstractEdge f, HashMap<SimpleEdge, Integer> segmentHash) {

        if (f instanceof SimpleEdge) {
            return segmentHash.get(f).toString();
        }
        SuperEdge se= (SuperEdge) f;
        if (se.isPend())
            return getEID(se.getEdges()[0], segmentHash)+ ","+ getEID(se.getEdges()[1], segmentHash);

        // else
        StringBuilder sb= new StringBuilder(getEID(se.getEdges()[0], segmentHash));
        for (int i = 1; i < se.getEdges().length; i++) {
            String sep= (se.getEdges()[i- 1].getHead()== se.getEdges()[i].getTail()? ".": ":");
            sb.append(sep+ getEID(se.getEdges()[i], segmentHash));
        }
        return sb.toString();
    }

}
