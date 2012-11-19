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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;

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
	static boolean debug= false;    // TODO replace by JUnit tests

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
     * @param readLen minimum length of a read
     * @param realReads number of reads observed after annotation mapping
     */
	private GraphLPsolver(AnnotationMapper aMapper, int readLen, int realReads) {
		this.aMapper= aMapper;
		this.readLen= readLen;
		this.nrMappingsObs= realReads;
	}

    /**
     * Extended constructor, with read and insert size attributes.
     * @param aMapper splicing graph with annotation mapped reads
     * @param readLen minimum length of a read
     * @param insertMinMax set of {minimum,maximum} of the observed insert size
     * @param realReads number of reads observed after annotation mapping
     * @param considerBothStrands flag for consideration of reads mapping in anti-/sense
     * @param pairedEnd flag for paired-end reads
     */
	public GraphLPsolver(AnnotationMapper aMapper, int readLen, int[] insertMinMax, int realReads, 
			boolean considerBothStrands, boolean pairedEnd) {
		this(aMapper, readLen, realReads);
		this.costSplitWC= considerBothStrands;
		this.pairedEnd= pairedEnd;
		this.insertMinMax= insertMinMax;
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
                int[] area = e.getFrac(aTt, getReadLen(), dir);
                int reads = m.get(area[0], area[1], tlen, dir);
                int sum = m.getSum(dir);
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
     * Returns the (minimal) length of read mappings.
     */
    public int getReadLen() {
		return readLen;
	}

    /**
     * Overwrites the (minimal) length of read mappings.
     */
    public void setReadLen(int readLen) {
		this.readLen = readLen;
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
			float rpk= aMapper.getNrMappingsMapped()* 1000f/ avgLen;
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

        if (lpOutFName == null && fileLPdir != null) {
            lpOutFName = fileLPdir + File.separator
                    + aMapper.trpts[0].getGene().getLocusID().replace(":", "_")
                    + SFX_LPOUT;
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
            res= getLPsolve().solve();
        } catch (LpSolveException e) {
            throw new RuntimeException(e);
        }
        return res;
    }


    /**
     * Algorithm to set up system of linear equations
     * and call the LP solver.
     */
	public strictfp void run() {

        debug= false;
        // TODO init time management
		//long t0= System.currentTimeMillis();

		// initialize LP program
        setConstraints(true, null);
        getLPsolve();
        constraintCtr= 0;
        HashMap<String, Integer> tMap= setConstraints(false, null);

        // solve
        int ret= solve(getLPoutFileName());
        if (ret!= 0)
            debug= false;

		// append additional debug info
		if (debug || getLPoutFileName() != null) {
			getLPsolve().printLp();
			getLPsolve().printObjective();
			getLPsolve().printSolution(1);

			// additional stream only afterwards
			try {
				PrintStream p= new PrintStream(new FileOutputStream(getLPoutFileName(), true));
                setConstraints(true, p);
            } catch (Exception e) {
                Log.error("[FATAL] failed to set lp output to:\n\t"+ getLPoutFileName(), e);
			}
		}
			
		// get transcription expression levels		
		trptExprHash= getResult(tMap);
		//normalizeBack2LocusExpr(trptExprHash);
		getLPsolve().deleteLp();	// closes file outFName
		
		// output debug info
		if (debug) {
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
        for (Transcript trpt : trpts) {
            int tlen = trpt.getExonicLength();
            UniversalMatrix m = profile.getMatrix(tlen);
            double f = m.getNfactor(0.2d);
            double x = trptExprHash.get(trpt.getTranscriptID());
            x *= f;
            assert (!Double.isNaN(x));
            trptExprHash.put(trpt.getTranscriptID(), x);
        }

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
                if ((!pairedEnd) && edge.length() >= getReadLen()) {

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
     * @param p stream to output additional reporting information for introspection
     */
	private void getConstraints(AbstractEdge e, long[] sig, IntVector v,
			HashMap<AbstractEdge, IntVector> mapE, boolean sense, boolean count, PrintStream p) {

		// for the edge itself
		IntVector w= mapE.get(e);
		if (count) 
			++constraintCtr;
		else {
			v.add(++constraintCtr);	// for transcript fraction
			if (w== null)
				w= new IntVector();
			w.add(constraintCtr);	// edge consistency
		}
		mapE.put(e, w);
		if (count&& p!= null)
			p.println(constraintCtr);

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
			if (count) {
				++constraintCtr;
				if (p!= null)
					p.println(se+"\t"+(sense?"sense":"asense")+"\t"+constraintCtr);
			} else {
				v.add(++constraintCtr);	// for transcript fraction
			}
			w= mapE.get(se);
			if (!count) {
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
				if (count) {
					++constraintCtr;
					if (p!= null)
						p.println(se2+"\t"+(sense?"sense":"asense")+"\t"+constraintCtr);
				} else {
					v.add(++constraintCtr);	// tx
				}
				w= mapE.get(se2);
				if (!count) {
					if (w== null)
						w= new IntVector();
					w.add(constraintCtr); // for edge consistency
				}
				mapE.put(se2, w);
			}
		}
		
		
	}
	
	
	HashMap<String, Double> mapCCheck= null;

	/**
     * Iterates the constraints for all edges and, counts them (<code>count</code> is <code>true</code>) or adds them to
     * the system of linear equations (<code>count</code> is <code>false</code>). Also the cost weights in the
     * objective function are set. Finally, a hash that maps transcript IDs to their corresponding expression value
     * constraint indices is provided.
     *
     * @param count flag, if <code>count</code> is <code>0</code> only counting of indices is performed. Otherwise,
     *              restrictions on the respective contraints are added to the linear program
     * @param p stream to output additional debug info for LP reports
     * @return hash that maps transcript IDs to the contraint indices corresponding to their expression levels
	 */
	public HashMap<String, Integer> setConstraints(boolean count, PrintStream p) {
		
		// transcript constraint variables
		Transcript[] trpts= aMapper.trpts;
		IntVector v= null, w= null, u= null;
		HashMap<String, Integer> tMap= null;
		if (count)
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
			v= new IntVector();	// indices for transcript/part 
			w= new IntVector();	// indices for cost function
			u= new IntVector();	// observation, bases for cost function
		}
		
		// iterate edges
		AbstractEdge[] edges= aMapper.getExonicEdgesInGenomicOrder();
		if (!count) {
			mapCCheck= new HashMap<String, Double>();
            for (Transcript trpt : trpts) mapCCheck.put(trpt.getTranscriptID(), 0d);
		}
        for (AbstractEdge e : edges) {
            if (!e.isExonic())
                continue;

            // the base edge
            Transcript[] tt = aMapper.decodeTset(e.getTranscripts());
            // sense/anti
            for (int sa = 0; sa < 2; ++sa) {

                HashMap<AbstractEdge, IntVector> mapE = new HashMap<AbstractEdge, IntVector>();    // BUG?

                for (Transcript aTt : tt) {
                    if (!count) {
                        v.removeAll();
                    } else if (p != null)
                        p.print(e + "\t" + (sa == 0 ? "sense" : "asense") + "\t" + aTt + "\t");

                    long[] sig = aMapper.encodeTset(aTt);
                    getConstraints(e, sig, v, mapE, sa == 0, count, p);

                    // add transcript constraint
                    if (!count) {
                        int[] idx = new int[v.length + 1]; // obs parts+ tx frac
                        System.arraycopy(v.vector, 0, idx, 0, v.length);
                        idx[idx.length - 1] = tMap.get(aTt.getTranscriptID());
                        double[] val = new double[idx.length];
                        Arrays.fill(val, 1d);
                        int tlen = aTt.getExonicLength();
                        UniversalMatrix m = profile.getMatrix(tlen);
                        double f = m.getFrac(
                                aTt.getExonicPosition(e.getDelimitingPos(true)),
                                aTt.getExonicPosition(e.getDelimitingPos(false)),
                                tlen,
                                sa == 0 ? Constants.DIR_FORWARD : Constants.DIR_BACKWARD);

                        assert (!(Double.isInfinite(f) || Double.isNaN(f)));
                        mapCCheck.put(aTt.getTranscriptID(),
                                mapCCheck.get(aTt.getTranscriptID()) + f);
                        val[val.length - 1] = -f;
                        if (debug && !count) {
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
                        addConstraintToLp(idx, val, LpSolve.EQ, 0);
                        ++restrNr;
                    }
                } // iterate transcripts

                // add edge constraints
                AbstractEdge[] ee = new AbstractEdge[mapE.size()];
                mapE.keySet().toArray(ee);

                // total obs
/*				int sumObs= 0;
				for (int j = 0; j < ee.length; j++) {
					Edge f= ee[j];
					boolean paird= (f instanceof SuperEdge)&& ((SuperEdge) f).isPend();
					int nr= ((paird|| sa== 0)? f.getReadNr(): f.getRevReadNr());
					sumObs+= nr;
				}
*/
                for (AbstractEdge f : ee) {
                    boolean paird = (f instanceof SuperEdge) && ((SuperEdge) f).isPend();

                    int nr = ((paird || sa == 0) ? ((MappingsInterface) f).getMappings().getReadNr()
                            : ((MappingsInterface) f).getMappings().getRevReadNr());
                    v = mapE.remove(f);
                    if (count)
                        constraintCtr += 2;
                    else {
                        int[] idx = new int[v.length + 2];    // +/-
                        System.arraycopy(v.vector, 0, idx, 0, v.length);
                        int c = ++constraintCtr;
                        // plus on not paired edges at 0-cost, it substracts from obs
                        if (paird || !pairedEnd) {
                            w.add(c);
                            u.add(nr);
                        }
                        idx[idx.length - 2] = c;
                        // plus has to be limited, it substracts
                        int lim = (paird || !pairedEnd) ? Math.max(nr - 1, 0) : nr;
                        try {
                            getLPsolve().setUpbo(constraintCtr, lim);
                        } catch (LpSolveException e1) {
                            e1.printStackTrace();
                        }

                        c = ++constraintCtr;
                        // do not limit adding, even with f= 100 unsolvable systems
                        // adding reads always costs, also on single edges
                        w.add(c);
                        u.add(nr);
                        idx[idx.length - 1] = c;
                        double[] val = new double[idx.length];
                        Arrays.fill(val, 1d);
                        val[val.length - 1] = -1d;
                        if (debug && !count) {
                            StringBuilder sb = new StringBuilder(f.toString());
                            sb.append(": ");
                            for (int k = 0; k < idx.length; k++) {
                                sb.append(val[k] == 1 ? "+C" : "-C");
                                sb.append(idx[k] + " ");
                            }
                            sb.append("= ").append(nr);
                            Log.debug(sb.toString());
                        }
                        addConstraintToLp(idx, val, LpSolve.EQ, nr);
                        ++restrNr;
                    }
                }
            }
        } // end all edges
		
		if (count)
			return null;
		
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
	
}
