package fbi.genome.sequencing.rnaseq.reconstruction;

import fbi.genome.model.SpliceSite;
import fbi.genome.model.Transcript;
import fbi.genome.model.commons.DoubleVector;
import fbi.genome.model.commons.IntVector;
import fbi.genome.model.constants.Constants;
import fbi.genome.model.splicegraph.Edge;
import fbi.genome.model.splicegraph.Graph;
import fbi.genome.model.splicegraph.SuperEdge;
import lpsolve.LpSolve;
import lpsolve.LpSolveException;

import java.io.*;
import java.util.*;
//import gphase.solexa.simulation.one.RandomReadSimulator;

/**
 * ok, this version gives freedom to the expectation model.
 * @author micha
 *
 */
public class GraphLPsolver2 implements ReadStatCalculator {
	
	static final String[] COND_SYMBOLS= new String[] {"", " <= ", " >= ", " = "}; 
	public final static String[] COSTS_NAMES= new String[] {"LIN", "LOG"};
	public final static byte COSTS_LINEAR= 0, COSTS_LOG= 1;
	public final static String SFX_LPOUT= ".lp";	 
	public final static double DBL_MOST_LITTLE_VAL= 0.000000001, MY_INFINITY= 100;
	static boolean debug= false;

	
	Graph g= null;
	LpSolve lpSolve= null;
	Hashtable<Object,int[]> constraintHash= null;	// for synchronizing different accesses to same constraints
	public int constraintCtr= 0;
	int restrNr= 0;	// for restrictions
	int readLen= 0;
	int[] insertLen= null;
	double valObjFunc= 0d, valDeltaPos= 0d, valDeltaNeg= 0d;
	double[] result= null;	// the primal solution
	HashMap<Object,Double> trptExprHash= null;
	Profile profile= null;
	boolean flow= true;
	int nrMappingsObs= 0;
	
	public GraphLPsolver2(Graph aGraph, int readLen, int realReads) {
		this.g= aGraph;
		this.readLen= readLen;
		this.nrMappingsObs= realReads;
	}
	public GraphLPsolver2(Graph aGraph, int readLen, int[] insertMinMax, int realReads, 
			boolean considerBothStrands, boolean pairedEnd) {
		this(aGraph, readLen, realReads);
		this.costSplitWC= considerBothStrands;
		this.pairedEnd= pairedEnd;
		this.insertMinMax= insertMinMax;
	}
	
	public void setProfile(Profile profile) {
		this.profile= profile;
	}
	
	
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
		for (int i = 0; i < ids.length; i++) {
			Object oo= map.get(ids[i]);
			int[] c= getConstraintHash().get(oo);
			p.print(ids[i]+"\tC"+c[0]);
			for (int j = 1; j < c.length; j++) 
				p.print("\tC"+c[j]);
			if (oo instanceof Edge) {
				Edge e= (Edge) oo;
				int a= e.getReadNr();
				//int b= e.getPossReadNr();
				p.print("\t"+Integer.toString(a));
				//p.print("\t"+Integer.toString(b));
				//p.print("\t"+Double.toString(((double) a)/((double) b)));
				Transcript[] t= g.decodeTset(e.getTranscripts());
				for (int j = 0; j < t.length; j++) 
					p.print("\t"+t[j].getTranscriptID());				
			} else if (oo instanceof Transcript) {
				Transcript t= (Transcript) oo;
				p.print("\t"+t.getExonicLength());
			}
			p.print("\n");
		}
	}
	
	double nFactor= Double.NaN;

	public double getNFactor() {
		if (Double.isNaN(nFactor)|| true) {
			double fictReads= 0;
			Iterator iter= getTrptExprHash().keySet().iterator();
			while (iter.hasNext()) {
				Object o= iter.next();
				if (!(o instanceof String))
					continue;
				fictReads+= getTrptExprHash().get(o);
			}

			if (fictReads< nrMappingsObs)
				++nrUnderPredicted;
			else
				++nrOverPredicted;
			if (fictReads== 0^ nrMappingsObs== 0)
				System.currentTimeMillis(); // TODO can happen, when reads are where none expected
			if (fictReads> 0.000001) {	// TODO: avoid large scaling; was 0, avoid NaN; 0.5 too large
				nFactor= nrMappingsObs/ fictReads;
			} else 
				nFactor= 1d;
		}	
		return nFactor;
	}
	public double getNFactor(Transcript t) {
		
/*		long[] sig= g.encodeTset(new Transcript[] {t});
		Vector<Edge> v= g.getEdges(sig, pairedEnd);
		TProfile supa= getSuperProfileMap().get(t.getTranscriptID());
		double sum= 0d;
		for (int i = 0; i < v.size(); i++) {
			Edge e= v.elementAt(i);
			sum+= supa.getAreaFrac(e.getFrac(t, readLen), readLen, 
					TProfile.DIR_FORWARD);	// TODO i==1?TProfile.DIR_BACKWARD:
		}
*/
		
		double real= getTrptExprHash().get(t.getTranscriptID());
		
		double tNorm= getNFactor()* real;	// getNFactor()/ sum
		return tNorm;
	}
	
	/**
	 * @deprecated see getNFactor()
	 * @param trptExprHash2
	 */
	private void normalizeBack2LocusExpr(HashMap<Object, Double> trptExprHash2) {
		double fictReads= 0;
		Iterator iter= trptExprHash2.keySet().iterator();
		while (iter.hasNext()) {
			Object o= iter.next();
			if (!(o instanceof String))
				continue;
			fictReads+= trptExprHash2.get(o);
		}

		double nFactor= 1d;
		if (fictReads> 0)	// avoid NaN
			nFactor= nrMappingsObs/ fictReads;
		
		Object[] keys= trptExprHash2.keySet().toArray();
		for (int i = 0; i < keys.length; i++) {
			double fictVal= trptExprHash2.remove(keys[i]);
			double realVal= fictVal* nFactor;
			trptExprHash2.put(keys[i], realVal);
		}
		
	}

	private int[] getModelSize() {

		int cols= 0, rows= 0;
		Iterator iter= getConstraintHash().keySet().iterator();
		while (iter.hasNext()) {
			Object o= iter.next();
			cols+= ((int[]) getConstraintHash().get(o)).length;
		
			if (o instanceof Edge)
				++rows;
			else if (o instanceof Region)
				++rows;
				
		}
		++rows; //upper stabi

		int[] dim= new int[] {rows, cols};
		return dim;
	}

	private double[] createArray(int constIdx) {
		double[] a= new double[getConstraintHash().size()+1];
		for (int i = 0; i < a.length; i++) 
			a[i]=(i== constIdx)?1d:0d;
		return a;
	}

	/**
	 * create an array where certain indices are set to certain values
	 * 
	 * @param constIdx indices of constraints
	 * @param constVal values of constraints
	 * @return
	 */
	private double[] createArray(int[] constIdx, double[] constVal) {
		
		double[] a= new double[constraintCtr+ 1];
		for (int i = 0; i < a.length; i++) 
			a[i]= 0d;
		for (int i = 0; i < constIdx.length; i++) 
			a[constIdx[i]]= constVal[i];		
		return a;
	}
	
	private double[] createArray(int[] constIdx) {
		
		double[] a= new double[constraintCtr+ 1];
		for (int i = 0; i < a.length; i++) 
			a[i]= 0d;
		for (int i = 0; i < constIdx.length; i++) 
			a[constIdx[i]]= 1d;		
		return a;
	}
	
	boolean costSplitWC= false, costUseMultimap= false, pairedEnd= false;	
	int ret= 0;
	public byte costModel= COSTS_LOG, costSplit= 1;
	public float[] costBounds= null;
	boolean writeFile= false;
	File fileLPdir= null, fileLPinput= null;
	IntVector costIdx;
	DoubleVector costVal;
	int[] insertMinMax; 

	
	/**
	 * @param e
	 * @return
	 */
	int[] getConstraintIDs(Edge e) {
		
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
			allTrptIdx = new int[g.trpts.length];
			for (int i = 0; i < g.trpts.length; i++) 
				allTrptIdx[i]= getConstraintHash().get(g.trpts[i])[0];
			Arrays.sort(allTrptIdx);
		}

		return allTrptIdx;
	}
	
	/**
	 * model size must have been determined
	 * constraints[] organization:
	 * +splitW1 ... +splitWN -splitW1 ... -splitW2 (+multiW)
	 * +splitC1 ... +splitCN -splitC1 ... -splitC2 (+multiC)
	 * 
	 * offset(+split.X,-split.X)= costSplit
	 * 
	 * @param e
	 */
	private int[] reuseIdx;
	private double[] reuseVal;
	void setConstraints(Edge e) {
		
		int[] a= getConstraintHash().get(e);
		assert(a!= null);
		int[] t= getAllTrptIdx();
		int restrEdgeConstr= costSplitWC? a.length/ 2: a.length;
		if (reuseIdx== null) {
			reuseIdx= new int[restrEdgeConstr+ t.length];
			for (int i = 0; i < t.length; i++) 
				reuseIdx[reuseIdx.length- t.length+ i]= t[i];	// c_t, invariant
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
		Transcript[] tt= g.decodeTset(e.getTranscripts());
		for (int i = 0; i < nrRestr; i++, offset+= restrEdgeConstr) {

			byte dir= i== 0? Constants.DIR_FORWARD: Constants.DIR_BACKWARD;
			
			// set individual edge constraints
			System.arraycopy(a, offset, reuseIdx, 0, restrEdgeConstr);
			
			// fill in transcript constraints
			double totVal= 0d;
			for (int j = 0; j < tt.length; j++) {
				int tlen= tt[j].getExonicLength();
				UniversalMatrix m= getMatrixMap().get(tt[j].getTranscriptID());
				int[] area= e.getFrac(tt[j], readLen, dir);
				int reads= m.get(area[0], area[1], tlen, dir);
				int sum= m.getSum(dir);
				double val= reads/ (double) sum;
				totVal+= val;
				if (val< 0|| Double.isNaN(val)|| Double.isInfinite(val)) { 
					System.err.println("invalid val "+Double.toString(val)+" "+g.trpts[0]+" "+e);
				}
				int tOffset= getConstraintHash().get(tt[j])[0]- t[0];	// gotta be consecutive
				reuseVal[restrEdgeConstr+ tOffset]= val;
			}
			if (totVal== 0) {
				System.err.println("edge with 0 expectation! "+g.trpts[0]+" "+e);
				if (1== 1)
					continue; // TODO !!! we throw it out ?!! brutal...
				
			}
			
			
			// set costs, bounds and cross-links for constraints
			if (costUseMultimap) {				
				costIdx.add(reuseIdx[reuseIdx.length- 1]); 
				costVal.add(0d);	// multi
			}
			setConstraintsCostsConstant(a, restrEdgeConstr, e); 
			
			// create restrictions
			if (writeFile)
				writeRowLP(reuseIdx, reuseVal, LpSolve.EQ, i==1?e.getRevReadNr():e.getReadNr());
			else
				addConstraintToLp(reuseIdx, reuseVal, LpSolve.EQ, i==1?e.getRevReadNr():e.getReadNr());
			++restrNr;
			
		}	// restriction iterator: i
		
	}
		
	/**
	 * a[] organization:
	 * +splitW1 ... +splitWN -splitW1 ... -splitW2 (+multiW)
	 * +splitC1 ... +splitCN -splitC1 ... -splitC2 (+multiC)
	 * @param a
	 * @param baseSize size of W/C cluster
	 * @param e
	 */
	private void setConstraintsCostsConstant(int[] a, int baseSize, Edge e) {
		
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
						long obs= i== 0? e.getReadNr(): e.getRevReadNr();
						double obs1= obs== 0? 1: obs;
						double div= obs1/ costSplit;
						double maxSub= div, maxAdd= div;
						if(!Double.isNaN(costBounds[0])) {
							//maxSub= Math.max(0, (obs1/ costBounds[0])- ((costSplit-1)*div));
							maxSub= Math.max(0, (obs1* costBounds[0])- ((costSplit-1)*div));
							getLPsolve().setUpbo(a[i*baseSize+j], maxSub);	// ub of substracting values, +1 in matrix
						}
						if (!Double.isNaN(costBounds[1])) {
							maxAdd= Math.max(0, (obs1*costBounds[1])- ((costSplit-1)*div));
							;//getLPsolve().setUpbo(a[i+ j+ costSplit], maxAdd);		// ub of adding reads, -1 in matrix
						}
					} catch (Exception ex) {
						ex.printStackTrace();
					}
				
			}
		}
		
	}
	
	/**
	 * 
	 * @param dd	return value, array size is elementar block (either plus or minus)
	 * @param obs
	 * @return
	 * @deprecated debug
	 */
	private void setConstraintsCostsLogarithmical(int[] a, int baseLen, Edge e) {
						
		long obs= e.getReadNr();
		double obs1= (obs==0)?1:obs;
		double incr= obs1/ costSplit;
		
		double lastX= 0, lastY= 0d, m= 0;
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
							new double[] {-max, 1d}, LpSolve.LE, 0d);	// connect bool
					++restrNr;
					getLPsolve().setInt(a[mIdx+costOffset-1], true);
					getLPsolve().setUpbo(a[mIdx+costOffset-1], 1);
					addConstraintToLp(
							new int[]{a[mIdx+costOffset-1],a[mIdx+costOffset]}, 
							new double[] {-max, 1d}, LpSolve.LE, 0d);	// connect bool
					++restrNr;
				} catch (LpSolveException exx) {
					exx.printStackTrace();
				}
			}

				
			lastX= newX;
			lastY= newY;
			
		}	// iterate basesize (wo multi)
						
		
	}
		
	public int getReadLen() {
		return readLen;
	}

	public void setReadLen(int readLen) {
		this.readLen = readLen;
	}

	public HashMap<Object, Double> getTrptExprHash() {
		return trptExprHash;
	}

	public double getValObjFunc() {
		return valObjFunc;
	}

	public double getValDeltaNeg() {
		return valDeltaNeg;
	}

	public double getValDeltaPos() {
		return valDeltaPos;
	}

	/**
	 * @deprecated 
	 * @param e
	 * @param t
	 * @return
	 */
	private double calcEdgeError(Edge e, Transcript t) {
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

	public int[] getInsertLen() {
		return insertLen;
	}

	public void setInsertLen(int[] insertLen) {
		this.insertLen = insertLen;
	}

	HashMap<String, UniversalMatrix> matrixMap; 
	HashMap<String, UniversalMatrix> getMatrixMap() {
		if (matrixMap == null) {
			matrixMap = new HashMap<String, UniversalMatrix>();
			int lenSum= 0;
			for (int i = 0; i < g.trpts.length; i++) 
				lenSum+= g.trpts[i].getExonicLength();
			float avgLen= lenSum/ g.trpts.length;
			float rpk= g.getMappedReads()* 1000f/ avgLen;
			for (int i = 0; i < g.trpts.length; i++) {
				int tlen= g.trpts[i].getExonicLength();
				UniversalMatrix m= profile.getMatrix(tlen, rpk);
				matrixMap.put(g.trpts[i].getTranscriptID(), m);
			}
		}

		return matrixMap;
	}
	
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

	
	public double getReads(Vector<Edge> v, byte dir, long[] sig, boolean normalized) {
		return getReadsAvg(v, dir, g, sig, false, normalized);
	}
	
	public static double[] bounds2rel(int[] bounds, int len) {
		double[] dd= new double[bounds.length];
		for (int i = 0; i < dd.length; i++) 
			dd[i]= bounds[i]/ (double) len;
		return dd;
	}
	
	static int nrUnderPredicted= 0, nrOverPredicted= 0;
	private double fracs= 0;
	public double getReadsAvg(Vector<Edge> v, byte dir, Graph g, long[] sig, boolean excl, boolean normalized) {
		
		double reads= 0, fracSum= 0;
		for (int i = 0; i < v.size(); i++) {
			Edge e= v.elementAt(i);
			if (!e.isExonic())
				continue;

			long[] trpts= e.getTranscripts();
			long[] inter= Graph.intersect(trpts,sig);
			if (Graph.isNull(inter)|| (excl&& !Graph.equalSet(trpts, sig)))
				continue;	// here, and for superedges !!!

			fracs= 0d;
			if (pairedEnd) {
				for (int j = 0; v.elementAt(i).getSuperEdges()!= null&& 
						j < v.elementAt(i).getSuperEdges().size(); j++) {
					SuperEdge se= v.elementAt(i).getSuperEdges().elementAt(j);
					if (!se.isPend())
						continue;
					int cnt= 0;
					for (int k = 0; k < se.getEdges().length; k++) 
						if (se.getEdges()[k]== e)
							++cnt; 
					reads+= getReadsAvgCalc(se, dir, sig, excl, cnt, normalized);
				}
			} else 
				reads+= getReadsAvgCalc(e, dir, sig, excl, 1, normalized);
			
			fracSum+= fracs;
		}
		
		//assert(fracs== 1);	// only for transcripts
		return reads;
	}
	

	double getReadsAvgCalc(Edge e, byte dir, long[] sig, boolean excl, int cnt, boolean normalized) {
		
		long[] trpts= e.getTranscripts();
		long[] inter= Graph.intersect(trpts,sig);
		if (Graph.isNull(inter)|| (excl&& !Graph.equalSet(trpts, sig)))
			return 0d;
		Transcript[] t= g.decodeTset(inter);

		// (de?)normalize from transcript coverage
		int length= -1;

		double reads= 0; // fracs= 0;;
		for (int j = 0; j < t.length; j++) {
			
			if (t[j].getExonicLength()< readLen)
				continue;
			
			int[] bounds= e.getFrac(t[j], readLen);
			if (length== -1)
				length= bounds[1]- bounds[0]+ 1;
			else
				try{assert(length== (bounds[1]- bounds[0]+ 1));}catch(AssertionError err){
					if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
						System.err.println("[ASSERTION FAILED] bounds differ in GraphLPSolver.getReadsAvgCalc()\n\t" +
								(bounds[1]- bounds[0]+ 1)+" <> "+ length);
				}
			
			double pred= getTrptExprHash().get(t[j].getTranscriptID());
			double nf= getNFactor();

			if (normalized) {
				double cov= (cnt* pred* nf)/ t[j].getExonicLength();
				double freq= cov* length;
				reads+= freq;
			} else {
				TSuperProfile supa= null; // getSuperProfileMap().get(t[j].getTranscriptID());
				double[] relBounds= bounds2rel(bounds, t[j].getExonicLength()- readLen);
				double frac= supa.getAreaFrac(g, t[j], relBounds, readLen, insertMinMax, dir);
				fracs+= frac;
				reads+= (cnt* pred* nf* frac);
	//			if (t[j].getTranscriptID().endsWith("RA")) {
	//				System.out.print((cnt* pred* nf* frac)+"\t"+(cnt*frac)+"\t"+frac+"\t"+e);
	//				SuperEdge se= (SuperEdge) e;
	//				for (int i = 0; i < se.getEdges().length; i++) 
	//					System.out.print("\t"+ se.getEdges()[i]);
	//				System.out.println();
	//			}
	//			for (int i = 0; i < se.getEdges().length; i++) 
	//				System.out.print("\t"+ se.getEdges()[i]);
	//			System.out.println();

			//System.out.println(t[j].getTranscriptID()+" "+(pred* nf* frac)+" "+(pred* frac)+" "+pred+" "+frac);
			}
		}
		
		return reads;
	}

	public double getReadsAvg_old(Vector<Edge> v, byte dir, Graph g, long[] sigExcl) {
		double sum= 0;
		for (int i = 0; i < v.size(); i++) {
			if (sigExcl!= null&& !Graph.isNull(Graph.intersect(v.elementAt(i).getTranscripts(), sigExcl)))
				continue;

			double partsum= 0;
			double sf= g.decodeCount(v.elementAt(i).getTranscripts());
			if (dir>= 0)
				partsum+= v.elementAt(i).getReadNr();
			if (dir<= 0)
				partsum+= v.elementAt(i).getRevReadNr();
			int[]  a= getConstraintHash().get(v.elementAt(i));
			if (a== null)
				return 0;
			if (costSplitWC)
				assert(a.length% 2== 0);
			int lower= costSplitWC? (dir>= 0? 0: a.length/ 2): 0;
			int upper= costSplitWC? (dir>= 0? a.length/ 2: a.length): a.length;
			assert((upper- lower)% 2== 0);
			int half= (upper- lower)/ 2;
			for (int j = lower; j < upper; j++) 
				partsum+= (j< half? -1: 1)* result[restrNr+ a[j]];
			sum+= (partsum* getNFactor(null))/ sf;
		}
		
		if (sum< 0.0000000001)
			sum= 0;	// TODO only negatives?
		return sum;
	}

	
//	public static final int LE = 1;
//	public static final int GE = 2;
//	public static final int EQ = 3;
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
			buffy.write(COND_SYMBOLS[cond]+Double.toString(rhs)+";\n");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	BufferedWriter lpWriter;
	BufferedWriter getLPWriter() {		
		if (lpWriter == null) {
			try {
				File dir= new File("I:\\solexa\\simulation\\lp");
				fileLPinput= File.createTempFile(g.trpts[0].getTranscriptID(), ".lp", dir);
				lpWriter = new BufferedWriter(new FileWriter(fileLPinput));
			} catch (Exception e) {
				e.printStackTrace();
			}			
		}

		return lpWriter;
	}

	public boolean isFlow() {
		return flow;
	}

	public void setFlow(boolean flow) {
		this.flow = flow;
	}

	private double getMaxPlus(Edge e, Transcript t) {
		
		double x= g.getMaxFlux(readLen);
		int[] a= e.getFrac(t, readLen);
		int len= a[1]- a[0]+ 1;
		x*= len;
		x-= e.getReadNr()+ e.getRevReadNr();
		assert(x> 0);
		return x;
	}
		
	PrintStream p= null;
	public strictfp void run() {
		
//		if (costModel== COSTS_LINEAR&& costBounds== null)
//			costBounds= new int[] {100,100};
		
		long t0= System.currentTimeMillis();
		setConstraints(true);
		getLPsolve();
		constraintCtr= 0;
		debug= false;		
/*		if (g.trpts[0].getTranscriptID().startsWith("NM_001168507")||
				g.trpts[0].getTranscriptID().startsWith("NM_173453"))
			debug= true;
*/			

		// set up program
		HashMap<String, Integer> tMap= setConstraints(false);
		
	
//		for (int i = 0; i < g.trpts.length; i++) {
//			checkPercent(g.trpts[i]);
//		}
		
		String tmpOutFName= null;
		if (fileLPdir!= null) {
			tmpOutFName= fileLPdir+ File.separator
				+ g.trpts[0].getGene().getGeneID().replace(":", "_")
				+ SFX_LPOUT;
		}
		try {
			//getLPsolve().printLp();
			if (tmpOutFName== null) 
				getLPsolve().setVerbose(LpSolve.CRITICAL);	//shut up ! IMPORTANT, SEVERE, CRITICAL
			else
				getLPsolve().setOutputfile(tmpOutFName);
			
			t0= System.currentTimeMillis();
			// getLPsolve().setScaling(LpSolve.SCALE_CURTISREID);
			ret= getLPsolve().solve();
/*			 NOMEMORY (-2)  	Out of memory
			   OPTIMAL (0) 	An optimal solution was obtained
			SUBOPTIMAL (1) 	The model is sub-optimal. Only happens if there are integer variables and there is already an integer solution found. The solution is not guaranteed the most optimal one.

			 							* A timeout occured (set via set_timeout or with the -timeout option in lp_solve)
			 							* set_break_at_first was called so that the first found integer solution is found (-f option in lp_solve)
			 							* set_break_at_value was called so that when integer solution is found that is better than the specified value that it stops (-o option in lp_solve)
			 							* set_mip_gap was called (-g/-ga/-gr options in lp_solve) to specify a MIP gap
			 							* An abort function is installed (put_abortfunc) and this function returned TRUE
			 							* At some point not enough memory could not be allocated 

			INFEASIBLE (2) 	The model is infeasible
			UNBOUNDED (3) 	The model is unbounded
			DEGENERATE (4) 	The model is degenerative
			NUMFAILURE (5) 	Numerical failure encountered
			USERABORT (6) 	The abort routine returned TRUE. See put_abortfunc
			TIMEOUT (7) 	A timeout occurred. A timeout was set via set_timeout
			PRESOLVED (9) 	The model could be solved by presolve. This can only happen if presolve is active via set_presolve
			PROCFAIL (10) 	The B&B routine failed
			PROCBREAK (11) 	The B&B was stopped because of a break-at-first (see set_break_at_first) or a break-at-value (see set_break_at_value)
			FEASFOUND (12) 	A feasible B&B solution was found
			NOFEASFOUND (13) 	No feasible B&B solution found
*/
			if (ret!= 0) 
				debug= true;
			
		} catch (LpSolveException e) {
			e.printStackTrace();
		}

		// write out
		if (debug|| tmpOutFName!= null) {
			getLPsolve().printLp();
			getLPsolve().printObjective();
			getLPsolve().printSolution(1);
			
			// additional stream only afterwards
			try {
				p= new PrintStream(new FileOutputStream(tmpOutFName, true));				

			} catch (Exception e) {
				//e.printStackTrace();
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("[FATAL] failed to set lp output to:\n\t"+ tmpOutFName);
			}
			
			setConstraints(true);
		}
			
		// get transcription expression levels		
		trptExprHash= getResult(tMap);
		//normalizeBack2LocusExpr(trptExprHash);
		getLPsolve().deleteLp();	// closes file outFName
		
		if (debug) {
			System.err.flush(); // doesnt work
			p= System.out;
			for (int i = 0; i < g.trpts.length; i++) {
				; //getAllPercent(g.trpts[i]);
			}
		}
		
		// output additionally
		if (p!= null) {
			Iterator<Object> idIter= trptExprHash.keySet().iterator();
			while(idIter.hasNext()) {
				Object o= idIter.next();
				if (!(o instanceof String))
					continue;
				String id= (String) o;
				p.println(id+" "+trptExprHash.get(id));
			}
			p.println();
			p.println("Settings:");
			p.println("paired-end\t"+pairedEnd);
			if (costBounds!= null) {
				if (!Double.isNaN(costBounds[0]))
					p.print("cost boundaries:\tlower /"+costBounds[0]);
				if (!Double.isNaN(costBounds[1]))
					p.print(" upper *"+costBounds[1]);
				p.println();
			}
			p.print("costfunc\t");
			if (costModel== COSTS_LINEAR)
				p.println("linear");
			else if(costModel== COSTS_LOG)
				p.println("log");
			printConstraintHash(p);
			
			p.println();
			p.println("Transcripts:");
			for (int i = 0; i < g.trpts.length; i++) {
				p.print(g.trpts[i].getTranscriptID()+"\t");
				SpliceSite[] ss= g.trpts[i].getSpliceSitesAll();
				for (int j = 0; j < ss.length; j++) 
					p.print(ss[j].toString());
				p.println();
			}
		}
		
		// close file
		if (p!= null)
			p.flush();
		if (tmpOutFName!= null) 
			p.close();
		
//		System.err.println("solved "+g.trpts[0].getTranscriptID()+": "+g.trpts.length+" trpts, "+constraintCtr+" constr, "+restrNr+" restr"
//				+ ((System.currentTimeMillis()- t0)/ 1000)+ " sec.");
		if (debug)
			System.currentTimeMillis();
	}

	private void printConstraintHash(PrintStream p) {
		Iterator iter= getConstraintHash().keySet().iterator();
		while (iter.hasNext()) {
			Object o= iter.next();
			if (!(o instanceof Edge))
				continue;
			p.print(o);
			p.print(":\t\t");
			int[] c= getConstraintHash().get(o);
			for (int i = 0; i < c.length; i++) 
				p.print("C"+Integer.toString(c[i])+" ");
			p.println();
		}
		
		iter= getConstraintHash().keySet().iterator();
		while (iter.hasNext()) {
			Object o= iter.next();
			if (!(o instanceof Transcript))
				continue;
			p.print(o);
			p.print(":\t\t");
			int[] c= getConstraintHash().get(o);
			for (int i = 0; i < c.length; i++) 
				p.print("C"+Integer.toString(c[i])+" ");
			p.println();
		}
	}

	private HashMap<Object, Double> getResult() {
		result= new double[1+ restrNr+ constraintCtr];
		try {
			getLPsolve().getPrimalSolution(result);
		} catch (LpSolveException e1) {
			e1.printStackTrace();
		}
		valObjFunc= result[0];	
		trptExprHash= new HashMap<Object,Double>(result.length, 1f);
		
		// transcripts
 		Transcript[] trpts= g.trpts;
		for (int i = 0; i < trpts.length; i++) { 
			int[] c= getConstraintHash().get(trpts[i]);
			double x= result[restrNr+c[0]];
			if (Double.isInfinite(x)|| Double.isNaN(x))
				System.err.println("Num.error");
			trptExprHash.put(trpts[i].getTranscriptID(), x);
		}
		
		// edges
		Object[] keys= constraintHash.keySet().toArray();
		for (int i = 0; i < keys.length; i++) {
			if (!(keys[i] instanceof Edge))
				continue;
			Edge e= (Edge) keys[i];
			int[] c= constraintHash.get(e);
			trptExprHash.put(e, e.getReadNr()+e.getRevReadNr()+result[restrNr+ c[0]]- result[restrNr+ c[1]]);
		}
			
		return trptExprHash;
	}
	
	private HashMap<Object, Double> getResult(HashMap<String, Integer> tMap) {
		
		result= new double[1+ restrNr+ constraintCtr];
		try {
			getLPsolve().getPrimalSolution(result);
		} catch (LpSolveException e1) {
			e1.printStackTrace();
		}
		valObjFunc= result[0];	
		trptExprHash= new HashMap<Object,Double>();
		
		// transcripts
 		Transcript[] trpts= g.trpts;
 		double sum= 0;
		for (int i = 0; i < trpts.length; i++) { 
			int c= tMap.get(trpts[i].getTranscriptID());
			if (Double.isNaN(c))
				System.currentTimeMillis();
			double x= result[restrNr+ c];
			double tot= mapCCheck.get(trpts[i].getTranscriptID());
			tot/= 2d;
			x/= tot;
			sum+= x;
			if (Double.isNaN(x))
				System.currentTimeMillis();
			trptExprHash.put(trpts[i].getTranscriptID(), x);
		}
		
		// nfac
		double nfac= nrMappingsObs/ (double) sum;
		for (int i = 0; i < trpts.length; i++) {
			double x= trptExprHash.get(trpts[i].getTranscriptID());
			x*= nfac;
			if (Double.isNaN(x))
				System.currentTimeMillis();
			trptExprHash.put(trpts[i].getTranscriptID(), x);
		}
		
		// locus norm
		for (int i = 0; i < trpts.length; i++) {
			int tlen= trpts[i].getExonicLength();
			UniversalMatrix m= profile.getMatrix(tlen);
			double f= m.getNfactor(0.2d);
			double x= trptExprHash.get(trpts[i].getTranscriptID());
			x*= f;
			if (Double.isNaN(x))
				System.currentTimeMillis();
			trptExprHash.put(trpts[i].getTranscriptID(), x);
		}
		
		
		return trptExprHash;
	}

	void addRestrictions() {
				
		costIdx= new IntVector();
		costVal= new DoubleVector();
		Iterator iter= getConstraintHash().keySet().iterator();
		while (iter.hasNext()) {
			Object o= iter.next();
			if (!(o instanceof Edge))
				continue;	// skip transcripts

			Edge e= (Edge) o;
			assert(e.length()!= 0);	// can happen aparently, if substracting from neigbor pos
			
			setConstraints(e);
		}
		
	}
	
	void setObjectiveFunction() {
		// set objective function
		double[] a= createArray(costIdx.toIntArray(),costVal.toDoubleArray());
		try {
			getLPsolve().setObjFn(a);
			getLPsolve().setMinim();
		} catch (LpSolveException e) {
			e.printStackTrace();
		}
	}

	public static int calcReadsInTranscript(int len, int readLen) {
		int cnt= (len+readLen-1)-(2*(readLen-1));	// expected reads, normalized to transcript length
		return cnt;
	}

	public boolean isPairedEnd() {
		return pairedEnd;
	}

	public void setPairedEnd(boolean pairedEnd) {
		this.pairedEnd = pairedEnd;
	}

	public void setCostSplit(byte costSplit) {
		this.costSplit = costSplit;
	}

	/**
	 * #(Edge,Transcript) x (int[],int[][])
	 * watch out with hash function of edges
	 * @return
	 */
	public Hashtable<Object,int[]> getConstraintHash() {
		if (constraintHash == null) {
	
			constraintHash= new Hashtable<Object,int[]>();
			
			// edges
			Edge[] edges= g.getExonicEdgesInGenomicOrder();
			for (int i = 0; i < edges.length; i++) {
				if (!edges[i].isExonic())
					continue;
				// add the edge itself
				if ((!pairedEnd)&& edges[i].length()>= readLen) {
					
					constraintHash.put(edges[i], getConstraintIDs(edges[i]));
				}
				for (int j = 0; edges[i].getSuperEdges()!= null&& j < edges[i].getSuperEdges().size(); j++) {
					SuperEdge se= edges[i].getSuperEdges().elementAt(j);
					if (se.getEdges()[0]!= edges[i])
						continue;	
					assert(!constraintHash.containsKey(se)); // paired ends are iterated twice, ej maybe more
					if (se.isPend()) {
						if (!pairedEnd) 							
							continue;	// should not incur, PEs never added
						// else add
					} else {
						if (pairedEnd) {	// paired-end of EJ
							for (int k = 0; se.getSuperEdges()!= null&& k < se.getSuperEdges().size(); k++) {
								SuperEdge se2= se.getSuperEdges().elementAt(k);
								assert(se2.isPend());
								if (se2.getEdges()[0]!= se)
									continue;
								assert(!constraintHash.containsKey(se2));
								constraintHash.put(se2, getConstraintIDs(se2));
							}
							continue;
						} // else add 
					}
					
					// add EJ, ggf PE
					constraintHash.put(se, getConstraintIDs(se));
					
				}
			}
			
			// transcripts
			//System.err.println("hash "+constraintHash.size());
			Transcript[] trpts= g.trpts;
			for (int i = 0; i < trpts.length; i++) 
				constraintHash.put(trpts[i],new int[] {++constraintCtr});
	
		}
	
		return constraintHash;
	}

	public void setFileLPdir(File dir) {
		this.fileLPdir = dir;
	}
	public double getAllExpectedFracs(Graph g, HashMap<String, TSuperProfile> supaMap, long[] partition, Edge[] edges, int readLen) {
		Transcript[] t= g.decodeTset(partition);
		double val= 0d;
		for (int i = 0; i < t.length; i++) {
			TSuperProfile profile= supaMap.get(t[i].getTranscriptID());
			int[] coord= SuperEdge.getPEfrac(t[i], readLen, edges);
			val+= profile.getAreaFrac(g, t[i],
					bounds2rel(coord, t[i].getExonicLength()- readLen), 
					readLen, insertMinMax, Constants.DIR_BOTH);
		}
	
		return val;
	}
	
	
	private void getConstraints(Edge e, long[] sig, IntVector v,
			HashMap<Edge, IntVector> mapE, boolean sense, boolean count) {
		
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
			if (Graph.isNull(Graph.intersect(se.getTranscripts(), sig)))
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
				if (Graph.isNull(Graph.intersect(se2.getTranscripts(), sig)))
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
	 * #(Edge,Transcript) x (int[],int[][])
	 * watch out with hash function of edges
	 * @return
	 */
	public HashMap<String, Integer> setConstraints(boolean count) {
		
		//boolean debug= false;
//		if (g.trpts[0].getTranscriptID().equals("ENST00000323441")) {
//			debug= true;
//		}
		
		
		// transcript constraint variables
		Transcript[] trpts= g.trpts;
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
			for (int i = 0; i < trpts.length; i++) { 
				tMap.put(trpts[i].getTranscriptID(), ++constraintCtr);
				if (p!= null)
					p.println(trpts[i].getTranscriptID()+"\t"+constraintCtr);
			}
			v= new IntVector();	// indices for transcript/part 
			w= new IntVector();	// indices for cost function
			u= new IntVector();	// observation, bases for cost function
		}
		
		// iterate edges
		Edge[] edges= g.getExonicEdgesInGenomicOrder();
		if (!count) {
			mapCCheck= new HashMap<String, Double>();
			for (int i = 0; i < trpts.length; i++) 
				mapCCheck.put(trpts[i].getTranscriptID(), 0d);
		}
		for (int i = 0; i < edges.length; i++) {
			Edge e= edges[i];
			if (!e.isExonic())
				continue;
			
			// the base edge
			Transcript[] tt= g.decodeTset(e.getTranscripts());
			// sense/anti
			for (int sa= 0; sa< 2; ++sa) {
				
				HashMap<Edge, IntVector> mapE= new HashMap<Edge, IntVector>();	// BUG?
				
				for (int x = 0; x < tt.length; ++x) {
					if (!count) {
						v.removeAll();
					} else if (p!= null)
						p.print(e+"\t"+(sa==0?"sense":"asense")+"\t"+tt[x]+"\t");

					long[] sig= g.encodeTset(tt[x]);
					getConstraints(e, sig, v, mapE, sa== 0, count);
					
					// add transcript constraint
					if (!count) {
						int[] idx= new int[v.length+ 1]; // obs parts+ tx frac
						System.arraycopy(v.vector, 0, idx, 0, v.length);
						idx[idx.length- 1]= tMap.get(tt[x].getTranscriptID());
						double[] val= new double[idx.length];
						Arrays.fill(val, 1d);
						int tlen= tt[x].getExonicLength();
						UniversalMatrix m= profile.getMatrix(tlen);
						double f= m.getFrac(
									tt[x].getExonicPosition(e.getFrac(true)),
									tt[x].getExonicPosition(e.getFrac(false)),
									tlen,
									sa== 0?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
						//System.err.println(f);
						if (Double.isInfinite(f)|| Double.isNaN(f)) {
							System.err.println("infinite value");
							f= m.getFrac(
									tt[x].getExonicPosition(e.getFrac(true)),
									tt[x].getExonicPosition(e.getFrac(false)),
									tlen,
									sa== 0?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
						}
						mapCCheck.put(tt[x].getTranscriptID(),
								mapCCheck.get(tt[x].getTranscriptID())+ f);
						val[val.length- 1]= -f;
						if (debug&& !count) {
							StringBuilder sb= new StringBuilder(e.toString());
							sb.append(": ");
							for (int k = 0; k < idx.length; k++) {
								sb.append(val[k]>0?"+":"");
								sb.append(val[k]%1==0?((int) val[k]):val[k]);
								sb.append("C");
								sb.append(idx[k]+" ");
							}
							sb.append("= 0");
							System.out.println(sb);
						}
						addConstraintToLp(idx, val, LpSolve.EQ, 0);
						++restrNr;
					}
				} // iterate transcripts
				
				// add edge constraints
				Edge[] ee= new Edge[mapE.size()];
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
				for (int j = 0; j < ee.length; j++) {
					Edge f= ee[j];
					
//					if (f.toString().equals("-1611848--1611703^-1611848--1609293]PE"))
//						System.currentTimeMillis();
					boolean paird= (f instanceof SuperEdge)&& ((SuperEdge) f).isPend();
					int nr= ((paird|| sa== 0)? f.getReadNr(): f.getRevReadNr());
					v= mapE.remove(f);
					if (count)
						constraintCtr+= 2;
					else {
						int[] idx= new int[v.length+ 2];	// +/-
						System.arraycopy(v.vector, 0, idx, 0, v.length);
						int c= ++constraintCtr;
						// plus on not paired edges at 0-cost, it substracts from obs
						if (paird|| !pairedEnd) {
							w.add(c);
							u.add(nr);
						}
						idx[idx.length- 2]= c;
						// plus has to be limited, it substracts
						int lim= (paird|| !pairedEnd)? Math.max(nr- 1, 0): nr; 
						try {
							getLPsolve().setUpbo(constraintCtr, lim);
						} catch (LpSolveException e1) {
							e1.printStackTrace();
						}
							
						c= ++constraintCtr;
						// do not limit adding, even with f= 100 unsolvable systems
						// adding reads always costs, also on single edges
						w.add(c);
						u.add(nr);
						idx[idx.length- 1]= c;
						double[] val= new double[idx.length];
						Arrays.fill(val, 1d);
						val[val.length- 1]= -1d;
						if (debug&& !count) {
							StringBuilder sb= new StringBuilder(f.toString());
							sb.append(": ");
							for (int k = 0; k < idx.length; k++) {
								sb.append(val[k]==1?"+C":"-C");
								sb.append(idx[k]+" ");
							}
							sb.append("= "+nr);
							System.out.println(sb);
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
//		double[] a= createArray(w.toIntArray());	// linear costs
		double min= Math.exp(-nrMappingsObs);
		double[] costs= new double[w.size()];
		for (int i = 0; i < costs.length; i++) {
			double x= w.get(i);
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
		
		// TODO consistency check
		Object[] oo= mapCCheck.keySet().toArray();
		for (int i = 0; i < oo.length; i++) {
			double val= mapCCheck.get(oo[i]);
			//System.err.println("check "+ val);
			if (Math.abs(2d- val)> 0.2d)
				System.err.println("Fraction inconsistency "+ oo[i]+"\t"+val);
		}
		
		
		return tMap;
	}
	
}
