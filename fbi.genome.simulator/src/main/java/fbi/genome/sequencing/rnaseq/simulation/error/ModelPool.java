package fbi.genome.sequencing.rnaseq.simulation.error;

import commons.ByteArrayCharSequence;
import commons.Progressable;

import fbi.genome.io.GEMobject;
import fbi.genome.io.ThreadedBufferedByteArrayStream;
import fbi.genome.model.commons.IntVector;
import fbi.genome.model.constants.Constants;
import fbi.genome.model.constants.PrintstreamProgressable;
import fbi.genome.sequencing.rnaseq.simulation.FluxSimulatorSettings;

import java.io.*;
import java.util.Arrays;
import java.util.Random;
import java.util.Vector;
import java.util.regex.Matcher;

public class ModelPool {

	public static final byte[] QLEVEL_ILLUMINA= new byte[] {-40,40}, QLEVEL_CNT_PHRED= new byte[] {-100,100}; 
	public static byte[] qualLevels= QLEVEL_ILLUMINA;
	public static final String ERR_FIELD_TAG= "#", ERR_ID= "MODEL", QUAL_ID= "QUALITIES";
	
	Vector<ErrorModel> modelV= null;
	CrossTalkTable xTalkTable= null;
	float tQual;
	long totalCases;
	int readLength= -1;
	
	public static boolean parsingStopped= false;
	
	public static void main(String[] args) {
		
		long t0= System.currentTimeMillis();
		ModelPool myPool= null;
		
		try {
			myPool= readlFromGEM(
					new File("c:\\workspace\\Genome\\resources\\formats\\contamination_4_qualities.0.map"), 
					33, new PrintstreamProgressable(System.err));
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.err.println("found "+myPool.modelV.size()+" models.");
		//System.err.println("\t"+myPool.getTotalErrors()+" errors.");
		
		//myPool.writeModel("c:\\workspace\\Genome\\resources\\formats\\contamination_4_qualities.0.modeltest");
		
		System.err.println("[TICTAC] took "+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
	}
	
	
	double[] goodFellas;
	public void addGoodGuys(int start, int extension,
			//ByteArrayCharSequence bseq,
			ByteArrayCharSequence cs) {
		for (int i = 0; xTalkTable.hasQualities()&& i < extension; i++) {
//			incr(goodFellas, 
//					cs.byteAt(start+ i) 
//					- (ModelPool.qualLevels== ModelPool.QLEVEL_ILLUMINA?64:33)
//					- ModelPool.qualLevels[0],
//				 goodCases);	// +1-1
			byte dbg= cs.byteAt(start+ i);
			char db2= cs.charAt(start+ i);
			int p= cs.byteAt(start+ i) 
					- (ModelPool.qualLevels[1]== ModelPool.QLEVEL_ILLUMINA[1]?64:33)
					- ModelPool.qualLevels[0];
			++goodFellas[p];
		}
	}
	
	public ModelPool(int readLength, boolean withQualities) {
		this.readLength= readLength;
		modelV= new Vector<ErrorModel>();
		xTalkTable= new CrossTalkTable(withQualities);
		goodFellas= new double[ModelPool.qualLevels[1]
		            - ModelPool.qualLevels[0]+ 1];
	}
	
	public long[] getPEMstarts(int readLen) {
		long[] a= new long[readLen];
		int max= -1;
		for (int i = 0; i < modelV.size(); ++i) {
			if (!(modelV.elementAt(i) instanceof PositionErrorModel))
				continue;
			int p= ((PositionErrorModel) modelV.elementAt(i)).getStart();
			if (p>= a.length) {
				long[] b= new long[p+1];
				System.arraycopy(a, 0, b, 0, a.length);
				a= b;
			}
			++a[p];
			max= Math.max(max, p);
		}
		
/*		if((max+1)!= a.length) {
			long[] b= new long[max+1];
			System.arraycopy(a, 0, b, 0, max+1);
			a= b;
		}
*/		
		
		return a;
	}
	
	
	public static ModelPool readlFromGEM(File f, float tQual, Progressable prog) throws Exception {
		try {
			ModelPool babeAtThePool= null;
			
			ThreadedBufferedByteArrayStream buffy= new ThreadedBufferedByteArrayStream(100000,f,true);
			ByteArrayCharSequence cs= new ByteArrayCharSequence(10000);
			GEMobject go= new GEMobject();
			int[] err= new int[3];
			byte[] tmpError= null;
			long bytesRead= 0, totBytes= f.length();
			int perc= 0;
			if (prog!= null) {
				prog.setString("parsing error model from GEM alignment ");
				prog.setValue(0);
			}
			long readCtr= 0;
			for (ByteArrayCharSequence line = buffy.readLine(cs); (!parsingStopped)&& line.end!= 0; line = buffy.readLine(cs)) {
				bytesRead+= line.length()+ 1;	// fsep
				++readCtr;
				if (bytesRead*10d/ totBytes> perc) {
					if (prog!= null)
						prog.progress();
					//System.err.print("*");
					++perc;
				}
				
				if (go.reuse(line)== null)
					continue;
				if (babeAtThePool== null) {
					babeAtThePool= new ModelPool(go.getReadlength(), go.hasQualities());
					babeAtThePool.tQual= tQual;
				} else {
					if (go.getReadlength()!= babeAtThePool.readLength)
						throw new RuntimeException("Inconsistent readlength "+go.getReadlength()
								+" vs "+babeAtThePool.readLength);
				}
				
				ByteArrayCharSequence quals= null;
				if (go.hasQualities())
					quals= go.getQualities();
				if (tmpError== null|| tmpError.length< go.getReadlength())
					tmpError= new byte[go.getReadlength()];
				else for (int i = 0; i < tmpError.length; i++) 
					tmpError[i]= 0;
					
				ByteArrayCharSequence ms= go.getMatch(-1);
				Matcher m= null;
				while (true) {	// mismatches 
					m= go.getMismatch(ms, m, err);
					if (m== null)
						break;
					else {
						int q= Integer.MIN_VALUE;
						if (quals!= null)
							q= getQual(quals.byteAt(err[0]- 1));
						babeAtThePool.addError(err, null, q);	// for xtalk
						tmpError[err[0]-1]= 1;
					}
				}

				// add problems
				for (int i = 0; i < tmpError.length; i++) {
					if (tmpError[i]!=1)
						continue;
					int j= 0;
					for (j = i+1; j < tmpError.length&& 
						(tmpError[j]!=0||
						(quals!= null&& getQual(quals.byteAt(j))< tQual)); j++);
					
					babeAtThePool.add(i, j-1, quals); 
					i= j;
				}
				
				// add good guys
				for (int i = 0; i < tmpError.length; i++) {
					if (tmpError[i]!=0|| (quals!= null&& getQual(quals.byteAt(i))< tQual))
						continue;
					int j= 0;
					for (j = i+1; j < tmpError.length&& 
						tmpError[j]!=1&&
						(quals!= null&& getQual(quals.byteAt(j))>= tQual); j++);
					
					babeAtThePool.addGoodGuys(i, j- i, quals); 
					i= j;
				}
//				if (!Graph.isNull(tmpError)) {
//					int cnt= 0, last= -1;
//					for (int i = 0; (cnt< quals.length())&& (i < tmpError.length); ++i) {
//						long tmp= tmpError[i];
//						for (int j = 0; (cnt< quals.length()) 
//								&&(j < 64); ++j, ++cnt) {
//							if ((tmp& 1)== 1) 
//								last= (last< 0)?cnt:last;
//							else { 
//								if (last>= 0) {
//									add(last, cnt-1, quals);
//									last= -1;
//								}
//							}
//							tmp= tmp>>1;
//						}
//					}
//				}
				
			}
			
			if (parsingStopped) {
				buffy.setStop(true);				
				if (prog!= null)
					prog.setString("stopped.");
				parsingStopped= false;
				return null;
			}
			
			babeAtThePool.init(readCtr);
			
			if (prog!= null)
				prog.setString("parsed "+readCtr+" reads.");
			return babeAtThePool;
			
		} catch (Exception e) {
			if (prog!= null)
				prog.setString("failed.");
			e.printStackTrace();
		}
		
		return null;				
	}
	
	public void write(File f, Progressable prog) {
		
		try {
			if (prog!= null) {
				prog.setString("Writing model ");
				prog.setMinimum(0);
				prog.setMaximum(10);
				prog.setValue(0);
			}
			BufferedWriter buffy= new BufferedWriter(new FileWriter(f));
			buffy.write(this.toString());
			if (prog!= null)
				prog.setValue(10);
			buffy.flush();
			buffy.close();
		} catch (Exception e) {
			if (prog!= null)
				prog.setString("Error writing model");
			e.printStackTrace();
		}
	}
	
	public static ModelPool read(File f, Progressable prog, FluxSimulatorSettings settings) {
		
		try {
			if (prog!= null) {
				prog.setString("Reading error model ");
				prog.setMinimum(0);
				prog.setMaximum(10);
				prog.setValue(0);
			}
			
			BufferedReader buffy= new BufferedReader(new FileReader(f));
			ModelPool babe= null;
			for (String s= buffy.readLine(); s!= null;) {
				if (s.length()== 0) {
					s= buffy.readLine();
					continue;
				}
				s= s.trim();
				if (s.startsWith(ERR_FIELD_TAG)) {
					if (s.charAt(ERR_FIELD_TAG.length())!= ' ') 
						s= ERR_FIELD_TAG+ " "+ s.substring(ERR_FIELD_TAG.length());
					String[] ss= s.split("\\s");
					ss[1]= ss[1].toUpperCase();
					if (ss[1].equals(ERR_ID)) {
						boolean qual= false;
						if (ss.length< 4)
							throw new RuntimeException(ERR_ID+" must contain 4 tokens: # "+ERR_ID+" readlength instances");
						if (ss.length> 4)
							qual= true;
						int len= Integer.parseInt(ss[2]);
						if (len!= settings.getReadLength()) {
							settings.setReadLength(len);
						}
						babe= new ModelPool(len, qual);
						babe.totalCases= Long.parseLong(ss[3]);
						
						if (qual) {
							qualLevels= new byte[] {Byte.parseByte(ss[4]), Byte.parseByte(ss[5])};
							babe.tQual= Float.parseFloat(ss[6]);
							
							s= buffy.readLine();
							ss= s.split("\t");
							for (int i = 0; i < ss.length; i++) 
								babe.goodFellas[i]= Double.parseDouble(ss[i]);
						}
						s= buffy.readLine();
						
					} else if (babe== null) {
						throw new RuntimeException(ERR_ID+" has to be first in file!");

					// Xtalk table
					} else if (ss[1].equals(CrossTalkTable.ERR_ID)) {
						ss[2]= ss[2].toUpperCase();
						int p1= Arrays.binarySearch(CrossTalkTable.SYMBOLS, ss[2].charAt(0));
						while ((s= buffy.readLine())!= null&& (!s.startsWith(ERR_FIELD_TAG))) {
							if (s.length()== 0)
								continue;
							s= s.trim();
							ss= s.split("\\s+");
							int p2= 0, x= 0;
							if (babe.xTalkTable.hasQualities()) {
								p2= Integer.parseInt(ss[0]);
								x= 1;
							}
							double sumP= 0d, p= 0d;
							for (int i = x; i < ss.length; i++) {
								p= Double.parseDouble(ss[i]);
								sumP+= p;
								if (babe.xTalkTable.hasQualities())
									((double[][][])babe.xTalkTable.table)[p1][p2- ModelPool.qualLevels[0]][i-x]= p;
								else
									((double[][])babe.xTalkTable.table)[p1][i-x]= p;
							}
							int last= (int) Math.round(p);
							int sum= (int) Math.round(sumP);
							// check for probabilities, 
							// border case 0,0,..,1 can be both, 
							// but is intrinsically CDF
							if (sum== 1&& last!= 1) {	
								// convert table to cdf
								sumP= 0d; p= 0d;
								for (int i = x; i < ss.length; i++) {
									if (babe.xTalkTable.hasQualities()) {
										p= ((double[][][])babe.xTalkTable.table)[p1][p2- ModelPool.qualLevels[0]][i-x];
										p+= sumP;
										((double[][][])babe.xTalkTable.table)[p1][p2- ModelPool.qualLevels[0]][i-x]= p;
										sumP= p;
									} else {
										p= ((double[][])babe.xTalkTable.table)[p1][i-x];
										p+= sumP;
										((double[][])babe.xTalkTable.table)[p1][i-x]= p;
										sumP= p;
									}
								}
								assert(Math.round(sumP)== 1);
							}
						}
						
					} else if (ss[1].equals(PositionErrorModel.ERR_ID)) {
						int start= Integer.parseInt(ss[2]), ext= Integer.parseInt(ss[3]);
						if (babe!= null) {
							int len= babe.readLength;
							if (start+ ext> len)
								throw new RuntimeException("Positional error out of range: "+start+ " + "+ext +" > "+ len);
						}
						PositionErrorModel pem= new PositionErrorModel(
								babe.getXTalkTable().hasQualities(),
								start, ext);
						pem.setBaseProbability(Double.parseDouble(ss[4]));
						while ((s= buffy.readLine())!= null&& (!s.startsWith(ERR_FIELD_TAG))) {
							if (s.length()== 0)
								continue;
							ss= s.split("\\s");
							int pos= Integer.parseInt(ss[0]);
							for (int i = 1; i < ss.length; i++) 
								pem.quals[pos- start][i-1]= Double.parseDouble(ss[i]);							
						}
						babe.modelV.add(pem);
					}
				} else
					s= buffy.readLine();
			}
			
			if (prog!= null) {
				prog.finish();
			}
			buffy.close();
			return babe;
			
		} catch (Exception e) {
			if (prog!= null)
				prog.setString("Error reading model");
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				e.printStackTrace();
		}
		return null;
	}
	
	void add(int first, int last, ByteArrayCharSequence quals) {
		int j, delta= last- first+ 1;
		for (j = 0; j < modelV.size(); j++) {
			if (modelV.elementAt(j) instanceof PositionErrorModel) { 
				if (((PositionErrorModel) modelV.elementAt(j)).add(
						first, delta, quals))
					break;
			}
		}
		
		if (j== modelV.size()) {
			PositionErrorModel pm= new PositionErrorModel(quals!= null, first, delta);
			pm.add(first, delta, quals);
			modelV.add(pm);
		}
	}

	void addError(int[] err, long[] tmpError, int qual) {
//		if (tmpError== null|| err[0]>= tmpError.length*64)
//			tmpError= extend(tmpError, err[0]+1);
//		int p= err[0]/ 64;
//		int b= err[0]% 64;
//		tmpError[p]|= (1<<(b-1)); 
		xTalkTable.add((char) err[1],(char) err[2], qual);
	}

	public static long[] extend(long[] a, int newLen) {
		long[] nu= new long[newLen];
		for (int i = 0; i < a.length; i++) 
			nu[i]= a[i];
		return nu;
	}

	public static void initZero(long[] a) {
		for (int i = 0; a!= null&& i < a.length; i++) 
			a[i]= 0l;
	}
	
	@Override
	public String toString() {
		StringBuilder sb= new StringBuilder(ERR_FIELD_TAG+" "+ERR_ID+" "
				+ Long.toString(totalCases)+ " "
				+ Integer.toString(readLength));
				
		if (hasQualities())
			sb.append(" "+ Integer.toString(qualLevels[0])+" "+Integer.toString(qualLevels[1])
				+" "+ Float.toString(tQual));
		
		sb.append(" "+ "\n");

		for (int i = 0; i < goodFellas.length; i++) {
			sb.append(Double.toString(goodFellas[i]));
			sb.append("\t");
		}
		sb.deleteCharAt(sb.length()- 1);
		sb.append("\n");
		
		sb.append("\n");
		sb.append(xTalkTable.toString());

		sb.append("\n");
		for (int i = 0; i < modelV.size(); i++) {
			sb.append(modelV.elementAt(i).toString());
			sb.append("\n");
		}
		
		return sb.toString();
	}


	public long[] getPEMends(int readLen) {
		
		long[] a= new long[readLen];
		int max= -1;
		for (int i = 0; i < modelV.size(); ++i) {
			if (!(modelV.elementAt(i) instanceof PositionErrorModel))
				continue;
			int p= ((PositionErrorModel) modelV.elementAt(i)).getStart()
					+ ((PositionErrorModel) modelV.elementAt(i)).extension- 1;
			if (p>= a.length) {
				long[] b= new long[p+1];
				System.arraycopy(a, 0, b, 0, a.length);
				a= b;
			}
			++a[p];
			max= Math.max(max, p);
		}
		
/*		if((max+1)!= a.length) {
			long[] b= new long[max+1];
			System.arraycopy(a, 0, b, 0, max+1);
			a= b;
		}
*/		
		
		return a;
	}
	
	public long[] getPEMoverlap(int readLen) {
		
		long[] a= new long[readLen];
		for (int i = 0; i < a.length; i++) 
			a[i]= 0;
		int max= -1;
		for (int i = 0; i < modelV.size(); ++i) {
			if (!(modelV.elementAt(i) instanceof PositionErrorModel))
				continue;
			int p= ((PositionErrorModel) modelV.elementAt(i)).getStart();
			int q= ((PositionErrorModel) modelV.elementAt(i)).extension;			
			if (q>= a.length) {
				long[] b= new long[p+1];
				System.arraycopy(a, 0, b, 0, a.length);
				a= b;
			}
			for (int j = 0; j < q; j++) 
				++a[p+ j];
			max= Math.max(max, p+ q- 1);
		}
		
/*		if((max+1)!= a.length) {
			long[] b= new long[max+1];
			System.arraycopy(a, 0, b, 0, max+1);
			a= b;
		}
*/		
		
		return a;
	}


	public CrossTalkTable getXTalkTable() {
		return xTalkTable;
	}
	
	public int[] getLengthDistr() {
		int[] len= new int[getMaxLen()];
		for (int i = 0; i < modelV.size(); i++) {
			if (!(modelV.elementAt(i) instanceof PositionErrorModel))
				continue;
			++len[((PositionErrorModel) modelV.elementAt(i)).extension-1];
		}
		return len;
	}
	
	public boolean hasQualities() {
		if (xTalkTable== null)
			return false;
		return xTalkTable.hasQualities();
	}
	
	public void init(long totalReads) {

		totalCases= totalReads;
		
		freq2cumuP(goodFellas);

		xTalkTable.fill();
		if (xTalkTable.hasQualities()) {
			double[][][] a= (double[][][]) xTalkTable.getTable();
			for (int i = 0; i < a.length; i++) 
				for (int j = 0; j < a[i].length; j++) 
					freq2cumuP(a[i][j]);
		} else {
			double[][] a= (double[][]) xTalkTable.getTable();
			for (int i = 0; i < a.length; i++) 
				freq2cumuP(a[i]);
		}
		
		
		for (int i = 0; i < modelV.size(); i++) {
			if (!(modelV.elementAt(i) instanceof PositionErrorModel))
				continue;
			PositionErrorModel pem= (PositionErrorModel) modelV.elementAt(i);
			for (int j = 0; pem.quals!= null&& j < pem.quals.length; j++) 
				freq2cumuP(pem.quals[j]);
			pem.setBaseProbability(
					pem.getBaseProbability()/ (double) totalReads);
		}
	}
	
	public static void freq2cumuP(double[] a) {
		
		if (a== null)
			return;
		
		double sum= 0;
		for (int i = 0; i < a.length; i++) 
			sum+= a[i];
		if (sum== 0)
			for (int i = 0; i < a.length; i++) {
				a[i]= (1d/ a.length)+ ((i> 0)?a[i- 1]:0);
			}
		else
			for (int i = 0; i < a.length; i++) {
				a[i]= (a[i]/ sum)+ ((i> 0)?a[i- 1]:0);
			}
		assert(a[a.length- 1]== 1d);
	}

	Random rndModelChooser= new Random(), rndMutator= new Random();
	IntVector idxModelV= new IntVector(2,2);
	public void apply(char[] seq, byte[] quals) {
				
		// check model(s)
		idxModelV.reset(); 
		for (int i = 0; i < modelV.size(); i++) {
			double r= rndModelChooser.nextDouble();
			if (r< modelV.elementAt(i).getBaseProbability()) {
				idxModelV.add(i);
			}
		}
		
		// 1. fill with good quals
		for (int i = 0; hasQualities()&& i < quals.length; i++) {
			double r= rndMutator.nextDouble();
			int p= Arrays.binarySearch(goodFellas, r);
			p= (p<0)?-p-1:p;
			quals[i]= (byte) (p+ qualLevels[0]);
		}
		
		// 2. downgrade
		for (int i = 0; hasQualities()&& i < idxModelV.size(); i++) {
			modelV.elementAt(i).apply(quals);
		}
		
		// 3. mutate seq
		if (hasQualities()) {
			for (int i = 0; i < seq.length; i++) {
				double pe= Math.pow(10, quals[i]/ -10d);
				double r= rndMutator.nextDouble();
				if (r< pe) {
					seq[i]= xTalkTable.mutate(seq[i], quals[i]);
				}
			}
		} else {
			for (int i = 0; i < idxModelV.size(); i++) {
				PositionErrorModel model= (PositionErrorModel) modelV.get(idxModelV.get(i));
				for (int j = 0; j < model.extension; j++) {
					seq[j+ model.start]= xTalkTable.mutate(seq[j], Byte.MIN_VALUE); 
				}
			}
		}
	}

	public int getMaxLen() {
		int len= 0;
		for (int i = 0; i < modelV.size(); i++) {
			if (!(modelV.elementAt(i) instanceof PositionErrorModel))
				continue;
			len= Math.max(((PositionErrorModel) modelV.elementAt(i)).extension, len);
		}
		return len;
	}


	public Vector<ErrorModel> getModelV() {
		return modelV;
	}

	public double[] getGoodFellas() {
		return goodFellas;
	}

	public static int getQual(byte ascii) {
		return (ascii
				- (ModelPool.qualLevels[1]== ModelPool.QLEVEL_ILLUMINA[1]?64:33));
				//- (byte) ModelPool.qualLevels[0]
	}
	public static byte getAscii(int val) {
		return (byte) (val
				+ (ModelPool.qualLevels[1]== ModelPool.QLEVEL_ILLUMINA[1]?64:33));
	}
	
	public static void incr(double[] ds, int pos, long cases) {
			double last= 0;
			for (int i = 0; i < ds.length; ++i) {
				double curr= (ds[i]- last)* cases;
				if (i== pos)
					++curr;
	//			double chk= (last/ (cases+ 1));
				last= ds[i];
				ds[i]= (curr/ (cases+ 1))+ (i>0?ds[i-1]:0);
				if (i>0&& ds[i]< ds[i-1])
					System.currentTimeMillis();
			}
			assert(ds[ds.length-1]== 1d);
		}

	public long getTotalCases() {
		return totalCases;
	}

	public int getReadLength() {
		return readLength;
	}

	public void apply(ByteArrayCharSequence cs, int seqStart) {
				
		// check model(s) 
		idxModelV.reset(); 
		for (int i = 0; i < modelV.size(); i++) {
			double r= rndModelChooser.nextDouble();
			if (r< modelV.elementAt(i).getBaseProbability()) {
				idxModelV.add(i);
			}
		}

		// generate qualities (FASTQ)
		int seqEnd= cs.end, len= seqEnd- seqStart;
		if (hasQualities()) {

			cs.append("\n+\n");
			seqEnd+= 3;
			cs.ensureLength(cs.end, len);	// for qualities
			byte[] a= cs.a;
			
			// 1. fill with good quals
			while (cs.end < seqEnd+ len) {
				double r= rndMutator.nextDouble();
				int p= Arrays.binarySearch(goodFellas, r);
				p= (p<0)?-p-1:p;
				a[cs.end++]= getAscii(p+ qualLevels[0]);
			}

			// 2. downgrade
			for (int i = 0; i < idxModelV.size(); i++) {
				modelV.elementAt(idxModelV.get(i)).apply(a, seqEnd, cs.end);
			}
			
			// 3. mutate seq
			for (int i = 0; i < len; ++i) {
				double pe= Math.pow(10, a[seqEnd+ i]/ -10d);
				double r= rndMutator.nextDouble();
				if (r< pe) {
					a[seqStart+ i]= 
						(byte) xTalkTable.mutate((char) a[seqStart+ i], a[seqEnd+ i]);
				}
			}
		} else {
			//  just mutate sequence
			byte[] a= cs.a;
			for (int i = 0; i < idxModelV.size(); i++) {
				PositionErrorModel model= (PositionErrorModel) modelV.get(idxModelV.get(i));
				int modelStart= model.start;
				for (int j = 0; j < model.extension; j++) {
					a[seqStart+ modelStart+ j]= 
						(byte) xTalkTable.mutate((char) a[seqStart+ modelStart+ j], Byte.MIN_VALUE);	// ? min 
				}
			}
		}
	}

}