package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.StringUtils;
import fbi.commons.file.FileHelper;
import fbi.commons.io.IOHandler;
import fbi.commons.io.IOHandlerFactory;
import fbi.commons.thread.StoppableRunnable;
import fbi.commons.thread.SyncIOHandler2;
import fbi.commons.thread.ThreadedQWriter;
import fbi.genome.io.BufferedByteArrayReader;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.Gene;
import fbi.genome.model.Graph;
import fbi.genome.model.Transcript;
import fbi.genome.model.commons.DoubleVector;
import fbi.genome.model.constants.Constants;

import genome.sequencing.rnaseq.simulation.distributions.AbstractDistribution;
import genome.sequencing.rnaseq.simulation.distributions.EmpiricalDistribution;
import genome.sequencing.rnaseq.simulation.distributions.NormalDistribution;

import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.special.Gamma;

import java.io.*;
import java.util.*;

public class Fragmenter implements StoppableRunnable {
	public static final String FRG_FILE_TAB= Constants.TAB;
	public static final byte  MODE_NOT_INITED= -1, MODE_WRITE_INITIAL= 0, MODE_FRAG_EZ= 1, MODE_NEBU= 2, MODE_FRAG= 3, MODE_RT= 4, MODE_FILT_REJ= 5, MODE_FILT_ACC= 6, MODE_FILT_MH= 7;
	public static final double DBL_LOW_ONE= 0.999d; 	//1d- Double.MIN_VALUE;
	private static final String PROCESSOR_TAG= "Processor-";
	private static int processorID= 0;
	private static char CHAR_TAB= '\t';
	private static byte BYTE_NL= (byte) '\n';
	
	public static double[] toCDF(double[] b) {
		
		double sum= 0f;
		for (int i = 0; i < b.length; ++i) 
			sum+= (double) b[i];
		double cumu= 0d;
		for (int i = 0; i < b.length; ++i) {
			cumu+= b[i]/ (double) sum;
			b[i]= (float) cumu;
		}
		b[b.length- 1]= 1f;	// correct for rounding errors
		
		return b;
	}

	
	class Processor extends Thread {
		
		ByteArrayCharSequence cs= null;
		boolean end= false;
		long interval= 10; 	// default 10ms
		Object lock= null;
		byte mode= MODE_NOT_INITED;
		int origLen;
		CharSequence locID;
		FileOutputStream fos;
		
		/**
		 * Lengths of fragments obtained by tokenizing the molecule.
		 */
		int[] fragments= null;
		
		/**
		 * Number of fragments in <code>fragments[]</code>.
		 */
		int fragmentNb= 0;
		
		private Processor() {
			super(PROCESSOR_TAG+ processorID++);
		}
		
		Processor(byte mode, boolean threaded) {
			this();
			this.mode= mode;
			if (mode== MODE_WRITE_INITIAL)
				this.cs= new ByteArrayCharSequence("1\t1\tname");
			if (threaded)
				lock= new Object();
		}
		
		@Override
		public void run() {
			if (lock!= null&& isAlive()) {
				// sync all, else deadlocks
				synchronized (lock) {
					while((!end)&& !Fragmenter.this.stop) {
						if (((mode== MODE_WRITE_INITIAL&& locID== null)
								|| (mode!= MODE_WRITE_INITIAL&& cs== null))
								&& (!Fragmenter.this.stop)&& (!end)) 
							try {
								lock.wait();
							} catch (InterruptedException e) {
								; // Exceptions slow
							}
						if ((mode== MODE_WRITE_INITIAL&& locID!= null)|| (mode!= MODE_WRITE_INITIAL&& cs!= null))
							process();
					}
				}
			} else
				process();
		}
		
		public boolean ready() {
			if (mode== MODE_WRITE_INITIAL)
				return (locID== null);
			return (cs== null);
		}
		
		void process() {
			
			if (mode== MODE_WRITE_INITIAL) {
				processInitial();
				rw.writeLine(cs, fos);
				locID= null;
				return;
			}
			
			int start= cs.getTokenInt(0);
			int end= cs.getTokenInt(1);
			assert(start<= end);
			int len= end- start+ 1;
//			if (processFragDiss(len))
//				return;
			if (len<= settings.getFragNBlambda())
				++tgtMols;

			ByteArrayCharSequence id= cs.getToken(2);
			addFragCount(id, 1l);
			tstCurrMols+= 1;

			if (mode== Fragmenter.MODE_FILT_REJ) 
				processFilterRejection(start, end, len, id, filterDist, true);
			else if (mode== Fragmenter.MODE_FILT_ACC) 
				processFilterRejection(start, end, len, id, filterDist, false);
			else if (mode== Fragmenter.MODE_FILT_MH) 
				processFilterMCMC(start, end, len, id, originalDist, filterDist);
			else if (mode== Fragmenter.MODE_NEBU) 
				processNebu(true, start, end, len, id);
			else if (mode== Fragmenter.MODE_FRAG) 
				processFrag(false, start, end, len, id);
				//processFragUniform(false, start, end, len, id);
			else if (mode== Fragmenter.MODE_RT) {
				//rw.writeLine(cs, 1);
				processRT(start, end, len, id);
			}
				
//			if (Fragmenter.this.isStop())
//				close();
			// cleanup
			cs= null;
		}
	
		
		void processFilter_new(int start, int end, int len, ByteArrayCharSequence id) {
			
			double p= filtFac* Math.exp(-Math.pow(len- filtMu, 2)/ (2d* filtSigSquare));
			p*= 10;
			if (p> 0) {
				double r= rndGel.nextDouble();
				if (r<= p) {
					rw.writeLine(cs, fos);
					if (Fragmenter.this.plotter!= null) 
						Fragmenter.this.plotter.plot(start, end, -1, id);
				} else
					--newMols;
			} else
				--newMols;
				

		}
		
		void processNebu_v10(int start, int end, int len, ByteArrayCharSequence id) {
			int bp= (int) nextGaussianDouble(rndNebu, start+1, end);	//	Math.ceil() f* slow 
			int len2= Math.min(bp- start, end- bp+ 1);	// shorter
			int len1=(bp-start==len2)?end-bp+1:bp-start;	// longer
//			double val= Math.exp((-1d)* Math.pow(len2- settings.lambda, 2)/ 
//					Math.pow(settings.sigma* settings.lambda, 2));	// either len or len2
//			double val= Math.exp((-1d)* Math.pow(Math.max(len- settings.lambda,0), 2)/ 
//					Math.pow(settings.lambda, 2));	// TODO either len or len2
//			double val= Math.exp((-1d)* (len/ settings.lambda));	
//			double val= Math.exp((-1d)* Math.pow(Math.max(0,len- settings.getLambda()),2));	// gives only 0 or 1 :s
			
			double val= Math.exp((-1d)* Math.pow(len2/ settings.fragNBlambda,2));	// TODO either len or len2
			val= 1- val;	// 0 f. len -> settings.lamda
			double r= rndBreak.nextDouble();
				
			if (r> val) {	// does not break (!)
				++lenNb;
				increaseFragCount(id);
				
				// plot
				if (Fragmenter.this.plotter!= null)
					Fragmenter.this.plotter.plot(start, end, -1, id);
				if (len>= Fragmenter.this.settings.getFiltMin()&&
						len<= Fragmenter.this.settings.getFiltMax())
					++Fragmenter.this.tgtMols;
				
				// write
				rw.writeLine(cs, fos);
				
			} else {
				++newMols;
				lenSum+= len1;
				lenSum+= len2;
				lenNb+= 2;
				maxLen= Math.max(maxLen,len1);
				maxLen= Math.max(maxLen,len2);
				minLen= Math.min(minLen,len1);
				minLen= Math.min(minLen,len2);
				if (len1>= Fragmenter.this.settings.getFiltMin()&&
						len1<= Fragmenter.this.settings.getFiltMax())
					++Fragmenter.this.tgtMols;
				if (len2>= Fragmenter.this.settings.getFiltMin()&&
						len2<= Fragmenter.this.settings.getFiltMax())
					++Fragmenter.this.tgtMols;

				// plot
				if (Fragmenter.this.plotter!= null) {
					Fragmenter.this.plotter.plot(start, bp-1, -1, id);
					Fragmenter.this.plotter.plot(bp, end, -1, id);
				}
				
				// write
				ByteArrayCharSequence cs2= cs.cloneCurrentSeq();
				assert(bp-1>= start);
				cs2.replace(1, bp- 1);
				rw.writeLine(cs2, fos);	// getNebuLine(start,bp-1,id)+"\n"
				assert(bp<= end);
				cs.replace(0, bp);	// no longer needed
				rw.writeLine(cs, fos);	// getNebuLine(bp,end,id)+"\n"
				// id is invalid now !
			}
				
			
		}
		
		void processFrag_old(int start, int end, int len, ByteArrayCharSequence id) {
            double k= 3.4;	// 4
            
            // weibull random sampling ?
            // mean= la* (ln(2))^(1/k)
            // double la= (len/2d)/ Math.pow(Math.log(2), (1d/k)); 
            // start+ Math.round(sampleWeibull(rndFrag, la,k))
            
            int bp= -1; 
            while (bp> end|| bp< start)
                bp=start+ rndBP.nextInt(end-start+1);	

            // quatsch, f. hydrolyse nur abh. von Gesamtlaenge
            int len2= Math.min(bp- start, end- bp+ 1);
            int len1=(bp-start==len2)?end-bp+1:bp-start;
            // double val= Math.exp((-1d)* Math.pow(Math.max(0,len- settings.lambda),2));    // TODO either len or len2
            double val= 1d;
            if (len2> settings.getFragNBlambda())	// 
                val= Math.pow(len2- settings.getFragNBlambda(), -2);    // 0 f. len -> settings.lamda 
            val= 1-val;
            double r= rndBreak.nextDouble();
                
            if (r> val) {    // does not break (!)
            	lenSum+= len;                            	
            	++lenNb;
				maxLen= Math.max(maxLen,len);
				minLen= Math.min(minLen,len);
				rw.writeLine(cs, fos);
				increaseFragCount(id);
				if (Fragmenter.this.plotter!= null)
					Fragmenter.this.plotter.plot(start, end, -1, id);
				if (len>= Fragmenter.this.settings.getFiltMin()&&
						len<= Fragmenter.this.settings.getFiltMax())
					++Fragmenter.this.tgtMols;
			} else {
				++newMols;
				//lengthV.add(end- bp+ 1);
				//lengthV.add(bp- start);
				lenSum+= len1;
				lenSum+= len2;
				lenNb+= 2;
				maxLen= Math.max(maxLen,len1);
				maxLen= Math.max(maxLen,len2);
				minLen= Math.min(minLen,len1);
				minLen= Math.min(minLen,len2);

				// plot
				if (len1>= Fragmenter.this.settings.getFiltMin()&&
						len1<= Fragmenter.this.settings.getFiltMax())
					++Fragmenter.this.tgtMols;
				if (len2>= Fragmenter.this.settings.getFiltMin()&&
						len2<= Fragmenter.this.settings.getFiltMax())
					++Fragmenter.this.tgtMols;
				//setFragCount(token[2], 2l);
				if (Fragmenter.this.plotter!= null) {
					Fragmenter.this.plotter.plot(start, bp-1, -1, id);
					Fragmenter.this.plotter.plot(bp, end, -1, id);
				}

				// write
//					writer.write(getNebuLine(start,bp-1,id)+"\n");
//					writer.write(getNebuLine(bp,end,id)+"\n");
				ByteArrayCharSequence cs2= cs.cloneCurrentSeq();
				cs2.replace(1, bp- 1);
				//writer.write(cs2.a, cs2.start, cs2.end);	// getNebuLine(start,bp-1,id)+"\n"
				rw.writeLine(cs2, fos);
				cs.replace(0, bp);	// no longer needed
				//writer.write(cs.a, cs.start, cs.end);	// getNebuLine(bp,end,id)+"\n"
				rw.writeLine(cs, fos);
				// id is invalid now
			}						
			
			// so nicht
//			IntVector vec= new IntVector();
//			int bp= start;
//			while (bp< end) {
//				int i= (int) Math.round(rdi.nextUniform(settings.getFiltMin(), settings.getFiltMax()));
//				vec.addElement(i);
//				bp+= i;
//			}
//			bp-= vec.elementAt(vec.size()-1);
//			vec.removeElement(vec.size()-1);
//			
//			try {
//				if (vec.size()== 0) {
//					writer.write(s+"\n");
//					if (Fragmenter.this.plotter!= null)
//						Fragmenter.this.plotter.plot(start, end, token[2]);
//					if (len>= Fragmenter.this.settings.filtMin&&
//							len<= Fragmenter.this.settings.filtMax)
//						++Fragmenter.this.tgtMols;
//				} else {
//					--newMols;
//					int offs= (int) Math.round((len- bp)/ 2d);
//					for (int i = 0; i < vec.size(); offs+= vec.elementAt(i), i++) {
//						++newMols;
//						int newEnd= offs+vec.elementAt(i)-1;
//						writer.write(getNebuLine(offs,newEnd,token[2])+"\n");
//						if (Fragmenter.this.plotter!= null) 
//							Fragmenter.this.plotter.plot(offs, newEnd, token[2]);
//					}
//					Fragmenter.this.tgtMols+= vec.size();
//				}
//			} catch (Exception e) {
//				e.printStackTrace();
//			}
			
		
		}
		
		private int[] where= new int[10], starts= new int[10];
		void processRT(int start, int end, int len, ByteArrayCharSequence id) {

			if (len< 6)
				processFragNot(start, end, len, id);
			
			CharSequence s= mapTxSeq.get(id);
			s= s.subSequence(start, end);	// TODO why end +1 ??
			double[] fwCDF= getRTCDF(s, true), rvCDF= getRTCDF(s, false);
			
			// choose new 3' end
			int howmany= 1+ (int) Math.round(rndHowMany.nextPoisson(len/ 100d));	// TODO Poisson..;
			howmany= 100;
			if (howmany== 0)
				return;
			int new3Prime= -1;
			for (int i = 0; i < howmany; i++) {
				int p= Arrays.binarySearch(rvCDF, rtRndWhere.nextDouble());
				p= (p>= 0? p: -(p+1)); 
				new3Prime= Math.max(start+ p, new3Prime);
			}
			
			// choose new 5' end
//			howmany= (int) Math.round(rndHowMany.nextPoisson(((new3Prime- start)+1)/ 50d));	// TODO Poisson..;
//			if (howmany== 0)
//				return;
			int new5Prime= new3Prime- 1;
			for (int i = 0; i < howmany; i++) {
				double r= rtRndWhere.nextDouble();
				r= 1- r;
				int p= Arrays.binarySearch(fwCDF, r);
				p= (p>= 0? p: -(p+1)); 
				new5Prime= Math.min(start+ p, new5Prime);
			}
			
			int nuLen= new3Prime- new5Prime+ 1;
			updateMedian(nuLen);
			cs.replace(0, new5Prime);
			cs.replace(1, new3Prime);
			cumuLen+= nuLen;
			incLinesWrote();
			++newMols;
			rw.writeLine(cs, fos);	// id is invalid now
			if (nuLen<= settings.getFragNBlambda())
				++tgtMols;
			if (len<= settings.getFragNBlambda())
				--tgtMols;

		}

		/**
		 * @deprecated
		 * @param s
		 * @param b
		 * @return
		 */
		private double[] getRTCDF(CharSequence s, boolean b) {
			
			double[] a= new double[s.length()];
			Arrays.fill(a, 0d);
			HashMap<CharSequence, double[]> map= null; // b? mapSense: mapAsense;			
			for (int i = 0; i < a.length- (mapMotifLen+ mapKmerLen); ++i) {
				double p= 1;
				int x= (b? mapMotifFirst: mapMotifLen- mapMotifFirst+ 1);
				
				for (int j = 0; j < mapMotifLen; j++) {
					CharSequence t= s.subSequence(i+ j, i+ j+ mapKmerLen);
					p*= map.get(t)[j]; // map.get(t)[j]; // Math.log(map.get(t)[j]);
				}
				
				if (p> 20)
					System.currentTimeMillis();
				a[i+ x]= p;
			}
			
			return Fragmenter.toCDF(a);
		}

		RandomDataImpl rndImpl= new RandomDataImpl(); // DEBUG
		void processFrag_curr(int start, int end, int len, ByteArrayCharSequence id) {
            
			double lambda= settings.getFragNBlambda();
/*			if (len<= lambda) {	// does not break
				processFragNot(start, end, len, id);
				return;
			}
*/			
			
			
			// first check whether frag gets broken
			//double lfac= len/ (double) lambda;	// len > lambda
			//double val= 1d- Math.exp((lambda- len)/ lambda);	// len> lambda, P(b)
			double val= 
				// linear
				// lambda/ len;
				//(len- lambda)/ lambda;
				//(len- lambda)/ (2* lambda);
				
				// polynomial
				//Math.pow((len- lambda)/ lambda, 2d);	// len> lambda, P(b)
			
			
				// exponential
				//1d- Math.exp((lambda- len)/ lambda);	// len> lambda, P(b)
				//1d- Math.exp(-Math.pow(len/ lambda, 2d));	// len> lambda, P(b)
				//1d- (1d/ Math.pow(len- lambda, 2d));	// len> lambda, P(b)
				//1d- Math.pow(lambda/ len, 2);
				//1d- Math.exp(-Math.pow(len/ lambda, 2));
				//1d- (1d/ len);
				1d;
			
				//1d- Math.pow(Math.exp(-(len/ lambda)), 2d);	// len> lambda, P(b)
				
			if (val< 1d) {
	            double r= rndBreak.nextDouble();
	            if (r> val) {   // does not break (!)
	            	processFragNot(start, end, len, id);
	            	return;
	            }
			}
            
            // else
			int howMany= (int) rndHowMany.nextPoisson(Math.pow(len/ lambda, 2d)); // len > lambda, nr events
			if (howMany== 0) {
//				howMany= 1;
				processFragNot(start, end, len, id);
				return;
			}
			if (howMany>= where.length)
				howMany= where.length- 1;
			//howMany= 1;	// DEBUG
			//double mid= len/ 2d;
			for (int i = 0; i < howMany; i++) {
				where[i]= 
//				where[i]= start+ (int) rndHowMany.nextGaussian(len/ 2d, len/ 4d);
//				if (where[i]< start|| where[i]> end)
//					--i;
				//start+ rndFrag.nextInt(len);	// bp, len exclusive
					
				//start+ (len/ 2);
				(int) nextGaussianDouble(rndBP, start, end);
//				start+ (int) (rndImpl.nextGaussian(mid, mid/ 2));
//				while (where[i]< start|| where[i]> end)
//					where[i]= start+ (int) (rndImpl.nextGaussian(mid, mid/ 2));
			}
			// DEBUG
			int len1= where[0]- start;
			int len2= end- where[0]+ 1;
			if (len1< len2) {
				int h= len1;
				len1= len2;
				len2= h;
			}
			
/*			if (len2<= lambda) {	// does not break
				processFragNot(start, end, len, id);
				return;
			}
*/			
/*			double val= 
				//1d- Math.exp(-Math.pow(len2/(double) lambda, 2d));
				 //(len2/ lambda);	// len> lambda, P(b)
				1d- (1/ (1+ Math.exp(-Math.pow(len2/(double) lambda, 2d))));
			double r= rndBreak.nextDouble();
			if (r> val) {
				processFragNot(start, end, len, id);
				return;
			}
*/			
			
			if (howMany> 1) 
				Arrays.sort(where, 0, howMany);
			int lastVal= start;
			if (howMany> 1)
				id= id.cloneCurrentSeq();	// gets invalid with replace
			int nowStart= start, nowEnd= end;
			
			for (int i = 0; i < howMany; i++) {
				if (where[i]== lastVal)
					continue;	// 0-length frag, 2x hydrolized at same point
				++newMols;
				int nuLen= where[i]- lastVal;	// (lastVal- 1)+ 1
				lenSum+= nuLen;
				++lenNb;
				if (nuLen>= Fragmenter.this.settings.getFiltMin()&&
						nuLen<= Fragmenter.this.settings.getFiltMax())
					++Fragmenter.this.tgtMols;
				
				// plot
				int nuEnd= where[i]- 1;	// floor(bp)
				if (Fragmenter.this.plotter!= null) 
					Fragmenter.this.plotter.plot(lastVal, nuEnd, -1, id);

				// write
				if (lastVal!= nowStart) {
					cs.replace(0, lastVal);
					nowStart= lastVal;
				}
				if (nuEnd!= nowEnd) {
					cs.replace(1, nuEnd);
					nowEnd= nuEnd;
				}
				assert(lastVal<= nuEnd);
				rw.writeLine(cs, fos);	// id is invalid now
				lastVal= where[i];
			}
			// last one, always-even if lastVal== end
			if (howMany> 0) {
				if (where[howMany- 1]!= nowStart) 
					cs.replace(0, where[howMany- 1]);
				if (nowEnd!= end) 
					cs.replace(1, end);
				assert(where[howMany- 1]<= end);
				rw.writeLine(cs, fos);	// id is invalid now
			}
		}

		private void processFragNot(int start, int end, int len,
				ByteArrayCharSequence id) {
        	lenSum+= len;                            	
        	++lenNb;
			maxLen= Math.max(maxLen,len);
			minLen= Math.min(minLen,len);
			incLinesWrote();
			rw.writeLine(cs, fos);
			//increaseFragCount(id);	// 20101215 deactivated, counts for pro-file
			if (Fragmenter.this.plotter!= null)
				Fragmenter.this.plotter.plot(start, end, -1, id);
		}

		private boolean multiBreaks= false; 
		/**
		 * @deprecated not in use
		 * @param start
		 * @param end
		 * @param len
		 * @param id
		 */
		void processNebu(int start, int end, int len, ByteArrayCharSequence id) {
		            
					// fragment breaks dependent on bp
					double lambda= settings.getFragNBlambda()/ 2d;
					
					int howMany= 1;
					if (multiBreaks) {
						howMany= (int) rndHowMany.nextPoisson(len/ lambda); // len > lambda, nr events
						if (howMany== 0)
							howMany= 1;	// assure one
						else if (howMany>= where.length)
							howMany= where.length- 1;	// lim mem
					}
					for (int i = 0; i < howMany; i++) {
						where[i]= 
						(int) nextGaussianDouble(rndNebu, start, end);
						if (i> 1&& where[i]- where[i- 1]+ 1< minLen)
							minLen= where[i]- where[i- 1]+ 1;
					}
					
					// break condition
					int minLen= Math.min(where[0]- start+ 1, end- where[0]+ 1);
					double val= 
						//1d- Math.exp(-Math.pow(minLen/(double) lambda, 2d));
						1d- Math.exp(-Math.pow(minLen/ lambda, 2d));	// len> lambda, P(b)
					double r= rndBreak.nextDouble();
					if (r> val) {
						processFragNot(start, end, len, id);
						return;
					}
					
					
					if (howMany> 1) 
						Arrays.sort(where, 0, howMany);
					int lastVal= start;
					if (howMany> 1)
						id= id.cloneCurrentSeq();	// gets invalid with replace
					int nowStart= start, nowEnd= end;
					addFragCount(id, (long) (howMany));
					for (int i = 0; i < howMany; i++) {
						if (where[i]== lastVal)
							continue;	// 0-length frag, 2x hydrolized at same point
						++newMols;
						int nuLen= where[i]- lastVal;	// (lastVal- 1)+ 1
						lenSum+= nuLen;
						++lenNb;
						if (nuLen>= Fragmenter.this.settings.getFiltMin()&&
								nuLen<= Fragmenter.this.settings.getFiltMax())
							++Fragmenter.this.tgtMols;
						
						// plot
						int nuEnd= where[i]- 1;	// floor(bp)
						if (Fragmenter.this.plotter!= null) 
							Fragmenter.this.plotter.plot(lastVal, nuEnd, -1, id);
		
						// write
						if (lastVal!= nowStart) {
							cs.replace(0, lastVal);
							nowStart= lastVal;
						}
						if (nuEnd!= nowEnd) {
							cs.replace(1, nuEnd);
							nowEnd= nuEnd;
						}
						assert(lastVal<= nuEnd);
						rw.writeLine(cs, fos);	// id is invalid now
						lastVal= where[i];
					}
					// last one, always-even if lastVal== end
					if (where[howMany- 1]!= nowStart) 
						cs.replace(0, where[howMany- 1]);
					if (nowEnd!= end) 
						cs.replace(1, end);
					assert(where[howMany- 1]<= end);
					rw.writeLine(cs, fos);	// id is invalid now
					if (end- where[howMany- 1]+ 1<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
				}

		void processFragUniform_save(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
			
			double lambda= settings.getFragNBlambda();
			double val= 
				Math.log10(len/ (lambda/ 5d));	//len> lambda must be able to not break
				//Math.log10(len/ (lambda/ 2d));
				
			if (val< 1d) {
		        double r= rndBreak.nextDouble();
		        if (r> val) {   // does not break (!)
		        	processFragNot(start, end, len, id);
		        	return;
		        }
			} 
			
			double r= rndBP.nextDouble();
			int bp= (int) (len* r);
		
			++newMols;
			if (bp<= Fragmenter.this.settings.getFragNBlambda())
				++Fragmenter.this.tgtMols;
			cs.replace(1, start+ bp);
			incLinesWrote();
			rw.writeLine(cs, fos);	// id is invalid now
				
			++newMols;
			if (len- bp<= Fragmenter.this.settings.getFragNBlambda())
				++Fragmenter.this.tgtMols;
			cs.replace(0, start+ bp);
			cs.replace(1, end);
			incLinesWrote();
			rw.writeLine(cs, fos);	// id is invalid now
		}
		
		private int selectRTSite(ByteArrayCharSequence id, int start, int end, boolean firstStrand) {
		
			float[] a= (firstStrand? mapPWMasense.get(id): mapPWMsense.get(id));
			float[] b= new float[end- start+ 1];
			//int len= Math.min(b.length, a.length- start);
			System.arraycopy(a, start, b, 0, b.length);
			//if (b.length> len)
			//	Arrays.fill(b, len, b.length- len, (byte) 'A');
//			for (int i = 0; i < b.length; i++) {
//				b[i]= (float) Math.exp(10* b[i]);
//			}
			b= toCDF(b);

			int howmany= Math.max(10, b.length/ 60); 	// hexamers
			int pos= firstStrand? -1: b.length; 
			for (int i = 0; i < 1; i++) {
				float r= rtRndWhere.nextFloat();
				int p= Arrays.binarySearch(b, r);
				if (p< 0) 
					p= -(p+ 1);
				pos= firstStrand? Math.max(pos, p): Math.min(pos, p);
			}
			
			return pos;
		}

		private float[] toCDF(float[] b) {
			
			double sum= 0f;
			for (int i = 0; i < b.length; ++i) 
				sum+= (double) b[i];
			double cumu= 0d;
			for (int i = 0; i < b.length; ++i) {
				cumu+= b[i]/ (double) sum;
				b[i]= (float) cumu;
			}
			b[b.length- 1]= 1f;	// correct for rounding errors, TODO output warning if d> thr?
			
			return b;
		}

		private boolean processFragDiss(int nuLen) {
			if (true)
				return false;
			double term= Math.pow(nuLen/ (settings.getFragNBlambda()/ 5d), 2d);	// lambda/ 5d ?
			double val= 
				// Math.exp(-Math.pow(Math.log10(nuLen),2d));
				Math.exp(-term);
			double r= rndDisappear.nextDouble();
			return (r<= val);
		}

		
		double dnorm(double x, double mu) {
			
			double a1= Math.pow(x- mu, 2d);
			
			return Math.exp(-a1/ 2d)/ Math.sqrt(2* Math.PI);
			
		}

		int lastLen= -1;
		double lastP= Double.NaN, lastQ= Double.NaN;
		/**
		 * see <code>http://personal.strath.ac.uk/gary.koop/extra_material_on_metropolis.pdf</code>,
		 * <code>http://www.ps.uci.edu/~markm/numerical_methods/Metropolis%96Hastings%20algorithm%20-%20Wikipedia,%20the%20free%20encyclopedia.pdf</code>
		 * @param start
		 * @param end
		 * @param len
		 * @param id
		 */
		void processFilterMCMC(int start, int end, int len, ByteArrayCharSequence id, 
				AbstractDistribution dGenerate, AbstractDistribution[] dProposal) {

			// first value always accepted to init algorithm (but not output)
			double p= 0d;
			for (int i = 0; i < dProposal.length; i++) {
				p+= dProposal[i].getP(len);
			}
			if (lastLen< 0) {
				lastLen= len;
				lastP= p;
				return;
			}
			
			// Metropolis/Hastings/Ema
			double a1= p/ lastP;
			double a2= dGenerate.getP(lastLen, len)/ dGenerate.getP(len, lastLen);
			double a= a1* a2;
			
			// accept 
			if (a>= 1|| rndGel.nextDouble()<= a) {
				lastLen= len;
				lastP= p;
				incLinesWrote();
				rw.writeLine(cs, fos);
				if (Fragmenter.this.plotter!= null) 
					Fragmenter.this.plotter.plot(start, end, len, id);
				addFragCount(id, 1l);
			} else
				--newMols;

		}

		void processInitial() {
			
			// transcript variation
			int start= 0;
			int end= origLen- 1;
			if (!Double.isNaN(settings.getTssMean())) {
				start= origLen;
				while (start>= Math.min(100, origLen))
					start= (int) Math.round(rndTSS.nextExponential(Math.min(settings.getTssMean(),origLen/4)));	// exp mean= 25: exceeds bounds, nextGaussian(1,100))-100;
				double r= rndPlusMinus.nextDouble();
				if (r< 0.5)
					start= -start;
			}
			if (!(Double.isNaN(settings.getPolyAshape())|| Double.isNaN(settings.getPolyAscale()))) {
				int pAtail= 301;
				while (pAtail> 300)
					pAtail= (int) Math.round(sampleWeibull(rndPA, settings.getPolyAscale(), settings.getPolyAshape()));	// 300d, 2d
				end= origLen+ pAtail; 
			}
			if (end< origLen) {
				if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
					System.err.println("[ERROR] end < length in Fragmenter$Processor.processInitial(): "+ end+ "<"+ origLen);
			}

			int newLen= end- start+ 1;
			assert(start< end);
			lenSum+= newLen;
			cumuLen+= newLen;
			++lenNb;
			maxLen= Math.max(maxLen,newLen);
			minLen= Math.min(minLen,newLen);
			
			
			//String line= start+ TAB+ end+ TAB+ settings.getProfiler().getIds()[i]+"\n";
			//ByteArrayCharSequence ccs= cs.cloneCurrentSeq();
			cs.replace(0, start);
			cs.replace(1, end);
			cs.replace(2, locID);
			
			// RFU: chr_tid separated
/*			int p= 0;
			while(locID.charAt(++p)!= ':');
			cs.end= cs.p2;
			cs.ensureLength(cs.end, p);
			byte[] a= cs.a;
			for (int j = 0; j < p; j++) 
				a[cs.end++]= (byte) locID.charAt(j); 
			cs.append(BYTE_SEP_LC_TX);
			cs.append(locID);
*/			
			

		}

		public FileOutputStream getFos() {
			return fos;
		}

		public void setFos(FileOutputStream fos) {
			this.fos = fos;
		}

		public void close() {
    		end= true;
    		while (isAlive())
        		interrupt();
	    		try {
					join();
				} catch (InterruptedException e) {
					; // :)
				}
    	}

		void processFrag_110108(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
				            
							double lambda= nebu? settings.getFragNBlambda():settings.getFragNBlambda();	// nebu / 2d
		
							if (!nebu) {
								// first check whether frag gets broken
								//double val= 
									// linear
									// lambda/ len;
									//(len- lambda)/ lambda;
									//(len- lambda)/ (2* lambda);
									
									// polynomial
									//Math.pow((len- lambda)/ lambda, 2d);	// len> lambda, P(b)
								
								
									// exponential
									//1d- Math.exp((lambda- len)/ lambda);	// len> lambda, P(b)
									//1d- Math.exp(-Math.pow(len/ lambda, 2d));	// len> lambda, P(b)
									//1d- Math.pow(Math.exp(-(len/ lambda)), 2d);	// len> lambda, P(b)
									
									//1d- Math.pow(len/ lambda, -0.5d);
									//Math.log(1d+ Math.pow(len/ lambda, 0.5));
									//Math.sqrt(Math.log(1d+len/ lambda));
									//Math.log(Math.sqrt(1d+ len/lambda));
									//1d- Math.exp(-len);
									//1d- (1d/ Math.pow(len, len/lambda));
									
									//1d- (lambda/ Math.pow(len, 1));	// friday
									
		//							double val= 1d;
		//							if (len< lambda)
		//								val= 1d/ Math.log10(lambda- len);
								
								double val= 
									Math.log10(len/ (lambda/ 5d));	//len> lambda must be able to not break
									//Math.log10(len/ (lambda/ 2d));
									
								if (val< 1d) {
						            double r= rndBreak.nextDouble();
						            if (r> val) {   // does not break (!)
						            	processFragNot(start, end, len, id);
						            	return;
						            }
								} 
									
							}
							
				            
				            // else
							int howMany= (len> 1?1:0);
							if (multiBreaks) {
								howMany= (int) rndHowMany.nextPoisson(len/ lambda); // nr events
								if (howMany== 0&& len> 1)
									howMany= 1;
								else if (howMany>= where.length)
									howMany= where.length- 1;
							}
							if (howMany== 0) {
								processFragNot(start, end, len, id);
								return;
							}
							int minLen= Integer.MAX_VALUE, maxLen= Integer.MIN_VALUE;
							for (int i = 0; i < howMany; i++) {
								where[i]= start- 1;
								while (where[i]< start|| where[i]> end) {
									
									// before 110108
									where[i]= //(int) nextGaussianDouble(rndBP, start, end);	// bp
										(int) (nebu? 
										//nextGaussianDouble(rndNebu, start, end):
										//start+ rdiNebuBP.nextGaussian(len/ 2d, Math.sqrt(len)):	// bp
										start+ rdiNebuBP.nextGaussian(len/ 2d, len/ 4d):	// bp
										
										//start+ nextDoubleExp(rndBP)* len);
											// something that falls exponential, either weibull or gauss/normal
										//start+ (int) (sampleWeibull(rndBP, 0.5, 5)*  len));
										
										nextGaussianDouble(rndBP, start, end));
										//start+ 5+ ((int) ((rndBP.nextGaussian()* (len- 10))/ 2d)));
										//start+ rndFrag.nextInt(len);	// bp, len exclusive
		
									
								}
								int nuLen= where[i]- (i== 0? start: where[i- 1])+ 1;
								if (nuLen< minLen)
									minLen= nuLen;
								if (nuLen> maxLen)
									maxLen= nuLen;
							}
							int lastLen= end- where[howMany- 1]+ 1;
							if (lastLen< minLen)
								minLen= lastLen;
							if (lastLen> maxLen)
								maxLen= lastLen;
		
							// nebu break condition, after bp has been choosen
							if (nebu) {
								double val= 
									//1d- Math.exp(-Math.pow(minLen/(double) lambda, 2d));
									// nebu
									//1d- Math.exp(-Math.pow(minLen/ lambda, 2d));	// len> lambda, P(b)
									1d- Math.exp(-(minLen/ lambda));	// len> lambda, P(b)
								
								
									//1d- (lambda/ (2* Math.pow(maxLen, 2)));
									//1d- ((lambda/ 2d)/ maxLen);
									
									//1d- Math.pow(maxLen, -(len/ lambda));
									//Math.log10(lambda)/ Math.log10(len);
									//Math.exp(-Math.pow(maxLen/ lambda, 2d));
									//Math.pow()
								double r= rndBreak.nextDouble();
								if (r> val) {
									processFragNot(start, end, len, id);
									return;
								}
							}
		
							
							if (howMany> 1) 
								Arrays.sort(where, 0, howMany);
							int lastVal= start;
							if (howMany> 1)
								id= id.cloneCurrentSeq();	// gets invalid with replace
							int nowStart= start, nowEnd= end;
							addFragCount(id, (long) (howMany));
							tstNewMols+= howMany;
							for (int i = 0; i < howMany; i++) {
								if (where[i]== lastVal)
									continue;	// 0-length frag, 2x hydrolized at same point
								int nuLen= where[i]- lastVal;	// (lastVal- 1)+ 1
								if (processFragDiss(nuLen))
									continue;
		
								++newMols;
								lenSum+= nuLen;
								++lenNb;
								if (nuLen<= Fragmenter.this.settings.getFragNBlambda())
									++Fragmenter.this.tgtMols;
								
								// plot
								int nuEnd= where[i]- 1;	// floor(bp)
								if (Fragmenter.this.plotter!= null) 
									Fragmenter.this.plotter.plot(lastVal, nuEnd, -1, id);
				
								// write
								if (lastVal!= nowStart) {
									cs.replace(0, lastVal);
									nowStart= lastVal;
								}
								if (nuEnd!= nowEnd) {
									cs.replace(1, nuEnd);
									nowEnd= nuEnd;
								}
								assert(lastVal<= nuEnd);
								incLinesWrote();
								rw.writeLine(cs, fos);	// id is invalid now
								lastVal= where[i];
							}
							// last one, always-even if lastVal== end
							int nuLen= end- where[howMany- 1]+ 1;
							if (processFragDiss(nuLen))
								return;
							if (where[howMany- 1]!= nowStart) 
								cs.replace(0, where[howMany- 1]);
							if (nowEnd!= end) 
								cs.replace(1, end);
		//					if (where[howMany- 1]> end)
		//						System.currentTimeMillis();
							assert(where[howMany- 1]<= end);
							incLinesWrote();
							rw.writeLine(cs, fos);	// id is invalid now
							if (nuLen<= settings.getFragNBlambda())
								++tgtMols;
							if (len<= settings.getFragNBlambda())
								--tgtMols;
						}

		void processRT_PWM(int start, int end, int len, ByteArrayCharSequence id) {
			
			//int[] myWhere= where, myStarts= starts;
			--newMols;	// substract the original one	
			int right= selectRTSite(id, start, end, true);
			if (right< 0)
				return;
			// foff of poly?
			int left= selectRTSite(id, start, start+ right, false);
			if (left< 0|| left>= right)
				return;
			
			// write ! may be iterated multiple times due to howmany
			id= id.cloneCurrentSeq();	// id invalid after replace operation !
			++newMols;
			//addFragCount(id, (long) 1);
			cs.replace(0, start+ left);
			cs.replace(1, start+ right);
			incLinesWrote();
			rw.writeLine(cs, fos);
			
		}
		
		void processRT_Bernoulli(int start, int end, int len, ByteArrayCharSequence id) {
		
				//int[] myWhere= where, myStarts= starts;
				--newMols;	// substract the original one	
				int right= selectRTSite(id, start, end, true);
				if (right< 0)
					return;
				// foff of poly?
				int left= selectRTSite(id, start, start+ right, false);
				if (left< 0)
					return;
				
				// write ! may be iterated multiple times due to howmany
				id= id.cloneCurrentSeq();	// id invalid after replace operation !
				++newMols;
				//addFragCount(id, (long) 1);
				cs.replace(0, start+ left);
				cs.replace(1, start+ right);
				incLinesWrote();
				rw.writeLine(cs, fos);
				
		}

		private void determineBPos(ByteArrayCharSequence id, int start, int end, int[] where, int howMany, boolean nebu) {
			
			float[] b= null;
			if (mapPWMsense!= null) {
				float[] a= mapPWMsense.get(id);
				b= new float[end- start+ 1];
				//int len= Math.min(b.length, a.length- start);
				System.arraycopy(a, start, b, 0, b.length);
				//if (b.length> len)
				//	Arrays.fill(b, len, b.length- len, (byte) 'A');
			}
			
			int len= end- start+ 1;
			for (int i = 0; i < howMany; i++) {
				where[i]= start- 1;
				while (where[i]< start|| where[i]> end) {
					
					if (b== null) {
						where[i]= //(int) nextGaussianDouble(rndBP, start, end);	// bp
							(int) 
							(nebu? 
							(start+ rdiNebuBP.nextGaussian(len/ 2d, len/ 4d))	// bp
							
							:nextGaussianDouble(rndBP, start, end))
							;
					} else {
						if (settings.getBpDistr()== FluxSimulatorSettings.DISTR_NORMAL) {							
							for (int j = 0; j < b.length; j++) {
								double x= ((j* TWICE_CUT_OFF_GAUSSIAN_VAL)/ len)- CUT_OFF_GAUSSIAN_VAL;
								b[j]= (float) (NORM_FACTOR* Math.pow(E_POW_12, x* x));
							}
						} else if (settings.getBpDistr()== FluxSimulatorSettings.DISTR_GAUSS) {
							double mu= len/ settings.getBpDistrPar()[0], 
									si= len/ settings.getBpDistrPar()[1];
							for (int j = 0; j < b.length; j++) {
								double d= j- mu;
								b[j]= (float) (NORM_FACTOR* (1d/ si)* Math.exp((d*d)/(-2* si* si)));
							}
						}
						b= toCDF(b);
						float r= rndBP.nextFloat();
						int p= Arrays.binarySearch(b, r);
						if (p< 0) 
							p= -(p+ 1);
						where[i]= start+ p;
					}
					
					
				}
			}
		}

		private float[] toPDF(float[] b) {
			
			double sum= 0f;
			for (int i = 0; i < b.length; ++i) 
				sum+= (double) b[i];
			for (int i = 0; i < b.length; ++i) 
				b[i]/= sum;
			
			return b;
		}

		private int selectRTSiteBernoulli(ByteArrayCharSequence id, int start, int end, boolean firstStrand) {
		
			float[] a= (firstStrand? mapPWMasense.get(id): mapPWMsense.get(id));
			float[] b= new float[end- start+ 1];
			//int len= Math.min(b.length, a.length- start);
			System.arraycopy(a, start, b, 0, b.length);
			//if (b.length> len)
			//	Arrays.fill(b, len, b.length- len, (byte) 'A');
			b= toPDF(b);
			
			for (int i = firstStrand? b.length- 1: 0; 
				(firstStrand&& i>= 0)|| ((!firstStrand)&& i< b.length); 
				i+= (firstStrand? -1: 1)) {
				double r= rtRndWhere.nextDouble();
				if (r< b[i])
					return i;
			}
			
			return -1;
		}

		void processFragUniform(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {

			int bp= -1;
			if (mapPWMsense== null) {
				double r= rndBP.nextDouble();
				bp= (int) (len* r);
			} else {
				if (mapPWMsense!= null) {
					float[] a= mapPWMsense.get(id);
					float[] b= new float[end- start+ 1];
					System.arraycopy(a, start, b, 0, b.length);
					float r= rndBP.nextFloat();
					bp= Arrays.binarySearch(b, r);
					if (bp< 0) 
						bp= -(bp+ 1);

				}
			}
			
			
//			int bp= (int) nextGaussianDouble(rndBP, start, end)- start;
			int minLen= Math.min(bp, len- bp);
			
			double f= minLen/ 15d;
			double e= Math.exp(-(f- 1));
			double ee= Math.pow(e, 0.1d);	// not perfect, but ok
			double pe= 1d- ee;	// len> lambda, P(b)
			
//			double pp= f< 1? 0: Math.pow(f- 1, 1d);
			
//			double d= minLen- 200d;
//			double dp= Math.pow(d, 2d);
//			double pd= 1d- (1d/ Math.exp(dp));
//			double pd2= 1d- (1d/ Math.exp(dp/ 200d));
			
	        double r= rndBreak.nextDouble();
	        if (r> pe) {   // does not break (!)
	        	processFragNot(start, end, len, id);
	        	return;
	        }
			
			++newMols;
			if (bp<= Fragmenter.this.settings.getFragNBlambda())
				++Fragmenter.this.tgtMols;
			cs.replace(1, start+ bp);
			incLinesWrote();
			rw.writeLine(cs, fos);	// id is invalid now
				
			++newMols;
			if (len- bp<= Fragmenter.this.settings.getFragNBlambda())
				++Fragmenter.this.tgtMols;
			cs.replace(0, start+ bp);
			cs.replace(1, end);
			incLinesWrote();
			rw.writeLine(cs, fos);	// id is invalid now
		}

		void processFilter_110210(int start, int end, int len, ByteArrayCharSequence id) {
					//RandomDataImpl rndGel= new RandomDataImpl();
					//rndGel.nextGaussian(mu, sigma);
					//rndGel.nextPoisson(mean)
					//if (len>= getDeltaFiltMin()&& len<= getDeltaFiltMax()) {
						int segLen= getSegregatedLength(len, true);
						if (segLen>= Fragmenter.this.settings.getFiltMin()&&
								segLen<= Fragmenter.this.settings.getFiltMax()) {
							incLinesWrote();
							rw.writeLine(cs, fos);
							if (Fragmenter.this.plotter!= null) {
								Fragmenter.this.plotter.plot(start, end, segLen, id);
			//					if (!Fragmenter.this.plotter.plot(start, end, segLen, id))
			//						Fragmenter.this.plotter.plot(start, end, segLen, id);
							}
							
							addFragCount(id, 1l);
						} else
							--newMols;
		//			} else
		//				--newMols;
					// increaseFragCount(token[2]);
						
		
				}

		void processFragPaolo(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
            
			// Paolos version
			int nr= (int) Math.round(rdiNebuBP.nextUniform(1, 5));
			nr= 7;
			
			for (int i = 0; i < nr; i++) {
			
				if (len< 2)
					continue;
				
				// calc cdf
				double[] cdf= new double[len- 1];
				try {
//					PrintWriter pp= new PrintWriter("test_cdf");
					for (int j = 0; j < cdf.length; ++j) {
						double p= getWeibullProb(j, 200, 2);
						p+= (j== 0? 0: cdf[j- 1]);
						cdf[j]= p;
//						pp.println(j+"\t"+cdf[j]);
					}
//					pp.flush();
//					pp.close();
					System.currentTimeMillis();
				} catch (Exception e) {
					e.printStackTrace();
				}
				
				// "extract" length (Palo's language)
				double r= rndBP.nextDouble()* cdf[cdf.length- 1];
				int bp= Arrays.binarySearch(cdf, r);
				if (bp< 0)
					bp= -(bp+ 1);
				//bp= 200;
				//bp= Math.max(bp, len- bp);
				
				if (bp<= 0|| bp+ 1>= len) {
					System.currentTimeMillis();
					continue;
				}
				//boolean ff= rndBreak.nextBoolean();
				double ff= 1d- rndBreak.nextDouble();
				if (ff>= 0.5) {
					int nuLen= bp- start+ 1;
					updateMedian(nuLen);
					cs.replace(0, start);
					cs.replace(1, start+ bp);
					cumuLen+= nuLen;
					incLinesWrote();
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
					start= start+ bp+ 1;
					len= end- start+ 1;
					
				} else {
					int nuLen= end- bp;
					updateMedian(nuLen);
					cs.replace(0, end- bp);
					cs.replace(1, end);
					cumuLen+= nuLen;
					incLinesWrote();
					++newMols;
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
					end= end- bp- 1;
					len= end- start+ 1;
				}
				
				if (len== 0)
					System.currentTimeMillis();
			}
			
		}				
			

		
		
		void processFragBreakNDeg(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
				            
			double lambda= nebu? settings.getFragNBlambda():settings.getFragNBlambda();	// nebu / 2d

			if (len<= 10) {
				cumuLen+= len;
				updateMedian(len);
				processFragNot(start, end, len, id);
				return;
			}
				
			
			int bp= -1;
			while (bp< start|| bp >= end)
				//bp= (int) nextGaussianDouble(rndBP, start, end);
				bp= (int) rdiNebuBP.nextGaussian(start+ (len/ 2d), 70);
/*			double pp= rndBP.nextDouble();
			
			bp= (int) ((len/2d)+ (Math.sqrt(Math.abs(pp- 0.5d))* (Math.sqrt(0.5)* len/ 2d)));
			double pp2= (0.5d- (-2d* (pp- 0.5d)* (pp- 0.5d)));
			bp= (int) (len* pp2);
			bp= (pp> 0.5? len- bp: bp);
*/			
			// bp= (int) (pp* len); // uniform
			int min= Math.min(bp- start+ 1, end- bp);
			int max= len- min;
//			bp+= start;

			double p= 1d- Math.exp(-Math.pow(len/ 700d, 2d));
			if (p< 1d) {
	            double r= rndBreak.nextDouble();
	            if (r> p) {   // does not break (!)
	            	bp= end;
	            }	            
			} else 
				bp= end;
			
			--newMols;
			if(bp> len/ 2d) {
				if (bp< end) {
					processFragDeg(start, bp, id, 2);
					processFragDeg(bp+1, end, id, 3);
				} else
					processFragDeg(start, bp, id, 0);
			} else {
				if (bp< end) {
					processFragDeg(start, bp, id, 3);
					processFragDeg(bp+1, end, id, 1);
				} else
					processFragDeg(start, bp, id, 0);
			}
			
		}

		private void processFragDeg(int start, int end, ByteArrayCharSequence id, int mod) {
        	
			int len= end- start+ 1;
			double p= 1- (250d/ len);
			double r= rndBreak.nextDouble();
			if (mod== 3|| r> p) {
				updateMedian(len);
				cs.replace(0, start);
				cs.replace(1, end);
				cumuLen+= len;
				incLinesWrote();
				rw.writeLine(cs, fos);	// id is invalid now
				if (len<= settings.getFragNBlambda())
					++tgtMols;
				if (len<= settings.getFragNBlambda())
					--tgtMols;
				return;
			}

			
			// hydrolyze
			if (len< 0)
				System.currentTimeMillis();
			int deg= len;
			double mu= Math.min(200d, len/ 2d);
        	while(deg>= len- 1)
        		//deg= 1+ (int) rdiNebuBP.nextExponential(Math.min(200d,len/2d));	//
        		deg= (int) rdiNebuBP.nextUniform(11, Math.min(250, len));
        		//deg= rdiNebuBP.nextGaussian(mu, mu/2);
        	
        	double ff= -1d;
        	if (mod== 1)
        		ff= 0;
        	else if (mod== 2)
        		ff= 1;
        	else if (mod== 0) 
        		ff= rndBreak.nextDouble();
        	
        	len-= (int) (len- deg);
        	if (ff< 0.5d) {
        		processFragWrite(start, start+ deg- 1, deg);
    			processFragWrite(start+ deg, end, len);
        	} else { 
        		processFragWrite(start, end- deg, len);
    			processFragWrite(end- deg+ 1, end, deg);
        	}
        	
		}
		
		private void processFragWrite(int start, int end, int len) {
			
			if (start> end)
				System.currentTimeMillis();
			updateMedian(len);
			cs.replace(0, start);
			cs.replace(1, end);
			cumuLen+= len;
			incLinesWrote();
			++newMols;
			rw.writeLine(cs, fos);	// id is invalid now
			if (len<= settings.getFragNBlambda())
				++tgtMols;
			if (len<= settings.getFragNBlambda())
				--tgtMols;

		}

		private void updateMedian(int nuLen) {
			
			int[] a= med;
			if (med[0]== 0) {
				med[0]= nuLen;
				for (int i = 0; i < a.length- 1; i++) {
					if (med[i]> med[i+ 1]) {
						int h= med[i];
						med[i]= med[i+ 1];
						med[i+ 1]= h;
					}
				}
				return;
			}
			
			int p= Arrays.binarySearch(med, nuLen);
			if (p>= 0)
				return;

			p= -(p+ 1);			
			if (p== 0|| p>= med.length)
				return;
			if (p+ 1< med.length)
				System.arraycopy(med, p, med, p+ 1, med.length- p- 1);
			med[p]= nuLen;
		}

		void processFragDeg(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
						            
					double lambda= nebu? settings.getFragNBlambda():settings.getFragNBlambda();	// nebu / 2d
		
					int bp= -1;
					while (bp< start|| bp >= end)
						bp= (int) nextGaussianDouble(rndBP, start, end);
						//bp= (int) rdiNebuBP.nextGaussian(start+ (len/ 2d), Math.min(100, len/4d));
		/*			double pp= rndBP.nextDouble();
					
					bp= (int) ((len/2d)+ (Math.sqrt(Math.abs(pp- 0.5d))* (Math.sqrt(0.5)* len/ 2d)));
					double pp2= (0.5d- (-2d* (pp- 0.5d)* (pp- 0.5d)));
					bp= (int) (len* pp2);
					bp= (pp> 0.5? len- bp: bp);
		*/			
					// bp= (int) (pp* len); // uniform
					int min= Math.min(bp- start+ 1, end- bp);
					int max= len- min;
		//			bp+= start;
		
					double p= 1d- Math.exp(-Math.pow(min/ 400d, 2d));
					if (p< 1d) {
			            double r= rndBreak.nextDouble();
			            if (r> p) {   // does not break (!)	            	
			            	if (r< Math.pow(len/ 200d, 3d)) {
				            	double deg= len;
				            	while(deg>= len- 1)
				            		deg= rdiNebuBP.nextExponential(Math.min(50, len/10));
				            	int nuLen= (int) (len- deg);
				            	double ff= rndBreak.nextDouble();
			        			updateMedian(nuLen);
			        			cumuLen+= nuLen;
				            	if (ff< 0.5d) {
				        			cs.replace(0, start+ (int) deg);
				        			cs.replace(1, end);
				            	} else {
				        			cs.replace(0, start);
				        			cs.replace(1, end- (int) deg);
				            	}
			        			incLinesWrote();
			        			rw.writeLine(cs, fos);	// id is invalid now
				            	return;
			            	} else {
			            		processFragNot(start, end, len, id);
			            		return;
			            	}
			            }
			            
					}
					
					int nuLen= bp- start+ 1;
					updateMedian(nuLen);
					cs.replace(0, start);
					cs.replace(1, bp);
					cumuLen+= nuLen;
					incLinesWrote();
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
					
					nuLen= end- bp;
					updateMedian(nuLen);
					cs.replace(0, bp+ 1);
					cs.replace(1, end);
					cumuLen+= nuLen;
					incLinesWrote();
					++newMols;
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
						
					
				}

		void processFrag110214(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
						            
					double lambda= nebu? settings.getFragNBlambda():settings.getFragNBlambda();	// nebu / 2d
		
		        	double deg= len;
		        	while(deg>= len- 1)
		        		deg= rdiNebuBP.nextExponential(10* Math.pow(len, 1d/2d));
		        	double ff= rndBreak.nextDouble();
		        	if (ff< 0.5d) 
		    			start+= (int) deg;
		        	else 
		    			end-= (int) deg;
		        	int nuLen= (int) (len- deg);
		        	
					int bp= -1;
					while (bp< start|| bp >= end)
						bp= (int) nextGaussianDouble(rndBP, start, end);
						//bp= (int) rdiNebuBP.nextGaussian(start+ (nuLen/ 2d), Math.min(nuLen/4d, 100));
		/*			double pp= rndBP.nextDouble();
					
					bp= (int) ((len/2d)+ (Math.sqrt(Math.abs(pp- 0.5d))* (Math.sqrt(0.5)* len/ 2d)));
					double pp2= (0.5d- (-2d* (pp- 0.5d)* (pp- 0.5d)));
					bp= (int) (len* pp2);
					bp= (pp> 0.5? len- bp: bp);
		*/			
					// bp= (int) (pp* len); // uniform
					int min= Math.min(bp- start+ 1, end- bp);
					int max= len- min;
		//			bp+= start;
		
					double p= 1d- Math.exp(-Math.pow(min/ 200d, 2d));
					if (p< 1d) {
			            double r= rndBreak.nextDouble();
			            if (r> p) {   // does not break (!)	            	
			    			updateMedian(nuLen);
			    			cumuLen+= nuLen;
			    			cs.replace(0, start);
			    			cs.replace(1, end);
			            	processFragNot(start, end, nuLen, id);
			            	return;
			            }
			            
					}
					
					nuLen= bp- start+ 1;
					updateMedian(nuLen);
					cs.replace(0, start);
					cs.replace(1, bp);
					cumuLen+= nuLen;
					incLinesWrote();
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
					
					nuLen= end- bp;
					updateMedian(nuLen);
					cs.replace(0, bp+ 1);
					cs.replace(1, end);
					cumuLen+= nuLen;
					incLinesWrote();
					++newMols;
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
						
					
				}

		void processFrag110228_newPaper(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
		            
					//double lambda= nebu? settings.getLambda():settings.getLambda();	// nebu / 2d

					// x ~ C (d-d_0)^delta
					
					int bp=
						//(int) sampleWeibull(rndBP, scale, shape);
						(int) Math.round((len- 2)* rndBP.nextDouble());	// avoid border pos
					
					int minLen= Math.min(bp, len- bp);
					double pb= Math.pow(Math.max((minLen- 100d) / (double) lastMaxLen, 0), 3d);
					double r= rndBreak.nextDouble();
					if (r> pb) {	// does not break
						updateMedian(len);
		    			cumuLen+= len;
		    			maxLen= Math.max(maxLen, len);
		    			cs.replace(0, start);
		    			cs.replace(1, end);
		            	processFragNot(start, end, len, id);
						return;
					}
						
					
		//			double pb= 
		//				1d- getWeibullProb(len, 300d, 2);
						//1d- (200d/ len);
			
					if (bp+ 1> len- 1|| bp< 0)
						System.currentTimeMillis();

					bp+= start;
					int nuLen= bp- start+ 1;
					updateMedian(nuLen);
					if (start< 0|| start> bp)
						System.currentTimeMillis();
					cs.replace(0, start);
					cs.replace(1, bp);
					cumuLen+= nuLen;
					maxLen= Math.max(nuLen, maxLen);
					incLinesWrote();
					//++newMols;	// one is the currMol
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
					
					nuLen= end- bp;
					updateMedian(nuLen);
					if (bp+ 1< 0|| bp+ 1> end)
						System.currentTimeMillis();
					cs.replace(0, bp+ 1);
					cs.replace(1, end);
					cumuLen+= nuLen;
					maxLen= Math.max(nuLen, maxLen);
					incLinesWrote();
					++newMols;
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
					
				}

		void processRT110217(int start, int end, int len, ByteArrayCharSequence id) {
		
					int howmany= 0;
					if (Fragmenter.this.settings.getRtMode().equals(FluxSimulatorSettings.PAR_RT_MODE_RANDOM)) {
						double poissonMean= (len* 2d)/ (settings.getRTmaxLen()+ settings.getRTminLen()); 
							// (len/ (double) medLen);
						//poissonMean*= 1000;
						if (poissonMean> 0) {
							try {
								howmany= (int) Math.round(rndHowMany.nextPoisson(poissonMean));	// TODO Poisson..
								//System.err.println("OK: mean "+ poissonMean+ ", howmany "+ howmany);
							} catch (Exception e) {
								System.err.println("mean "+ poissonMean+ ", howmany "+ howmany);
								e.printStackTrace();
							}
						}
						if (howmany== 0)
							++howmany;	// guarantee one priming event
						if (howmany> where.length)
							howmany= where.length;
						for (int i = 0; i < howmany; i++) {
							// add weighted biases here 
							/*if (mapTxPosition!= null) {						
								double[] fragCDF= getCDF(mapTxPosition.get(id), start, end);
								double r= rtRndWhere.nextDouble();
								int p= Arrays.binarySearch(fragCDF, r);
								if (p< 0)
									p= -(p+1);
								where[i]= p;
							} else
							*/
								where[i]= 
								//start+ (int) (rtRndWhere.nextDouble()* len);	// skip Math.round() and -1
								//end- (int) Math.exp(-(rtRndWhere.nextDouble()* len));
								end- (int) Math.sqrt(rtRndWhere.nextDouble()* len);
						}
						Arrays.sort(where,0,howmany);
						
					} else if (Fragmenter.this.settings.getRtMode().equals(FluxSimulatorSettings.PAR_RT_MODE_POLY_DT)) {
						howmany= 1;
						int truLen= 0;
		//				int p= ByteArrayCharSequence.indexOf(id, (char) BYTE_SEP_LC_TX, 0, id.length());
		//				ByteArrayCharSequence tid= id.subSequence(p, id.length());
						truLen= settings.getProfiler().getLength(id);	// tid
						int lenPA= end- truLen;
						// rtRndWhere.nextDouble()* len
						if (lenPA< 0) // TODO kill after chr@tid introduced
							where[0]= end;
						else 
							where[0]= truLen+ (int) (rtRndWhere.nextDouble()* lenPA)- 1;	// Math.round() f* slow
		//				if (where[0]< start)
		//					System.currentTimeMillis();
					}
		
					// fall-off
					int lastLo= Integer.MAX_VALUE;
					int lastInt= -1;
					for (int i = howmany-1; i >= 0; --i) {
						
						if (where[i]>= lastLo) {
							double r= fiftyFiftyRnd.nextDouble();
							int dist= Math.min(where[i]- start, where[lastInt]-where[i]);								
							// 500= length that is unlikely to be unwound
							double val=  Math.exp((-1d)*dist/ settings.getRTminLen());
							if (r< val) {	// can unwind, settings.getMaxRTLen()
								where[i]= Integer.MIN_VALUE;
								continue;
							} else 
								starts[lastInt]= where[i];
						}
						
						int diff= settings.getRTmaxLen()- settings.getRTminLen();
						int foffPos= settings.getRTminLen();
						if (diff!= 0)
							foffPos+= rndPolPt.nextInt(diff); 
						
						//foffPos= Math.min(where[i]- start+ 1,foffPos);	
						starts[i]= Math.max(where[i]- foffPos+ 1, start);
						lastLo= starts[i];
						lastInt= i;
					}
					
					//int[] myWhere= where, myStarts= starts;
					--newMols;	// substract the original one	
					
					id= id.cloneCurrentSeq();	// id invalid after replace operation !
					addFragCount(id, (long) (howmany));
					int nowStart= start, nowEnd= end;
					for (int i = 0; i < howmany; i++) {
						if (where[i]== Integer.MIN_VALUE)// || where[i]- starts[i]< settings.getMinRTLen())	// || where[i]- starts[i]< settings.getMinRTLen()
							continue;
						++newMols;
						int nuStart= starts[i], nuEnd= where[i]; 
		//				if (nuStart< start|| nuEnd> end)
		//					System.currentTimeMillis();
						assert(nuStart>= start&& nuEnd<= end);
		//				int len1= where[i]- starts[i]+ 1;
						
						// plot
						if (Fragmenter.this.plotter!= null)
							Fragmenter.this.plotter.plot(nuStart, nuEnd, 0, id);	
						
		//				if (nuStart> nuEnd)
		//					System.currentTimeMillis();
						assert(nuStart<= nuEnd);
						
						// write ! may be iterated multiple times due to howmany
						if (nuStart!= nowStart) {
							cs.replace(0, nuStart);
							nowStart= nuStart;
						}
						if (nuEnd!= nowEnd) {
							cs.replace(1, nuEnd);
							nowEnd= nuEnd;
						}
						incLinesWrote();
						rw.writeLine(cs, fos);
					}
				}

		void processFrag110223_Weibull(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
		            
					// 110223_dynamicWeibull
			
					//double lambda= nebu? settings.getLambda():settings.getLambda();	// nebu / 2d
					
					if (len<= 1) {	// does not break
						updateMedian(len);
		    			cumuLen+= len;
		    			cs.replace(0, start);
		    			cs.replace(1, end);
		            	processFragNot(start, end, len, id);
						return;
					}
			
			
					// dynamic weibull
					double scale= 200d; 	//avgLength/ 2;
					double shape= 2d; // 5d- roundCtr;
					
					int bp= 
						//(int) Math.round((len- 2)* rndBP.nextDouble());	// avoid border pos
						(int) sampleWeibull(rndBP, scale, shape);
					
		//			double pb= 
		//				1d- getWeibullProb(len, 300d, 2);
						//1d- (200d/ len);
					
					if (bp+ 1> len- 1|| bp< 0) {
						System.currentTimeMillis();
						updateMedian(len);
		    			cumuLen+= len;
		    			cs.replace(0, start);
		    			cs.replace(1, end);
		            	processFragNot(start, end, len, id);
						return;
					}
					
/*					double p0= getWeibullProb(len, scale, shape);
					double p1= getWeibullProb(bp+ 1, scale, shape);
					double p2= getWeibullProb(len-bp, scale, shape);
					double pb= Math.min(p1, p2)/ p0;
					double r= rndBreak.nextDouble();
					if (r> pb) {	// does not break
						updateMedian(len);
		    			cumuLen+= len;
		    			cs.replace(0, start);
		    			cs.replace(1, end);
		            	processFragNot(start, end, len, id);
						return;
					}
*/			
					double ff= rndBreak.nextDouble();
					if (ff>= 0.5)
						bp= len- bp- 2;
					
					bp+= start;
					int nuLen= bp- start+ 1;
					updateMedian(nuLen);
					if (start< 0|| start> bp)
						System.currentTimeMillis();
					cs.replace(0, start);
					cs.replace(1, bp);
					cumuLen+= nuLen;
					incLinesWrote();
					//++newMols;	// one is the currMol
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
					
					nuLen= end- bp;
					updateMedian(nuLen);
					if (bp+ 1< 0|| bp+ 1> end)
						System.currentTimeMillis();
					cs.replace(0, bp+ 1);
					cs.replace(1, end);
					cumuLen+= nuLen;
					incLinesWrote();
					++newMols;
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
					
				}

		void processFrag110225delta(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
		            
					//double lambda= nebu? settings.getLambda():settings.getLambda();	// nebu / 2d
		
					// x ~ C (d-d_0)^delta
					double pb= Math.pow(Math.max(len- 100d, 0), 3d) / (double) lastMaxLen;
					
					
					double r= rndBreak.nextDouble();
					if (r> pb) {	// does not break
						updateMedian(len);
		    			cumuLen+= len;
		    			maxLen= Math.max(maxLen, len);
		    			cs.replace(0, start);
		    			cs.replace(1, end);
		            	processFragNot(start, end, len, id);
						return;
					}
					
					int bp= 
						(int) Math.round((len- 2)* rndBP.nextDouble());	// avoid border pos
						//(int) sampleWeibull(rndBP, scale, shape);
					
		//			double pb= 
		//				1d- getWeibullProb(len, 300d, 2);
						//1d- (200d/ len);
			
					if (bp+ 1> len- 1|| bp< 0)
						System.currentTimeMillis();
		
					bp+= start;
					int nuLen= bp- start+ 1;
					updateMedian(nuLen);
					if (start< 0|| start> bp)
						System.currentTimeMillis();
					cs.replace(0, start);
					cs.replace(1, bp);
					cumuLen+= nuLen;
					maxLen= Math.max(nuLen, maxLen);
					incLinesWrote();
					//++newMols;	// one is the currMol
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
					
					nuLen= end- bp;
					updateMedian(nuLen);
					if (bp+ 1< 0|| bp+ 1> end)
						System.currentTimeMillis();
					cs.replace(0, bp+ 1);
					cs.replace(1, end);
					cumuLen+= nuLen;
					maxLen= Math.max(nuLen, maxLen);
					incLinesWrote();
					++newMols;
					rw.writeLine(cs, fos);	// id is invalid now
					if (nuLen<= settings.getFragNBlambda())
						++tgtMols;
					if (len<= settings.getFragNBlambda())
						--tgtMols;
					
				}

		void processFrag110211(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
				            
							double lambda= nebu? settings.getFragNBlambda():settings.getFragNBlambda();	// nebu / 2d
		
							if (!nebu) {
								// first check whether frag gets broken
								//double val= 
									// linear
									// lambda/ len;
									//(len- lambda)/ lambda;
									//(len- lambda)/ (2* lambda);
									
									// polynomial
									//Math.pow((len- lambda)/ lambda, 2d);	// len> lambda, P(b)
								
								
									// exponential
									//1d- Math.exp((lambda- len)/ lambda);	// len> lambda, P(b)
									//1d- Math.exp(-Math.pow(len/ lambda, 2d));	// len> lambda, P(b)
									//1d- Math.pow(Math.exp(-(len/ lambda)), 2d);	// len> lambda, P(b)
									
									//1d- Math.pow(len/ lambda, -0.5d);
									//Math.log(1d+ Math.pow(len/ lambda, 0.5));
									//Math.sqrt(Math.log(1d+len/ lambda));
									//Math.log(Math.sqrt(1d+ len/lambda));
									//1d- Math.exp(-len);
									//1d- (1d/ Math.pow(len, len/lambda));
									
									//1d- (lambda/ Math.pow(len, 1));	// friday
									
		//							double val= 1d;
		//							if (len< lambda)
		//								val= 1d/ Math.log10(lambda- len);
								
								double val= 
									Math.log10(len/ (lambda));	//len> lambda must be able to not break
									//Math.log10(len/ (lambda/ 2d));
									
								if (val< 1d) {
						            double r= rndBreak.nextDouble();
						            if (r> val) {   // does not break (!)
						            	processFragNot(start, end, len, id);
						            	return;
						            }
								} 
									
							}
							
				            
				            // else
							int howMany= (len> 1?1:0);
							if (multiBreaks) {
								howMany= (int) rndHowMany.nextPoisson(len/ lambda); // nr events
								if (howMany== 0&& len> 1)
									howMany= 1;
								else if (howMany>= where.length)
									howMany= where.length- 1;
							}
							if (howMany== 0) {
								processFragNot(start, end, len, id);
								return;
							}
							
							determineBPos(id, start, end, where, howMany, nebu);
							
							// determine len_min, len_max
							int minLen= Integer.MAX_VALUE, maxLen= Integer.MIN_VALUE;
							for (int i = 0; i < howMany; i++) {
								int nuLen= where[i]- (i== 0? start: where[i- 1])+ 1;
								if (nuLen< minLen)
									minLen= nuLen;
								if (nuLen> maxLen)
									maxLen= nuLen;
							}
							int lastLen= end- where[howMany- 1]+ 1;
							if (lastLen< minLen)
								minLen= lastLen;
							if (lastLen> maxLen)
								maxLen= lastLen;
		
							// nebu break condition, after bp has been choosen
							if (nebu) {
								double val= 
									//1d- Math.exp(-Math.pow(minLen/(double) lambda, 2d));
									// nebu
									//1d- Math.exp(-Math.pow(minLen/ lambda, 2d));	// len> lambda, P(b)
									1d- Math.exp(-(minLen/ lambda));	// len> lambda, P(b)
								
								
									//1d- (lambda/ (2* Math.pow(maxLen, 2)));
									//1d- ((lambda/ 2d)/ maxLen);
									
									//1d- Math.pow(maxLen, -(len/ lambda));
									//Math.log10(lambda)/ Math.log10(len);
									//Math.exp(-Math.pow(maxLen/ lambda, 2d));
									//Math.pow()
								double r= rndBreak.nextDouble();
								if (r> val) {
									processFragNot(start, end, len, id);
									return;
								}
							}
		
							
							if (howMany> 1) 
								Arrays.sort(where, 0, howMany);
							int lastVal= start;
							if (howMany> 1)
								id= id.cloneCurrentSeq();	// gets invalid with replace
							int nowStart= start, nowEnd= end;
							addFragCount(id, (long) (howMany));
							tstNewMols+= howMany;
							for (int i = 0; i < howMany; i++) {
								if (where[i]== lastVal)
									continue;	// 0-length frag, 2x hydrolized at same point
								int nuLen= where[i]- lastVal;	// (lastVal- 1)+ 1
								if (processFragDiss(nuLen))
									continue;
		
								++newMols;
								lenSum+= nuLen;
								++lenNb;
								if (nuLen<= Fragmenter.this.settings.getFragNBlambda())
									++Fragmenter.this.tgtMols;
								
								// plot
								int nuEnd= where[i]- 1;	// floor(bp)
								if (Fragmenter.this.plotter!= null) 
									Fragmenter.this.plotter.plot(lastVal, nuEnd, -1, id);
				
								// write
								if (lastVal!= nowStart) {
									cs.replace(0, lastVal);
									nowStart= lastVal;
								}
								if (nuEnd!= nowEnd) {
									cs.replace(1, nuEnd);
									nowEnd= nuEnd;
								}
								cumuLen+= nowEnd- nowStart+ 1;
								assert(lastVal<= nuEnd);
								incLinesWrote();
								rw.writeLine(cs, fos);	// id is invalid now
								lastVal= where[i];
							}
							// last one, always-even if lastVal== end
							int nuLen= end- where[howMany- 1]+ 1;
							if (processFragDiss(nuLen))
								return;
							if (where[howMany- 1]!= nowStart) 
								cs.replace(0, where[howMany- 1]);
							if (nowEnd!= end) 
								cs.replace(1, end);
		//					if (where[howMany- 1]> end)
		//						System.currentTimeMillis();
							assert(where[howMany- 1]<= end);
							cumuLen+= nuLen;
							incLinesWrote();
							rw.writeLine(cs, fos);	// id is invalid now
							if (nuLen<= settings.getFragNBlambda())
								++tgtMols;
							if (len<= settings.getFragNBlambda())
								--tgtMols;
						}

		boolean out1= true, out2= true, out3= true;
		
		/**
		 * Implements a uniform random fragmentation model according to 
		 * Tenchov and Yanev.
		 * 
		 * @param nebu
		 * @param start
		 * @param end
		 * @param len
		 * @param id
		 */
		void processFrag(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
				            
			// get parameters
			double d0= settings.getFragURd0();
			assert(d0>= 1); // problem with integer breakpoints, when fragment size << 1 !
			double delta= getFragURdelta(len);  
			double eta= getFragUReta();
			//eta= 170;
			
			// [2a] dmax= eta((delta- 1)/delta)^(1/ delta)
			double dmax= eta* Math.pow((delta- 1)/ delta,1/ delta);
			// expectancy E of the Weibull distribution W(delta,eta)
			// [3a] E= d0+ eta*gamma(1/delta+ 1)
			double E= d0+ eta* Math.exp(Gamma.logGamma(1d+ 1d/delta));
			double D2= Math.pow(eta, 2)* 
							(Math.exp(Gamma.logGamma((2d/ delta)+ 1d))- 
							Math.pow(Math.exp(Gamma.logGamma((1d/ delta)+ 1d)), 2d));
			
			// DEBUG
			if (out1&& len< 1000) {
				System.err.println(len+"\t"+delta+"\t"+E+"\t"+D2+"\t"+dmax);
				out1= false;
				System.currentTimeMillis();
			} else if (out2&& len> 1200&& len< 1500) {
				System.err.println(len+"\t"+delta+"\t"+E+"\t"+D2+"\t"+dmax);
				out2= false;
				System.currentTimeMillis();
			} else if (out3&& len> 2000) {
				System.err.println(len+"\t"+delta+"\t"+E+"\t"+D2+"\t"+dmax);
				out3= false;
				System.currentTimeMillis();
			}
			
			// determine n, the number of fragments (i.e. (n-1) breakpoints)
			double nn= ((double) len)/ E;			
			int n= (int) nn;
				// (int) Math.round(nn+ (rndBreak.nextGaussian()));	 
				// (int) rndImpl.nextPoisson(nn);		// too variable
			double nr= nn- n;
			double r= rndBreak.nextDouble();
			if (r<= nr)
				++n;

			// molecule does not break
		    if (n<= 1|| len<= 1|| (n*d0)>= len) {	
				updateMedian(len);
    			cumuLen+= len;
    			cs.replace(0, start);
    			cs.replace(1, end);
            	processFragNot(start, end, len, id);
				return;
			}
		    
			// uniformly cut (n-1) times unit space
		    double [] x = new double[n];
		    for(int i= 0; i< (n-1); ++i) 
		      x[i]= rndBreak.nextDouble();	  
		    x[x.length- 1]= 1;	// last breakpoint is end
			
		    // get fragment lengths (in unit space)
		    Arrays.sort(x);
		    for (int i= (n- 1);i> 0; --i)
		      x[i]-= x[i-1];

		    // compute c, transform to molecule space
		    float sum = 0;
		    for(int i= 0; i< x.length; i++) {
		      sum+=Math.pow(x[i], 1/delta);
		    }  
		    double c = Math.pow((len - n* d0)/ sum, -delta); 	
		    for (int i= 0; i< n; i++)
		      x[i]=  d0 + Math.pow(x[i]/ c, (1/delta));

		    double dsum= 0;
		    for (int i = 0; i < n; i++) {
		    	int nuStart= start+ (int) Math.ceil(dsum); 
		    	dsum+= x[i];
				int nuEnd= start+ (i== (n-1)? len: (int) Math.floor(dsum));
		    	//double frac= dsum/ len;
				int nuLen= (nuEnd- nuStart)+ 1;
				if (nuLen< 0)
					System.currentTimeMillis();
				cs.replace(0, nuStart);
				cs.replace(1, nuEnd);
				cumuLen+= nuLen;
				updateMedian(nuLen);
				incLinesWrote();
				rw.writeLine(cs, fos);	// id is invalid now
			}
		    assert(Math.round(dsum)== len);
		}

		void processFrag_110316_LAST(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
		            
					// iterative WEibull
			
					//double lambda= nebu? settings.getLambda():settings.getLambda();	// nebu / 2d
					
					if (len<= 1) {	// does not break
						updateMedian(len);
		    			cumuLen+= len;
		    			cs.replace(0, start);
		    			cs.replace(1, end);
		            	processFragNot(start, end, len, id);
						return;
					}
			
			
					// dynamic weibull
					double scale= 200d; 	//avgLength/ 2;
					double shape= 1.5d; // 5d- roundCtr;
					
					int bp= 
						//(int) Math.round((len- 2)* rndBP.nextDouble());	// avoid border pos
						1+ (int) sampleWeibull(rndBP, scale, shape);
					
		//			double pb= 
		//				1d- getWeibullProb(len, 300d, 2);
						//1d- (200d/ len);
		
					while (bp< len- 1) {
						
						double ff= rndBreak.nextDouble();
						if (ff>= 0.5) {
							cs.replace(0, end- bp+ 1);
							cs.replace(1, end);
							incLinesWrote();
							rw.writeLine(cs, fos);	// id is invalid now
							end-= bp;
						} else {
							cs.replace(0, start);
							cs.replace(1, start+ bp- 1);
							incLinesWrote();
							rw.writeLine(cs, fos);	// id is invalid now
							start+= bp;
						}
						
						// update
						++newMols;
						len-= bp;
						cumuLen+= bp;
						updateMedian(bp);
						bp= 1+ (int) sampleWeibull(rndBP, scale, shape);
					}
					cs.replace(0, start);
					cs.replace(1, end);
					cumuLen+= len;
					updateMedian(len);
					incLinesWrote();
					rw.writeLine(cs, fos);	// id is invalid now
		
				}

		void processFragEma1(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
				            
			if (len<= 1) {	// does not break
				updateMedian(len);
				cumuLen+= len;
				cs.replace(0, start);
				cs.replace(1, end);
		    	processFragNot(start, end, len, id);
				return;
			}
		
			// Ema's Model
			float c= 1f, delta= 3f, d0= 200;
			int n= Math.round(len/ d0);	// TODO Poisson
			
		    // uniformly cut in hyperspace
		    float vol =  c* (float) Math.pow(len,delta);
		    float [] fragments = new float [n+1], x = new float [n];
		    for(int i= 0; i< (n-1); ++i) 
		      x[i]= rndBreak.nextFloat()* vol;  
		    x[x.length- 1]= vol;	// last breakpoint is end
		    
		    Arrays.sort(x);	// transform to fragment lengths (in hyperspace)
		    for (int i= (n- 1);i> 0; --i)
		      x[i]-= x[i-1];
		
		    // transform fragment lengths to linear space
		    for (int i= 0; i< n; i++)
		      x[i]=  d0 + (float) Math.pow(x[i]/c,(1/delta));
		    double sum= 0;
		    for (int i = 0; i < n; i++) {
		    	sum+= x[i];
		    	double frac= sum/ len;
		    	int lo= (int) Math.floor(frac);
		    	start= (i== 0? start: end+ 1); // TODO overflow proof ?!
				end= start+ (i== (n-1)? len: lo);
				
				int nuLen= (end- start)+ 1;
				cs.replace(0, start);
				cs.replace(1, end);
				cumuLen+= nuLen;
				updateMedian(nuLen);
				incLinesWrote();
				rw.writeLine(cs, fos);	// id is invalid now
			}
		}

		void processFragIter(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
						            
					// Ema's Model
					float delta= 3f, d0= 1;
					int n= (int) Math.round(len/ 190d);	// TODO Poisson
		//			if (n== 0)
		//				n= 1;
					
					if (n> 0) {
						n= (int) rndHowMany.nextPoisson(n);
						n= 2;
					} else
						n= 0;
					
					//n= (int) Math.ceil(len/ 190d)+ 1;
					
					//n= (int) rndHowMany.nextPoisson(len/ 250d);
					//n= (int) Math.round(2+ rndHowMany.nextUniform(0, 1));
				    float [] x = new float [n];
		
				    if (n== 0|| len<= 1) {	// does not break
						updateMedian(len);
		    			cumuLen+= len;
		    			cs.replace(0, start);
		    			cs.replace(1, end);
		            	processFragNot(start, end, len, id);
						return;
					}
			
					
					
					//int l = molecule.length();
				    //System.err.printf("molecule length:%d\n",l);
				    //float[] uni_cuts = uniform_cuts (n);
		
					// uniformly cut in hyperspace
				    //float vol =  c* (float) Math.pow(len,delta);
				    for(int i= 0; i< (n-1); ++i) 
				      x[i]= rndBreak.nextFloat();	// * vol  
				    x[x.length- 1]= 1;	// last breakpoint is end
				    Arrays.sort(x);	// transform to fragment lengths (in hyperspace)
				    for (int i= (n- 1);i> 0; --i)
				      x[i]-= x[i-1];
		
				    // compute c, derive eta
				    float sum = 0;
				    for(int i= 0; i< x.length; i++) {
				      sum+=Math.pow(x[i], 1/delta);
				    }  
				    float c = (len - n* d0)/ sum; 	// n+ 1 ???
				    c= (float) Math.pow(c, -delta);
				    float eta = (float)Math.pow((1/(c*n)),(1/delta));
				    
		
				    // transform fragment lengths to linear space
				    for (int i= 0; i< n; i++)	// TODO (n- 1)???
				      x[i]=  d0 + (float) Math.pow(x[i]/ c, (1/delta));
		
				    double dsum= 0;
				    for (int i = 0; i < n; i++) {
				    	int nuStart= start+ (int) Math.ceil(dsum); 
				    	dsum+= x[i];
						int nuEnd= start+ (i== (n-1)? len: (int) Math.floor(dsum));
				    	double frac= dsum/ len;
						
						int nuLen= (nuEnd- nuStart)+ 1;
						cs.replace(0, nuStart);
						cs.replace(1, nuEnd);
						cumuLen+= nuLen;
						updateMedian(nuLen);
						incLinesWrote();
						rw.writeLine(cs, fos);	// id is invalid now
					}
				}

		void processNebu(boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {

			// parameters:
			// pb: M, lambda, len
			// bp: Sigma, length
			double lambda= Fragmenter.this.settings.getFragNBlambda();
			double sigma= Fragmenter.this.settings.getFragNBsigma();
			double M= Fragmenter.this.settings.getFragNBm();
			int recDepth= Fragmenter.this.nebuRecursionDepth;
			double C= Fragmenter.this.nebuC;
			if (fragments== null)
				fragments= new int[(int) Math.pow(2, recDepth)];
			Arrays.fill(fragments, -1);
			fragments[0]= len;
			fragmentNb= 1;
			
			// now break it!
			for (int i = 0; i < recDepth; ++i) {
				for (int j = 0; fragments[j]> 0; ++j) {

					// breakpoint location
					// N(length/2,sigma)= (N(0,1)*sigma*(length/2)+ length/2
					int L= fragments[j]; 
					double rr= rndBP.nextGaussian();
					int bp= (int) ((rr* sigma
							* ((L-1)/2d))+ (L-1)/2d);	// bp index [0;L[						
					//bp= (int) rdiNebuBP.nextGaussian(len/ 2d, len/ 4d);
					bp= (int) nextGaussianDouble(rndBP, 0, L- 1);
					
					// breaking probability (pb)
					// pb= 1- exp^(-((x-C)/lambda)^M)
					int minL= (bp< (L- bp)? bp+1 : L-bp-1);
					double pb= minL< C? 0d: 1d- Math.exp(-Math.pow((minL- C)/lambda, M));
					double r= rndBreak.nextDouble();
					if (r> pb)
						continue;
					
					// fragment j breaks
					int rest= fragmentNb- j- 1;
					if (rest> 0) {
						assert(j+ 2+ rest<= fragments.length);
						System.arraycopy(fragments, j+ 1, fragments, j+ 2, rest);
					}
					assert(bp+1> 0&& L- bp- 1> 0);
					fragments[j]= bp+ 1;
					fragments[j+1]= L- bp- 1;	// L- (bp+1)
					++fragmentNb;
					++j;
				}
				
			}
			
			// write result to disk
			for (int j = 0; fragments[j]> 0; ++j) {
				cumuLen+= fragments[j];
				cs.replace(0, start);
				start+= fragments[j];
				cs.replace(1, start- 1);
				rw.writeLine(cs, fos);	// id is invalid now
				
				// update stats
				incLinesWrote();
				if (fragments[j]<= settings.getFragNBlambda())
					++tgtMols;
				if (len<= settings.getFragNBlambda())
					--tgtMols;
			}
			assert(start== len);
			
		}

		
		void processFilterRejection(int start, int end, int len, ByteArrayCharSequence id, 
				AbstractDistribution[] d, boolean probDistr) {
			
			// get (possibly cumulative) probability for length being in result
			double plen= 0d;
			for (int i = 0; i < d.length; i++) {
				double p= (probDistr? d[i].getP(len): d[i].getRelFreq(len));
				plen+= d[i].getWeight()* p;
			}
			
			// Bernoulli trial
			if (plen> 1|| rndGel.nextDouble()< plen) {
				incLinesWrote();
				rw.writeLine(cs, fos);
				if (Fragmenter.this.plotter!= null) 
					Fragmenter.this.plotter.plot(start, end, len, id);
				addFragCount(id, 1l);
			} else
				--newMols;
		}

		/**
		 * see <code>http://personal.strath.ac.uk/gary.koop/extra_material_on_metropolis.pdf</code>,
		 * <code>http://www.ps.uci.edu/~markm/numerical_methods/Metropolis%96Hastings%20algorithm%20-%20Wikipedia,%20the%20free%20encyclopedia.pdf</code>
		 * @param start
		 * @param end
		 * @param len
		 * @param id
		 */
		void processFilterHastings(int start, int end, int len, ByteArrayCharSequence id) {
			
			// Metropolis/Hastings/Ema
			
			if (len< gelSizeMin|| len> gelSizeMax) {
				--newMols;
				return;
			}
			
			// first value always accepted (?)
			int pos= getBin(len), pos2= getOriBin(len);
			if (lastLen< 0) {
				lastLen= len;
				lastP= gelProb[pos];
				
				lastQ= oriProb[pos2];
				return;
			}
			
			double p= gelProb[pos],
					q= oriProb[pos2];
			
			double a1= p/ lastP;
			//double a2= dnorm(lastLen, len)/ dnorm(len, lastLen);
			double a2= lastQ/ q;
			double a= a1* a2;
			
			// accept 
			if (a>= 1|| rndGel.nextDouble()<= a) {
				lastLen= len;
				lastP= p;
				incLinesWrote();
				rw.writeLine(cs, fos);
				if (Fragmenter.this.plotter!= null) 
					Fragmenter.this.plotter.plot(start, end, len, id);
				addFragCount(id, 1l);
			} else
				--newMols;
		
		}
	}
	
	
	// 2.85 - <0.5% 
	private static final double CUT_OFF_GAUSSIAN_VAL= 2.85f, TWICE_CUT_OFF_GAUSSIAN_VAL= 2* CUT_OFF_GAUSSIAN_VAL, NORM_FACTOR= 1d/ Math.sqrt(2d* Math.PI), E_POW_12= Math.exp(-1d/ 2d);
	
	public static final String PAR_SAMPLE_REJECTION= "RJ", PAR_SAMPLE_ACCEPTANCE= "AC", PAR_SAMPLE_MH= "MH";
	public static final char FILTER_DISTRIBUTION_NORMAL= 'N', FILTER_DISTRIBUTION_UNIFORM= 'U', FILTER_DISTRIBUTION_WEIBULL= 'W';
	
	/**
	 * Parses command line string and returns corresponding subsampling mode identifier.
	 * @param s string that describes sampling method
	 * @return sampling mode identifier
	 */
	public static byte parseFilterSampling(String s) {
		s= s.toUpperCase();
		if (s.equals(PAR_SAMPLE_REJECTION))
			return MODE_FILT_REJ;
		if (s.equals(PAR_SAMPLE_ACCEPTANCE))
			return MODE_FILT_ACC;
		if (s.equals(PAR_SAMPLE_MH))
			return MODE_FILT_MH;
		return MODE_NOT_INITED;
	}
	
	/**
	 * Parses the distribution argument, either a path to a file with an empirical   
	 * description of the function, or a set of analytically described functions. 
	 * @param s string describing the distribution
	 * @return a (set of) distributions as initialized from the argument
	 */
	public static AbstractDistribution[] parseFilterDistribution(String s, double min, double max, int nrBins, boolean fragFile) {
		 
		File f= new File(s);
		if (f.exists())
			try {
				return new AbstractDistribution[] {EmpiricalDistribution.create(f, min, max, nrBins, fragFile)};
			} catch (Exception e) {
				e.printStackTrace();
				return null;
			}
		
		// no file
		String[] ss= s.split("+");
		AbstractDistribution[] d= new AbstractDistribution[ss.length];
		for (int i = 0; i < ss.length; i++) {
			d[i]= parseFilterDistributionFunction(ss[i]);
		}
		
		return d;
	}
	
	
	
	private static AbstractDistribution parseFilterDistributionFunction(String s) {
		
		// find parameter block delimiters
		int p= s.indexOf('('), q= s.indexOf(')');
		assert(p> 0&& q>p);
		
		// parse parameters
		String[] ss= s.substring(p+ 1, q).split(",");
		double[] par= new double[ss.length];
		for (int i = 0; i < par.length; i++) 
			par[i]= Double.parseDouble(ss[i]);
		
		// nature of function
		AbstractDistribution d= null;
		if (s.charAt(p- 1)== FILTER_DISTRIBUTION_NORMAL)
			d= (par.length== 1? new NormalDistribution(par[0]): new NormalDistribution(par[0],par[1]));
		else if (s.charAt(p- 1)== FILTER_DISTRIBUTION_UNIFORM)
			; // TODO
		else if (s.charAt(p- 1)== FILTER_DISTRIBUTION_WEIBULL)
			; // TODO
		
		// weight (in sums of functions)
		if (p> 1) {
			double f= Double.parseDouble(s.substring(0, p- 1));
			d.setWeight(f); 
		}
		
		return d;
	}

	/**
	 * min<= r <= max
	 * @param random
	 * @param min
	 * @param max
	 * @return
	 */
	//public static RandomDataImpl rndINebu= new RandomDataImpl();
	public static strictfp double nextGaussianDouble(Random random, double min, double max) {
		double rdm= 3d;	// gaussian value, stddev 1
		while (rdm< -CUT_OFF_GAUSSIAN_VAL|| rdm> CUT_OFF_GAUSSIAN_VAL)
			rdm= random.nextGaussian();
			//rdm= rndINebu.nextGaussian(0d, 1d);
		double mid= ((double) min)+ (max- min)/ 2f;
		double realValue= mid+ (rdm* (max-mid)/ CUT_OFF_GAUSSIAN_VAL);	// 0..CUT_OFF_GAUSSIAN = mid..max
		assert(realValue>= min&& realValue<= max);
		return realValue;
	}
	
	public static void main(String[] args) {
		try {
			Random rnd= new Random();
			PrintWriter p= new PrintWriter("testWeibull.txt");
			for (int i = 0; i < 1000; i++) {
				double x= sampleWeibull(rnd, 200d, 2);
				p.println(x);
			}
			p.flush();
			p.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * 
	 * @param rnd
	 * @param lamda scale
	 * @param k shape
	 * @return
	 */
	public static double sampleWeibull(Random rnd, double lamda, double k) {
		double ret= lamda* Math.pow((-1d* Math.log(1d- rnd.nextDouble())), (1d/k));
		return ret;
	}
	
	

	
	/**
	 * 
	 * @param lamda scale
	 * @param k shape
	 * @return
	 */
	public static double getWeibullProb(double x, double lamda, double k) {
		double ret= (k/ lamda)* Math.pow(x/ lamda, k-1)* Math.exp(-Math.pow(x/ lamda, k));
			//lamda* Math.pow((-1d* Math.log(1d- rnd.nextDouble())), (1d/k));
		return ret;
	}

	
	boolean stop= false;
	long molInit, molRT, molFrag, molFilt;
	FluxSimulatorSettings settings;
	Hashtable<CharSequence,Long> mapFrags;
	ProgressablePlotter plotter;
	HashMap<CharSequence, double[]> mapWeightSense= null, mapWeightAsense= null;
	
	/**
	 * @deprecated
	 */
	HashMap<CharSequence, float[]> mapPWMsense= null, mapPWMasense= null;	
	PWM pwmSense, pwmAsense;
	HashMap<CharSequence, CharSequence> mapTxSeq= null;
	int mapKmerLen= -1, mapMotifLen= -1, mapMotifFirst= -1;
	byte mode= MODE_NOT_INITED;
	
	// gel bins
	static int GEL_NB_BINS_LENGTH= 100, GEL_NB_BINS_GC= 10;
	int gelSizeMin= Integer.MAX_VALUE, gelSizeMax= -1, gelSizeBin= -1;
	double gelGCMin= Double.MAX_VALUE, gelGCMax= -1d, gelGCbin= -1d;
	double[] gelProb= null;
	
	protected int getBin(int len) {
		int x= (int) ((GEL_NB_BINS_LENGTH- 1)* ((len- gelSizeMin)/ (double) (gelSizeMax- gelSizeMin)));
		if (gelGCbin> 0) 
			x*= GEL_NB_BINS_GC;
		return x;
	}
	
	protected int getBin(double gc) {
		return (int) ((gc- gelGCMin)/ gelGCbin);
	}
	
	protected int getBin(int len, double gc) {
		
		int x= 0, y= 0;
		assert(gelSizeBin> 0);
		x= getBin(len);
		if (gelGCbin> 0) 
			y= getBin(gc);
		
		return (x+ y);
	}
	
	
	
	protected void initBins(File f, double[] pp, int min, int max) {
		try {

			HashMap<Integer, DoubleVector> map2D= new HashMap<Integer, DoubleVector>();
			BufferedReader buffy= new BufferedReader(new FileReader((f)));
			long total= 0; String[] ss;
			for (int i = 0; i < pp.length; i++) 
				pp[i]= 0d;
			for (String s= null; (s= buffy.readLine())!= null;++total) {
				ss= s.split("\\s");
				int pos1= Integer.parseInt(ss[0]);
				int pos2= Integer.parseInt(ss[1]);
				int len= pos2- pos1+ 1;
				if (len< min|| len> max)
					continue;
				int x= getOriBin(len);				
				++pp[x];
				++total;
			}
			buffy.close();
			for (int i = 0; i < pp.length; i++) 
				pp[i]/= total;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private int getOriBin(int len) {
		int x= (int) ((oriProb.length- 1)* ((len- oriMin)/ (double) (oriMax- oriMin)));
		return x;
	}

	public Fragmenter(FluxSimulatorSettings settings) {
		this.settings= settings;
	}

	public Fragmenter(FluxSimulatorSettings settings, byte mode) {
		this(settings);
		this.mode= mode;
	}

	static RandomDataImpl rndImplGelLen= new RandomDataImpl();
	public static final int GEL_RESOLUTION_LIMIT= 100000;	// max frag that moves
	public static int getSegregatedLength(int len, boolean approx) {
/*		int nuLen= (int) rndImplGelLen.nextGaussian(len, 1000d/ len);	// -(len-1)
		if (1== 1)
			return nuLen;
*/			
		if (len< GEL_RESOLUTION_LIMIT) {
			int passes= 100; // GEL_RESOLUTION_LIMIT- len;
			double rPasses= 0d;
			if (passes>= 1000|| approx) 
				rPasses= rndImplGelLen.nextGaussian(passes, Math.sqrt(passes));
			else 
				rPasses= rndImplGelLen.nextPoisson(passes);

			return (int) (len* rPasses/ passes);
		}
		return len;
	}
	
	int deltaFiltMin= -1, deltaFiltMax= -1;
	public int getDeltaFiltMin() {
		if (deltaFiltMin< 0) 
			initDeltaRange();

		return deltaFiltMin;
	}
	public int getDeltaFiltMax() {
		if (deltaFiltMax< 0) 
			initDeltaRange();

		return deltaFiltMax;
	}
	private void initDeltaRange() {
		if (settings.isFilter()) {
			int irange= (int) ((settings.getFiltMax()- settings.getFiltMin()+ 1)* 1.5d);
			deltaFiltMin= Math.max(settings.getFiltMin()- irange, 0);
			deltaFiltMax= settings.getFiltMax()+ irange;
		}
	}
	
	public void run() {
		
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
			System.err.println("[LIBRARY] creating the cDNA libary");
		
		// TODO adapt init() and check whether all parameters are ok
//		if (!init())
//			throw new RuntimeException();
		
		if (settings.getFrgFile().exists()&& checkFrgFile())
			return;
		
		// REUSE frg file
		//if (settings.getFrgFile().exists()) {
			//FileHelper.move(settings.getFrgFile(), tmpFile, null);
		//} else {
			tmpFile= writeInitialFile();
			if (tmpFile== null)
				return; // FluxSimulator.exit(-1);
		//}

			
		// do it
		if (settings.isFragB4RT()) {
			rndBreak= new Random();
			if (settings.isFragment()) {
				byte mode= MODE_NOT_INITED;
				if (settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_NB_DEPRECATED)
						|| settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_NB)) {
					mode= MODE_NEBU;
					rdiNebuBP= new RandomDataImpl();
				} else if (settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_UR_DEPRECATED)
						|| settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_UR)) {
					mode= MODE_FRAG;
					rdiNebuBP= new RandomDataImpl();
					rndBP= new Random();
				} 
				if (isStop()|| !process(mode)) {
					stopProcessors();
					return;
				}
			}
			if (settings.isRT()) {
				initRTpar(settings.getRTminGC(), settings.getRTmaxGC());
				if (settings.getRTmotif()!= null) { 
					getMapTxSeq();
					//getWeights(settings.getRTmotif());
					try {
						pwmSense= PWM.create2(settings.getRTmotif());
						for (int i = 0; i < 100; i++) 
							pwmSense.multiply();
						pwmSense.makePDF();
						mapWeightSense= getMapWeight(getMapTxSeq(), pwmSense);
						pwmSense.invert();
						pwmAsense= pwmSense;
						mapWeightAsense= getMapWeight(getMapTxSeq(), pwmAsense);
					} catch (Exception e) {
						e.printStackTrace();
					}

				}
				if (!process(MODE_RT)) {
					return;
				}
			}
			
		// start with RT
		} else {
			if (settings.isRT()) {
				initRTpar(settings.getRTminGC(), settings.getRTmaxGC());
				if (settings.getRTmotif()!= null) {
					getMapTxSeq();
					getWeights(settings.getRTmotif());
				}
				if (!process(MODE_RT)) {
					return;
				}
			}
			rndBreak= new Random();
			if (settings.isFragment()) {
				byte mode= MODE_NOT_INITED;
				if (settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_NB_DEPRECATED)
						|| settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_NB)) {
					mode= MODE_NEBU;
					rdiNebuBP= new RandomDataImpl();
				} else if (settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_UR_DEPRECATED)
						|| settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_UR)) {
					mode= MODE_FRAG;
					rndBP= new Random();
				} else if (settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_EZ)) {
					mode= MODE_FRAG_EZ;
					rnd1= new Random();
					rnd2= new Random();
					rnd3= new Random();
					try {
						getMapTxSeq();
						pwmSense= PWM.create2(settings.getFragEZmotif());
						mapWeightSense= getMapWeight(getMapTxSeq(), pwmSense);
						pwmSense.invert();
						pwmAsense= pwmSense;
						mapWeightAsense= getMapWeight(getMapTxSeq(), pwmAsense);
					} catch (Exception e) {
						e.printStackTrace();
					} 
				}
				if (isStop()|| !process(mode)) {
					stopProcessors();
					return;
				}
			}
		}
		
		
		// size selection
		if (settings.isFilter()) {
			
			//if (settings.getFileFilterDistr()!= null) {
			//	initGelBins(new File(settings.getFilterDistribution()));
			//}
			System.err.println("Initializing Selected Size distribution");
			//Log.progressStart("Initializing Selected Size distribution");
			filterDist= parseFilterDistribution(settings.getFilterDistribution(), Double.NaN, Double.NaN, GEL_NB_BINS_LENGTH, false);
			//Log.progressFinish();
			EmpiricalDistribution eDist= ((EmpiricalDistribution) filterDist[0]); 			
			
			if (filterDist[0] instanceof EmpiricalDistribution) {
				
				System.err.println("Initializing Current Size distribution");
				//Log.progressStart("Initializing Current Size distribution");	// TODO to be done during earlier steps
				originalDist= parseFilterDistribution(tmpFile.getAbsolutePath(), eDist.getMin(), eDist.getMax(),
						eDist.getBins().length, true)[0];
				// Log.progressFinish();
				
				eDist.normalizeToPrior((EmpiricalDistribution) originalDist);
			}
			
			byte mode= parseFilterSampling(settings.getFilterSampling());
			if ((!process(mode))) {
				stopProcessors();
				return;
			}
		}
		
		if (isStop()) {
			if (tmpFile!= null&& tmpFile.exists())
				tmpFile.delete();
			if (tmpWriteFile!= null&& tmpWriteFile.exists())
				tmpWriteFile.delete();
			stopProcessors();
			return;
		}
		if (!writeFinalFile())
			throw new RuntimeException();
		if (isStop()) {
			stopProcessors();
			return;
		}
		if (!FluxSimulatorSettings.appendProfile(settings, FluxSimulatorSettings.PRO_COL_NR_FRG, mapFrags))	// writeProFile()
			throw new RuntimeException();
		stopProcessors();
	}

	private GFFReader getGFFReader() {
		GFFReader gffReader = new GFFReader(settings.getRefFile().getAbsolutePath());	
		try {
			if (!gffReader.isApplicable()) {
				File refFile= gffReader.createSortedFile();
				if (refFile== null)
					return null;
				settings.setRefFile(new File(settings.getProFile().getParent()+ File.separator+ refFile.getName()));
				if (!refFile.equals(settings.getRefFile())) {
					if (!FileHelper.move(refFile, settings.getRefFile(), null))
						settings.setRefFile(refFile);
				}
				gffReader= new GFFReader(settings.getRefFile().getAbsolutePath());
			}
			gffReader.setSilent(true);
			gffReader.setStars(true);
			
		} catch (Exception e) {
			return null;
		}
	
		return gffReader;
		
	}
	
	
	private HashMap<CharSequence, double[]> getMapWeight(HashMap<CharSequence, CharSequence> mapSeq, PWM pwm) {
		HashMap<CharSequence, double[]> map= new HashMap<CharSequence, double[]>(mapSeq.size(), 1f);
		Iterator<CharSequence> iter= mapSeq.keySet().iterator();
		while (iter.hasNext()) { 
			CharSequence id= iter.next();
			CharSequence seq= getMapTxSeq().get(id);
			double[] a= new double[seq.length()];
			for (int p = 0; p < a.length; ++p) {
				double pb= pwm.apply(seq, p);
				a[p]= pb;
				assert(!(Double.isNaN(pb)|| Double.isInfinite(pb)));
			}
			map.put(id, a);
		}
		return map;
	}
	private HashMap<CharSequence, CharSequence> getMapTxSeq() {

		if (mapTxSeq == null) {
			mapTxSeq= new HashMap<CharSequence, CharSequence>(10000);
			GFFReader reader= getGFFReader();		
			if (settings.genDir!= null)
				Graph.overrideSequenceDirPath= settings.genDir.getAbsolutePath();
			try {
				reader.read();
				for (Gene[] g; (!stop)&& (g= reader.getGenes())!= null; reader.read()) {				
					for (int i = 0; (!stop)&& i < g.length; i++) {
						for (int j = 0; (!stop)&& j < g[i].getTranscripts().length; j++) {
							Transcript t= g[i].getTranscripts()[j];
							String s= t.getSplicedSequence();	// TODO check chr pre-loading
							ByteArrayCharSequence combID= new ByteArrayCharSequence(g[i].getGeneID());
							combID.append((byte) FluxSimulatorSettings.SEP_LOC_TID);
							combID.append(t.getTranscriptID());
							mapTxSeq.put(combID, s);
						}
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		return mapTxSeq;
	}
	
	private void getWeights(File filePWM) {
		try {
			pwmSense= PWM.create2(filePWM);
			mapWeightSense= getMapWeight(getMapTxSeq(), pwmSense);
			pwmSense.invert();
			pwmAsense= pwmSense;
			mapWeightAsense= getMapWeight(getMapTxSeq(), pwmAsense);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	

	
	private boolean writeFinalFile() {
		
		Log.progressStart("Copying results");

		if (FileHelper.move(tmpFile, settings.getFrgFile(), null)) 
			return true;
		return false;
	}

	private boolean checkFrgFile() {
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			System.err.print("\tchecking file "+settings.getFrgFile().getAbsolutePath()+" ");
			System.err.flush();
		}
		
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
			System.err.println(" OK");
		
		return false;
	}
	

	boolean init() {
		
		if (settings.getProFile()== null) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				System.err.println("\t[OOPS] no input for nebulizing");
			return false;
		} else {
			int nbTx= 0;
			try {
				FileHelper.countLines(settings.getProFile().getCanonicalPath());
			} catch (IOException e) {
				if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
					e.printStackTrace();
				return false;
			}
			mapFrags= new Hashtable<CharSequence,Long>(nbTx, 1f);
		}
		
		if (settings.getTmpDir()== null) {
			settings.setTmpDir(new File(System.getProperty("java.io.tmpdir")));
		}
		if ((!settings.getTmpDir().exists())|| (!settings.getTmpDir().canWrite())) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\t[UHUPS] there is something wrong with the scratch directory:");
				try {
					System.err.println("\t"+settings.getTmpDir().getCanonicalPath());
				} catch (IOException e) {
					if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS) 
						e.printStackTrace();
				}
			}
			return false;
		} else {
			try {
				tmpFile= File.createTempFile(FluxSimulatorSettings.TMP_PFX, FluxSimulatorSettings.TMP_SFX);
				tmpWriteFile= File.createTempFile(FluxSimulatorSettings.TMP_PFX, FluxSimulatorSettings.TMP_SFX);
			} catch (IOException e) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("\t[NOFILE] problems creating temporary files");
				if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS) 
					e.printStackTrace();
				return false;
			}

		}
		
		if (	settings.getFragNBlambda()< 0|| Double.isNaN(settings.getFragNBlambda())
				|| (settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_NB_DEPRECATED)&& (settings.getThold()< 0|| Double.isNaN(settings.getThold())))
				// || settings.getSigma()< 0
				//|| (!settings.getProFile().getParentFile().canWrite()) // killed: NPException when run on file without leading path
				) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.print("\t[TOOLITTLE] there are parameters missing: ");
//				if (settings.getSigma()< 0)
//					System.err.print(settings.PAR_FRAG_SIGMA+" ");
				if (settings.getFragNBlambda()< 0)
					System.err.print(settings.PAR_FRAG_NB_LAMBDA+" ");
				if (settings.getThold()< 0)
					System.err.print(settings.PAR_FRAG_NB_THOLD+" ");
				System.err.println("\n");
			}
			return false;
			
		} else {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\t"+settings.PAR_FRAG_NB_LAMBDA+"\t"+settings.getFragNBlambda());
				// System.err.println("\t"+settings.PAR_FRAG_SIGMA+"\t"+settings.getSigma());
				System.err.println("\t"+settings.PAR_FRAG_NB_THOLD+"\t"+settings.getThold());
				System.err.println("\t"+settings.PAR_FRAG_SUBSTRATE_DEPRECATED+"\t"+settings.isFragB4RT());
				System.err.println("\t"+settings.PAR_RT_MIN+"\t"+settings.getRTminLen());
				System.err.println("\t"+settings.PAR_RT_MAX+"\t"+settings.getRTmaxLen());
				System.err.println("\t"+settings.PAR_RT_MODE+"\t"+settings.getRtMode());
				System.err.println("\t"+settings.PAR_FILT+"\t"+settings.isFilter());
				System.err.println("\t"+settings.PAR_FILT_MIN+"\t"+settings.getFiltMin());
				System.err.println("\t"+settings.PAR_FILT_MAX+"\t"+settings.getFiltMax());

				try {
					System.err.println("\t"+FluxSimulatorSettings.PAR_PRO_FNAME+"\t"+settings.getProFile().getCanonicalPath());
					System.err.println("\t"+FluxSimulatorSettings.PAR_TMP_FNAME+"\t"+settings.getTmpDir().getCanonicalPath());
				} catch (IOException e) {
					;	// :)
				}
				System.err.println();
			}
		}
		
		return true;
	}
	
	Random rndBreak; 
	Random rndBP= new Random();
	Random rndNebu= new Random();
	Random rndDisappear= new Random();
	Random rndGel= new Random();
	Random rnd1= new Random(), rnd2= new Random(), rnd3= new Random();
	
	//GaussianRndThread rndNebu;	
	RandomDataImpl rdiNebuBP;
	long newMols= 0, currMols= 0, tgtMols= 0,
		tstNewMols= 0, tstCurrMols= 0;	// un-synchronized, but should be ok -- are estimates anyway
	
	/**
	 * Parameter to adapt breaking probability distribution to 0.5 for (1.5*lambda).
	 */
	double nebuC= -1;
	
	/**
	 * Number of Iterations for recursive nebulization.
	 */
	int nebuRecursionDepth= -1;
	
	File tmpFile, tmpWriteFile;
	
	// tmporary variables re-used (nebulization, RT)
	int[] index1= null, index2= null; 

	
	boolean rtranscribe() {
		
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.print("\treverse transcribing +");
		
		return process(MODE_RT);	
	}
	
	boolean nebulize() {
		
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			System.err.print("\tnebulizing +");
			System.err.println("lamda "+settings.getFragNBlambda()
					//+", sigma "+settings.getSigma()
					+", thold "+settings.getThold());
		}
		
		return process(MODE_NEBU);
	}

	boolean filter() {
		
		return process(MODE_FILT_REJ);
	}

	
	
	ThreadedQWriter qwriter= null;
	
	
	RandomDataImpl rndHowMany= new RandomDataImpl();
	
	Random rtRndWhere= new Random();
	Random rndPolPt= new Random(), fiftyFiftyRnd= new Random();
	long lenSum;
	int lenNb;
	public static final int RT_MIN_FRAGMENT_LENGTH= 200;
	private static final String TAB= "\t";
	/**
	 * @deprecated
	 * @return
	 */
	private ThreadedQWriter getQWriter() {
/*		ThreadedQWriter qwriter= new ThreadedQWriter(tmpWriteFile) {
			
			@Override
			public void writeAll() {
				while ((!q.isEmpty())) {	// DEADLOCKS: && (!Fragmenter.this.isStop()) because flush() waits for q to become empty
					String s= (String) q.poll();
					if (s== null)
						continue;
					++currMols;
					//String[] token= s.split("\\s");	
					StringTokenizer toki= new StringTokenizer(s, TAB);
					//int start= Integer.parseInt(token[0]), end= Integer.parseInt(token[1]);
					int start= Integer.parseInt(toki.nextToken()), end= Integer.parseInt(toki.nextToken()); 
					//assert(start< end);
					if (start>= end) {
						System.currentTimeMillis();
						continue;
					}
					int len= end- start+ 1;
					String id= toki.nextToken(); // token[2];
					
					if (Fragmenter.this.mode== Fragmenter.MODE_FILT) {
						
						Random rnd= new Random();
						RandomDataImpl rndGel= new RandomDataImpl();
						//rndGel.nextGaussian(mu, sigma);
						//rndGel.nextPoisson(mean)
						if (len>= Fragmenter.this.settings.getFiltMin()&&
								len<= Fragmenter.this.settings.getFiltMax())
							try {
								writer.write(s+"\n");
								if (Fragmenter.this.plotter!= null) 
									Fragmenter.this.plotter.plot(start, end, id);
								setFragCount(id, 1l);
							} catch (IOException e) {
								e.printStackTrace();
							}
						else
							--newMols;
						// increaseFragCount(token[2]);
							
							
					} else if (Fragmenter.this.mode== Fragmenter.MODE_NEBU) {
						// (int) Math.ceil(rdiNebu.nextGaussian(start+ len/2d, settings.sigma* settings.lambda))
						int bp= (int) nextGaussianDouble(rndNebu, start+1, end);	//	Math.ceil() f* slow 
//						while (bp>= end|| bp<= start+1)
//							bp= (int) Math.ceil(rdiNebu.nextGaussian(start+ len/2d, settings.sigma* settings.lambda));
						int len2= Math.min(bp- start, end- bp+ 1);
						int len1=(bp-start==len2)?end-bp+1:bp-start;
//						double val= Math.exp((-1d)* Math.pow(len2- settings.lambda, 2)/ 
//								Math.pow(settings.sigma* settings.lambda, 2));	// TODO either len or len2
//						double val= Math.exp((-1d)* Math.pow(Math.max(len- settings.lambda,0), 2)/ 
//								Math.pow(settings.lambda, 2));	// TODO either len or len2
//						double val= Math.exp((-1d)* (len/ settings.lambda));	// TODO either len or len2
						//Fragmenter.this.settings.thold= 0.05;
						
						double val= Math.exp((-1d)* Math.pow(Math.max(0,len- settings.getLambda()),2));	// TODO either len or len2
						//double val= Math.exp((-1d)* Math.pow(len2/ settings.lambda,2));	// TODO either len or len2
						val= 1- val;	// 0 f. len -> settings.lamda
						double r= rndBreak.nextDouble();
							
						try {
							if (r> val) {	// does not break (!)
								//lengthV.add(end- start+ 1);
								lenSum+= len;
								maxLen= Math.max(maxLen,len);
								minLen= Math.min(minLen,len);
								++lenNb;
								writer.write(s+"\n");
								increaseFragCount(id);
								if (Fragmenter.this.plotter!= null)
									Fragmenter.this.plotter.plot(start, end, id);
								if (len>= Fragmenter.this.settings.getFiltMin()&&
										len<= Fragmenter.this.settings.getFiltMax())
									++Fragmenter.this.tgtMols;
							} else {
								++newMols;
//								lengthV.add(end- bp+ 1);
//								lengthV.add(bp- start);
								
								lenSum+= len1;
								lenSum+= len2;
								lenNb+= 2;
								maxLen= Math.max(maxLen,len1);
								maxLen= Math.max(maxLen,len2);
								minLen= Math.min(minLen,len1);
								minLen= Math.min(minLen,len2);
								
								writer.write(getNebuLine(start,bp-1,id)+"\n");
								writer.write(getNebuLine(bp,end,id)+"\n");
								if (len1>= Fragmenter.this.settings.getFiltMin()&&
										len1<= Fragmenter.this.settings.getFiltMax())
									++Fragmenter.this.tgtMols;
								if (len2>= Fragmenter.this.settings.getFiltMin()&&
										len2<= Fragmenter.this.settings.getFiltMax())
									++Fragmenter.this.tgtMols;
								//setFragCount(token[2], 2l);
								if (Fragmenter.this.plotter!= null) {
									Fragmenter.this.plotter.plot(start, bp-1, id);
									Fragmenter.this.plotter.plot(bp, end, id);
								}								
							}						
						} catch (Exception e) {
							e.printStackTrace();
						}
						
						
					} else if (Fragmenter.this.mode== Fragmenter.MODE_FRAG) {

	                     double k= 3.4;	// 4
	                        
                        // mean= la* (ln(2))^(1/k)
                        double la= (len/2d)/ Math.pow(Math.log(2), (1d/k)); 
                        
                        // weibull random sampling
                        int bp= -1; 
                        while (bp> end|| bp< start)
                            bp=start+ rndFrag.nextInt(end-start+1);	// start+ Math.round(sampleWeibull(rndFrag, la,k))
                        
                        int len2= Math.min(bp- start, end- bp+ 1);
                        int len1=(bp-start==len2)?end-bp+1:bp-start;
                        // double val= Math.exp((-1d)* Math.pow(Math.max(0,len- settings.lambda),2));    // TODO either len or len2
                        double val= 1d;
                        if (len2> settings.getLambda())	// 
                            val= Math.pow(len2- settings.getLambda(), -2);    // 0 f. len -> settings.lamda 
                        val= 1-val;
                        double r= rndBreak.nextDouble();
                            
                        try {
                            if (r> val) {    // does not break (!)
								//lengthV.add(end- start+ 1);
                            	lenSum+= len;                            	
                            	++lenNb;
								maxLen= Math.max(maxLen,len);
								minLen= Math.min(minLen,len);
								writer.write(s+"\n");
								increaseFragCount(id);
								if (Fragmenter.this.plotter!= null)
									Fragmenter.this.plotter.plot(start, end, id);
								if (len>= Fragmenter.this.settings.getFiltMin()&&
										len<= Fragmenter.this.settings.getFiltMax())
									++Fragmenter.this.tgtMols;
							} else {
								if (len< 100)
									System.currentTimeMillis();
								++newMols;
								//lengthV.add(end- bp+ 1);
								//lengthV.add(bp- start);
								lenSum+= len1;
								lenSum+= len2;
								lenNb+= 2;
								maxLen= Math.max(maxLen,len1);
								maxLen= Math.max(maxLen,len2);
								minLen= Math.min(minLen,len1);
								minLen= Math.min(minLen,len2);
								writer.write(getNebuLine(start,bp-1,id)+"\n");
								writer.write(getNebuLine(bp,end,id)+"\n");
								if (len1>= Fragmenter.this.settings.getFiltMin()&&
										len1<= Fragmenter.this.settings.getFiltMax())
									++Fragmenter.this.tgtMols;
								if (len2>= Fragmenter.this.settings.getFiltMin()&&
										len2<= Fragmenter.this.settings.getFiltMax())
									++Fragmenter.this.tgtMols;
								//setFragCount(token[2], 2l);
								if (Fragmenter.this.plotter!= null) {
									Fragmenter.this.plotter.plot(start, bp-1, id);
									Fragmenter.this.plotter.plot(bp, end, id);
								}
							}						
						} catch (Exception e) {
							e.printStackTrace();
						}
						
						// so nicht
//						IntVector vec= new IntVector();
//						int bp= start;
//						while (bp< end) {
//							int i= (int) Math.round(rdi.nextUniform(settings.getFiltMin(), settings.getFiltMax()));
//							vec.addElement(i);
//							bp+= i;
//						}
//						bp-= vec.elementAt(vec.size()-1);
//						vec.removeElement(vec.size()-1);
//						
//						try {
//							if (vec.size()== 0) {
//								writer.write(s+"\n");
//								if (Fragmenter.this.plotter!= null)
//									Fragmenter.this.plotter.plot(start, end, token[2]);
//								if (len>= Fragmenter.this.settings.filtMin&&
//										len<= Fragmenter.this.settings.filtMax)
//									++Fragmenter.this.tgtMols;
//							} else {
//								--newMols;
//								int offs= (int) Math.round((len- bp)/ 2d);
//								for (int i = 0; i < vec.size(); offs+= vec.elementAt(i), i++) {
//									++newMols;
//									int newEnd= offs+vec.elementAt(i)-1;
//									writer.write(getNebuLine(offs,newEnd,token[2])+"\n");
//									if (Fragmenter.this.plotter!= null) 
//										Fragmenter.this.plotter.plot(offs, newEnd, token[2]);
//								}
//								Fragmenter.this.tgtMols+= vec.size();
//							}
//						} catch (Exception e) {
//							e.printStackTrace();
//						}
						
					
					} else if (Fragmenter.this.mode== Fragmenter.MODE_RT) {
						
						int howmany= 0;
						if (Fragmenter.this.settings.getRtMode().equals(FluxSimulatorSettings.PAR_RT_MODE_RANDOM)) {
							double poissonMean= (len* 2d)/ (settings.getMaxRTLen()+ settings.getMinRTLen()); 
								// (len/ (double) medLen);
							//poissonMean*= 1000;
							if (poissonMean> 0)
								howmany= (int) Math.round(rtRndHowMany.nextPoisson(poissonMean));	
							if (howmany== 0)
								++howmany;	// guarantee one priming event
							if (howmany> where.length)
								howmany= where.length;
							for (int i = 0; i < howmany; i++) 
								where[i]= start+ (int) Math.round(rtRndWhere.nextDouble()* len)- 1;
							Arrays.sort(where,0,howmany);
							
						} else if (Fragmenter.this.settings.getRtMode().equals(FluxSimulatorSettings.PAR_RT_MODE_POLY_DT)) {
							howmany= 1;
							int truLen= settings.getProfiler().getLength(id);
							int lenPA= end- truLen;
							// rtRndWhere.nextDouble()* len
							where[0]= truLen+ (int) (rtRndWhere.nextDouble()* lenPA)- 1;	// Math.round() f* slow
						}

						// fall-off
						int lastLo= Integer.MAX_VALUE;
						int lastInt= -1;
						for (int i = howmany-1; i >= 0; --i) {
							
							if (where[i]>= lastLo) {
								double r= fiftyFiftyRnd.nextDouble();
								int dist= Math.min(where[i]- start, where[lastInt]-where[i]);								
								// 500= length that is unlikely to be unwound
								double val=  Math.exp((-1d)*dist/ settings.getMinRTLen());
								if (r< val) {	// can unwind, settings.getMaxRTLen()
									where[i]= Integer.MIN_VALUE;
									continue;
								} else 
									starts[lastInt]= where[i];
							}
							
							int diff= settings.getMaxRTLen()- settings.getMinRTLen();
							int foffPos= settings.getMinRTLen();
							if (diff!= 0)
								foffPos+= rndPolPt.nextInt(diff); 
							
							//foffPos= Math.min(where[i]- start+ 1,foffPos);	
							starts[i]= Math.max(where[i]- foffPos+ 1, start);
							lastLo= starts[i];
							lastInt= i;
						}
						
						//int[] myWhere= where, myStarts= starts;
						--newMols;	// substract the original one						
						for (int i = 0; i < howmany; i++) {
							if (where[i]== Integer.MIN_VALUE)// || where[i]- starts[i]< settings.getMinRTLen())	// || where[i]- starts[i]< settings.getMinRTLen()
								continue;
							++newMols;
							try {
								int nuStart= starts[i], nuEnd= where[i]; 
								assert(nuStart>= start&& nuEnd<= end);
								writer.write(getNebuLine(nuStart,nuEnd,id)+"\n");
								int len1= where[i]- starts[i]+ 1;
								if (len1>= Fragmenter.this.settings.getFiltMin()&&
										len1<= Fragmenter.this.settings.getFiltMax())
									++Fragmenter.this.tgtMols;
								if (Fragmenter.this.plotter!= null)
									Fragmenter.this.plotter.plot(nuStart, nuEnd, id);							
							} catch (IOException e) {
								if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
									e.printStackTrace();
								if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
									System.err.println("\tError during writing RT file");
								return;
							}
							increaseFragCount(id);
						}						
						
					}
				}	// end while
				if (Fragmenter.this.isStop())
					close();
			}
			
		};
*/
		
		return qwriter;
	}
	
	public static final String MODE_RT_MESSAGE= "Reverse Transcription", MODE_NEBU_MESSAGE= "Nebulization", MODE_FILT_MESSAGE= "Segregating cDNA", MODE_FRAG_MESSAGE= "Fragmentation UR", MODE_FRAG_EZ_MESSAGE= "Enzymatic Digestion";
	
	/**
	 * Default median size after fragmentation, according to 2010 Illumina protocol.
	 */
	private static final double DEFAULT_MED_SIZE = 200;
	
	byte getMode() {
		if (processorPool== null|| processorPool.length== 0)
			return MODE_NOT_INITED;
		return processorPool[0].mode;
	}
	
	void setMode(byte mode) {
		if (processorPool== null|| processorPool.length== 0)
			return;
		processorPool[0].mode= mode;
	}
	
	long prevMols= 0;
	SyncIOHandler2 rw;
	/**
	 * @deprected no longer in use
	 */
	double filtMu, filtSigSquare, filtFac;
	long cumuLen;
	int[] med= new int[3];
	AbstractDistribution originalDist= null;
	AbstractDistribution[] filterDist= null;

	int roundCtr= 0;
	double avgLength= 0;
	private double fragUReta;
	double[] oriProb;
	int oriMin= -1, oriMax= -1;
	FileOutputStream fos;
	boolean process(byte mode) {
		
		if (tmpFile== null)
			return false;
		
		try {
			long t0= System.currentTimeMillis();
			double breakRatio= 1;
			Long longNull= new Long(0);
			double thrTgt= 0.95;
			double tgtFrac= 0d;
			cumuLen= 0;
			int[] allLength= settings.getProfiler().getLen();
			long[] allMos= settings.getProfiler().getMolecules();
			for (int i = 0; i < allLength.length; i++) {
				cumuLen+= allMos[i]* allLength[i];
				prevMols+= allMos[i];
			}
			
			if (mode== MODE_NEBU) 
				initNebulizationParameters(
						settings.getFragNBlambda(),
						settings.getFragNBm(),
						settings.getThold(),
						settings.getProfiler().getMaxMoleculeLength()
				);
			
			// TODO estimate required disk space
			
			Processor[] processors= getProcessorPool(mode, Math.min(settings.getMaxThreads(), 1));
			
			for (roundCtr= 0;(!isStop())&& roundCtr< 1; ++roundCtr, mode= getMode()) {
				
				avgLength= cumuLen/ (double) prevMols; 
				
				//lengthV= new IntVector();
				lastMaxLen= maxLen;
				maxLen= 0;
				cumuLen= 0;
				lenSum= 0l;
				lenNb= 0;
				minLen= 0; maxLen= 0;
				if (mapFrags!= null) 
					mapFrags.clear();	// 20101215 re-init
								
				String msg= null;
				if (mode== MODE_RT)
					msg= MODE_RT_MESSAGE;
				else if (mode== MODE_NEBU) {
					msg= MODE_NEBU_MESSAGE+"-"+Integer.toString(roundCtr+1);
				} else if (mode== MODE_FILT_REJ|| mode== MODE_FILT_ACC|| mode== MODE_FILT_MH) {
					msg= MODE_FILT_MESSAGE;
				} else if (mode== MODE_FRAG) {
					msg= MODE_FRAG_MESSAGE+"-"+Integer.toString(roundCtr+1);
				} else if (mode== MODE_FRAG_EZ)
					msg= MODE_FRAG_EZ_MESSAGE;
				
				if (plotter!= null) {
					plotter.reset(msg);
				}
				Log.progressStart(msg);

				Object[] keys= mapFrags.keySet().toArray();
				for (int i = 0; i < keys.length; i++) 
					mapFrags.put((CharSequence) keys[i], longNull);	// 20101205: changed from cast to string, ClassCastException for ByteArrayCharSequence
				
				currMols= 0; newMols= 0; tgtMols= 0;
				tstNewMols= 0; tstCurrMols= 0;
				cntLinesWr= 0;
				for (int i = 0; i < med.length; i++) 
					med[i]= 0;
				
				// IO 
//				BufferedReader buffy= new BufferedReader(new FileReader(tmpFile));
//				ThreadedQWriter qwriter= getQWriter();
//				qwriter.init();
//				qwriter.start();
//				ThreadedBufferedByteArrayStream buffy= 
//					new ThreadedBufferedByteArrayStream(10* 1024* 1024, fis, true, false);
				//BufferedOutputStream outStream
				//writer= new BufferedOutputStream(new FileOutputStream(tmpWriteFile), 1024* 1024);
				//ThreadedBufferedByteArrayStream inBacs= null;
				FileInputStream fis= new FileInputStream(tmpFile);
				fos= new FileOutputStream(tmpWriteFile);
				rw= new SyncIOHandler2(2);
				rw.addStream(fis, 10* 1024);
				rw.addStream(fos, 10* 1024);				
				for (int i = 0; i < processorPool.length; i++) 
					processorPool[i].setFos(fos);
				if (FluxSimulatorSettings.optDisk) 
					rw.start();
				
				ByteArrayCharSequence cs= new ByteArrayCharSequence(100);
				
				int perc= 0;
				long byteTot= tmpFile.length(), byteNow= 0l;
				//for (String s; (!isStop())&& (s= buffy.readLine())!= null;/*++currMols*/) {
				//for (buffy.readLine(cs); cs.end> 0; buffy.readLine(cs)) {
				while((!isStop())&& rw.readLine(fis, cs)> 0) {
					byteNow+= cs.length()+ 1;	// TODO fs length
					Log.progress(byteNow, byteTot);
					
					++currMols;
					
					cs.resetFind();
					int start= cs.getTokenInt(0);
					int end= cs.getTokenInt(1);
					assert(start<= end);
					int len= end- start+ 1;
//					if (processFragDiss(len))
//						return;
					if (len<= settings.getFragNBlambda())
						++tgtMols;

					ByteArrayCharSequence id= cs.getToken(2);
					addFragCount(id, 1l);
					tstCurrMols+= 1;

					if (mode== MODE_WRITE_INITIAL) {
						processInitial(cs, len, id);	// CHECK
						rw.writeLine(cs, fos);
						continue;
					}
					
					ByteArrayCharSequence ccs= cs.cloneCurrentSeq(); // TODO check whether needed
					if (mode== Fragmenter.MODE_FILT_REJ) 
						processFilterRejection(ccs, start, end, len, id, filterDist, true);
					else if (mode== Fragmenter.MODE_FILT_ACC) 
						processFilterRejection(ccs, start, end, len, id, filterDist, false);
					else if (mode== Fragmenter.MODE_FILT_MH) 
						processFilterMCMC(ccs, start, end, len, id, originalDist, filterDist);
					else if (mode== Fragmenter.MODE_NEBU) 
						processNebu(ccs, true, start, end, len, id);
					else if (mode== Fragmenter.MODE_FRAG) 
						processFrag(ccs, false, start, end, len, id);
						//processFragUniform(false, start, end, len, id);
					else if (mode== Fragmenter.MODE_FRAG_EZ) 
						processFragEnzyme(ccs, start, end, len, id);
					else if (mode== Fragmenter.MODE_RT) {
						//rw.writeLine(cs, 1);
						processRT(ccs, start, end, len, id);
					}	
					
					
				}
				
				rw.close();
				if (FluxSimulatorSettings.optDisk) 
					rw.join();
				
				//fis.close();
				//writer.flush();
				//writer.close();
				
				//buffy.close();
/*				if (isStop())
					qwriter.close();
				qwriter.flush();
				qwriter.close();
*/				
/*				if (mode== MODE_NEBU) {					
					rndNebu.setStop();
					rndNebu.interrupt();
				}
*/				
/*				while (qwriter.isAlive()) // while deadlocks
					try {
						qwriter.join();
					} catch (InterruptedException e) {
						qwriter.close(); // ?? Thread.currentThread().interrupt();
					}
*/
				//Distribution dist= new Distribution(lengthV.toIntArray());
				medLen= lenSum/ (double) lenNb;
				lastMaxLen= maxLen;
				lastMinLen= minLen;
				//lengthV= null;
				
				if (isStop())
					break;
				
				breakRatio= newMols/ (2d* currMols);
				tgtFrac= ((double) tgtMols)/ (currMols+ newMols);
/*				if (mode== MODE_NEBU|| mode== MODE_FRAG)
					System.out.println("bratio "+Double.toString(breakRatio)
							+ ", all "+ (currMols+newMols)
							+ ", in "+ currMols
							+ ", new "+ newMols
							+ ", tgt "+ tgtMols
							+ ", trat "+ tgtFrac
					);
*/					
				//DEBUG
/*								File save= new File("N:\\tmp\\round_"+roundCtr);
								if (save.exists())
									save.delete();									
								FileHelper.copy(tmpFile, save);
*/								
								
				boolean b= tmpFile.delete();
				if (!b)
					throw new IOException("Couldn't delete source");
				b= tmpWriteFile.renameTo(tmpFile);
				if (!b)
					throw new IOException("Couldn't move file");;
								
				Log.progressFinish(null, true);
				long total= 0;
				Iterator iter= mapFrags.keySet().iterator();
				while(iter.hasNext()) {
					Long val= mapFrags.get(iter.next());
					if (val!= null)
						total+= val;
				}
				if (Constants.verboseLevel>= Constants.VERBOSE_NORMAL) 
					System.err.println("\t\t"+(currMols+ newMols)+ " mol, "+total+": "+currMols+","+tstCurrMols+" "+newMols+","+tstNewMols);
				if (Constants.verboseLevel>= Constants.VERBOSE_NORMAL) {
					System.err.println("\t\t"+(currMols+ newMols)+ " mol: in "+ currMols+ ", new "+ newMols+ ", out "+ cntLinesWr);
					System.err.println("\t\tavg Len "+ (cumuLen/ (float) currMols)+ ", maxLen "+ maxLen);
					System.err.println("tgt frac "+tgtFrac);
				}
				if (plotter!= null) {
					plotter.paint();
					plotter.setMolTot(currMols+ newMols);
				}
				prevMols= currMols+ newMols;

			}

			if (!isStop()) {
				/*for (int i = 0; i < processors.length; i++) {
					if (processors[i]== null)
						continue;
					while (processors[i].isAlive()) {
						synchronized (processors[i].lock) { 
							if (processors[i].ready()) {
								processors[i].stop= true;
								//processors[i].interrupt();
								processors[i].lock.notifyAll();
							}
						}
						if (processors[i].stop)
							try {
								processors[i].join();
							} catch (InterruptedException e) {
								; // :)
							}
					}
				}
*/				
				if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
					System.err.println(" OK");
				System.gc();
				return true;
			}
			
		} catch (Exception e) {
			if (Constants.verboseLevel>= Constants.VERBOSE_SHUTUP)
				e.printStackTrace();
		}

		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			if (isStop())
				System.err.print(" FAILED");
			else
				System.err.print(" STOPPED");

		return false;

	}
	
	private void deconvolve(double[] oriProb, double[] gelProb) {
		assert(oriProb.length== gelProb.length);
		double sum= 0d, max= 0d;
		for (int i = 0; i < gelProb.length; i++) {
			gelProb[i]= (oriProb[i]== 0? 0: gelProb[i]/ oriProb[i]);
			sum+= gelProb[i];
			max= (gelProb[i]> max? gelProb[i]: max);
		}
		for (int i = 0; i < gelProb.length; i++) 
			gelProb[i]/= max;
	}

	double rtC= Double.NaN;
	public void initRTpar(double minGC, double maxGC) {
		rtC= minGC* (1.5d+ Math.log(0.5d));
		
	}


	private void initNebulizationParameters(double lambda, double M, double thold,int maxLen) {
		// determine C
		// C~ f(lambda,M), adjust that 1.5 lambda => pb= 0.5
		// C= lambda(1.5- (-ln(0.5))^(1/M))
		// e.g., C=486 f. M=9, lambda=900
		nebuC= lambda* (1.5d- Math.pow(-Math.log(0.5d), 1d/M));
		
		// expected recursion depth
		// p_t<= 0.1
		// tmax= ceil( log2((maxlen-C)/(lambda*(-ln(0.9)^(1/M))) )
		// e.g. len= 1500 -> ceil(0.53), len= 15k -> ceil(4.37)
		// maxLen= 10000; lambda= 900; nebuC= 486;
		nebuRecursionDepth= 
			(int) Math.ceil(Math.log10((maxLen- nebuC)/ (lambda* Math.pow(-Math.log(1d- thold), 1d/M)))/ Math.log10(2));

	}
	
	private double getFragURdelta(double len) {
		
		if (Double.isNaN(settings.getFragURdelta()))
			return Math.max(Math.log10(len), 1);
		return settings.getFragURdelta();
	}
	
	/**
	 * Provides eta ("intensity of fragmentation") of the uniform random fragmentation process;
	 * if no eta has been provided as input parameter, eta is optimized to provide the median 
	 * molecule length a value that corresponds to the median of the subsequently filtered 
	 * values, or constant <code>DEFAULT_MED_SIZE</code>.
	 * @return
	 */
	private double getFragUReta() {
		
		if (Double.isNaN(settings.getFragUReta())) {
			double medLen= settings.getProfiler().getMedMoleculeLength();
			double medDelta= getFragURdelta(medLen);
			double d0= settings.getFragURd0();
			double medFilt= DEFAULT_MED_SIZE; 
			if (settings.isFilter())
				medFilt= 170; //getFiltMedian();
			
			settings.setFragUReta((medFilt- d0)/ Math.exp(Gamma.logGamma(1d+ 1d/medDelta)));
		}

		return settings.getFragUReta();
	}


	Processor[] processorPool;
	/**
	 * @deprecated use Settings
	 * @return
	 */
	private boolean writeProFile() {
		try {
			System.err.print("\tupdating PROfile ");
			System.err.flush();
			BufferedReader buffy= new BufferedReader(new FileReader(settings.getProFile()));
			BufferedWriter wright= new BufferedWriter(new FileWriter(settings.getTmpDir()+File.separator+settings.getProFile().getName()+FluxSimulatorSettings.TMP_SFX));
			String[] token;
			long bytesRead= 0, bytesTotal= settings.getProFile().length();
			int perc= 0, lineCtr= 0;
			String nullStr= Double.toString(0d)+FluxSimulatorSettings.PRO_FILE_TAB+Long.toString(0);
			for (String s= null; (s= buffy.readLine())!= null;++lineCtr) {
				
				bytesRead+= s.length()+ FluxSimulatorSettings.PRO_FILE_CR.length();
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					perc= StringUtils.printPercentage(perc, bytesRead, bytesTotal, System.err);
				
				token= s.split(FluxSimulatorSettings.PRO_FILE_TAB);
				wright.write(s);
				wright.write(FluxSimulatorSettings.PRO_FILE_TAB);
				if (mapFrags.containsKey(token[0])) {
					long absCnt= mapFrags.get(token[0]);
					double relFreq= absCnt/ (double) currMols;
					wright.write(Double.toString(relFreq));
					wright.write(FluxSimulatorSettings.PRO_FILE_TAB);
					wright.write(Long.toString(absCnt));
				} else 
					wright.write(nullStr);
				
				wright.write(FluxSimulatorSettings.PRO_FILE_CR);
				if (lineCtr%1000== 0)
					wright.flush();
			}
			buffy.close();
			wright.flush();
			wright.close();
			
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println(" OK");
			return true;
			
		} catch (Exception e) {
			if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
				e.printStackTrace();
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println(" FAILED");
			return false;
		}
		
	}
	
	private void updateStats(int start, int end, int len, ByteArrayCharSequence id) {
		cumuLen+= len;
		maxLen= Math.max(maxLen, len);
	}
	
	private void increaseFragCount(CharSequence string) {
		addFragCount(string, 1l);
	}
	private synchronized void addFragCount(CharSequence string, Long i) {
//		if (1== 1)
//			return;
		if (mapFrags.containsKey(string)) {
			long otherVal= mapFrags.get(string)+ i;	// 20101215: whyever, dont use in put clause
													// see tstCurMol and tstNewMol divergence
			mapFrags.put(string, otherVal);
		} else
			mapFrags.put(string, i);
		
	}
	
	int cntLinesWr= 0;
	/**
	 * @deprecated
	 */
	private synchronized void incLinesWrote() {
		++cntLinesWr;
	}


	private int minLen, maxLen, lastMinLen, lastMaxLen;
	private double medLen= 0;
	public static final byte BYTE_SEP_LC_TX= '@';
	RandomDataImpl rndTSS;
	Random rndPA, rndPlusMinus;
	int cntMolInit= 0;
	private File writeInitialFile() {
		rndTSS= new RandomDataImpl();
		rndPA= new Random();
		rndPlusMinus= new Random();
		lenSum= 0;	// for RT priming length, here with median
		lenNb= 0;
		cntMolInit= 0;
		long t0= System.currentTimeMillis();
		try {
			//if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				//System.err.print("\tiniting fragments +");
			Log.progressStart("Initializing Fragmentation File");

			//BufferedReader buffy= new BufferedReader(new FileReader(settings.proFile));
			//String[] token;
			//int lineCtr= 0;
			//long bytesRead= 0, totalBytes= settings.proFile.length();
//			for (String s; (!stop)&& (s= buffy.readLine())!= null; ++lineCtr) {
//				bytesRead+= s.length()+ 1;	// line.separator
//				
//				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
//					// perc= MyFormatter.printPercentage(perc, bytesRead, totalBytes, System.err);
//					int frac= (int) (bytesRead* 10d/ totalBytes);
//					if (frac> perc) {
//						perc= frac;
//						Constants.progress.setValue(perc);
//					}
//				}
//				token= s.split("\\s");
//				int nb= Integer.parseInt(token[3]);
//				int len= Integer.parseInt(token[1]);
//				if (lineCtr% 1000== 0)
//					writer.flush();
//			}
//			buffy.close();
			
			FileOutputStream fos= new FileOutputStream(tmpFile);	
			rw= new SyncIOHandler2(1);
			rw.addStream(fos, 10* 1024);
//			if (FluxSimulatorSettings.optDisk)  
//				rw.start();
			Profiler profiler= settings.getProfiler();
			Processor[] processors= getProcessorPool(MODE_WRITE_INITIAL, Math.min(settings.getMaxThreads(), 4));
			for (int i = 0; i < processors.length; i++) 
				processors[i].setFos(fos);
			int perc= 0; molInit= 0; medLen= 0; minLen= Integer.MAX_VALUE; maxLen= Integer.MIN_VALUE;
			for (int i = 0; (!isStop())&& i < settings.getProfiler().getMolecules().length; i++) {
				
				Log.progress(i, settings.getProfiler().getMolecules().length);

				if (profiler.getLen()[i]>= Fragmenter.this.settings.getFiltMin()&&
						settings.getProfiler().getLen()[i]<= Fragmenter.this.settings.getFiltMax())
					++Fragmenter.this.tgtMols;
				if (plotter!= null) 
					plotter.addBase(settings.getProfiler().getCombinedID(i), 
							settings.getProfiler().getLen()[i], 
							(int) settings.getProfiler().getMolecules()[i]);
				
				int origLen= settings.getProfiler().getLen()[i];
				StringBuilder compID= new StringBuilder(profiler.getLocIDs()[i]);
				compID.append(FluxSimulatorSettings.SEP_LOC_TID);
				compID.append(profiler.getIds()[i]);
				String compIDstring= compID.toString();
//				if (tid== null) {
//					if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
//						System.err.println("[ERROR] TID is null in Fragmenter.writeInitialFile()");
//				}
				for (int x = 0; (!isStop())&& x < settings.getProfiler().getMolecules()[i]; x++) {

					++cntMolInit;
					if (processors.length== 1) {
						processors[0].origLen= origLen;
						processors[0].locID= compIDstring; //locID; 	// not in same line !
						processors[0].run();
					} else {
						boolean searching= true;
						while (searching && !stop) {
							for (int j = 0; (!stop)&& j < processors.length; ++j) {
								if (processors[j].ready()) {
									processors[j].origLen= origLen;
									processors[j].locID= compIDstring; //locID; 	// not in same line !
									// deadlocks :(
									synchronized (processors[j].lock) { 
										processors[j].lock.notifyAll();
										//System.err.println("notify "+ccs);
									}
	//								//processors[i].interrupt();
									searching= false;
									break;
								}
							}
						}
					}

					
						//getNebuLine(start,end,settings.getProfiler().getIds()[i])+"\n";
				}
				molInit+= settings.getProfiler().getMolecules()[i];
			}
			for (int i = 0; i < processors.length; i++) {
				processors[i].close();	// 20101205: combine control flow to Processor.stop
/*				processors[i].stop= true;
				processors[i].interrupt();
				processors[i].join();
*/				
				processors[i]= null;
			}
			this.processorPool= null;
			
			rw.close();
			if (FluxSimulatorSettings.optDisk) 
				rw.join();

			prevMols= molInit;
			if (plotter!= null)
				plotter.setMolTot(molInit);
			
			//System.out.println("targets: "+tgtMols);
			Log.progressFinish("OK", true);
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\t"+ cntMolInit+ " mol inited");
			}
			
			medLen= lenSum/ (double) lenNb;		
			lastMinLen= minLen;
			lastMaxLen= maxLen;
			return tmpFile;
		} catch (Exception e) {
			if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
				e.printStackTrace();
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println(" FAILED");
			return null;
		}		
		
	}

	public static String getNebuLine(int start, int end, String ancestor) {
		return getNebuLine(Integer.toString(start), Integer.toString(end), ancestor);
	}
	public static String getNebuLine(String start, String end, String ancestor) {
		return start+ FluxSimulatorSettings.PRO_FILE_TAB+ end+ FluxSimulatorSettings.PRO_FILE_TAB+ ancestor;
	}

	public boolean isStop() {
		return stop;
	}

	public boolean setStop(boolean stop) {
		if (stop)
			return setStop();
		this.stop = stop;
		return true;
	}

	public ProgressablePlotter getPlotter() {
		return plotter;
	}

	public void setPlotter(ProgressablePlotter plotter) {
		this.plotter = plotter;
	}
	
    /**
     * Returns an error message if something is broken or missing and null if everything is fine
     *
     * @return message error message or null
     */
	public String isReady() {
        if(settings == null) return "No Setting specified!";
        if(settings.getProfiler() == null) return "Profiler is not initialized!";
        if(settings.getProfiler().getMolecules() == null) return "Profiler has no molecules!";
        if(settings.getFrgFile() == null) return "No fragmentation file specified!";
		return null;
	}
	
	public boolean isFinished() {
		if (settings.getFrgFile()!= null&& settings.getFrgFile().exists()&& settings.getFrgFile().length()!= 0)
			return true;
		return false;
	}

	public boolean loadStats() {
		try {
			Log.progressStart("initializing library");
			molInit= 0;			
			for (int i = 0; (!isStop())&& i < settings.getProfiler().getMolecules().length; i++) 
				molInit+= settings.getProfiler().getMolecules()[i];
			if (plotter!= null)
				plotter.setMolTot(molInit);
			//BufferedReader buffy= new BufferedReader(new FileReader(settings.getFrgFile()));
			BufferedByteArrayReader buffy= new BufferedByteArrayReader();			
			FileInputStream fos= new FileInputStream(settings.getFrgFile());
			String[] token= null;
			//String s;
			long totBytes= settings.getFrgFile().length(), bytesRead= 0l;
			int perc= 0;
			ByteArrayCharSequence cs= new ByteArrayCharSequence(50);
			for (currMols= 0, newMols= 0; (!isStop())&& (buffy.readLine(fos, cs).length()> 0)/*(s= buffy.readLine())!= null*/; ++currMols) {
				cs.resetFind();
				bytesRead+= cs.length()+ 1;
				Log.progress(bytesRead, totBytes);
				
				if (plotter!= null) {
					int start= cs.getTokenInt(0); 	//Integer.parseInt(token[0]);
					int end= cs.getTokenInt(1); // Integer.parseInt(token[1]);
					plotter.plot(start, end, -1, cs.getToken(2)/*token[2]*/);
				}
			}			
			//buffy.close();
			fos.close();
			Log.progressFinish();
			if ((!isStop())&& plotter!= null) {
				plotter.paint();
				plotter.setMolTot(currMols);
			}
			return true;
		} catch (Exception e) {
			return false; // :)
		}
	}
	
	public long getMoleculeNr() {
		return currMols+ newMols;
	}

	private void stopProcessors() {
		for (int i = 0; processorPool!= null&& i < processorPool.length; i++) {
			if (processorPool[i]== null)
				continue;
			processorPool[i].close();	//20101205: combine control flow to Processor.stop
/*			processorPool[i].stop= true;
			processorPool[i].interrupt();
			while (processorPool[i].isAlive())
				try {
					processorPool[i].join();
				} catch (InterruptedException e) {
					; // :)
				}
*/
			processorPool[i]= null;
		}

	}
	
	public boolean setStop() {
		this.stop = true;
		return true;
	}

	public static int presetThreads= 4;
	private Processor[] getProcessorPool(byte mode, int threads) {		
		if (processorPool == null|| processorPool.length!= threads|| mode!= getMode()) {
			for (int i = 0; processorPool!= null&& i < processorPool.length; i++) {
				processorPool[i].close();
			}
			processorPool = new Processor[threads];
			for (int i = 0; i < processorPool.length; i++) { 
				processorPool[i]= new Processor(mode, threads> 1);
				if (threads> 1)
					processorPool[i].start();
			}			
		}
		
		return processorPool;
	}

		
	boolean process_110407(byte mode) {
			
			if (tmpFile== null)
				return false;
			
			try {
				long t0= System.currentTimeMillis();
				double breakRatio= 1;
				Long longNull= new Long(0);
				double thrTgt= 0.95;
				double tgtFrac= 0d;
				cumuLen= 0;
				int[] allLength= settings.getProfiler().getLen();
				long[] allMos= settings.getProfiler().getMolecules();
				for (int i = 0; i < allLength.length; i++) {
					cumuLen+= allMos[i]* allLength[i];
					prevMols+= allMos[i];
				}
				
				if (mode== MODE_NEBU) 
					initNebulizationParameters(
							settings.getFragNBlambda(),
							settings.getFragNBm(),
							settings.getThold(),
							settings.getProfiler().getMaxMoleculeLength()
					);
				
				// TODO estimate required disk space
				Processor[] processors= getProcessorPool(mode, Math.min(settings.getMaxThreads(), 1));			
				for (roundCtr= 0;(!isStop()) 
						&& (mode== MODE_FRAG&& roundCtr< 1)
							|| ((mode== MODE_FILT_REJ|| mode==MODE_RT)&& roundCtr< 1)
	//					((mode== MODE_NEBU&& breakRatio> settings.getThold())	// && breakRatio> settings.getThold()
	//							|| (mode== MODE_FRAG&& ((cumuLen== -1|| (cumuLen/ (float) currMols)> 400)))		// && roundCtr< 4, tgtFrac< thrTgt||	 
	//							|| ((mode== MODE_RT|| mode== MODE_FILT)&& roundCtr< 1))
						; ++roundCtr, mode= getMode()) {
					
					avgLength= cumuLen/ (double) prevMols; 
					
					//lengthV= new IntVector();
					lastMaxLen= maxLen;
					maxLen= 0;
					cumuLen= 0;
					lenSum= 0l;
					lenNb= 0;
					minLen= 0; maxLen= 0;
					if (mapFrags!= null) 
						mapFrags.clear();	// 20101215 re-init
									
					String msg= null;
					if (mode== MODE_RT)
						msg= MODE_RT_MESSAGE;
					else if (mode== MODE_NEBU) {
						msg= MODE_NEBU_MESSAGE+"-"+Integer.toString(roundCtr+1);
						//rndNebu= new GaussianRndThread(100);
						//rndNebu.start();
					} else if (mode== MODE_FILT_REJ) {
						msg= MODE_FILT_MESSAGE;
						filtSigSquare= ((settings.getFiltMax()- settings.getFiltMin())/ 2d);
						filtMu= settings.getFiltMin()+ filtSigSquare;
						filtFac= 1d/ Math.sqrt(2d* Math.PI* filtSigSquare);
						
					} else if (mode== MODE_FRAG)
						msg= MODE_FRAG_MESSAGE+"-"+Integer.toString(roundCtr+1);
					if (plotter!= null) {
						plotter.reset(msg);
					}
					Log.progressStart(msg);
	
					Object[] keys= mapFrags.keySet().toArray();
					for (int i = 0; i < keys.length; i++) 
						mapFrags.put((CharSequence) keys[i], longNull);	// 20101205: changed from cast to string, ClassCastException for ByteArrayCharSequence
					
					currMols= 0; newMols= 0; tgtMols= 0;
					tstNewMols= 0; tstCurrMols= 0;
					cntLinesWr= 0;
					for (int i = 0; i < med.length; i++) 
						med[i]= 0;
					
					// IO 
	//				BufferedReader buffy= new BufferedReader(new FileReader(tmpFile));
	//				ThreadedQWriter qwriter= getQWriter();
	//				qwriter.init();
	//				qwriter.start();
	//				ThreadedBufferedByteArrayStream buffy= 
	//					new ThreadedBufferedByteArrayStream(10* 1024* 1024, fis, true, false);
					//BufferedOutputStream outStream
					//writer= new BufferedOutputStream(new FileOutputStream(tmpWriteFile), 1024* 1024);
					//ThreadedBufferedByteArrayStream inBacs= null;
					FileInputStream fis= new FileInputStream(tmpFile);
					FileOutputStream fos= new FileOutputStream(tmpWriteFile);
					rw= new SyncIOHandler2(2);
					rw.addStream(fis, 10* 1024);
					rw.addStream(fos, 10* 1024);
					for (int i = 0; i < processorPool.length; i++) 
						processorPool[i].setFos(fos);
					if (FluxSimulatorSettings.optDisk) 
						rw.start();
					
					ByteArrayCharSequence cs= new ByteArrayCharSequence(100);
					
					int perc= 0;
					long byteTot= tmpFile.length(), byteNow= 0l;
					//for (String s; (!isStop())&& (s= buffy.readLine())!= null;/*++currMols*/) {
					//for (buffy.readLine(cs); cs.end> 0; buffy.readLine(cs)) {
					while((!isStop())&& rw.readLine(fis, cs)> 0) {
						byteNow+= cs.length()+ 1;	// TODO fs length
						Log.progress(byteNow, byteTot);
						++currMols;
						if (processors.length== 1) {
							ByteArrayCharSequence ccs= cs.cloneCurrentSeq();
							processors[0].cs= ccs;
							processors[0].run();
						} else {
							boolean searching= true;
							while (searching && !stop) {
								for (int i = 0; (!stop)&& i < processors.length; i++) {
									// 100330: took into sync block
		//							if (processors[i].ready()) {
		//								ByteArrayCharSequence ccs= cs.cloneCurrentSeq();
		//								processors[i].cs= ccs; 	// not in same line !
									// deadlocks :(
									synchronized (processors[i].lock) { 
										if (processors[i].ready()) {
											ByteArrayCharSequence ccs= cs.cloneCurrentSeq();
											processors[i].cs= ccs;
											processors[i].lock.notifyAll();
											searching= false;
										}
										//System.err.println("notify "+ccs);
									}
		//								//processors[i].interrupt();
									break;
								}
		/*						if (searching)
									try {
										Thread.currentThread().sleep(10);
									} catch (InterruptedException e) {
										; // :)
									}
								else
									break;
		*/												
							}
						}
					}
					
					rw.close();
					if (FluxSimulatorSettings.optDisk) 
						rw.join();
					
					//fis.close();
					//writer.flush();
					//writer.close();
					
					//buffy.close();
	/*				if (isStop())
						qwriter.close();
					qwriter.flush();
					qwriter.close();
	*/				
	/*				if (mode== MODE_NEBU) {					
						rndNebu.setStop();
						rndNebu.interrupt();
					}
	*/				
	/*				while (qwriter.isAlive()) // while deadlocks
						try {
							qwriter.join();
						} catch (InterruptedException e) {
							qwriter.close(); // ?? Thread.currentThread().interrupt();
						}
	*/
					//Distribution dist= new Distribution(lengthV.toIntArray());
					medLen= lenSum/ (double) lenNb;
					lastMaxLen= maxLen;
					lastMinLen= minLen;
					//lengthV= null;
					
					if (isStop())
						break;
					
					breakRatio= newMols/ (2d* currMols);
					tgtFrac= ((double) tgtMols)/ (currMols+ newMols);
	/*				if (mode== MODE_NEBU|| mode== MODE_FRAG)
						System.out.println("bratio "+Double.toString(breakRatio)
								+ ", all "+ (currMols+newMols)
								+ ", in "+ currMols
								+ ", new "+ newMols
								+ ", tgt "+ tgtMols
								+ ", trat "+ tgtFrac
						);
	*/					
					//DEBUG
	/*								File save= new File("N:\\tmp\\round_"+roundCtr);
									if (save.exists())
										save.delete();									
									FileHelper.copy(tmpFile, save);
	*/								
									
					boolean b= tmpFile.delete();
					if (!b)
						throw new IOException("Couldn't delete source");
					b= tmpWriteFile.renameTo(tmpFile);
					if (!b)
						throw new IOException("Couldn't move file");;
									
					Log.progressFinish(null, true);
					long total= 0;
					Iterator iter= mapFrags.keySet().iterator();
					while(iter.hasNext()) {
						Long val= mapFrags.get(iter.next());
						if (val!= null)
							total+= val;
					}
					if (Constants.verboseLevel>= Constants.VERBOSE_NORMAL) 
						System.err.println("\t\t"+(currMols+ newMols)+ " mol, "+total+": "+currMols+","+tstCurrMols+" "+newMols+","+tstNewMols);
					if (Constants.verboseLevel>= Constants.VERBOSE_NORMAL) {
						System.err.println("\t\t"+(currMols+ newMols)+ " mol: in "+ currMols+ ", new "+ newMols+ ", out "+ cntLinesWr);
						System.err.println("\t\tavg Len "+ (cumuLen/ (float) currMols)+ ", maxLen "+ maxLen);
						System.err.println("tgt frac "+tgtFrac);
					}
					if (plotter!= null) {
						plotter.paint();
						plotter.setMolTot(currMols+ newMols);
					}
					prevMols= currMols+ newMols;
	
				}
	
				if (!isStop()) {
					/*for (int i = 0; i < processors.length; i++) {
						if (processors[i]== null)
							continue;
						while (processors[i].isAlive()) {
							synchronized (processors[i].lock) { 
								if (processors[i].ready()) {
									processors[i].stop= true;
									//processors[i].interrupt();
									processors[i].lock.notifyAll();
								}
							}
							if (processors[i].stop)
								try {
									processors[i].join();
								} catch (InterruptedException e) {
									; // :)
								}
						}
					}
	*/				
					if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
						System.err.println(" OK");
					System.gc();
					return true;
				}
				
			} catch (Exception e) {
				if (Constants.verboseLevel>= Constants.VERBOSE_SHUTUP)
					e.printStackTrace();
			}
	
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				if (isStop())
					System.err.print(" FAILED");
				else
					System.err.print(" STOPPED");
	
			return false;
	
		}

	protected void initGelBins(File f) {
		try {
	
			HashMap<Integer, DoubleVector> map2D= new HashMap<Integer, DoubleVector>();
			BufferedReader buffy= new BufferedReader(new FileReader((f)));
			int total= 0; String[] ss;
			for (String s= null; (s= buffy.readLine())!= null;++total) {
				ss= s.split("\\s");
				int len= Integer.parseInt(ss[0]);
				if (len< gelSizeMin)
					gelSizeMin= len;
				if (len> gelSizeMax)
					gelSizeMax= len;
				DoubleVector vGC= map2D.get(len);
				if (vGC== null) {
					vGC= new DoubleVector(1);
					map2D.put(len, vGC);
				}
				if (ss.length> 1) {
					double gc= Double.parseDouble(ss[1]);
					if (gc< gelGCMin)
						gelGCMin= gc;
					if (gc> gelGCMax)
						gelGCMax= gc;
					vGC.add(gc);
				} else
					vGC.add(-1d);
			}
			buffy.close();
			GEL_NB_BINS_LENGTH= gelSizeMax- gelSizeMin+ 1;
			
			// init bin values
			gelSizeBin= (gelSizeMax- gelSizeMin)/ GEL_NB_BINS_LENGTH;
			if (gelGCMax> 0) 
				gelGCbin= (gelGCMax- gelGCMin)/ GEL_NB_BINS_GC;
			
			// binning 
			gelProb= new double[GEL_NB_BINS_LENGTH* (gelGCMax> 0? GEL_NB_BINS_GC: 1)];
			Iterator<Integer> iter= map2D.keySet().iterator();
			int debug43= 0;
			while(iter.hasNext()) {
				int len= iter.next();
				int p= getBin(len);
				DoubleVector vGC= map2D.get(len);
				double pLen= vGC.size()/ (double) total;
				gelProb[p]+= pLen; // come in unsorted
				if (p== 56)
					++debug43;
			}
			if (gelGCMax> 0) {
				iter= map2D.keySet().iterator();
				double frac= 1d/ total;
				while (iter.hasNext()) {
					int len= iter.next();
					int bin= getBin(len);
					double pLenFrac= gelProb[bin]* frac;
					gelProb[bin]= 0;
					DoubleVector vGC= map2D.get(len);
					for (int i = 0; i < vGC.size(); i++) {
						int bin2= getBin(vGC.get(i));
						gelProb[bin+ bin2]+= pLenFrac;
					} 
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	boolean collectOutputDistribution;
	
	boolean processNew(InputStream in, byte mode) {
			
			if (tmpFile== null)
				return false;
			
			try {
				long t0= System.currentTimeMillis();
				double breakRatio= 1;
				Long longNull= new Long(0);
				double thrTgt= 0.95;
				double tgtFrac= 0d;
				cumuLen= 0;
				int[] allLength= settings.getProfiler().getLen();
				long[] allMos= settings.getProfiler().getMolecules();
				for (int i = 0; i < allLength.length; i++) {
					cumuLen+= allMos[i]* allLength[i];
					prevMols+= allMos[i];
				}
				
				if (mode== MODE_NEBU) 
					initNebulizationParameters(
							settings.getFragNBlambda(),
							settings.getFragNBm(),
							settings.getThold(),
							settings.getProfiler().getMaxMoleculeLength()
					);
				
				// TODO estimate required disk space
				
				Processor[] processors= getProcessorPool(mode, Math.min(settings.getMaxThreads(), 1));
				
				for (roundCtr= 0;(!isStop())&& roundCtr< 1; ++roundCtr, mode= getMode()) {
					
					avgLength= cumuLen/ (double) prevMols; 
					
					//lengthV= new IntVector();
					lastMaxLen= maxLen;
					maxLen= 0;
					cumuLen= 0;
					lenSum= 0l;
					lenNb= 0;
					minLen= 0; maxLen= 0;
					if (mapFrags!= null) 
						mapFrags.clear();	// 20101215 re-init
									
					String msg= null;
					if (mode== MODE_RT)
						msg= MODE_RT_MESSAGE;
					else if (mode== MODE_NEBU) {
						msg= MODE_NEBU_MESSAGE+"-"+Integer.toString(roundCtr+1);
						//rndNebu= new GaussianRndThread(100);
						//rndNebu.start();
					} else if (mode== MODE_FILT_REJ) {
						msg= MODE_FILT_MESSAGE;
						filtSigSquare= ((settings.getFiltMax()- settings.getFiltMin())/ 2d);
						filtMu= settings.getFiltMin()+ filtSigSquare;
						filtFac= 1d/ Math.sqrt(2d* Math.PI* filtSigSquare);
						//oriMin= 15; oriMax= 1500;
						//oriProb= new double[100];
						oriMin= gelSizeMin; oriMax= gelSizeMax;
						oriProb= new double[gelProb.length];
						initBins(tmpFile, oriProb, oriMin, oriMax);
						deconvolve(oriProb, gelProb);
					} else if (mode== MODE_FRAG)
						msg= MODE_FRAG_MESSAGE+"-"+Integer.toString(roundCtr+1);
					if (plotter!= null) {
						plotter.reset(msg);
					}
					Log.progressStart(msg);
	
					Object[] keys= mapFrags.keySet().toArray();
					for (int i = 0; i < keys.length; i++) 
						mapFrags.put((CharSequence) keys[i], longNull);	// 20101205: changed from cast to string, ClassCastException for ByteArrayCharSequence
					
					currMols= 0; newMols= 0; tgtMols= 0;
					tstNewMols= 0; tstCurrMols= 0;
					cntLinesWr= 0;
					for (int i = 0; i < med.length; i++) 
						med[i]= 0;
					
					// IO 
	//				BufferedReader buffy= new BufferedReader(new FileReader(tmpFile));
	//				ThreadedQWriter qwriter= getQWriter();
	//				qwriter.init();
	//				qwriter.start();
	//				ThreadedBufferedByteArrayStream buffy= 
	//					new ThreadedBufferedByteArrayStream(10* 1024* 1024, fis, true, false);
					//BufferedOutputStream outStream
					//writer= new BufferedOutputStream(new FileOutputStream(tmpWriteFile), 1024* 1024);
					//ThreadedBufferedByteArrayStream inBacs= null;
					FileInputStream fis= new FileInputStream(tmpFile);
					FileOutputStream fos= new FileOutputStream(tmpWriteFile);
					rw= new SyncIOHandler2(2);
					rw.addStream(fis, 10* 1024);
					rw.addStream(fos, 10* 1024);
					for (int i = 0; i < processorPool.length; i++) 
						processorPool[i].setFos(fos);
					if (FluxSimulatorSettings.optDisk) 
						rw.start();
					
					ByteArrayCharSequence cs= new ByteArrayCharSequence(100);
					
					int perc= 0;
					long byteTot= tmpFile.length(), byteNow= 0l;
					//for (String s; (!isStop())&& (s= buffy.readLine())!= null;/*++currMols*/) {
					//for (buffy.readLine(cs); cs.end> 0; buffy.readLine(cs)) {
					while((!isStop())&& rw.readLine(fis, cs)> 0) {
						byteNow+= cs.length()+ 1;	// TODO fs length
						Log.progress(byteNow, byteTot);
						
						//qwriter.add(s);
						////qwriter.writeAll();
						++currMols;
						if (processors.length== 1) {
							ByteArrayCharSequence ccs= cs.cloneCurrentSeq();
							processors[0].cs= ccs;
							processors[0].run();
						} else {
							boolean searching= true;
							while (searching && !stop) {
								for (int i = 0; (!stop)&& i < processors.length; i++) {
									// 100330: took into sync block
		//							if (processors[i].ready()) {
		//								ByteArrayCharSequence ccs= cs.cloneCurrentSeq();
		//								processors[i].cs= ccs; 	// not in same line !
									// deadlocks :(
									synchronized (processors[i].lock) { 
										if (processors[i].ready()) {
											ByteArrayCharSequence ccs= cs.cloneCurrentSeq();
											processors[i].cs= ccs;
											processors[i].lock.notifyAll();
											searching= false;
										}
										//System.err.println("notify "+ccs);
									}
		//								//processors[i].interrupt();
									break;
								}
		/*						if (searching)
									try {
										Thread.currentThread().sleep(10);
									} catch (InterruptedException e) {
										; // :)
									}
								else
									break;
		*/												
							}
						}
					}
					
					rw.close();
					if (FluxSimulatorSettings.optDisk) 
						rw.join();
					
					//fis.close();
					//writer.flush();
					//writer.close();
					
					//buffy.close();
	/*				if (isStop())
						qwriter.close();
					qwriter.flush();
					qwriter.close();
	*/				
	/*				if (mode== MODE_NEBU) {					
						rndNebu.setStop();
						rndNebu.interrupt();
					}
	*/				
	/*				while (qwriter.isAlive()) // while deadlocks
						try {
							qwriter.join();
						} catch (InterruptedException e) {
							qwriter.close(); // ?? Thread.currentThread().interrupt();
						}
	*/
					//Distribution dist= new Distribution(lengthV.toIntArray());
					medLen= lenSum/ (double) lenNb;
					lastMaxLen= maxLen;
					lastMinLen= minLen;
					//lengthV= null;
					
					if (isStop())
						break;
					
					breakRatio= newMols/ (2d* currMols);
					tgtFrac= ((double) tgtMols)/ (currMols+ newMols);
	/*				if (mode== MODE_NEBU|| mode== MODE_FRAG)
						System.out.println("bratio "+Double.toString(breakRatio)
								+ ", all "+ (currMols+newMols)
								+ ", in "+ currMols
								+ ", new "+ newMols
								+ ", tgt "+ tgtMols
								+ ", trat "+ tgtFrac
						);
	*/					
					//DEBUG
	/*								File save= new File("N:\\tmp\\round_"+roundCtr);
									if (save.exists())
										save.delete();									
									FileHelper.copy(tmpFile, save);
	*/								
									
					boolean b= tmpFile.delete();
					if (!b)
						throw new IOException("Couldn't delete source");
					b= tmpWriteFile.renameTo(tmpFile);
					if (!b)
						throw new IOException("Couldn't move file");;
									
					Log.progressFinish(null, true);
					long total= 0;
					Iterator iter= mapFrags.keySet().iterator();
					while(iter.hasNext()) {
						Long val= mapFrags.get(iter.next());
						if (val!= null)
							total+= val;
					}
					if (Constants.verboseLevel>= Constants.VERBOSE_NORMAL) 
						System.err.println("\t\t"+(currMols+ newMols)+ " mol, "+total+": "+currMols+","+tstCurrMols+" "+newMols+","+tstNewMols);
					if (Constants.verboseLevel>= Constants.VERBOSE_NORMAL) {
						System.err.println("\t\t"+(currMols+ newMols)+ " mol: in "+ currMols+ ", new "+ newMols+ ", out "+ cntLinesWr);
						System.err.println("\t\tavg Len "+ (cumuLen/ (float) currMols)+ ", maxLen "+ maxLen);
						System.err.println("tgt frac "+tgtFrac);
					}
					if (plotter!= null) {
						plotter.paint();
						plotter.setMolTot(currMols+ newMols);
					}
					prevMols= currMols+ newMols;
	
				}
	
				if (!isStop()) {
					/*for (int i = 0; i < processors.length; i++) {
						if (processors[i]== null)
							continue;
						while (processors[i].isAlive()) {
							synchronized (processors[i].lock) { 
								if (processors[i].ready()) {
									processors[i].stop= true;
									//processors[i].interrupt();
									processors[i].lock.notifyAll();
								}
							}
							if (processors[i].stop)
								try {
									processors[i].join();
								} catch (InterruptedException e) {
									; // :)
								}
						}
					}
	*/				
					if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
						System.err.println(" OK");
					System.gc();
					return true;
				}
				
			} catch (Exception e) {
				if (Constants.verboseLevel>= Constants.VERBOSE_SHUTUP)
					e.printStackTrace();
			}
	
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				if (isStop())
					System.err.print(" FAILED");
				else
					System.err.print(" STOPPED");
	
			return false;
	
		}

	void processRT(ByteArrayCharSequence cs, int start, int end, int len, ByteArrayCharSequence id) {
			
			if (index1== null) 
				index1= new int[getRTeventNr(settings.getProfiler().getMaxMoleculeLength())];
			double[] wSense= null, wAsense= null;
			double n1= Double.NaN, n2= Double.NaN;
			int txLen= -1, howmany= -1;
			if (settings.getRTmotif()== null) {
				if (settings.getRtMode().equals(FluxSimulatorSettings.PAR_RT_MODE_POLY_DT)) {
					txLen= settings.getProfiler().getLength(id);
					if (end< txLen- 1)	// 0-based coordinates
						return;
					else
						howmany= getRTeventNr(txLen- end);
				} else 
					howmany= getRTeventNr(len);
			} else {
				howmany= getRTeventNr(len);
				wAsense= mapWeightAsense.get(id);
				n2= toCDF(wAsense, start, end);
				wSense= mapWeightSense.get(id);
				n1= toCDF(wSense, start, end);
			}
			if (howmany== 0)
				return;
			//processFragNot(start, end, len, id);

			// choose new 3' end(s)
			for (int i = 0; i < howmany; i++) {
				int p;
				if (wAsense== null) {
					if (settings.getRtMode().equals(FluxSimulatorSettings.PAR_RT_MODE_POLY_DT)) {
						p= end+ (int) Math.floor(rtRndWhere.nextDouble()* (txLen- end));
					} else {
						assert(settings.getRtMode().equals(FluxSimulatorSettings.PAR_RT_MODE_RANDOM));
						p= start+ (int) Math.round(rtRndWhere.nextDouble()* (len- 1));
					}
				} else {
					p= Arrays.binarySearch(wAsense, start, end, rtRndWhere.nextDouble());
					p= (p>= 0? p: -(p+1));
					// TODO put a minimum threshold on weight (thermodynamics affinity)?!
					++p; // anti-sense matrix, cut 3' of 0-position
				} 
				//new3Prime= Math.max(start+ p, new3Prime);
				index1[i]= p;
			}
			
			// extend first strand
			// resolve conflicts
			// choose new 5' (second strand synthesis) 
			// output
			Arrays.sort(index1,0,howmany);
			for (int i = howmany- 1; i >= 0; --i) {
				if (index1[i]< 0)
					continue; // got displaced
				int ext= settings.getRTminLen()+ (int) rnd2.nextDouble()* (settings.getRTmaxLen()- settings.getRTminLen());
				int to= Math.max(index1[i]- ext, start);
				// check GC
				double gc= getGCcontent(id, to, index1[i]- 1);	// bp is first nt of next fragment
				double pg= gc< rtC? 0d: 1d- Math.exp(- (gc- rtC)/settings.getRTminGC());
				Math.exp((-1d)*(gc- settings.getRTminGC())/ settings.getRTminGC());
				if (pg> 1|| rnd1.nextDouble()> pg) {
					index1[i]= -1;
					continue;
				}
				
				// resolve displacements
				boolean displaced= false;
				for (int j= i-1; j>= 0; --j) {
					if (index1[j]< 0)
						continue;
					if (index1[j]< to)
						break;
					double f;
					if (wAsense== null) {	// displacement function of distance
						int dist= Math.min(index1[i]- index1[j], index1[j]- start);								
						f=  Math.exp((-1d)*dist/ settings.getRTminLen());
					} else {	// displacement is a function of motifs
						double mi= wAsense[index1[i]- 1]- (index1[i]== 1? 0: wAsense[index1[i]- 2]);	// (-1) as cut after asense 0-pos
						double mj= wAsense[index1[j]- 1]- (index1[j]== 1? 0: wAsense[index1[j]- 2]);
						f= (mj== 0? 2: mi/ mj); // smaller -> more likely displaced
					}
					if (f> 1|| rnd3.nextDouble()<= f)
						index1[j]= -1; 	// displaced other pol
					else {
						displaced= true;
						break;
					}
				}
				if (displaced)
					continue;

				// choose 5'-end for second strand synthesis
				int howmany2= 5;	// roll 5 values
				int to2= Math.min(to+ 50, index1[i]- 1);	// within 50nt closest to 5'
				int new5Prime= index1[i];
				for (int j = 0; j < howmany2; j++) {
					int p;
					double r= rtRndWhere.nextDouble();
					if (wAsense== null)
						p= to+ (int) Math.floor(r* (to2- to));
					else {
						r= wSense[to]+ (r* (wSense[to2]- wSense[to]));
						p= Arrays.binarySearch(wSense, to, to2, r);
						p= (p>= 0? p: -(p+1));
					}
					new5Prime= Math.min(p, new5Prime);
				}
				
				// write it
				int nuLen= index1[i]- new5Prime+ 1;
				// updateMedian(nuLen);
				cs.replace(0, new5Prime);
				cs.replace(1, index1[i]);
				cumuLen+= nuLen;
				incLinesWrote();
				++newMols;
				rw.writeLine(cs, fos);	// id is invalid now
				if (nuLen<= settings.getFragNBlambda())
					++tgtMols;
				if (len<= settings.getFragNBlambda())
					--tgtMols;

			}
			
			// re-normalize
			if (settings.getRTmotif()!= null) {
				toPDF(wAsense, start, end, n2);
				toPDF(wSense, start, end, n1);				
			}
				
	}

	/**
	 * @param id	ID of the sequence
	 * @param i start index (included)
	 * @param j end index (included)
	 * @return
	 */
	private double getGCcontent(ByteArrayCharSequence id, int i, int j) {
		
		CharSequence seq= getMapTxSeq().get(id);
		int g= 0, n= 0;
		for (int k = i; k <= j; ++k) {
			char c= seq.charAt(k);
			if (c== 'G'|| c=='C')
				++g;
			else if (c!= 'N')
				++n;
		}
		
		return (g/ (double) (n+ g));
	}

	private int getRTeventNr(int len) {
		
		int howmany= (int) Math.ceil(len/ (double) 100);	//  settings.getRTPrimerPerNucleotides() 
//		howmany= (int) Math.round(rndHowMany.nextPoisson(((new3Prime- start)+1)/ 50d));	// TODO Poisson..;
		return howmany;
	}


	int[] tmpFragEnzyme= null;
	/**
	 * Implements enzymatic cleavage of the fragment according to a provided motif (PWM).
	 * 
	 * @param start
	 * @param end
	 * @param len
	 * @param id
	 */
	void processFragEnzyme(ByteArrayCharSequence cs, int start, int end, int len, ByteArrayCharSequence id) {
		
		assert(mapWeightSense!= null&& mapWeightAsense!= null);
		double[] wsense= mapWeightSense.get(id), wasense= mapWeightAsense.get(id);
		
		// short-hand value to get number of enzyme attacks as function of length
		// could be a parameter
		int n= (int) Math.sqrt(len);	
		int k= 0;
		int[] pos= (tmpFragEnzyme== null? new int[n]: tmpFragEnzyme);
		
		// get breakpoints
		double norm= Fragmenter.toCDF(wsense, start, end);
		for (int i = 0; i < n; i++) {
			
			// localize potential breakpoint
			boolean sense= rnd1.nextFloat()< 0.5;
			double[] a= (sense? wsense: wasense);
			int bp= Arrays.binarySearch(a, start, end+ 1, rnd2.nextDouble());
			bp= (bp>= 0? bp: -(bp+1));
			if (bp>= end)
				continue;			
			
			// lookup in array with bp positions
			double pb= (bp== 0? a[bp]: a[bp]- a[bp- 1]);	// breaking probability
			// ok, the pwm recognizes the base where it is cut
			// but has no concept about 5'/3'. If cutting always
			// consistently 5' of 0-position, the bp in sense
			// marks the start of the next fragment, when motif
			// is recognized in anti-sense, we got to add 1.
			// (e.g., palindrome ACGT cut at G is found at 
			// position 1166 in sense, and position 1165 in 
			// anti-sense)
			bp= (sense? bp: bp+1);	// set bp consistently to 1st nt of next fragment
			int idx= Arrays.binarySearch(pos, 0, k, bp);
			if (idx>= 0)
				continue;
			else
				idx= -(idx+1);
			
			// breaking probability
			pb/= (sense? pwmSense.getMaximumP(): pwmSense.getMaximumP()); // should be the same if only transposed
			if (rnd3.nextDouble()<= pb) {
				System.arraycopy(pos, idx, pos, idx+ 1, (k- idx)+ 1);
				pos[idx]= bp;
				++k;
			}
		}
		Fragmenter.toPDF(wsense, start, end, norm);
		
		// break
		int last= start;
		for (int i = 0; i < k; i++) {
			int now= start+ pos[i]- 1;	// end
			int nuLen= now- last;
			if (nuLen== 0)
				continue;
			//updateMedian(nuLen);
			++newMols;
			cumuLen+= nuLen;
			cs.replace(0, last);
			cs.replace(1, now);
			rw.writeLine(cs, fos);	
			++cntLinesWr;
			last= now+ 1;	// start
		}
		int now= end;
		int nuLen= now- last+ 1;
		//updateMedian(nuLen);
		cumuLen+= nuLen;
		cs.replace(0, last);
		cs.replace(1, now);
		rw.writeLine(cs, fos);	
		++cntLinesWr;
	}

	/**
	 * 
	 * @param id
	 * @param start
	 * @param end
	 * @param sense
	 * @return a value representing the sum of the values that have been turned into a CDF
	 */
	public static double toCDF(double[] a, int start, int end) {

		for (int i = start+ 1; i <= end; i++) 
			a[i]+= a[i- 1];
		double sum= a[end];
		if (sum== 0d)
			return sum;
		for (int i = start; i <= end; i++){
			a[i]/= sum;
			assert(!(Double.isNaN(a[i])|| Double.isInfinite(a[i])|| a[i]< 0|| a[i]> 1));
		}
		
		return sum;
	}

	public static void toPDF(double[] a, int start, int end, double sum) {
		if (sum== 0)
			return;
		for (int i = start; i <= end; i++) 
			a[i]*= sum;
		for (int i = end; i > start; --i) {
			a[i]-= a[i-1];
			assert(!(Double.isNaN(a[i])|| Double.isInfinite(a[i])|| a[i]< 0/*|| a[i]> 1*/));
		}
	}

	boolean out1= false, out2= false, out3= false; //debug for frag
	/**
	 * Implements a uniform random fragmentation model according to 
	 * Tenchov and Yanev.
	 * 
	 * @param nebu
	 * @param start
	 * @param end
	 * @param len
	 * @param id
	 */
	void processFrag(ByteArrayCharSequence cs, boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
			            
		// get parameters
		double d0= settings.getFragURd0();
		assert(d0>= 1); // problem with integer breakpoints, when fragment size << 1 !
		double delta= getFragURdelta(len);  
		double eta= getFragUReta();
		//eta= 170;
		
		// [2a] dmax= eta((delta- 1)/delta)^(1/ delta)
		double dmax= eta* Math.pow((delta- 1)/ delta,1/ delta);
		// expectancy E of the Weibull distribution W(delta,eta)
		// [3a] E= d0+ eta*gamma(1/delta+ 1)
		double E= d0+ eta* Math.exp(Gamma.logGamma(1d+ 1d/delta));
		double D2= Math.pow(eta, 2)* 
						(Math.exp(Gamma.logGamma((2d/ delta)+ 1d))- 
						Math.pow(Math.exp(Gamma.logGamma((1d/ delta)+ 1d)), 2d));
		
		// DEBUG
/*		if (out1&& len< 1000) {
			System.err.println(len+"\t"+delta+"\t"+E+"\t"+D2+"\t"+dmax);
			out1= false;
			System.currentTimeMillis();
		} else if (out2&& len> 1200&& len< 1500) {
			System.err.println(len+"\t"+delta+"\t"+E+"\t"+D2+"\t"+dmax);
			out2= false;
			System.currentTimeMillis();
		} else if (out3&& len> 2000) {
			System.err.println(len+"\t"+delta+"\t"+E+"\t"+D2+"\t"+dmax);
			out3= false;
			System.currentTimeMillis();
		}
*/		
		
		// determine n, the number of fragments (i.e. (n-1) breakpoints)
		double nn= ((double) len)/ E;			
		int n= (int) nn;
			// (int) Math.round(nn+ (rndBreak.nextGaussian()));	 
			// (int) rndImpl.nextPoisson(nn);		// too variable
		double nr= nn- n;
		double r= rndBreak.nextDouble();
		if (r<= nr)
			++n;
	
		// molecule does not break
	    if (n<= 1|| len<= 1|| (n*d0)>= len) {	
			updateMedian(len);
			cumuLen+= len;
			cs.replace(0, start);
			cs.replace(1, end);
	    	processFragNot(cs, start, end, len, id);
			return;
		}
	    
		// uniformly cut (n-1) times unit space
	    double [] x = new double[n];
	    for(int i= 0; i< (n-1); ++i) 
	      x[i]= rndBreak.nextDouble();	  
	    x[x.length- 1]= 1;	// last breakpoint is end
		
	    // get fragment lengths (in unit space)
	    Arrays.sort(x);
	    for (int i= (n- 1);i> 0; --i)
	      x[i]-= x[i-1];
	
	    // compute c, transform to molecule space
	    float sum = 0;
	    for(int i= 0; i< x.length; i++) {
	      sum+=Math.pow(x[i], 1/delta);
	    }  
	    double c = Math.pow((len - n* d0)/ sum, -delta); 	
	    for (int i= 0; i< n; i++)
	      x[i]=  d0 + Math.pow(x[i]/ c, (1/delta));
	
	    double dsum= 0;
	    for (int i = 0; i < n; i++) {
	    	int nuStart= start+ (int) Math.round(dsum); 
	    	dsum+= x[i];
			int nuEnd= start+ (int) Math.round(dsum)- 1;
	    	//double frac= dsum/ len;
			int nuLen= (nuEnd- nuStart)+ 1;
			if (nuLen< 0)
				System.currentTimeMillis();
			cs.replace(0, nuStart);
			cs.replace(1, nuEnd);
			cumuLen+= nuLen;
			updateMedian(nuLen);
			incLinesWrote();
			rw.writeLine(cs, fos);	// id is invalid now
		}
	    assert(Math.round(dsum)== len);
	}

	void processNebu(ByteArrayCharSequence cs, boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {
	
		if (nebuRecursionDepth< 1) {
			rw.writeLine(cs, fos);	// no breakage
			return;
		}
		
		// parameters:
		// pb: M, lambda, len
		// bp: Sigma, length
		double lambda= Fragmenter.this.settings.getFragNBlambda();
		double sigma= Fragmenter.this.settings.getFragNBsigma();
		double M= Fragmenter.this.settings.getFragNBm();
		int recDepth= Fragmenter.this.nebuRecursionDepth;
		double C= Fragmenter.this.nebuC;
		if (index1== null)
			index1= new int[(int) Math.pow(2, recDepth)];
		Arrays.fill(index1, -1);
		index1[0]= len;
		int fragmentNb= 1;
		
		// now break it!
		for (int i = 0; i < recDepth; ++i) {
			for (int j = 0; index1[j]> 0; ++j) {
	
				// breakpoint location
				// N(length/2,sigma)= (N(0,1)*sigma*(length/2)+ length/2
				int L= index1[j]; 
				//int bp= (int) rdiNebuBP.nextGaussian(len/ 2d, len/ 4d);
				//int bp= (int) ((rndBP.nextGaussian()* sigma
				//		* ((L-1)/2d))+ (L-1)/2d);	// bp index [0;L[						
				int bp= (int) nextGaussianDouble(rndBP, 0, L- 1);
				
				// breaking probability (pb)
				// pb= 1- exp^(-((x-C)/lambda)^M)
				int minL= (bp< (L- bp)? bp+1 : L-bp-1);
				double pb= minL< C? 0d: 1d- Math.exp(-Math.pow((minL- C)/lambda, M));
				double r= rndBreak.nextDouble();
				if (r> pb)
					continue;
				
				// fragment j breaks
				int rest= fragmentNb- j- 1;
				if (rest> 0) {
					assert(j+ 2+ rest<= index1.length);
					System.arraycopy(index1, j+ 1, index1, j+ 2, rest);
				}
				assert(bp+1> 0&& L- bp- 1> 0);
				index1[j]= bp+ 1;
				index1[j+1]= L- bp- 1;	// L- (bp+1)
				++fragmentNb;
				++j;
			}
			
		}
		
		// write result to disk
		for (int j = 0; index1[j]> 0; ++j) {
			cumuLen+= index1[j];
			cs.replace(0, start);
			start+= index1[j];
			cs.replace(1, start- 1);
			rw.writeLine(cs, fos);	// id is invalid now
			
			// update stats
			incLinesWrote();
			if (index1[j]<= settings.getFragNBlambda())
				++tgtMols;
			if (len<= settings.getFragNBlambda())
				--tgtMols;
		}
		assert(start== len);
		
	}

	int lastLen= -1; //MCMC
	double lastP= -1;
	/**
	 * see <code>http://personal.strath.ac.uk/gary.koop/extra_material_on_metropolis.pdf</code>,
	 * <code>http://www.ps.uci.edu/~markm/numerical_methods/Metropolis%96Hastings%20algorithm%20-%20Wikipedia,%20the%20free%20encyclopedia.pdf</code>
	 * @param start
	 * @param end
	 * @param len
	 * @param id
	 */
	void processFilterMCMC(ByteArrayCharSequence cs, int start, int end, int len, ByteArrayCharSequence id, 
			AbstractDistribution dGenerate, AbstractDistribution[] dProposal) {
	
		// first value always accepted to init algorithm (but not output)
		double p= 0d;
		for (int i = 0; i < dProposal.length; i++) {
			p+= dProposal[i].getP(len);
		}
		if (lastLen< 0) {
			lastLen= len;
			lastP= p;
			return;
		}
		
		// Metropolis/Hastings/Ema
		double a1= p/ lastP;
		double a2= dGenerate.getP(lastLen, len)/ dGenerate.getP(len, lastLen);
		double a= a1* a2;
		
		// accept 
		if (a>= 1|| rndGel.nextDouble()<= a) {
			lastLen= len;
			lastP= p;
			incLinesWrote();
			rw.writeLine(cs, fos);
			if (Fragmenter.this.plotter!= null) 
				Fragmenter.this.plotter.plot(start, end, len, id);
			addFragCount(id, 1l);
		} else
			--newMols;
	
	}

	void processFilterRejection(ByteArrayCharSequence cs, int start, int end, int len, ByteArrayCharSequence id, 
			AbstractDistribution[] d, boolean probDistr) {
		
		// get (possibly cumulative) probability for length being in result
		double plen= 0d;
		for (int i = 0; i < d.length; i++) {
			double p= (probDistr? d[i].getP(len): d[i].getRelFreq(len));
			plen+= d[i].getWeight()* p;
		}
		
		// Bernoulli trial
		if (plen> 1|| rndGel.nextDouble()< plen) {
			incLinesWrote();
			rw.writeLine(cs, fos);
			if (Fragmenter.this.plotter!= null) 
				Fragmenter.this.plotter.plot(start, end, len, id);
			addFragCount(id, 1l);
		} else
			--newMols;
	}

	private void processFragNot(ByteArrayCharSequence cs, int start, int end, int len,
			ByteArrayCharSequence id) {
		lenSum+= len;                            	
		++lenNb;
		maxLen= Math.max(maxLen,len);
		minLen= Math.min(minLen,len);
		incLinesWrote();
		rw.writeLine(cs, fos);
		//increaseFragCount(id);	// 20101215 deactivated, counts for pro-file
		if (Fragmenter.this.plotter!= null)
			Fragmenter.this.plotter.plot(start, end, -1, id);
	}

	private void updateMedian(int nuLen) {
		
		int[] a= med;
		if (med[0]== 0) {
			med[0]= nuLen;
			for (int i = 0; i < a.length- 1; i++) {
				if (med[i]> med[i+ 1]) {
					int h= med[i];
					med[i]= med[i+ 1];
					med[i+ 1]= h;
				}
			}
			return;
		}
		
		int p= Arrays.binarySearch(med, nuLen);
		if (p>= 0)
			return;
	
		p= -(p+ 1);			
		if (p== 0|| p>= med.length)
			return;
		if (p+ 1< med.length)
			System.arraycopy(med, p, med, p+ 1, med.length- p- 1);
		med[p]= nuLen;
	}

	void processInitial(ByteArrayCharSequence cs, int origLen, CharSequence locID) {
				
				// transcript variation
				int start= 0;
				int end= origLen- 1;
				if (!Double.isNaN(settings.getTssMean())) {
					start= origLen;
					while (start>= Math.min(100, origLen))
						start= (int) Math.round(rndTSS.nextExponential(Math.min(settings.getTssMean(),origLen/4)));	// exp mean= 25: exceeds bounds, nextGaussian(1,100))-100;
					double r= rndPlusMinus.nextDouble();
					if (r< 0.5)
						start= -start;
				}
				if (!(Double.isNaN(settings.getPolyAshape())|| Double.isNaN(settings.getPolyAscale()))) {
					int pAtail= 301;
					while (pAtail> 300)
						pAtail= (int) Math.round(sampleWeibull(rndPA, settings.getPolyAscale(), settings.getPolyAshape()));	// 300d, 2d
					end= origLen+ pAtail; 
				}
				if (end< origLen) {
					if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
						System.err.println("[ERROR] end < length in Fragmenter$Processor.processInitial(): "+ end+ "<"+ origLen);
				}
	
				int newLen= end- start+ 1;
				assert(start< end);
				lenSum+= newLen;
				cumuLen+= newLen;
				++lenNb;
				maxLen= Math.max(maxLen,newLen);
				minLen= Math.min(minLen,newLen);
				
				
				//String line= start+ TAB+ end+ TAB+ settings.getProfiler().getIds()[i]+"\n";
				//ByteArrayCharSequence ccs= cs.cloneCurrentSeq();
				cs.replace(0, start);
				cs.replace(1, end);
				cs.replace(2, locID);
				
				// RFU: chr_tid separated
	/*			int p= 0;
				while(locID.charAt(++p)!= ':');
				cs.end= cs.p2;
				cs.ensureLength(cs.end, p);
				byte[] a= cs.a;
				for (int j = 0; j < p; j++) 
					a[cs.end++]= (byte) locID.charAt(j); 
				cs.append(BYTE_SEP_LC_TX);
				cs.append(locID);
	*/			
				
	
			}
}
