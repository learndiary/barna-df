package fbi.genome.astalavista.statistics;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.sql.Time;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;
import java.util.Vector;

import net.maizegenetics.pal.statistics.FisherExact;
import fbi.genome.io.FileHelper;
import fbi.genome.io.gtf.GTFwrapper;
import fbi.genome.model.ASEvent;
import fbi.genome.model.commons.MyFile;
import fbi.genome.model.constants.Constants;
import fbi.genome.model.gff.GFFObject;

public class EventDavid extends Thread {

	final static String SFX_TBL= ".tbl", SFX_KAPPA= "_kappa";
	
	final static String DEFAULT_REDTABLE_PATH= "M:\\work\\claudia\\david_genbank\\test\\redtable_events.txt",
		DEFAULT_BACKTABLE_PATH="M:\\work\\claudia\\david_genbank\\test\\backtable_events.txt";
	
	static final String TERM_HEADER=
		"PVAL_CB\tCOUNT_Cp\tNOT_IN_GROUP\t%FRAC_C\tFCHANGE_CB\tCATEGORY\tTERM\t" +
		"count_Ca\tcount_Up\tcount_Ua\tcount_Bp\tcount_Ba\t" +
		"pVal_CUB\tpVal_UB\tpVal_CU\t"+
		"%frac_CU\tfc_CUB\t%FRAC_U\tfc_UB\n" +
		
		"f(Cp,Ca,Bp,Ba)\tCp\tUp-Cp\tCp/Ca\t(Cp*Ba)/(Ca*Bp)\tdbase/level\tterm\t" +
		"Ca\tUp\tUa\tBp\tBa\t" +
		"f(Cp,Ua,Bp,Ba)\tf(Up,Ua,Bp,Ba)\tf(Cp,Ca,Up,Ua)\t" + 
		"Cp/Ua\t(Cp*Ba)/(Ua*Bp)\tUp/Ua\t(Up*Ba)/(Ua*Bp)\n";
	


	
	class RedTableWrapper {
		
		final char[] symbols= new char[] {'0','1','2'};
		
		File f;
		int lineLength;
		
		public RedTableWrapper(String fPath) {
			f= new File(fPath);
		}
		
		private BufferedWriter buffy= null;
		/**
		 * @deprecated untested
		 * @param cbuf
		 * @param keepOpen
		 */
		public void appendLine(char[] cbuf, boolean keepOpen) {

			try {				
				if (buffy == null) 
					buffy = new BufferedWriter(new FileWriter(f, true));

				buffy.write(cbuf);
				
				if (!keepOpen) {
					buffy.flush();
					buffy.close();
				}
					
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		public void writeOutTable() {
			try {
				BufferedWriter buffy = new BufferedWriter(new FileWriter(f));
				
				int idLen= 0;
				for (int i = 0; i < rowIDs.length; i++) 
					idLen=(rowIDs[i].length()>idLen)?rowIDs[i].length():idLen;
				
				char[] cbuf= new char[idLen+1+redTable[0].length+1];
				cbuf[idLen]= '\t';
				cbuf[cbuf.length-1]= '\n';
				for (int i = 0; i < redTable.length; i++) {
					rowIDs[i].getChars(0,rowIDs[i].length(),cbuf,0);
					for (int j = rowIDs[i].length(); j < idLen; j++) 
						cbuf[j]= ' ';
					for (int j = 0; j < redTable[i].length; j++) 
						cbuf[idLen+1+j]= symbols[redTable[i][j]];
					buffy.write(cbuf);
				}
				
				// write "header"
				String tabStr= "\t";
				for (int i = 0; i < attrIDs.length; i++) {
					buffy.write(attrIDs[i]);
					if (i< attrIDs.length-1)
						buffy.write(tabStr);
				}
				buffy.write("\n");
				
				buffy.flush();
				buffy.close();
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		public int readInTable() {
			
			try {
				
				BufferedReader buffy= new BufferedReader(new FileReader(f));
				System.err.print("["+Constants.getTime()+"] Reading RedTable.. ");
				System.err.flush();
				
				// read always a pair of lines to find the last one!!
				String s= buffy.readLine(), s2= s;
				int idLen= -1, tabWidth= -1;
				for (idLen = 0; idLen < s2.length(); idLen++) 
					if (s2.charAt(idLen)== '\t')
						break;
				tabWidth= s2.length()- idLen- 1;
				
				Vector<byte[]> vRed= new Vector<byte[]>(100,100);
				Vector<String> vRowID= new Vector<String>(100,100);
				int perc= 0;
				long bytesRead= s2.length()+1, fileLen= f.length();
				byte[] row= new byte[tabWidth];
				while((s= buffy.readLine())!= null) {
					bytesRead+= s.length();
					if (bytesRead*10d/ fileLen> perc) {
						System.err.print("*");
						System.err.flush();
						++perc;
					}
					vRowID.add(s2.substring(0,idLen).trim());
					assert(s2.length()- idLen- 1== tabWidth);
					for (int i = idLen+1; i < s2.length(); i++)  {
						row[i-idLen-1]= (byte) (((byte) s2.charAt(i))- 48);
//						if (row[i- idLen- 1]!= 0)
//							System.currentTimeMillis();
					}
					vRed.add(row.clone());
					s2= s;
				}
				
				buffy.close();
				attrIDs= s2.split("\t");
				assert(attrIDs.length== row.length);
				idLists= new HashSet[attrGrpIDs.size()+1];	// last is gene
				for (int i = 0; i < attrIDs.length; i++) {
					String[] grpID= attrIDs[i].split(attrGrpIDseparator);
					assert(grpID.length== 2);
					
				}
				for (int i = 0; i < idLists.length-1; i++) 
					idLists[i]= new HashSet<Integer>();
				
				
				redTable= new byte[vRed.size()] [];
				rowIDs= new String[vRowID.size()];
				for (int i = 0; i < redTable.length; i++) {
					redTable[i]= vRed.elementAt(i);
//					for (int j = 0; j < row.length; j++) {
//						if (redTable[i][j]!= 0)
//							System.currentTimeMillis();
//					}
					rowIDs[i]= vRowID.elementAt(i);
				}				
				System.err.println(" "+redTable.length+"x"+redTable[0].length+" entries.\n");
				return 0;
				
			} catch (Exception e) {
				System.err.println("Redtable not found");	//e.printStackTrace();
				return -1;
			}
			
		}
		
		
		
		public void close() {
			if (buffy!= null)
				try {
					buffy.flush();
					buffy.close();
				} catch (Exception e) {
					; // :)
				}
		}
		
	}
	
	public static class IntegerCompound implements Comparable {
		public int i;
		public Object o;
		
		//@Override
		public int compareTo(Object o) {
			int j= ((IntegerCompound) o).i;
			if (i< j)
				return -1;
			if (j> i)
				return 1;
			return 0;
		}
	}
	/**
	 * from the page http://www.langsrud.com/fisher.htm
	 * @author micha
	 *
	 */
	// TODO maybe inefficient because all tailed tests are calculated
	public static class FisherExactTest {
		static final boolean debug= false;
		boolean console= false;
		
		double left, right, twotail;
		private int sn11, sn1_, sn_1, sn;
		private double sprob;
		
		public void reset() {
			left= 0; right= 0; twotail= 0;
			sn11= 0; sn1_= 0; sn_1= 0; sn= 0;
			sprob= 0;
		}
		
		public void fisherExact(int n11, int n12, int n21, int n22) {

			if (sprob!= 0)	// safety reset
				reset();
			
			if (console)
				System.out.println("       Fisher's Exact Test \nhttp://www.nr.no/~langsrud/fisher.htm"
						+ "\n------------------------------------------");
			
			if (n11 < 0)
				n11 *= -1;
			if (n12 < 0)
				n12 *= -1;
			if (n21 < 0)
				n21 *= -1;
			if (n22 < 0)
				n22 *= -1;

			int n1_ = n11 + n12;
			int n_1 = n11 + n21;
			int n = n11 + n12 + n21 + n22;
			double prob = exact(n11, n1_, n_1, n);

			if (twotail > 1)
				twotail = 1;
			if (console)
				System.out.println("\n" + " TABLE = [ " + n11 + " , " + n12 + " , "
						+ n21 + " , " + n22 + " ]" + "\n" + "Left   : p-value = "
						+ left + "\nRight  : p-value = " + right
						+ "\n2-Tail : p-value = " + twotail
						+ "\n------------------------------------------");
		}


		private double exact(int n11, int n1_, int n_1, int n) {
			if (debug)
				System.out.println("exact called " + n11 + ", " + n1_ + ", " + n_1
					+ ", " + n);
			double sleft, sright, sless, slarg;
			double p, prob;
			int i, j = 0;
			int max = n1_;
			if (n_1 < max)
				max = n_1;
			int min = n1_ + n_1 - n;
			if (min < 0)
				min = 0;
			if (min == max) {
				sless = 1;
				sright = 1;
				sleft = 1;
				slarg = 1;
				left= sless;
				right= slarg;
				twotail= sleft+ sright;
				return 1;
			}
			prob = hyper0(n11, n1_, n_1, n);
			sleft = 0;
			p = hyper(min);
			for (i = min + 1; p < 0.99999999 * prob; i++) {
				sleft += p;
				p = hyper(i);
			}
			i--;
			if (p < 1.00000001 * prob)
				sleft += p;
			else
				i--;
			sright = 0;
			p = hyper(max);
			for (j = max - 1; p < 0.99999999 * prob; j--) {
				sright += p;
				p = hyper(j);
				if (debug)
					System.out.println("j " + j + ", p " + p);
			}
			j++;
			if (debug)
				System.out.println("j " + j + "p " + p + ", prob" + prob);
			if (p < 1.00000001 * prob)
				sright += p;
			else
				j++;
			if (debug)
				System.out.println("i " + i + ", n11 " + n11 + ", j " + j);
			if (Math.abs(i - n11) < Math.abs(j - n11)) {
				sless = sleft;
				if (debug)
					System.out.println("sless1 now " + sless);
				slarg = 1 - sleft + prob;
			} else {
				sless = 1 - sright + prob;
				if (debug)
					System.out.println("sless1 now " + sless);
				slarg = sright;
			}
			if (debug)
				System.out.println("exact() returns " + prob);
			
			left= sless;
			right= slarg;
			twotail= sleft+ sright;
			return prob;
		}

		private double hyper(int n11) {
			return (hyper0(n11, 0, 0, 0));
		}

		private double hyper0(int n11i, int n1_i, int n_1i, int ni) {
			if (debug)
				System.out.println("hyper0() " + n11i + "," + n1_i + "," + n_1i + ","
					+ ni);
			if ((n1_i | n_1i | ni) == 0) {
				if (!(n11i % 10 == 0)) {
					if (n11i == sn11 + 1) {
						if (debug)
							System.out.println("choose A1");
						sprob *= ((double) (sn1_ - sn11) / (n11i))
								* (((double) sn_1 - sn11) / (n11i + sn - sn1_ - sn_1));
						sn11 = n11i;
						return sprob;
					}
					if (n11i == sn11 - 1) {
						if (debug)
							System.out.println("choose A2");
						sprob *= ((double) (sn11) / (sn1_ - n11i))
								* (((double) sn11 + sn - sn1_ - sn_1) / (sn_1 - n11i));
						sn11 = n11i;
						return sprob;
					}
				}
				sn11 = n11i;
			} else {
				if (debug)
					System.out.println("choose B");
				sn11 = n11i;
				sn1_ = n1_i;
				sn_1 = n_1i;
				sn = ni;
			}
			if (debug)
				System.out.println("hyper323() " + sn11 + "," + sn1_ + "," + sn_1 + ","
					+ sn);
			sprob = hyper_323(sn11, sn1_, sn_1, sn);
			if (debug)
				System.out.println("hyper323() returns " + sprob);
			return sprob;
		}

		private double lnfact(int n) {
			if (n <= 1)
				return (0);
			return (lngamm(n + 1));
		}

		private double lnbico(int n, int k) {
			return (lnfact(n) - lnfact(k) - lnfact(n - k));
		}

		private double hyper_323(int n11, int n1_, int n_1, int n) {
			if (debug)
				System.out.println(lnbico(n1_, n11) + "," + lnbico(n - n1_, n_1 - n11)
					+ "," + lnbico(n, n_1));
			return (Math.exp(lnbico(n1_, n11) + lnbico(n - n1_, n_1 - n11)
					- lnbico(n, n_1)));
		}

		private double lngamm(int z)
		// Reference: "Lanczos, C. 'A precision approximation
		// of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
		// Translation of Alan Miller's FORTRAN-implementation
		// See http://lib.stat.cmu.edu/apstat/245
		{
			double x = 0;
			x += 0.1659470187408462e-06 / (z + 7);
			x += 0.9934937113930748e-05 / (z + 6);
			x -= 0.1385710331296526 / (z + 5);
			x += 12.50734324009056 / (z + 4);
			x -= 176.6150291498386 / (z + 3);
			x += 771.3234287757674 / (z + 2);
			x -= 1259.139216722289 / (z + 1);
			x += 676.5203681218835 / (z);
			x += 0.9999999999995183;
			return (Math.log(x) - 5.58106146679532777 - z + (z - 0.5)
					* Math.log(z + 6.5));
		}

		public boolean isConsole() {
			return console;
		}

		public void setConsole(boolean console) {
			this.console = console;
		}

		public double getLeft() {
			return left;
		}

		public void setLeft(double left) {
			this.left = left;
		}

		public double getRight() {
			return right;
		}

		public void setRight(double right) {
			this.right = right;
		}

		public double getTwotail() {
			return twotail;
		}

		public void setTwotail(double twotail) {
			this.twotail = twotail;
		}

	}
	
	
	static byte[][] testLittle = new byte[][] { { 1, 1, 0, 0, 1, 0 },
			{ 1, 1, 0, 1, 1, 0 } };
	static byte[][] testBitMore = new byte[][] { { 1, 1, 0, 0, 1, 0 },
			{ 1, 1, 0, 1, 1, 0 }, { 1, 0, 0, 1, 1, 1 }, { 1, 1, 0, 0, 1, 1 },
			{ 0, 1, 1, 1, 1, 1 }, { 0, 0, 1, 1, 0, 1 }, { 0, 0, 1, 1, 0, 1 } };
	static float[][] testKappaTbl = new float[][] {
			{ 1, 1, 1, 0.35f, -0.5f, -0.5f, -0.5f, 0 },
			{ 1, 1, 1, 0.35f, -0.5f, -0.5f, -0.5f, 0 },
			{ 1, 1, 1, 0.35f, -0.5f, -0.5f, -0.5f, 0 },
			{ 0.35f, 0.35f, 0.35f, 1, 0.35f, 0.35f, 0.35f, -0.11f },
			{ -0.5f, -0.5f, -0.5f, 0.35f, 1, 1, 1, 0 },
			{ -0.5f, -0.5f, -0.5f, 0.35f, 1, 1, 1, 0 },
			{ -0.5f, -0.5f, -0.5f, 0.35f, 1, 1, 1, 0 },
			{ 0, 0, 0, -0.11f, 0, 0, 0, 1 } };
	
	
	float filter_min_association= 0.1f;
	
	/**
	 * Similarity Term Overlap (any value >=0; default = 4): 
	 * the minimum number of annotation terms overlapped between 
	 * two genes in order to be qualified for kappa calculation. 
	 * This parameter is to maintain necessary statistical power 
	 * to make kappa value more meaningful. The higher value, the 
	 * more meaningful the result is.
	 */
	int kappa_similarityOvl= 10;	// default 3 in 2008 webinterface
	
	/**
	 * Similarity Threshold (any value between 0 to 1; Default = 0.35): 
	 * the minimum kappa value to be considered biological significant. 
	 * The higher setting, the more genes will be put into unclustered 
	 * group, which lead to higher quality of functional classification 
	 * result with a fewer groups and a fewer gene members. Kappa value 
	 * 0.3 starts giving meaningful biology based on our genome-wide 
	 * distribution study. Anything below 0.3 have great chance to be noise.
	 */
	float kappa_similarityThr = 0.5f; // default 0.5 in 2008 webinterface
	
	float clustering_nbCloseThr = 0.5f;	// not in webinterface (!?!)

	/**
	 * Initial Group Members (any value >=2; default = 4): 
	 * the minimum gene number in a seeding group, which affects the 
	 * minimum size of each functional group in the final. In general, 
	 * the lower value attempts to include more genes in functional groups, 
	 * particularly generates a lot small size groups.
	 */
	int clustering_initialGroupMembershipThr = 50; // default 3 in 2008 webinterface

	/** 
	 * Final Group Members (any value >=2; default = 4): 
	 * the minimum gene number in one final group after �cleanup� procedure. 
	 * In general, the lower value attempts to include more genes in functional 
	 * groups, particularly generates a lot small size groups. It co-functions 
	 * with previous parameters to control the minimum size of functional groups. 
	 * If you are interested in functional groups containing only 2 or 3 genes, 
	 * you need to set it to a very low value. Otherwise, the small group will 
	 * not be displayed and will be put into the unclustered group.Final Group 
	 * Members (any value >=2; default = 4): the minimum gene number in one final 
	 * group after �cleanup� procedure. In general, the lower value attempts to 
	 * include more genes in functional groups, particularly generates a lot small 
	 * size groups. It co-functions with previous parameters to control the minimum 
	 * size of functional groups. If you are interested in functional groups 
	 * containing only 2 or 3 genes, you need to set it to a very low value. 
	 * Otherwise, the small group will not be displayed and will be put into the 
	 * unclustered group.
	 */
	int clustering_finalGroupMembershipThr= 50; // default 3 in 2008 webinterface
	
	/**
	 * Multi-linkage Threshold (any value between 0% to 100%; default = 50%): 
	 * It controls how seeding groups merge each other, i.e. two groups sharing 
	 * the same gene members over the percentage will become one group. The higher 
	 * percentage, in general, gives sharper separation i.e. it generates more 
	 * final functional groups with more tightly associated genes in each group. 
	 * In addition, changing the parameter does not contribute extra genes into 
	 * unclustered group.
	 */
	float clustering_multiLinkThr = 0.5f;	// default 0.5 in the 2008 webinterface

	class Term {
		Cluster cluster= null;
		double pVal= -1;	// default pVal, according to David CB
		int cpos, call, upos, uall, bpos, ball;
		int count;
		String category, term;
		Integer termID;
		public Term(int cpos, int call, int upos, int uall, int bpos, int ball, String catAndTerm) {
			this.cpos= cpos;
			this.call= call;
			this.upos= upos;
			this.uall= uall;
			this.bpos= bpos;
			this.ball= ball;
			String[] cas= catAndTerm.split(attrGrpIDseparator);
			assert(cas.length== 2);
			category= cas[0];
			termID= Integer.parseInt(cas[1]);	
		}
	
		public boolean isValid() {
			if (cpos>= Cluster.minElementCount)
				return true;
			return false;
		}
		
		public double getPValCU() {
//			FisherExactTest fisherTest= new FisherExactTest();
//			fisherTest.fisherExact(cpos-1,upos,call-cpos,uall-upos);	// EASE score
//			double pVal= fisherTest.getRight();
			pVal= EventDavid.this.getFisherTest().getRightTailedP(cpos-1,upos,call-cpos,uall-upos);
			return pVal;
		}
			
		public double getPValCB() {
			//FisherExactTest fisherTest= new FisherExactTest();
			//fisherTest.fisherExact(cpos-1,bpos,call-cpos,ball-bpos);	// EASE score
			//double pVal= fisherTest.getRight();
			//pVal= EventDavid.this.getFisherTest().getRightTailedP(cpos-1,bpos,call-cpos,ball-bpos);
			pVal= EventDavid.this.getFisherTest().getRightTailedP(cpos-1,call-cpos,bpos- cpos,ball- (bpos-cpos));

			return pVal;
		} // http://mbe.oxfordjournals.org/content/vol23/issue12/images/large/molbiolevolmsl111f01_lw.jpeg
		
		public double getPValDavid() {
			
			// a b
			// c d
			int a= upos- 1;
			int b= uall- a;
			int c= bpos- a;
			int d= ball- bpos- b;
			
			if (a< 0|| b< 0|| c< 0|| d< 0)
				System.currentTimeMillis();
			
			try {
				//FisherExact tst= new FisherExact(a+b+c+d);
				// outputs often 0 instead of 1
				// e.g., for (2123 , 5428 , 16153 , 5984119)
				//pVal= EventDavid.this.getFisherTest().getRightTailedP(a, b, c, d);
				//double pTest= EventDavid.this.getFisherTest().getRightTailedP(18 , 314 , 196 , 14431);
			} catch (Exception e) {
				System.currentTimeMillis();
			}
			// check.fisherExact(18 , 314 , 196 , 14431);
			// double chktst= check.getRight();
			getMyFisher().reset();
			getMyFisher().fisherExact(a, b, c, d);
			pVal= getMyFisher().getRight();
			
//			if (pVal== 0) {
//				float fu= (upos/ (float) uall);
//				float fb= (bpos/ (float) ball);
//				//if (fu> 0.1&& fb> 0.1)
//				System.err.println(fu+" vs "+fb+"("+a+","+b+","+c+","+d+")");
//			}
			return pVal;
		}
		
		FisherExactTest myFisher;
		private FisherExactTest getMyFisher() {

			if (myFisher == null) {
				myFisher = new FisherExactTest();
			}

			return myFisher;
		}

		public double getPValUB() {
//			FisherExactTest fisherTest= new FisherExactTest();
//			fisherTest.fisherExact(upos-1,bpos,uall-upos,ball-bpos);	// EASE score
//			double pVal= fisherTest.getRight();
			//pVal= EventDavid.this.g7030	etFisherTest().getRightTailedP(upos-1,uall-upos,ball-bpos,bpos);
			int a= upos- 1, b= uall-upos, c= (ball-(bpos-upos)), d= bpos-upos;
			pVal= EventDavid.this.getFisherTest().getRightTailedP(upos- 1, bpos-upos, uall-upos, (ball-(bpos-upos)));
			/*double x2= EventDavid.this.getFisherTest().getRightTailedP(a,b,d,c);
			double x3= EventDavid.this.getFisherTest().getRightTailedP(a,d,b,c);
			*/ // works
			/* double x1= EventDavid.this.getFisherTest().getRightTailedP(a,b,c,d);
			double x11= EventDavid.this.getFisherTest().getRightTailedP(b,a,c,d);
			double x12= EventDavid.this.getFisherTest().getRightTailedP(a,c,b,d);
			double x122= EventDavid.this.getFisherTest().getRightTailedP(b,c,a,d);
			double x13= EventDavid.this.getFisherTest().getRightTailedP(c,a,b,d);
			double x14= EventDavid.this.getFisherTest().getRightTailedP(c,b,a,d);
			double x22= EventDavid.this.getFisherTest().getRightTailedP(b,a,d,c);
			double x32= EventDavid.this.getFisherTest().getRightTailedP(b,d,a,c);
			double x4= EventDavid.this.getFisherTest().getRightTailedP(d,a,b,c);
			double x42= EventDavid.this.getFisherTest().getRightTailedP(d,b,a,c); */ // works not

			return pVal;
		}
		
		public double getPValCUB() {
//			FisherExactTest fisherTest= new FisherExactTest();
//			fisherTest.fisherExact(cpos-1,bpos,uall-cpos,ball-bpos);	// EASE score
//			double pVal= fisherTest.getRight();
			pVal= EventDavid.this.getFisherTest().getRightTailedP(cpos-1,uall-cpos,bpos,ball-bpos);
			return pVal;
		}
		
		public double getDefaultPVal() {
			if (pVal < 0) {
				pVal = getPValDavid();	// getPValCB();				
			}
			if (pVal<= 0) {
				System.currentTimeMillis();
//				getPValDavid();
//				getPValDavid();
			}

			return pVal;
		}
		
		@Override
		public String toString() {
			String s= 
				formatDouble(getDefaultPVal())+"\t"+cpos+"\t"+(upos-cpos)+"\t"+(cpos* 100f/call)+"\t"+(((float) cpos*ball)/(call*bpos))+"\t"
					+category+"\t"+term+"\t"+
				call+"\t"+upos+"\t"+uall+"\t"+bpos+"\t"+ball+"\t"+formatDouble(pVal)
				
				/*formatDouble(getPValCUB())+"\t"+formatDouble(getPValUB())+"\t"+formatDouble(getPValCU())+"\t"+
				(cpos*100f/uall)+"\t"+(((float) cpos*ball)/(uall*bpos))+"\t"+(upos*100f/uall)+"\t"+(((float)upos*ball)/(uall*bpos))*/
				+"\n";
			return s;
		}

		String formatDouble(double x) {
			return formatDouble(x, 7); // seven for heaven !
		}
		
		String formatDouble(double x, int pos) {
			String s= Double.toString(x);
			if (s.length()<= pos)
				return s;
			int epos= s.indexOf('E');
			if (epos>= 0) {
				int pfxPos= pos- (s.length()- epos);
				return (s.substring(0, pfxPos)+ s.substring(epos));	// could be rounded
			} 
			
			// else
			return s.substring(0, pos);			
		}
		
	}
	
	static class Cluster {
		static int clusterCtr= 0;
		static int minElementCount= 2; // lower bound how often ("how many genes") a term has to exhibit, there is a pb with "1" and the ease score.. 
		
		static ClusterByEScoreComparator defaultClusterByEScoreComparator= new ClusterByEScoreComparator();
		static class ClusterByEScoreComparator implements Comparator<Cluster> {
			//@Override
			public int compare(Cluster o1, Cluster o2) {
				if (o1.eScore< o2.eScore)
					return -1;
				if (o1.eScore> o2.eScore)
					return 1;
				return 0;
			}
		}
		
		int nr;
		double eScore;
		Term[] terms;
		HashMap<String,Vector<Integer>> mapReplace;
		String[] memberIDs= null;
		String geneTerms= null;
		
		public Cluster(double eScore, Term[] members) {
			this.nr= ++clusterCtr;
			this.eScore= eScore;
			
			this.terms= members;
			Arrays.sort(members, getDefaultTermByPValComparator());
			
			for (int i = 0; i < members.length; i++) {
				members[i].cluster= this;
				Vector<Integer> v= getMapReplace().get(members[i].category);
				if (v== null) {
					v= new Vector<Integer>();
					getMapReplace().put(members[i].category, v);
				}
				v.add(members[i].termID);
			}
		}
		
		public HashMap<String,Vector<Integer>> getMapReplace() {
			if (mapReplace == null) {
				mapReplace = new HashMap<String,Vector<Integer>>();				
			}

			return mapReplace;
		}
		
		public int getValidMembers() {
			int ctr= 0;
			for (int i = 0; i < terms.length; i++) {
				if (terms[i].isValid())
					++ctr;
			}
			return ctr;
		}
		
		public void replace(HashMap<String,TermReplacement[]> replacementMap) {
			for (int i = 0; i < terms.length; i++) {
				TermReplacement[] rep= replacementMap.get(terms[i].category);
				TermReplacement dummy= new TermReplacement(terms[i].termID, null);
				int p= Arrays.binarySearch(rep, dummy, TermReplacement.defaultTermReplacementByIDComparator);
				terms[i].term= rep[p].term;
			}
		}
				
		@Override
		public String toString() {
			int validNr= getValidMembers();
			StringBuffer sb= new StringBuffer("CLUSTER_ID "+nr+", eScore "+eScore+
					", members "+ memberIDs.length+
					", contributing terms "+validNr+" each joining >="+
					minElementCount+" members (in total "+terms.length+" terms in cluster)\n");
			if (geneTerms!= null) {
				sb.append("GENES:\n");
				sb.append(geneTerms);
			}
			sb.append("\nMEMBERS: ");
			for (int i = 0; i < memberIDs.length; i++) {
				sb.append(memberIDs[i]+", ");
			}
			sb.append("\nTERMS:\n");
			sb.append(TERM_HEADER);
			for (int i = 0; i < terms.length; i++) 
				if (terms[i].isValid())
					sb.append(terms[i].toString());
			return sb.toString();
		}
	}
	
	static class TermReplacement {
		
		static class TermReplacementByIDComparator implements Comparator<TermReplacement> {
			//@Override
			public int compare(TermReplacement o1, TermReplacement o2) {
				if (o1.termID< o2.termID)
					return -1;
				if (o1.termID> o2.termID)
					return 1;
				return 0;
			}
		}
		
		static TermReplacementByIDComparator defaultTermReplacementByIDComparator= new TermReplacementByIDComparator();
		
		int termID;
		String term;
		public TermReplacement(int termID) {
			this.termID= termID;
		}	
			
		public TermReplacement(int termID, String term) {
			this(termID);
			this.term= term;
		}
	}
	
	
	static class TermByPValComparator implements Comparator<Term> {
		//@Override
		public int compare(Term o1, Term o2) {
			if (o1.getDefaultPVal()< o2.getDefaultPVal())
				return -1;	
			if (o1.getDefaultPVal()> o2.getDefaultPVal())
				return 1;
			return 0;
		}
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		EventDavid dave= parseArguments(args);
		dave.run();
	}
	
	static EventDavid parseArguments(String[] args) {
		
		File inFile= null, bgFile= null, outFile= null, attrFile= null, annFile= null;
		HashMap<String, Float> attrKeys= null;
		String rowIDkey= null;
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-in")) {
				inFile= new File(args[++i]);
				continue;
			}
			if (args[i].equals("-out")) {
				outFile= new File(args[++i]);
				continue;
			}
			if (args[i].equals("-bg")) {
				bgFile= new File(args[++i]);
				continue;
			}
			if (args[i].equals("-ann")) {
				annFile= new File(args[++i]);
				continue;
			}
			if (args[i].equals("-attr")) {
				attrFile= new File(args[++i]);
				continue;
			}
			if (args[i].equals("-key")) {
				rowIDkey= args[++i];
				continue;
			}
			
		}
		
		if (rowIDkey== null)
			rowIDkey="GENBANK_ACCESSION";
		attrKeys= readAttrKeys(attrFile);
		
		if (inFile== null|| annFile== null|| attrKeys== null|| rowIDkey== null) {
			System.err.println("[FATAL] missing crucial parameters.");
			return null;
		}
		if (outFile== null) {
			String in= inFile.getName();
			String bg= "all";
			if (bgFile!= null)
				bg= bgFile.getName();
			String path= inFile.getParentFile().getAbsolutePath();
			outFile= new File(path+File.separator+in+"_"+bg);
		}
		
		EventDavid dave= new EventDavid(inFile, bgFile, annFile, outFile);
		dave.setAttrGrpIDs(attrKeys);
		dave.setRowIDkey(rowIDkey);
		return dave;
	}
	

	private static HashMap<String, Float> readAttrKeys(File attrFile) {
		
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(attrFile));
			String s;
			HashMap<String, Float> attrKeys= new HashMap<String, Float>();
			while ((s= buffy.readLine())!= null) {
				String[] ss= s.split("\t");
				if(ss.length==1) 
					attrKeys.put(ss[0], 1f);
				else
					attrKeys.put(ss[0], Float.parseFloat(ss[1]));
			}
			buffy.close();
			return attrKeys;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	static void testKappa() {
		output2D(testBitMore, System.out);
		EventDavid myKap = new EventDavid();
		// int[][][] conting= myKap.initContingencyTable(testLittle);
		// output2D(conting, System.out);
		float[][] kappa = myKap.getKappa(testBitMore);
		output2D(kappa, System.out);
	}

	static void testCluster() {
		EventDavid kap = new EventDavid();
		PreCluster[] clusters = kap.cluster(testKappaTbl);
		System.out.println("found " + clusters.length + " clusters:");
		for (int i = 0; i < clusters.length; i++) {
			for (int j = 0; j < clusters[i].a.length; j++) {
				System.out.print(clusters[i].a[j]);
				System.out.print("\t");
			}
			System.out.println();
		}
	}

	// could be binary encoded to save space
	byte[][] redTable; // (non?)redundant, 1= annotated, 0= not annotated, -1= not at all annotated in category
	
	File inFile, refFile, bgFile, annFile, outFile;
	String[] attrIDs, rowIDs;	// row/col IDs of redTable
	String rowIDkey;			// key value that identifies row
	HashMap<String, Float> attrGrpIDs;	// attrGrpIDs that are used
	HashMap<String, Integer> mapCountBg;	// maps attributeGrpID/attribute -> count
	HashMap<String, String> mapAttr2Grp;	// maps attribute to group ID
	HashSet[] idLists;	// last is geneID (?! will see)
	
	static final String attrGrpIDseparator= "@";
	
	private EventDavid() {
	}

	public EventDavid(File inFile, File bgFile, File annFile, File outFile) {
		this.inFile = inFile;
		this.bgFile = bgFile;
		this.annFile= annFile;
		this.outFile= outFile;
	}

	/**
	 * runs in |rows|^2 * |columns| time (check for faster non-naive algo!!) and
	 * (4*|rows|^2) memory (output size, does not replace input; 2 fields of
	 * inner array maybe redundant).
	 * 
	 * @param rowFeatureTable
	 * @return row-row similarites in 4 fields, [0]= both 0, [1]= row 0 and col
	 *         1, [2] row 1 and col 0, [3] both 1 has to be int, otherwise limit
	 *         to 255 common members of feature!
	 * @deprecated
	 */
	int[][][] getContingencyTable(byte[][] rowFeatureTable) {

		int[][][] sim = new int[rowFeatureTable.length][][];
		for (int i = 0; i < rowFeatureTable.length; i++) {
			sim[i] = new int[rowFeatureTable.length][];
			for (int j = 0; j < sim.length; j++) {
				sim[i][j] = new int[4];
				for (int k = 0; k < sim[i][j].length; k++)
					sim[i][j][k] = 0;
			}

			for (int j = 0; j < rowFeatureTable[i].length; j++) {
				assert (rowFeatureTable[i][j] == 1 || rowFeatureTable[i][j] == 0);

				for (int k = 0; k < rowFeatureTable.length; k++) {
					if (k == i)
						continue;

					int idx = (rowFeatureTable[i][j] > 0) ? ((rowFeatureTable[k][j] > 0) ? 3
							: 2)
							: ((rowFeatureTable[k][j] > 0) ? 1 : 0);
					++sim[i][k][idx];
				}
			}
		}

		return sim;
	}

	/**
	 * @deprecated
	 * @param simTable
	 * @return
	 */
	float[][] getKappa(int[][][] simTable) {
		int sumTot = 0;
		int[][] rowTotal = new int[simTable.length][], colTotal = new int[simTable.length][];
		for (int i = 0; i < simTable.length; i++) {
			colTotal[i] = new int[2];
			rowTotal[i] = new int[2];
			for (int j = 0; j < colTotal.length; j++) {
				colTotal[i][j] = 0;
				rowTotal[i][j] = 0;
			}
		}

		for (int i = 0; i < simTable.length; i++) {
			for (int j = 0; j < simTable[i].length; j++) {
				for (int k = 0; k < simTable[i][j].length; k++) {
					if (k == 0 || k == 3)
						colTotal[j][k % 2] += simTable[i][j][k];
					rowTotal[j][k / 2] += simTable[i][j][k];
					if (j > i) // only one half of the matrix, symetric
						sumTot += simTable[i][j][k];
				}
			}
		}

		// names lt. file:///D:/projects/claudia/david/help_kappa.htm
		float[][] kappa = new float[simTable.length][];
		for (int i = 0; i < kappa.length; i++) {
			kappa[i] = new float[simTable.length];
			for (int j = 0; j < kappa[i].length; j++) {
				if (i == j)
					continue;
				float O_ab = (simTable[i][j][0] + simTable[i][j][3])
						/ ((float) sumTot); // TODO check whether it is really
											// sum of all, not only gene
				float A_ab = ((rowTotal[i][0] * colTotal[i][0]) + (rowTotal[i][1] * colTotal[i][1]))
						/ ((float) sumTot * sumTot); // latter can be precalc
														// for speedup
				float K_ab = (O_ab - A_ab) / (1f - A_ab);
				kappa[i][j] = K_ab;
			}
		}

		return kappa;
	}

	private int[][] tmpConTable = new int[3][3];	// is reused
	int[][] getConTable(byte[] rowA, byte[] rowB) {
		assert (rowA.length == rowB.length);
		for (int i = 0; i < tmpConTable.length; i++)
			for (int j = 0; j < tmpConTable[i].length; j++)
				tmpConTable[i][j] = 0;
		for (int i = 0; i < rowB.length; i++) {
			if (rowA[i]==2|| rowB[i]==2)	// not annotated
				continue;
			++tmpConTable[rowA[i]][rowB[i]];
			++tmpConTable[rowA[i]][2];
			++tmpConTable[2][rowB[i]];
		}
		return tmpConTable;
	}

	/**
	 * filters terms off that are weakly associated with the foreground or 
	 * background.
	 * @param redTbl
	 * @param mapCountBg
	 */
	void preFilter(byte[][] redTbl, HashMap<String, Integer> mapCountBg) {
		
		System.err.print("["+Constants.getTime()+"] pre-filtering ");
		int perc= 0, eTerms= 0, eCat= 0, kTerms= 0;
		for (int i = 0; i < redTbl[0].length; i++) {

			if (i* 10d/ redTbl.length> perc) {
				System.err.print("*");
				System.err.flush();
				++perc;
			}
			
			String attr= attrIDs[i];
			if (attr.indexOf(attrGrpIDseparator)< 0)
				continue;	// skip categories
			
			String[] grpAttID= attr.split(attrGrpIDseparator);
			int bpos= mapCountBg.get(attr);	// background positive
			int ball= mapCountBg.get(grpAttID[0]);	// background all
			
//			if (ball== 0)
//				continue;
			
			float fb= bpos/ (float) ball, fu= 0;
			int upos= 0, uall= 0;

			for (int j = 0; j < redTbl.length; j++) {
				
				if (redTbl[j][i]== 1) 	
					++upos;
				if (redTable[j][i]< 2) 	// analogously, if (redTable[j][i]== 0|| redTable[j][i]== 1) {
					++uall;
			}
			
			if (uall!= 0)
				fu= upos/ (float) uall;
//			else
//				continue;
			
			// eliminate term
			double pVal= 0;
			if (upos> 0&& uall> 0&& bpos> 0) {
				Term t= new Term(0,0,upos,uall,bpos,ball,"1@1");
				pVal= t.getPValDavid();
			}
			if (pVal== 0) {// || fu< filter_min_association) {	
				
				++eTerms;
				
				--bpos;
				if (bpos== 0) 
					mapCountBg.remove(attr);
				--ball;
				if (ball== 0)
					mapCountBg.remove(grpAttID[0]);
				else
					mapCountBg.put(grpAttID[0], ball);
				--upos; 
				
				int eTermCat= 0;
				for (int j = 0; j < redTbl.length; j++) {
					if (redTbl[j][i]!= 1)
						continue;
					// count how often this term is annotated in the category
					int cat= 0;
					for (int k = 0; k < attrIDs.length; k++) {
						if (k== i)
							continue;
						int p= attrIDs[k].indexOf(attrGrpIDseparator);
						if (p<= 0)
							continue;
						if (attrIDs[k].substring(0,p).equals(grpAttID[0])&& redTbl[j][k]== 1) 
							++cat;
					}
					if (cat== 0)
						redTbl[j][i]= 0;
					else {	
						++eTermCat;
						redTbl[j][i]= 2;
						for (int k = 0; k < attrIDs.length; k++) {
							if (k== i)
								continue;
							int p= attrIDs[k].indexOf(attrGrpIDseparator);
							if (p<= 0)
								continue;
							if (attrIDs[k].substring(0,p).equals(grpAttID[0])) 
								redTbl[j][k]= 2;
						}
					}
				}
				
				if (ball== 0|| eTermCat> 0)
					++eCat;
				
			} else
				++kTerms;
		}
		
		System.err.println("\n\teliminated "+eTerms+" terms and "+eCat+" categories.");
		System.err.println("\tkept "+kTerms+" terms.");
	}
	
	
	float[][] getKappa(byte[][] depTbl) {

		float[][] kapTbl = null;
		if ((kapTbl= readKappa())!= null)
			return kapTbl;
		
		for (int i = 0; i < depTbl.length; i++) {
			for (int j = 0; j < depTbl[i].length; j++) {
				if (depTbl[i][j]!= 0)
					System.currentTimeMillis();
			}
		}
		
		kapTbl = new float[depTbl.length][];
		System.err.print("["+Constants.getTime()+"] Calculating kappa table ");
		for (int i = 0, perc= 0; i < depTbl.length; i++) {
			if (i* 10d/ depTbl.length> perc) {
				System.err.print("*");
				System.err.flush();
				++perc;
			}
			for (int j = i + 1; j < depTbl.length; j++) {
				
				// count term overlap between members
				int cnt= 0;
				for (int x = 0; x < depTbl[i].length; x++) {
					if (depTbl[i][x]== 1&& depTbl[j][x]== 1)
						++cnt;
				}
				
				float K_ab = -1;
				if (cnt> kappa_similarityOvl) {
					int[][] conTbl = getConTable(depTbl[i], depTbl[j]);
					int sumTot = 0;
					for (int k = 0; k < conTbl.length; k++)
						sumTot += conTbl[k][2]; // can also be row scores
	
					float O_ab = (conTbl[0][0] + conTbl[1][1]) / ((float) sumTot);
					float A_ab = ((conTbl[0][2] * conTbl[2][0]) + (conTbl[1][2] * conTbl[2][1]))
							/ ((float) sumTot * sumTot); // latter can be precalc
															// for speedup
					K_ab = 1;
					if (A_ab!= 1)
						K_ab= (O_ab - A_ab) / (1f - A_ab); // K_ab gets 1 if everything is aggreeing or disagreeing perfectly.
	//				else								// consequently also O_ab is one in this case and there is "0/0"
	//					System.currentTimeMillis();
				}
				// else it stays at -1 which should be less than kappa_similarity threshold s.t. they dont get clustered
				
				kapTbl[i] = (kapTbl[i] == null) ? new float[depTbl.length]
						: kapTbl[i];
				kapTbl[j] = (kapTbl[j] == null) ? new float[depTbl.length]
						: kapTbl[j];				
				kapTbl[i][j] = kapTbl[j][i] = K_ab;	// mirror
			}
		}

		saveKappa(kapTbl);
		
		return kapTbl;
	}

	private float[][] readKappa() {
		
		Vector<float[]> kapV= new Vector<float[]>();
		int chk= -1;
		System.err.print("["+Constants.getTime()+"] Reading kappa table ");
		try {
			String fName= getKappaFName();
			File f= new File(fName);
			if (!f.exists())
				return null;
			BufferedReader buffy= new BufferedReader(new FileReader(fName));
			int ctr= 0, perc= 0;
			long bytes= 0, total= new File(fName).length();
			for (String s = null; (s= buffy.readLine())!= null;++ctr) {
				bytes+= s.length();
				if (bytes* 10d/ total> perc) {
					System.err.print("*");
					System.err.flush();
					++perc;
				}
				String[] ss= s.split("\t");
				if (chk< 0)
					chk= ss.length;
				else if(chk!= ss.length)
					System.err.println("\n\t[ERROR] line "+ctr+" does not have "+chk+" elements.");
				float[] ff= new float[ss.length];
				for (int i = 0; i < ff.length; i++) 
					ff[i]= Float.parseFloat(ss[i]);
				kapV.add(ff);
			}
			buffy.close();
			
			float[][] kapTbl= new float[kapV.size()][];
			for (int i = 0; i < kapTbl.length; i++) 
				kapTbl[i]= kapV.elementAt(i);
			
			return kapTbl;
			
		} catch (Exception e) {
			System.err.println();
			e.printStackTrace();
		}
		
		return null;
	}

	private String getKappaFName() {
		String fname= inFile.getParent()
				+ File.separator
				+ FileHelper.getFileNameWithoutExtension(inFile.getAbsolutePath())+
				SFX_KAPPA+ SFX_TBL; 
		return fname;
	}
	
	private void saveKappa(float[][] kapTbl) {
		
		try {
			BufferedWriter writer= new BufferedWriter(new FileWriter(getKappaFName()));
			for (int i = 0; i < kapTbl.length; i++) {
				StringBuilder sb= new StringBuilder();
				for (int j = 0; j < kapTbl[i].length; j++) {
					sb.append(Float.toString(kapTbl[i][j]));
					sb.append("\t");
				}
				sb.deleteCharAt(sb.length()- 1);
				sb.append("\n");
				writer.write(sb.toString());
			}
			writer.flush();
			writer.close();
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	static void output2D(int[][][] tab, PrintStream p) {
		for (int i = 0; i < tab.length; i++) {
			for (int j = 0; j < tab[i].length; j++) {
				for (int k = 0; k < tab[i][j].length; k++) {
					p.print(Integer.toString(tab[i][j][k]));
					if (k < tab[i][j].length - 1)
						p.print(",");
				}
				p.print("\t");
			}
			p.println();
		}
	}

	static void output2D(byte[][] tab, PrintStream p) {
		for (int i = 0; i < tab.length; i++) {
			for (int j = 0; j < tab[i].length; j++) {
				p.print(Byte.toString(tab[i][j]));
				p.print("\t");
			}
			p.println();
		}
	}

	static void output2D(float[][] tab, PrintStream p) {
		for (int i = 0; i < tab.length; i++) {
			for (int j = 0; j < tab[i].length; j++) {
				p.print(Float.toString(tab[i][j]));
				p.print("\t");
			}
			p.println();
		}
	}

	/**
	 * Create qualified initial seeding groups: Each gene could form a initial
	 * seeding group (initial seeds) as long as it has close relationships (e.g.
	 * kappa >=0.35) with more than > 2 other members (i.e. 'initial group
	 * membership' threshold in DAVID interface) . In order to control the
	 * quality of the seeding groups, the qualified seeding groups (qualified
	 * seeds) need to meet the second condition, i.e. majority (>50%) of members
	 * in the seed should have close relationships (e.g. kappa >= 0.35) each
	 * other.
	 * 
	 * @param kappaTbl
	 * @return
	 */
	HashSet<PreCluster> getInitialSeeds(float[][] kappaTbl) {

		HashSet<PreCluster> seedMap = new HashSet<PreCluster>();
		int[] tmpMembers = new int[kappaTbl.length];
		for (int i = 0; i < kappaTbl.length; i++) {
			int tmpMemberCtr = 0;
			tmpMembers[tmpMemberCtr++] = i;

			// collect members according to close relationship threshold
			for (int j = 0; j < kappaTbl.length; j++) {
				if (i == j)
					continue;
				if (kappaTbl[i][j] >= kappa_similarityThr)
					tmpMembers[tmpMemberCtr++] = j;
			}

			// membership nr threshold
			if (tmpMemberCtr <= clustering_initialGroupMembershipThr)
				continue;

			// check for quality of qualified seeds
			int ctr = 0;
			for (int j = 1; j < tmpMemberCtr; j++) {
				for (int k = (j + 1); k < tmpMemberCtr; k++) {
					if (kappaTbl[tmpMembers[j]][tmpMembers[k]] >= kappa_similarityThr)
						++ctr;
				}
			}
			int totRelations = (tmpMemberCtr - 1) * (tmpMemberCtr - 2) / 2;
			float closeRatio = ((float) ctr) / totRelations;
			if (closeRatio <= clustering_nbCloseThr)
				continue;

			// create new seed
			int[] seed = new int[tmpMemberCtr];
			System.arraycopy(tmpMembers, 0, seed, 0, tmpMemberCtr);
			Arrays.sort(seed); // sorting can eventually be speed up -> move [0]
								// to correct pos, finshed.

			PreCluster a = new PreCluster(seed);
			seedMap.add(a);
		}

		return seedMap;
		// int[][] seeds= new int[seedMap.size()][];
		// int ctr= 0;
		// Iterator<int[]> iter= seedMap.iterator();
		// while (iter.hasNext())
		// seeds[ctr++]= iter.next();
		// return seeds;
	}

	PreCluster[] cluster(float[][] kappaTbl) {
		// IntArray.normal= false;
		HashSet<PreCluster> clusterSet = getInitialSeeds(kappaTbl);
		// IntArray.normal= true;
		PreCluster[] tmpClusterSet = new PreCluster[clusterSet.size()];
		boolean smaller = true;
		while (smaller) {
			int before = clusterSet.size();
			tmpClusterSet = clusterSet.toArray(tmpClusterSet);
			Outer: for (int i = 0; i < tmpClusterSet.length
					&& tmpClusterSet[i] != null; i++) {
				for (int j = i + 1; j < tmpClusterSet.length
						&& tmpClusterSet[j] != null; j++) {
					PreCluster shorter = (tmpClusterSet[i].a.length < tmpClusterSet[j].a.length) ? tmpClusterSet[i]
							: tmpClusterSet[j];
					PreCluster longer = (shorter == tmpClusterSet[j]) ? tmpClusterSet[i]
							: tmpClusterSet[j];
					int idCtr = 0, lo = 0;
					for (int k = 0; k < shorter.a.length; k++) {
//						int p = Arrays.binarySearch(longer.a, lo,
//								longer.a.length, shorter.a[k]);
						int p = Arrays.binarySearch(longer.a, shorter.a[k]);
						if (p >= 0) {
							++idCtr;
							lo = (p + 1);
						} else
							lo = -(p + 1);
					}

					int totGrp = shorter.a.length + longer.a.length - idCtr;
					float linkRatio = idCtr / ((float) totGrp);
					if (linkRatio > clustering_multiLinkThr) {
						int[] newGrp = new int[totGrp];
						int shortC = 0, longC = 0;
						for (int k = 0; k < newGrp.length; k++) {
							if (shortC == shorter.a.length) {
								newGrp[k] = longer.a[longC++];
								continue;
							}
							if (longC == longer.a.length) {
								newGrp[k] = shorter.a[shortC++];
								continue;
							}
							if (shorter.a[shortC] == longer.a[longC]) {
								newGrp[k] = shorter.a[shortC++];
								++longC;
							} else
								newGrp[k] = (shorter.a[shortC] < longer.a[longC]) ? shorter.a[shortC++]
										: longer.a[longC++];
						}
						clusterSet.remove(shorter);
						clusterSet.remove(longer);
						PreCluster newA = new PreCluster(newGrp);
						// IntArray.normal= false;
						clusterSet.add(newA);
						// IntArray.normal= true;
						break Outer;
					}
				}
			}

			int after = clusterSet.size();
			smaller = (after < before);

		}

		// final "clean-up" procedure, remove clusters with not enough genes
		PreCluster[] clusters = new PreCluster[clusterSet.size()];
		clusterSet.toArray(clusters);
		for (int i = 0; i < clusters.length; i++) {
			if (clusters[i].a.length< clustering_finalGroupMembershipThr)
				clusterSet.remove(clusters[i]);
		}
				
		clusters = new PreCluster[clusterSet.size()];
		clusterSet.toArray(clusters);
		
		System.err.println(" found "+ clusters.length+ " clusters.");
		return clusters;
	}

	// for outliers, members that participate in multiple clusters..
	static public int[] getClusterAssociations(int len, PreCluster[] clusters) {
		int[] data = new int[len];
		for (int i = 0; i < data.length; i++)
			data[i] = 0;
		for (int i = 0; i < clusters.length; i++) {
			for (int j = 0; j < clusters[i].a.length; j++) {
				++data[clusters[i].a[j]];
			}
		}
		return data;
	}
	
	public static TermByPValComparator getDefaultTermByPValComparator() {
		return defaultTermByPValComparator;
	}


	/**
	 * @deprecated
	 * @param tmpClusterSet
	 */
	void cluster(int[][] tmpClusterSet) {
		for (int i = 0; tmpClusterSet[i] != null; i++) {
			for (int j = i + 1; tmpClusterSet[j] != null; j++) {
				int[] shorter = (tmpClusterSet[i].length < tmpClusterSet[j].length) ? tmpClusterSet[i]
						: tmpClusterSet[j];
				int[] longer = (shorter == tmpClusterSet[j]) ? tmpClusterSet[i]
						: tmpClusterSet[j];
				int idCtr = 0, lo = 0;
				for (int k = 0; k < shorter.length; k++) {
//					int p = Arrays.binarySearch(longer, lo, longer.length,
//							shorter[k]);
					int p = Arrays.binarySearch(longer, shorter[k]);
					if (p >= 0) {
						++idCtr;
						lo = (p + 1);
					} else
						lo = -(p + 1);
				}

				int totGrp = shorter.length + longer.length - (2 * idCtr);
				float linkRatio = idCtr / ((float) totGrp);
				if (linkRatio > clustering_multiLinkThr) {

				}
			}
		}
	}

	static void testFisher() {
		int[] data = new int[] { 2, 40, 297, 29960 };
		int[] oldData = new int[4];
		FisherExactTest fisherman = new FisherExactTest();
		fisherman.console= true;
		fisherman.fisherExact(3,1,1,3); // right~0.24 (webpage defaults)
		fisherman.fisherExact(2, 40, 297, 29960);	// right~0.06	(David example)
		fisherman.fisherExact(18, 314, 196, 14431);	// right~1.4E-6 (David example)
	}

	
	// private FisherExactTest fisherTest= new FisherExactTest();
	private FisherExact fisherTest= null;
	private FisherExact getFisherTest() {
		if (fisherTest == null) {
			int bmax= 0;
			for (int i = 0; i < attrIDs.length; i++) {
				int bpos= mapCountBg.get(attrIDs[i]);	// background positive
				String[] grpAttID= attrIDs[i].split(attrGrpIDseparator);
				int ball= mapCountBg.get(grpAttID[0]);	// background all
				bmax= Math.max(bpos+ ball, bmax);	
			}
			
			int cmax= 0;
			for (int i = 0; i < redTable[0].length; i++) {
				for (int j = 0; j < redTable.length; j++) {
					int call= 0, cpos= 0;
					if (redTable[j][i]>= 0) {
						if (redTable[j][i]== 1) 	
							++cpos;
						if (redTable[j][i]< 2) 
							++call;
					}
					cmax= Math.max(cmax, cpos+call);
				}
			}

			fisherTest = new FisherExact(bmax+ cmax);
			
		}

		return fisherTest;
	}

	
	private int clusterCtr= 0;
	/**
	 * creates a cluster from a precluster
	 * @param preCluster
	 * @return
	 */
	Cluster createCluster(PreCluster preCluster) {
		
		// get features
		int[] present= new int[redTable[0].length];	// count annotated terms in cluster
		for (int i = 0; i < present.length; i++) 
			present[i]= 0;
		int ctr= 0, ctr2= 0;
		String[] members= new String[preCluster.a.length];
		for (int i = 0; i < preCluster.a.length; i++) {	// iterate over MEMBERS (genes, transcripts, events..)
			members[i]= rowIDs[preCluster.a[i]];
			for (int j = 0; j < redTable[preCluster.a[i]].length; j++) {
				present[j]+= (redTable[preCluster.a[i]][j]== 1)?1:0;	// only matches
				if (i== preCluster.a.length-1&& present[j]> 0) {	// >= minCount, if neglecting members without statistics
					++ctr;
					if (present[j]>= Cluster.minElementCount)
						++ctr2;
				}
			}
		}
		
		// EASE scores
		Vector<Term> elementV= new Vector<Term>(ctr);	// we want all members in the cluster..
		ctr= 0;
		for (int i = 0; i < present.length; i++) {	// iterate over TERMS
			if (present[i]== 0)	// < minCount, if neglecting members without statistics
				continue;
			int cPtr= 0;
			int cpos= 0, call= 0, upos= 0, uall= 0; // cluster positive/all, user positive/all
			for (int j = 0; j < redTable.length; j++) {
				
				if (redTable[j][i]>= 0) {
					if (redTable[j][i]== 1) {	
						++upos;
						if (cPtr< preCluster.a.length&& j== preCluster.a[cPtr]) 
							++cpos;
					} 
					if (redTable[j][i]< 2) {	// analogously, if (redTable[j][i]== 0|| redTable[j][i]== 1) {
						++uall;
						if (cPtr< preCluster.a.length&& j== preCluster.a[cPtr]) {
							++call;
							++cPtr;
						}
					} else { // 2= unannotated
						if (cPtr< preCluster.a.length&& j== preCluster.a[cPtr]) 
							++cPtr;
					}
				}
			}
			assert(cpos== present[i]);	// that should hold
			//String attr= attrIDs[i];
			int bpos= mapCountBg.get(attrIDs[i]);	// background positive
			String[] grpAttID= attrIDs[i].split(attrGrpIDseparator);
			int ball= mapCountBg.get(grpAttID[0]);	// background all

			// elements[ctr++]=
			//if ((upos/(float) uall)> 0.1&& (bpos/ (float) ball)> 0.1)
			Term t= new Term(cpos, call, upos, uall, bpos, ball, attrIDs[i]);
			if (t.getPValDavid()> 0)
				elementV.add(t);
			else
				System.err.println("pVal="+t.pVal+" ("+upos+","+uall+","+bpos+","+ball+")");
		}
		
		Cluster c= null;
		if (elementV.size()> 0) {
		
			Term[] elements= new Term[elementV.size()];
			for (int i = 0; i < elements.length; i++) 
				elements[i]= elementV.elementAt(i);

			// enrichment score
			// ok, for small pVals change to the sum of the logs
	//		double geoMean= 1;
	//		for (int i = 0; i < pVals.length; i++) 
	//			geoMean*= pVals[i];
	//		geoMean= Math.pow(geoMean, 1d/ctr);
	//		double eScore= -Math.log10(geoMean);
			double eScore= 0, radix= 1d/ctr2;	// .. but we only count the ones with pvalues for the escore
			for (int i = 0; i < elements.length; i++) {
				double pVal= elements[i].getDefaultPVal();
				if (!Double.isNaN(pVal))
					eScore+= -Math.log10(Math.pow(pVal,radix));
			}
			c= new Cluster(eScore, elements);
			c.memberIDs= members;
		
		} else
			System.currentTimeMillis();
		
		return c;
	}
	
	Cluster[] createClusters(PreCluster[] preClusters) {
		try {
			//BufferedWriter writer= new BufferedWriter(new FileWriter(outFile));
			//PrintStream p= new PrintStream(new BufferedOutputStream(new FileOutputStream(outFile)));
			Vector<Cluster> clusterV= new Vector<Cluster>(preClusters.length);
			System.err.print("["+new Time(System.currentTimeMillis())+"] creating clusters ");
			System.err.flush();
			for (int i = 0, perc= 0; i < preClusters.length; i++) {
				if (i* 10d/ preClusters.length> perc) {
					System.err.print("*");
					System.err.flush();
					++perc;
				}
				Cluster tmpCluster= createCluster(preClusters[i]);
				if (tmpCluster!= null)
					clusterV.add(tmpCluster);
			}
			
			Cluster[] clusters= new Cluster[clusterV.size()];
			for (int i = 0; i < clusters.length; i++) 
				clusters[i]= clusterV.elementAt(i);
			
			System.err.println(" "+ clusters.length+" clusters created.\n");
			//p.flush();
			//p.close();
			return clusters;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	
	void getInput() {
		getInput(false);
		getInput(true);
	}
	
	File redFile= null;
	void readBackground() {
		// try to get whole table
		File f= new File(DEFAULT_BACKTABLE_PATH);
		if (f.exists()) {
			System.err.print("["+Constants.getTime()+"] Reading background data ");
			System.err.flush();
			try {
				BufferedReader buffy= new BufferedReader(new FileReader(f));
				mapCountBg= new HashMap<String, Integer>();	
				String s;
				int perc= 0;
				long len= f.length(), bread= 0;
				while((s= buffy.readLine())!= null) {
					bread+= s.length()+1;
					if ((bread*10d/ len)> perc) {
						System.err.print("*");
						System.err.flush();
						++perc;
					}
					String[] ss= s.split("\t");		// id x count
					assert(ss.length== 2);
					int x= -1;
					try {
						x= Integer.parseInt(ss[1]);
					} catch (Exception e) {
						; // :)
					}
					// TODO check consistency here with foreground..
					mapCountBg.put(ss[0],x);
				}
				System.err.println(" "+mapCountBg.size()+" terms.\n");
				buffy.close();
				return;
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// try to get IDs
		HashSet map= new HashSet();	// reuse later for Int GeneIDs
		if (bgFile != null) {
			System.err.println("[INFO] Reading background IDs.. ");
			try {
				BufferedReader buffy= new BufferedReader(new FileReader(bgFile));
				String s;
				long t0= System.currentTimeMillis();
				while ((s= buffy.readLine())!= null) {
					map.add(s.trim());
				}
				System.err.println("[YIHAA] found "+map.size()+" unique transcript IDs, "
						+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		if (map.size()== 0)
			map= null;	// take all IDs

		assert(attrGrpIDs!= null);	// TODO sync with order of reading foreground..

		readAnnotation(map, false, true);
				
		// write out table
		try {
			BufferedWriter writer= new BufferedWriter(new FileWriter(f));
			Iterator<String> iter= mapCountBg.keySet().iterator();
			String s= null;
			while (iter.hasNext()) {
				s= iter.next();
				writer.write(s);
				writer.write("\t");
				s= Integer.toString(mapCountBg.get(s));
				writer.write(s);
				writer.write("\n");
			}
			writer.flush();
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	void readBackground_save() {
		
		try {		
			
				// try to get whole table
			File f= new File(DEFAULT_BACKTABLE_PATH);
			if (f.exists()) {
				System.err.print("Reading background data.. ");
				System.err.flush();
				BufferedReader buffy= new BufferedReader(new FileReader(f));
				mapCountBg= new HashMap<String, Integer>();
				String s;
				int perc= 0;
				long len= f.length(), bread= 0;
				while((s= buffy.readLine())!= null) {
					bread+= s.length()+1;
					if ((bread*10d/ len)> perc) {
						System.err.print("*");
						System.err.flush();
						++perc;
					}
					String[] ss= s.split("\t");		// id x count
					assert(ss.length== 2);
					int x= -1;
					try {
						x= Integer.parseInt(ss[1]);
					} catch (Exception e) {
						; // :)
					}
					// TODO check consistency here with foreground..
					mapCountBg.put(ss[0],x);
				}
				System.err.println();
				buffy.close();
				return;
			}
			
				// try to get IDs
			HashSet map= new HashSet();	// reuse later for Int GeneIDs
			if (bgFile != null) {
				System.err.println("[INFO] Reading background IDs.. ");
				BufferedReader buffy= new BufferedReader(new FileReader(bgFile));
				String s;
				long t0= System.currentTimeMillis();
				while ((s= buffy.readLine())!= null) {
					map.add(s.trim());
				}
				System.err.println("[YIHAA] found "+map.size()+" unique transcript IDs, "
						+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
			}
			if (map.size()== 0)
				map= null;	// take all IDs
			
//			if (useIDs== null)
//				readAttributes();
			assert(attrGrpIDs!= null);	// TODO sync with order of reading foreground..
			
			
			// get background counts
			System.err.print("[ANDNOW] count from background.. ");
			BufferedReader buffy= new BufferedReader(new FileReader(annFile));

			String s= buffy.readLine();	// header
			IntegerCompound[] useIDs= readAnnotationHeader(s);			
			
			int rowCtr= 1, perc= 0;
			long t0= System.currentTimeMillis(), bcount= 0, flen= annFile.length();
			int cntPlus= 0, cntMinus= 0, redRowCtr= 0;
			String inFileDelim= "\t";
			mapCountBg= new HashMap<String, Integer>();
			while((s= buffy.readLine())!= null) {
				++rowCtr;
				bcount+= s.length()+1;
				if (bcount*10f/ flen> perc) {
					System.err.print("*");
					System.err.flush();
					++perc;
				}
				
				StringTokenizer toki= new StringTokenizer(s, inFileDelim, true);

				int chk= 0;
				while(toki.hasMoreElements())
					if (toki.nextElement().equals(inFileDelim))
						++chk;
				if (chk!= 107) {
					System.err.println("\n\t[OHLALA] line "+rowCtr+" with "+chk+"<>107 tokens, skipped.");
					continue;
				}
				
				toki= new StringTokenizer(s, inFileDelim, true);
				s= toki.nextToken();	// transcriptID, must be there
				if (map!= null&& (!map.contains(s)))
					continue;
				s= toki.nextToken();
				s= toki.nextToken();	// geneID, ..or not
				if (s.equals(inFileDelim)) {
					continue;
					//noNameGeneIDCtr++;
				} else {
					s= toki.nextToken();	// delim
				}
				
				int ctr= 2, ctrIdx= 0;	// we are in field index[2] of the table
				s= toki.nextToken();
				while (ctrIdx< useIDs.length) {
					if (s.equals(inFileDelim)) {
						try {
							s= toki.nextToken();
						} catch (NoSuchElementException e) {
							System.err.println(rowCtr);
							e.printStackTrace();
						}
						if (ctr== useIDs[ctrIdx].i)
							++ctrIdx;
						++ctr;
						continue;
					}
					String sTmp= s;
					while (ctr< useIDs[ctrIdx].i) {
						s= toki.nextToken();
						if (s.equals(inFileDelim)) 
							++ctr;
					}
					s= sTmp;
					
					// read the info of the field
					if (s.trim().length()== 0) {	// not annotated
						; // :)
					} else {	// annotated						
						String[] ss= s.split(",\\s+");	// all int ids for this field
						for (int i = 0; i < ss.length; i++) {
							String id= useIDs[ctrIdx].o+attrGrpIDseparator+ss[i];
							int p= Arrays.binarySearch(attrIDs, id);
							if (mapCountBg.get(id)== null) 
								mapCountBg.put(id, 1);
							else
								mapCountBg.put(id, mapCountBg.get(id)+1);
							
							++cntPlus;
						}
						// count up for the group the annotated transcripts
						if (mapCountBg.get(useIDs[ctrIdx].o)== null)
							mapCountBg.put((String) useIDs[ctrIdx].o, 1);
						else
							mapCountBg.put((String) useIDs[ctrIdx].o, 
									mapCountBg.get((String) useIDs[ctrIdx].o)+ 1);

					}
					++ctrIdx;
					s= toki.nextToken();
				}
				// TODO this is a lazy safety check, shouldnt consume so much time
				if (s.equals(inFileDelim))
					++ctr;
				while(toki.hasMoreElements()) {
					s= toki.nextToken();
					if (s.equals(inFileDelim))
						++ctr;
				}
				assert(ctr==107);	// field 108= index 107
				++redRowCtr;
			}
			System.err.println();
			long tabSize= ((long) rowIDs.length)* attrIDs.length;
			System.err.println("[YES!] redTable ("+rowIDs.length+" x "+attrIDs.length+", "
					+MyFile.humanReadableSize(tabSize)+") inited, "
					+cntPlus+" pos, "+(tabSize-cntPlus-cntMinus)+" neg, "+cntMinus+" not annotated, "
					+((System.currentTimeMillis()- t0)/1000)+" sec.");


			// write out table
			BufferedWriter writer= new BufferedWriter(new FileWriter(f));
			Iterator<String> iter= mapCountBg.keySet().iterator();
			while (iter.hasNext()) {
				s= iter.next();
				writer.write(s);
				writer.write("\t");
				s= Integer.toString(mapCountBg.get(s));
				writer.write(s);
				writer.write("\n");
			}
			writer.flush();
			writer.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	
	
	private void readAttributes() {
		try { 
			System.err.println("["+Constants.getTime()+"] Reading annotation ");
			BufferedReader buffy= new BufferedReader(new FileReader(annFile));
			String s= buffy.readLine();	// header
			IntegerCompound[] useIDs= new IntegerCompound[attrGrpIDs.size()];
			String inFileDelim= "\t";
			StringTokenizer toki= new StringTokenizer(s, inFileDelim);
			int ctr= 0, ctrIdx= 0;
			while(toki.hasMoreElements()) {
				String f= toki.nextToken();
				if (attrGrpIDs.get(f)!= null) {
					useIDs[ctrIdx]= new IntegerCompound();
					useIDs[ctrIdx].o= f;
					useIDs[ctrIdx++].i= ctr;
					if (ctrIdx== useIDs.length)
						break;
				}
				++ctr;
			}
			System.err.println(" "+attrGrpIDs.size()+" attribute categories.\n");
			assert(ctrIdx== useIDs.length);
			Arrays.sort(useIDs);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	// read foreground first
	void getInput(boolean background) {
		try {
			assert(attrGrpIDs!= null);
			assert(rowIDkey!= null);
			
			Vector<HashMap<String,String>> vecAttr= null;
			if (background) {
				mapCountBg= new HashMap<String, Integer>();	// maps attributeGrpID/attribute -> count
			} else {
				vecAttr= new Vector<HashMap<String,String>>();
				assert(mapAttr2Grp== null);
				mapAttr2Grp= new HashMap<String, String>();	// maps attribute to group ID
			}
			
			File f= (background)?bgFile:inFile;
			GTFwrapper reader= new GTFwrapper(f.getAbsolutePath());
			reader.setReadGene(false);
			reader.setReadGTF(true);
			reader.setReadFeatures(new String[]{"as_event", "ds_event"});
			reader.read();
			for (GFFObject[] obj= reader.getGtfObj(); obj!= null; reader.read()) {
				for (int i = 0; i < obj.length; i++) {
					HashMap<String,String> map= obj[i].getAttributes();
					if (!background)
						vecAttr.add(map);
					Iterator<String> iter= map.keySet().iterator();
					while(iter.hasNext()) {
						String key= iter.next();
						if (attrGrpIDs.get(key)== null)
							continue;
						String val= map.get(key);
						if (background) {
							if (mapCountBg.containsKey(key))
								mapCountBg.put(key, mapCountBg.remove(key)+1);
							else
								mapCountBg.put(key, 1);
							if (mapCountBg.containsKey(val))
								mapCountBg.put(val, mapCountBg.remove(val)+1);
							else
								mapCountBg.put(val, 1);
							if(!mapAttr2Grp.containsKey(val))
								mapAttr2Grp.put(val, key);
						} else {
							if(!mapAttr2Grp.containsKey(val))
								mapAttr2Grp.put(val, key);
						}
					}
				}
			}
			
			if (!background) {
				redTable= new byte[vecAttr.size()][]; // (non?)redundant, 1= annotated, 0= not annotated, -1= not at all annotated in category
				rowIDs= new String[vecAttr.size()];
				attrIDs= new String[mapAttr2Grp.size()];
				Hashtable<String, Integer> tmpAttrIdx= new Hashtable<String,Integer>(attrIDs.length,1f);
				int ctr= 0;
				Iterator<String> iter= mapAttr2Grp.keySet().iterator();
				while(iter.hasNext()) {
					attrIDs[ctr]= iter.next();
					tmpAttrIdx.put(attrIDs[ctr], ctr);
					++ctr;
				}
				
				for (int i = 0; i < vecAttr.size(); i++) {
					redTable[i]= new byte[mapAttr2Grp.size()];
					for (int j = 0; j < redTable[i].length; j++) {
						String val= vecAttr.elementAt(i).get(mapAttr2Grp.get(attrIDs[j]));
						if (val== null)
							redTable[i][j]= -1;
						else if (val.equals(attrIDs[j]))
							redTable[i][j]= 1;
						else
							redTable[i][j]= 0;
					}
					
					rowIDs[i]= vecAttr.elementAt(i).get(rowIDkey);
				}
				
			}
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	boolean annotGenes= true;
	@Override
	public void run() {
		
		readEvents();		
		
		preFilter(redTable, mapCountBg);
		
		float[][] kappaTbl= getKappa(redTable);
		
		PreCluster[] preClusters= cluster(kappaTbl);
		Cluster[] clusters= createClusters(preClusters);
		
		clusters= getReplacementTerms(clusters);
		if (annotGenes)
			getGeneAnnotation(clusters);
		
		Arrays.sort(clusters, Cluster.defaultClusterByEScoreComparator);
		try {
			PrintStream p = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFile)));
			for (int i = clusters.length-1; i >= 0; --i) {
				p.println(clusters[i].toString());
			}		
			p.flush();
			p.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}
	
	private void getGeneAnnotation(Cluster[] clusters) {
		// get all events
		System.err.print("["+new Time(System.currentTimeMillis())+"] getting gene ids ");
		HashMap<String,Object> map= new HashMap<String,Object>();
		for (int i = 0; i < clusters.length; i++) {
			for (int j = 0; j < clusters[i].memberIDs.length; j++) {
				map.put(clusters[i].memberIDs[j], null);
			}
		}

		
		// get all genes of events
		try {
			GTFwrapper reader= new GTFwrapper(annFile.getAbsolutePath());
			reader.setReadGene(false);
			reader.setReadGTF(true);
			reader.setLimitGTFObs(100);
			reader.setSilent(false);
			reader.setReadFeatures(new String[] {"as_event", "ds_event"});	// GTFObject.FEATURE_ASEVENT
			GFFObject[] obj;
			for (reader.read(); (obj= reader.getGtfObj())!= null; reader.read()) {
				for (int i = 0; i < obj.length; i++) {
					String idKey= obj[i].getChromosome()+"_"+obj[i].getStart()+"_"+obj[i].getEnd()
								+"_"+obj[i].getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_SPLICE_CHAIN);
					if(!map.containsKey(idKey))
						continue;
					
					int[][] geneAnn= new int[3][];
					String s= obj[i].getAttribute("GENE_SYMBOL");
					if (s!= null) {
						String[] ss= s.split("[/,]");
						geneAnn[0]= new int[ss.length];
						for (int j = 0; j < ss.length; j++) 
							geneAnn[0][j]= Integer.parseInt(ss[j]);						
					}
					
					s= obj[i].getAttribute("Gene_Name");
					if (s!= null) {
						String[] ss= s.split("[/,]");
						geneAnn[1]= new int[ss.length];
						for (int j = 0; j < ss.length; j++) 
							geneAnn[1][j]= Integer.parseInt(ss[j]);						
					}
					
					s= obj[i].getAttribute("GENE_NAME");
					if (s!= null) {
						String[] ss= s.split("[/,]");
						geneAnn[2]= new int[ss.length];
						for (int j = 0; j < ss.length; j++) 
							geneAnn[2][j]= Integer.parseInt(ss[j]);						
					}

					map.put(idKey,geneAnn);
				}
			}
			System.err.println();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.err.println(" "+map.size()+" ids annotated with genes.\n");
		
		// make non-redundant
		HashSet<Integer> set0= new HashSet<Integer>(), set1= new HashSet<Integer>(), set2= new HashSet<Integer>();
		Iterator iter= map.values().iterator();
		while (iter.hasNext()) {	
			int[][] ano= (int[][]) iter.next();
			if (ano== null)
				continue;
			for (int i = 0; ano[0]!= null&& i < ano[0].length; i++) 
				set0.add(ano[0][i]);
			for (int i = 0; ano[1]!= null&& i < ano[1].length; i++) 
				set1.add(ano[1][i]);
			for (int i = 0; ano[2]!= null&& i < ano[2].length; i++) 
				set2.add(ano[2][i]);
		}
		Integer[] list0= new Integer[set0.size()], list1= new Integer[set1.size()], list2= new Integer[set2.size()];
		list0= set0.toArray(list0); list1= set1.toArray(list1); list2= set2.toArray(list2);
		set0= null; set1= null; set2= null;
		Arrays.sort(list0); Arrays.sort(list1); Arrays.sort(list2);
		String[] terms0= new String[list0.length], terms1= new String[list1.length], terms2= new String[list2.length];

		// get terms of genes
		try {
			File f= new File("M:\\work\\claudia\\david_genbank\\save\\39_GENE_SYMBOL.tbl");
			BufferedReader buffy= new BufferedReader(new FileReader(f));
			int ctr= 0;
			System.err.print("["+Constants.getTime()+"] GENE_SYMBOL ");
			System.err.flush();
			long fsize= f.length(), bread= 0;
			int perc= 0;
			for (String s= buffy.readLine();list0.length> 0&& s!= null;s= buffy.readLine()) {
				bread+= s.length()+1;
				if (bread*10f/fsize> perc) {
					System.err.print("*");
					System.err.flush();
					++perc;
				}
				//String[] ss= s.split("\\s+");
				int p= s.indexOf('\t');
				Integer id= new Integer(s.substring(0,p));	// ss[0]
				if (list0[ctr].equals(id)) {
					terms0[ctr]= s.substring(p+ 1); 	// ss[1];
					if (++ctr== list0.length)
						break;
				}
			}
			buffy.close();
			System.err.println(" "+ctr+" gene symbols.\n");
			
			f= new File("M:\\work\\claudia\\david_genbank\\save\\2_Gene_Name.tbl");
			buffy= new BufferedReader(new FileReader(f));
			ctr= 0; perc= 0;
			System.err.print("["+Constants.getTime()+"]Gene_Name ");
			System.err.flush();
			fsize= f.length(); bread= 0;
			for (String s= buffy.readLine();list1.length> 0&& s!= null;s= buffy.readLine()) {
				bread+= s.length()+1;
				if (bread*10f/fsize> perc) {
					System.err.print("*");
					System.err.flush();
					++perc;
				}
				//String[] ss= s.split("\\s+");
				int p= s.indexOf('\t');
				Integer id= new Integer(s.substring(0,p));	// ss[0]
				if (list1[ctr].equals(id)) {
					terms1[ctr]= s.substring(p+ 1); 	// ss[1];
					if (++ctr== list1.length)
						break;
				}
			}
			buffy.close();
			System.err.println(" "+ctr+" gene names.");

			f= new File("M:\\work\\claudia\\david_genbank\\save\\40_GENE_NAME.tbl");
			buffy= new BufferedReader(new FileReader(f));
			ctr= 0; perc= 0;
			System.err.print("["+Constants.getTime()+"]GENE_NAME ");
			System.err.flush();
			fsize= f.length(); bread= 0;
			for (String s= buffy.readLine();list2.length> 0&& s!= null;s= buffy.readLine()) {
				bread+= s.length()+1;
				if (bread*10f/fsize> perc) {
					System.err.print("*");
					System.err.flush();
					++perc;
				}
				// String[] ss= s.split("\\s+");
				int p= s.indexOf('\t');
				Integer id= new Integer(s.substring(0,p));	// ss[0]
				if (list2[ctr].equals(id)) {
					terms2[ctr]= s.substring(p+ 1); 	// ss[1];
					if (++ctr== list2.length)
						break;
				}
			}
			System.err.println(" "+ctr+" gene names.\n");
			buffy.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
	
		// put in event-map
		Object[] keys= map.keySet().toArray();
		for (int x = 0; x < keys.length; x++) {
			int[][] ano= (int[][]) map.remove(keys[x]);
			if (ano== null)
				continue;
			//String[] concat= new String[3];
			StringBuffer sb= new StringBuffer();
			for (int i = 0; ano[0]!= null&& i < ano[0].length; i++) {
				int p= Arrays.binarySearch(list0, ano[0][i]);
//				if (concat[0]== null)
//					concat[0]= terms0[p];
//				else
//					concat[0]+= " "+terms0[p];
				sb.append(" ");
				sb.append(terms0[p]);
			}
			for (int i = 0; ano[1]!= null&& i < ano[1].length; i++) {
				int p= Arrays.binarySearch(list1, ano[1][i]);
//				if (concat[1]== null)
//					concat[1]= terms1[p];
//				else
//					concat[1]+= " "+terms1[p];
				sb.append(" ");
				sb.append(terms1[p]);
			}
			for (int i = 0; ano[2]!= null&& i < ano[2].length; i++) {
				int p= Arrays.binarySearch(list2, ano[2][i]);
//				if (concat[2]== null)
//					concat[2]= terms2[p];
//				else
//					concat[2]+= " "+terms2[p];
				sb.append(" ");
				sb.append(terms2[p]);
			}
			map.put((String) keys[x],sb.toString().trim());
		}
		
		// lazily annotate clusters
		for (int i = 0; i < clusters.length; i++) {
			HashSet<String> clusterGenes= new HashSet<String>();
			for (int j = 0; j < clusters[i].memberIDs.length; j++) {
				clusterGenes.add((String) map.get(clusters[i].memberIDs[j]));
			}
			Object[] o= clusterGenes.toArray();
			StringBuffer sb= new StringBuffer();
			for (int j = 0; j < o.length; j++) {
				sb.append(o[j]);
				if (j< o.length-1)
					sb.append("\n");
			}

			clusters[i].geneTerms= sb.toString();
		}
	}

	void readEvents() {
		try {
			RedTableWrapper redStyle= new RedTableWrapper(DEFAULT_REDTABLE_PATH);
			int ok= redStyle.readInTable();
			if (ok== 0) {	// TODO read in backtable
				boolean ok2= readBacktable();
				if (ok2)
					return;
			}
			
			System.err.println("["+Constants.getTime()+"] Reading foreground IDs ");
			BufferedReader buffy= new BufferedReader(new FileReader(inFile));
			String s;
			refTIDMap= new HashSet<String>();	// reuse later for Int GeneIDs
			long t0= System.currentTimeMillis();
			while ((s= buffy.readLine())!= null) {
				refTIDMap.add(s.trim());
			}
			System.err.println(" "+refTIDMap.size()+" reference IDs, "
					+((System.currentTimeMillis()- t0)/ 1000)+" sec.\n");

			System.err.print("["+Constants.getTime()+"] reading foreground/background..");
			System.err.flush();
			HashMap<String,HashSet<Integer>> idsInUse= new HashMap<String,HashSet<Integer>>();
			Iterator<String> iter= attrGrpIDs.keySet().iterator();
			while(iter.hasNext()) {
				String grpID= iter.next();
				if (grpID.indexOf("domain")>= 0) {	// reference/variant specific
					idsInUse.put(grpID+"_R", new HashSet<Integer>());
					idsInUse.put(grpID+"_V", new HashSet<Integer>());
				} else
					idsInUse.put(grpID, new HashSet<Integer>());
			}
				
			// get foreground IDs
			Vector<String> rowIDsVec= new Vector<String>();
			mapCountBg= new HashMap<String, Integer>();
			
			readEventsAnnotation(true, rowIDsVec, idsInUse);
			
			rowIDs= new String[rowIDsVec.size()];
			rowIDs= rowIDsVec.toArray(rowIDs);
			rowIDsVec= null;
			Iterator<HashSet<Integer>> iter2= idsInUse.values().iterator();
			int sum= 0;
			while (iter2.hasNext()) {
				sum+= iter2.next().size();
			}			
			System.err.println("[INFO] "+sum+" different annotation items in foreground.");
			attrIDs= new String[sum];
			int ctr= 0;
			Iterator<String> iter4= idsInUse.keySet().iterator();
			while (iter4.hasNext()) {
				String gID= iter4.next();
				Iterator<Integer> iter3= idsInUse.get(gID).iterator();
				while (iter3.hasNext()) {
					String val= gID+attrGrpIDseparator+iter3.next();
					attrIDs[ctr++]= val;
				}
			}
			Arrays.sort(attrIDs);
			writeBacktable();
			System.gc();
				
			
			System.err.println("[INFO] Building red table "+rowIDs.length+" x "+sum);
			redTable= new byte[rowIDs.length][];
			for (int i = 0; i < rowIDs.length; i++) {
					redTable[i]= new byte[sum];
				for (int j = 0; j < redTable[i].length; j++) 
					redTable[i][j]= 0;	// .. not annotated
			}

			readEventsAnnotation(false, rowIDsVec, idsInUse);
			
			RedTableWrapper wrapper= new RedTableWrapper(DEFAULT_REDTABLE_PATH);
			wrapper.writeOutTable();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void readEventsAnnotation(boolean first, Vector<String> rowIDsVec, HashMap<String,HashSet<Integer>> idsInUse) {
		
		try {
			GTFwrapper reader= new GTFwrapper(annFile.getAbsolutePath());
			reader.setReadGene(false);
			reader.setReadGTF(true);
			reader.setLimitGTFObs(100);
			reader.setSilent(false);
			reader.setReadFeatures(new String[] {"as_event", "ds_event"});	// GTFObject.FEATURE_ASEVENT
			GFFObject[] obj;
			int rowCtr= 0, debug= 0;
			for (reader.read(); (obj= reader.getGtfObj())!= null; reader.read()) {
				++debug;
				for (int i = 0; i < obj.length; i++) {
					boolean foreGround= false;
					if (isForeground(obj[i])) {
						foreGround= true;
						String idKey= obj[i].getChromosome()+"_"+obj[i].getStart()+"_"+obj[i].getEnd()
							+"_"+obj[i].getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_SPLICE_CHAIN);
						if (first)
							rowIDsVec.add(idKey);
					} else {
						if (!first)
							continue;
					}
	
					
					Iterator<String> iter= attrGrpIDs.keySet().iterator();
					while(iter.hasNext()) {
						String grpID= iter.next();
						String value= obj[i].getAttribute(grpID);
						boolean noContent= true;
						if (value!= null) {
							if (grpID.indexOf("domain")>= 0) {	// reference/variant specific
								int refVar= Integer.parseInt(obj[i].getAttribute("ref"))-1;
								String[] valAB= value.split(",");
								for (int j = 0; j < valAB.length; j++) {
									String splitGID= null;
									if (j== refVar)
										splitGID= grpID+"_R";
									else
										splitGID= grpID+"_V";
									String[] comp= valAB[j].split("/");
									for (int k = 0; k < comp.length; k++) {	
										if(comp[k].length()== 0)
											continue;
										noContent= false;
										String compID= splitGID+attrGrpIDseparator+comp[k];
										if (first) {
											if (mapCountBg.get(compID)== null)
												mapCountBg.put(compID, 1);
											else 
												mapCountBg.put(compID, mapCountBg.get(compID)+1);
											if (mapCountBg.get(splitGID)== null)
												mapCountBg.put(splitGID, 1);
											else 
												mapCountBg.put(splitGID, mapCountBg.get(splitGID)+1);
											if (isForeground(obj[i]))
												idsInUse.get(splitGID).add(Integer.parseInt(comp[k]));
										} else {
											int p= Arrays.binarySearch(attrIDs, compID);
											redTable[rowCtr][p]= 1;
										}
									}
								}
							} else {
								String[] comp= value.split("[/,]");
								for (int k = 0; k < comp.length; k++) { 
									if(comp[k].length()== 0)
										continue;
									noContent= false;
									String compID= grpID+attrGrpIDseparator+comp[k];
									if (first) {
										if (mapCountBg.get(compID)== null)
											mapCountBg.put(compID, 1);
										else 
											mapCountBg.put(compID, mapCountBg.get(compID)+1);
										if (mapCountBg.get(grpID)== null)
											mapCountBg.put(grpID, 1);
										else 
											mapCountBg.put(grpID, mapCountBg.get(grpID)+1);
										if (isForeground(obj[i]))
											idsInUse.get(grpID).add(Integer.parseInt(comp[k]));
									} else {
										int p= Arrays.binarySearch(attrIDs, compID);
										redTable[rowCtr][p]= 1;										
									}
								}
							}
							
						} 
						
						if (noContent&& (!first)) {	// null or no value, unannotated
							String[] allGIDs= null;
							if (grpID.indexOf("domain")>= 0) 	// reference/variant specific
								allGIDs= new String[] {grpID+"_R", grpID+"_V"}; 
							else
								allGIDs= new String[] {grpID};
							for (int k = 0; k < allGIDs.length; k++) {
								Iterator<Integer> iterX= idsInUse.get(allGIDs[k]).iterator();
								while (iterX.hasNext()) {
									String key= allGIDs[k]+attrGrpIDseparator+iterX.next();
									int p= Arrays.binarySearch(attrIDs, key);
									redTable[rowCtr][p]= 2;
								}
							}
						}
						
					} // iterate all groups

					if (!first)
						++rowCtr;
					
				} // all objects
				
				//break;	// *** DEBUG ***

			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private void writeBacktable() {
		try {
			BufferedWriter writer= new BufferedWriter(new FileWriter(DEFAULT_BACKTABLE_PATH));
			Iterator<String> iter= mapCountBg.keySet().iterator();
			String s= null;
			while (iter.hasNext()) {
				s= iter.next();
				writer.write(s);
				writer.write("\t");
				s= Integer.toString(mapCountBg.get(s));
				writer.write(s);
				writer.write("\n");
			}
			writer.flush();
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private boolean readBacktable() {

		File f= new File(DEFAULT_BACKTABLE_PATH);
		if (f.exists()) {
			System.err.print("["+Constants.getTime()+"] Reading background data ");
			System.err.flush();
			try {
				BufferedReader buffy= new BufferedReader(new FileReader(f));
				mapCountBg= new HashMap<String, Integer>();
				String s;
				int perc= 0;
				long len= f.length(), bread= 0;
				while((s= buffy.readLine())!= null) {
					bread+= s.length()+1;
					if ((bread*10d/ len)> perc) {
						System.err.print("*");
						System.err.flush();
						++perc;
					}
					String[] ss= s.split("\t");		// id x count
					assert(ss.length== 2);
					int x= -1;
					try {
						x= Integer.parseInt(ss[1]);
					} catch (Exception e) {
						; // :)
					}
					// TODO check consistency here with foreground..
					mapCountBg.put(ss[0],x);
				}
				System.err.println(" "+mapCountBg.size()+" terms.\n");
				buffy.close();
				return true;
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return false;
	}

	HashSet<String> refTIDMap= null;
	boolean isForeground(GFFObject o) {
		String[] s= o.getTranscriptID().split("[,/]");
		for (int i = 0; i < s.length; i++) {
			int p= s[i].indexOf("_dup");
			if (p>= 0)
				s[i]= s[i].substring(0,p);
			if (refTIDMap.contains(s[i]))
				return true;
		}
		return false;
	}

	Cluster[] getReplacementTerms(Cluster[] clusters) {
		TermReplacement[] nullReplaceArray= new TermReplacement[0];
		
		HashMap<String,TermReplacement[]> globalMap= new HashMap<String,TermReplacement[]>();
		System.err.print("["+Constants.getTime()+"] getting replacement terms ");
		System.err.flush();
		int perc= 0;
		for (int i = 0; i < clusters.length; i++) {
			if (i* 10d/ clusters.length> perc) {
				System.err.print("*");
				System.err.flush();
				++perc;
			}
			Iterator<String> catIter= clusters[i].getMapReplace().keySet().iterator();
			TermReplacement[] savTerms= null;
			while (catIter.hasNext()) {
				String catID= catIter.next();
				Vector<Integer> v= clusters[i].getMapReplace().get(catID);
				TermReplacement[] terms= globalMap.get(catID);
				if (terms== null) {
					terms= new TermReplacement[v.size()];
					savTerms= nullReplaceArray;
				} else { 
					savTerms= terms;
					terms= new TermReplacement[terms.length+ v.size()];
					System.arraycopy(savTerms, 0, terms, 0, savTerms.length);
				}
				
				globalMap.put(catID, terms);
				int ctr= 0;
				for (int j = 0; j < v.size(); j++) {
					TermReplacement dummy= new TermReplacement(v.elementAt(j));
					int p= -1;
					if (savTerms!= null)
						p= Arrays.binarySearch(savTerms, dummy, TermReplacement.defaultTermReplacementByIDComparator);
					if (p< 0) // non-redundant
						terms[terms.length- v.size()+ ctr++]= dummy;
				}
				int offset= v.size()- ctr;
				if (offset!= 0) {
					savTerms= terms;
					terms= new TermReplacement[savTerms.length- offset];
					System.arraycopy(savTerms, 0, terms, 0, terms.length);
					globalMap.put(catID, terms);
				}
				
				Arrays.sort(terms, TermReplacement.defaultTermReplacementByIDComparator);
			}
		}
		
		Iterator<TermReplacement[]> iterTR= globalMap.values().iterator();
		while (iterTR.hasNext()) {
			TermReplacement[] tr= iterTR.next();
			Arrays.sort(tr, TermReplacement.defaultTermReplacementByIDComparator);
		}
		
		globalMap= getReplacementTerms(globalMap);
		
		for (int i = 0; i < clusters.length; i++) {
			clusters[i].replace(globalMap);
		}
		
		System.err.println(" "+globalMap.size()+" terms replaced.\n");
		
		return clusters;
	}

	public static String kBaseDir= "M:\\work\\claudia\\david_genbank\\save", SFX_kBase=".tbl";
	HashMap<String, TermReplacement[]> getReplacementTerms(
			HashMap<String, TermReplacement[]> globalMap) {
		
		try {
			String[] fNames= new File(kBaseDir).list();
			Iterator<String> catIter= globalMap.keySet().iterator();
			while(catIter.hasNext()) {
				String catID= catIter.next();
				String fileSfx= catID+SFX_kBase;
				if (catID.contains("domain"))
					fileSfx= "PFAM_DOMAINS_sorted.tbl";
				int x;
				for (x = 0; x < fNames.length; x++) 
					if (fNames[x].endsWith(fileSfx))
						break;
				BufferedReader buffy= new BufferedReader(new FileReader(kBaseDir+File.separator+fNames[x]));
				TermReplacement[] tReps= globalMap.get(catID);
				int ctr= 0, lineCtr= 0;
				String s;
				while (ctr< tReps.length&& (s= buffy.readLine())!= null) {
					++lineCtr;
					if (lineCtr== tReps[ctr].termID) {
						String[] ss= s.split("\t");
						if(ss.length!= 2) {
							System.err.println("\n\tBAD line "+s);
							System.err.println("\tin file "+kBaseDir+File.separator+fNames[x]);
							continue;
						}
						tReps[ctr].term= ss[1];
						++ctr;
					}
				}
				buffy.close();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return globalMap;
	}

	public String getRowIDkey() {
		return rowIDkey;
	}

	public void setRowIDkey(String rowIDkey) {
		this.rowIDkey = rowIDkey;
	}

	public HashMap<String, Float> getAttrGrpIDs() {
		return attrGrpIDs;
	}

	public void setAttrGrpIDs(HashMap<String, Float> attrGrpIDs) {
		this.attrGrpIDs = attrGrpIDs;
	}

	void readForeground() {
			
			try {			
				
				RedTableWrapper redStyle= new RedTableWrapper(DEFAULT_REDTABLE_PATH);
				int ok= redStyle.readInTable();
				if (ok== 0) 
					return;
				
				System.err.println("["+Constants.getTime()+"] Reading foreground IDs.. ");
				BufferedReader buffy= new BufferedReader(new FileReader(inFile));
				String s;
				HashSet map= new HashSet();	// reuse later for Int GeneIDs
				long t0= System.currentTimeMillis();
				while ((s= buffy.readLine())!= null) {
					map.add(s.trim());
				}
				System.err.println(" "+map.size()+" unique transcript IDs, "
						+((System.currentTimeMillis()- t0)/ 1000)+" sec.\n");
	
				System.err.print("[HEHO, LETS GO] prescanning..");
				int sum= readAnnotation(map, true, false);
								
				// INIT the things
				//redFile= new File("redTable");	// otherwise, 3GB memory.. background, foreground only 3MB			
				redTable= new byte[idLists[idLists.length-1].size()][];
				rowIDs= new String[redTable.length];
				Iterator iter= idLists[idLists.length-1].iterator();
				for (int i = 0; i < rowIDs.length; i++) {
					rowIDs[i]= (String) iter.next();
					redTable[i]= new byte[sum];
					for (int j = 0; j < redTable[i].length; j++) 
						redTable[i][j]= 0;	// .. not annotated
				}
				int ctr= 0;
				attrIDs= new String[sum];
				for (int i = 0; i < idLists.length-1; i++) {
					iter= idLists[i].iterator();
					while (iter.hasNext())
						attrIDs[ctr++]= useIDs[i].o+attrGrpIDseparator+iter.next();
				}
				Arrays.sort(attrIDs);
				
				// and now init the redTable..
				System.err.print("[ANDNOW] initing redtable.. ");
				readAnnotation(map, false, false);
				//assert(sum2== sum);

				redStyle.writeOutTable();
				
			} catch (Exception e) {
				e.printStackTrace();
			}

		}

	
	private void printUserDataStats() {
		// foreground stats
		Iterator<String> iter= attrGrpIDs.keySet().iterator();
		System.err.println("[INFO] annotations per category found in user data:");
		System.err.println("group\tterms\tpositives");
		while (iter.hasNext()) {
			String groupID= iter.next();
			int cntAttr= 0, cntPos= 0;
			for (int i = 0; i < attrIDs.length; i++) {
				if (attrIDs[i].startsWith(groupID)) {
					++cntAttr;	// attributes with that group;
					for (int j = 0; j < redTable.length; j++) 
						if (redTable[j][i]==1)
							++cntPos;
				}
			}
			
			System.err.println(groupID+"\t"+cntAttr+"\t"+cntPos);
		}
		System.err.println();
	}
	
	/**
	 * @deprecated uses string tokenizer
	 * @param map
	 * @param prescan
	 * @param background
	 * @return
	 */
	int readAnnotation_old(HashSet map, boolean prescan, boolean background) {
	
		try {
			if (background)
				mapCountBg= new HashMap<String, Integer>();
			
			BufferedReader buffy= new BufferedReader(new FileReader(annFile));
			String s= buffy.readLine();	// header

			if (useIDs== null)
				useIDs= readAnnotationHeader(s);
			
			if (!background) {
				idLists= new HashSet[useIDs.length+1];	// last is gene
				for (int i = 0; i < idLists.length-1; i++) 
					idLists[i]= new HashSet<Integer>();
				idLists[idLists.length-1]= new HashSet<String>();	// genes are not consistent, for instance gene "130, " in row 15084 contains some different annotation than in a line before
			}
			
			int rowCtr= 1, perc= 0, noNameGeneIDCtr= 1;;
			int cntPlus= 0, cntMinus= 0, redRowCtr= 0;
			long t0= System.currentTimeMillis(), bcount= 0, flen= annFile.length();
			while((s= buffy.readLine())!= null) {
				++rowCtr;
				bcount+= s.length()+1;
				if (bcount*10f/ flen> perc) {
					System.err.print("*");
					System.err.flush();
					++perc;
				}
				
				StringTokenizer toki= new StringTokenizer(s, inFileDelim, true);
	
				int chk= 0;
				while(toki.hasMoreElements())
					if (toki.nextElement().equals(inFileDelim))
						++chk;
				if (chk!= 107) {
					System.err.println("\n\t[OHLALA] line "+rowCtr+" with "+chk+"<>107 tokens, skipped.");
					continue;
				}
				
				toki= new StringTokenizer(s, inFileDelim, true);
				s= toki.nextToken();	// transcriptID, must be there
				if ((map!= null) && (!map.contains(s)))	// map only null when background/all
					continue;
				if (prescan)
					idLists[idLists.length-1].add(s);	// annotate transcripts
				s= toki.nextToken();
				s= toki.nextToken();	// geneID, ..or not
				if (s.equals(inFileDelim)) {
					continue;
					//noNameGeneIDCtr++;
				} else {
					s= toki.nextToken();	// delim
				}
				
				int ctr= 2;
				int ctrIdx= 0;	// we are in field index[2] of the table
				s= toki.nextToken();
				while (ctrIdx< useIDs.length) {
					if (s.equals(inFileDelim)) {
						try {
							s= toki.nextToken();
						} catch (NoSuchElementException e) {
							System.err.println(rowCtr);
							e.printStackTrace();
						}
						if (ctr== useIDs[ctrIdx].i)
							++ctrIdx;
						++ctr;
						continue;
					}
					//String sTmp= s;
					while (ctr< useIDs[ctrIdx].i) {
						s= toki.nextToken();
						if (s.equals(inFileDelim)) 
							++ctr;
					}
					//s= sTmp;
					s= toki.nextToken();	// this one should contain the value of the field if it is not empty
					
					// read the info of the field
					if (s.trim().length()== 0) {	// not annotated
						if ((!background)&& (!prescan)) {
							Iterator iter= idLists[ctrIdx].iterator();
							while (iter.hasNext()) {
								String id= useIDs[ctrIdx].o+attrGrpIDseparator+iter.next();
								int p= Arrays.binarySearch(attrIDs, id);
								redTable[redRowCtr][p]= 2;
								++cntMinus;
							}
						}
					} else {	// >0, annotated						
						String[] ss= s.split(",\\s+");
						for (int i = 0; i < ss.length; i++) {
							if (prescan&& (!background)) {
								Integer id= null;
								try {
									id= Integer.parseInt(ss[i]);
								} catch (Exception e) {
									System.err.println(rowCtr);
									e.printStackTrace();
								}
								idLists[ctrIdx].add(id);
							} else if (!background) {
								String id= useIDs[ctrIdx].o+attrGrpIDseparator+ss[i];
								int p= Arrays.binarySearch(attrIDs, id);
								redTable[redRowCtr][p]= 1;
							} else { // background
								String id= useIDs[ctrIdx].o+attrGrpIDseparator+ss[i];
								if (mapCountBg.get(id)== null) 
									mapCountBg.put(id, 1);
								else
									mapCountBg.put(id, mapCountBg.get(id)+1);
							}
							++cntPlus;
						}
						if (background) {	// count up for the group the annotated transcripts
							if (mapCountBg.get(useIDs[ctrIdx].o)== null)
								mapCountBg.put((String) useIDs[ctrIdx].o, 1);
							else
								mapCountBg.put((String) useIDs[ctrIdx].o, 
										mapCountBg.get((String) useIDs[ctrIdx].o)+ 1);
						}
						s= toki.nextToken();	// next field
					}
					++ctrIdx;					
				}
				// TODO this is a lazy safety check, shouldnt consume so much time
				if (s.equals(inFileDelim))
					++ctr;
				while(toki.hasMoreElements()) {
					s= toki.nextToken();
					if (s.equals(inFileDelim))
						++ctr;
				}
				assert(ctr==107);	// field 108= index 107
				++redRowCtr;
			}
			System.err.println();
			
			if (prescan) {
				int sum= 0;
				for (int i = 0; i < idLists.length-1; i++) 
					sum+= idLists[i].size();
				System.err.println("[PUUH] found "+idLists[idLists.length-1].size()+" genes (" +
						+(noNameGeneIDCtr-1)+") genes without name, "
						+sum+" unique attributes in total, "
						+((System.currentTimeMillis()- t0)/ 1000)+ "sec.");
				return sum;
			} else if (!background) {
				long tabSize= ((long) rowIDs.length)* attrIDs.length;
				System.err.println("[YES!] redTable ("+rowIDs.length+" x "+attrIDs.length+", "
						+MyFile.humanReadableSize(tabSize)+") inited, "
						+cntPlus+" pos, "+(tabSize-cntPlus-cntMinus)+" neg, "+cntMinus+" not annotated, "
						+((System.currentTimeMillis()- t0)/1000)+" sec.");
			}
			

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return -1;
	}

	String inFileDelim= "\t";
	IntegerCompound[] useIDs;

	static TermByPValComparator defaultTermByPValComparator= new TermByPValComparator();
	private IntegerCompound[] readAnnotationHeader(String s) {
		
		System.err.print("[INFO] Checking header.. ");
		useIDs= new IntegerCompound[attrGrpIDs.size()];	// maps group_ids to fields
		StringTokenizer toki= new StringTokenizer(s, inFileDelim);
		int ctr= 0, ctrIdx= 0;
		while(toki.hasMoreElements()) {
			String f= toki.nextToken();
			if (attrGrpIDs.get(f)!= null) {
				useIDs[ctrIdx]= new IntegerCompound();
				useIDs[ctrIdx].o= f;
				useIDs[ctrIdx++].i= ctr;
				if (ctrIdx== useIDs.length)
					break;
			}
			++ctr;
		}
		assert(ctrIdx== useIDs.length);
		Arrays.sort(useIDs);
		System.err.println("found "+attrGrpIDs.size()+" attribute categories to use.");
		
		return useIDs;
	}

	void readForeground_save() {
			
			try {			
				
				RedTableWrapper redStyle= new RedTableWrapper(DEFAULT_REDTABLE_PATH);
				int ok= redStyle.readInTable();
				if (ok== 0) 
					return;
				
				System.err.println("[INFO] Reading foreground IDs.. ");
				BufferedReader buffy= new BufferedReader(new FileReader(inFile));
				String s;
				HashSet map= new HashSet();	// reuse later for Int GeneIDs
				long t0= System.currentTimeMillis();
				while ((s= buffy.readLine())!= null) {
					map.add(s.trim());
				}
				System.err.println("[YIHAA] found "+map.size()+" unique transcript IDs, "
						+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
	
				buffy= new BufferedReader(new FileReader(annFile));
				s= buffy.readLine();	// header
				IntegerCompound[] useIDs= readAnnotationHeader(s);
				
				System.err.print("[HEHO, LETS GO] prescanning..");
				idLists= new HashSet[useIDs.length+1];	// last is gene
				for (int i = 0; i < idLists.length-1; i++) 
					idLists[i]= new HashSet<Integer>();
				idLists[idLists.length-1]= new HashSet<String>();	// genes are not consistent, for instance gene "130, " in row 15084 contains some different annotation than in a line before
				int rowCtr= 1, noNameGeneIDCtr= 1;
				t0= System.currentTimeMillis();
				int perc= 0;
				long bcount= 0, flen= annFile.length();
				while((s= buffy.readLine())!= null) {	
					++rowCtr;
					bcount+= s.length()+1;
					if (bcount*10f/ flen> perc) {
						System.err.print("*");
						System.err.flush();
						++perc;
					}
						
					StringTokenizer toki= new StringTokenizer(s, inFileDelim, true);
	
					int chk= 0;
					while(toki.hasMoreElements())
						if (toki.nextElement().equals(inFileDelim))
							++chk;
	
	//				if (rowCtr== 264218) {
	//					System.err.println(chk);
	//					System.currentTimeMillis();
	//				}
	
					if (chk!= 107) {
						System.err.println("\n\t[OHLALA] line "+rowCtr+" with "+chk+"<>107 tokens, skipped.");
						continue;
					}
					toki= new StringTokenizer(s, inFileDelim, true);
					s= toki.nextToken();	// transcriptID, must be there
					if (!map.contains(s))
						continue;
					idLists[idLists.length-1].add(s);	// annotate transcripts
					s= toki.nextToken();	
					s= toki.nextToken();	// geneID, ..or not
					if (s.equals(inFileDelim)) {
						continue;
	//					s= Integer.toString(noNameGeneIDCtr++);
	//					idLists[idLists.length-1].add(s);
					} else {
	//					if (idLists[idLists.length-1].contains(s))
	//						continue;	// gene already annotated
	//					idLists[idLists.length-1].add(s);
						s= toki.nextToken();	// delim
					}
					
					int ctr= 2, ctrIdx= 0;	// we are in field index[2] of the table
					s= toki.nextToken();
					while (ctrIdx< useIDs.length) {
						if (s.equals(inFileDelim)) {
							try {
								s= toki.nextToken();
							} catch (NoSuchElementException e) {
								System.err.println(rowCtr);
								e.printStackTrace();
							}
							if (ctr== useIDs[ctrIdx].i)
								++ctrIdx;
							++ctr;
							continue;
						}
						//String sTmp= s;
						while (ctr< useIDs[ctrIdx].i) {	// count up to field
							s= toki.nextToken();
							if (s.equals(inFileDelim)) 
								++ctr;
						}
						s= toki.nextToken();	// this one should contain the value of the field if it is not empty
						//s= sTmp;
						
						if (s.trim().length()> 0) {
							String[] ss= s.split(",\\s+");
							for (int i = 0; i < ss.length; i++) {
								Integer id= null;
								try {
									id= Integer.parseInt(ss[i]);
								} catch (Exception e) {
									System.err.println(rowCtr);
									e.printStackTrace();
								}
								idLists[ctrIdx].add(id);
							}
							s= toki.nextToken();	// next field
						}
						++ctrIdx;
					}
					// TODO this is a lazy safety check, shouldnt consume so much time
					if (s.equals(inFileDelim))
						++ctr;
					while(toki.hasMoreElements()) {
						s= toki.nextToken();
						if (s.equals(inFileDelim))
							++ctr;
					}
					assert(ctr==107);	// field 108= index 107
				}
				System.err.println();
				int sum= 0;
				for (int i = 0; i < idLists.length-1; i++) 
					sum+= idLists[i].size();
				System.err.println("[PUUH] found "+idLists[idLists.length-1].size()+" genes (" +
						+(noNameGeneIDCtr-1)+") genes without name, "
						+sum+" unique attributes in total, "
						+((System.currentTimeMillis()- t0)/ 1000)+ "sec.");
				
				
				
				// INIT the things
				//redFile= new File("redTable");	// otherwise, 3GB memory.. background, foreground only 3MB			
				redTable= new byte[idLists[idLists.length-1].size()][];
				rowIDs= new String[redTable.length];
				Iterator iter= idLists[idLists.length-1].iterator();
				for (int i = 0; i < rowIDs.length; i++) {
					rowIDs[i]= (String) iter.next();
					redTable[i]= new byte[sum];
					for (int j = 0; j < redTable[i].length; j++) 
						redTable[i][j]= 0;	// .. not annotated
				}
				int ctr= 0;
				attrIDs= new String[sum];
				for (int i = 0; i < idLists.length-1; i++) {
					iter= idLists[i].iterator();
					while (iter.hasNext())
						attrIDs[ctr++]= useIDs[i].o+attrGrpIDseparator+iter.next();
				}
				Arrays.sort(attrIDs);
				
				// and now init the redTable..
				System.err.print("[ANDNOW] initing redtable.. ");
				buffy= new BufferedReader(new FileReader(annFile));
				s= buffy.readLine();	// header
				rowCtr= 1;
				t0= System.currentTimeMillis();
				perc= 0;
				bcount= 0;
				int cntPlus= 0, cntMinus= 0, redRowCtr= 0;
				while((s= buffy.readLine())!= null) {
					++rowCtr;
					bcount+= s.length()+1;
					if (bcount*10f/ flen> perc) {
						System.err.print("*");
						System.err.flush();
						++perc;
					}
					
					StringTokenizer toki= new StringTokenizer(s, inFileDelim, true);
	
					int chk= 0;
					while(toki.hasMoreElements())
						if (toki.nextElement().equals(inFileDelim))
							++chk;
					if (chk!= 107) {
						System.err.println("\n\t[OHLALA] line "+rowCtr+" with "+chk+"<>107 tokens, skipped.");
						continue;
					}
					
					toki= new StringTokenizer(s, inFileDelim, true);
					s= toki.nextToken();	// transcriptID, must be there
					if (!map.contains(s))
						continue;
					s= toki.nextToken();
					s= toki.nextToken();	// geneID, ..or not
					if (s.equals(inFileDelim)) {
						continue;
						//noNameGeneIDCtr++;
					} else {
						s= toki.nextToken();	// delim
					}
					
					ctr= 2;
					int ctrIdx= 0;	// we are in field index[2] of the table
					s= toki.nextToken();
					while (ctrIdx< useIDs.length) {
						if (s.equals(inFileDelim)) {
							try {
								s= toki.nextToken();
							} catch (NoSuchElementException e) {
								System.err.println(rowCtr);
								e.printStackTrace();
							}
							if (ctr== useIDs[ctrIdx].i)
								++ctrIdx;
							++ctr;
							continue;
						}
						//String sTmp= s;
						while (ctr< useIDs[ctrIdx].i) {
							s= toki.nextToken();
							if (s.equals(inFileDelim)) 
								++ctr;
						}
						//s= sTmp;
						s= toki.nextToken();	// this one should contain the value of the field if it is not empty
						
						// read the info of the field
						if (s.trim().length()== 0) {	// not annotated
							iter= idLists[ctrIdx].iterator();
							while (iter.hasNext()) {
								String id= useIDs[ctrIdx].o+attrGrpIDseparator+iter.next();
								int p= Arrays.binarySearch(attrIDs, id);
								redTable[redRowCtr][p]= 2;
								++cntMinus;
							}
						} else {	// annotated						
							String[] ss= s.split(",\\s+");
							for (int i = 0; i < ss.length; i++) {
								String id= useIDs[ctrIdx].o+attrGrpIDseparator+ss[i];
								int p= Arrays.binarySearch(attrIDs, id);
								redTable[redRowCtr][p]= 1;
								++cntPlus;
							}
						}
						++ctrIdx;
						s= toki.nextToken();
					}
					// TODO this is a lazy safety check, shouldnt consume so much time
					if (s.equals(inFileDelim))
						++ctr;
					while(toki.hasMoreElements()) {
						s= toki.nextToken();
						if (s.equals(inFileDelim))
							++ctr;
					}
					assert(ctr==107);	// field 108= index 107
					++redRowCtr;
				}
				System.err.println();
				long tabSize= ((long) rowIDs.length)* attrIDs.length;
				System.err.println("[YES!] redTable ("+rowIDs.length+" x "+attrIDs.length+", "
						+MyFile.humanReadableSize(tabSize)+") inited, "
						+cntPlus+" pos, "+(tabSize-cntPlus-cntMinus)+" neg, "+cntMinus+" not annotated, "
						+((System.currentTimeMillis()- t0)/1000)+" sec.");
	
				redStyle.writeOutTable();
				
			} catch (Exception e) {
				e.printStackTrace();
			}
	
		}

	int readAnnotation(HashSet map, boolean prescan, boolean background) {
	
		try {
			if (background)
				mapCountBg= new HashMap<String, Integer>();
			
			BufferedReader buffy= new BufferedReader(new FileReader(annFile));
			String s= buffy.readLine();	// header
	
			if (useIDs== null)
				useIDs= readAnnotationHeader(s);
			
			if (!background) {
				idLists= new HashSet[useIDs.length+1];	// last is gene
				for (int i = 0; i < idLists.length-1; i++) 
					idLists[i]= new HashSet<Integer>();
				idLists[idLists.length-1]= new HashSet<String>();	// genes are not consistent, for instance gene "130, " in row 15084 contains some different annotation than in a line before
			}
			
			int rowCtr= 1, perc= 0, noNameGeneIDCtr= 1;;
			int cntPlus= 0, cntMinus= 0, redRowCtr= 0;
			long t0= System.currentTimeMillis(), bcount= 0, flen= annFile.length();
			while((s= buffy.readLine())!= null) {
				++rowCtr;
				bcount+= s.length()+1;
				if (bcount*10f/ flen> perc) {
					System.err.print("*");
					System.err.flush();
					++perc;
				}
				
				//StringTokenizer toki= new StringTokenizer(s, inFileDelim, true);
				String[] tokens= s.split("\t");
				int ep= s.length()-1;
				while (s.substring(ep,ep+1).equals(inFileDelim))
					--ep; 
				if (ep< s.length()-1) {
					String[] tokSave= tokens;
					tokens= new String[tokSave.length+ (s.length()-ep-1)];
					System.arraycopy(tokSave,0,tokens,0,tokSave.length);
					for (int i = tokSave.length; i < tokens.length; i++) {
						tokens[i]= "";
					}
				}
				assert(tokens.length== 108);
				
				if ((map!= null) && (!map.contains(tokens[0])))	// map only null when background/all
					continue;
				if (prescan)
					idLists[idLists.length-1].add(tokens[0]);	// annotate transcripts
				if (tokens[1].length()== 0) {	// geneID, ..or not
					continue;
					//noNameGeneIDCtr++;
				}
				
				for (int i = 0; i < useIDs.length; i++) {
					// read the info of the field
					if (tokens[useIDs[i].i].trim().length()== 0) {	// not annotated
						if ((!background)&& (!prescan)) {
							Iterator iter= idLists[i].iterator();
							while (iter.hasNext()) {
								String id= useIDs[i].o+attrGrpIDseparator+iter.next();
								int p= Arrays.binarySearch(attrIDs, id);
								redTable[redRowCtr][p]= 2;
								++cntMinus;
							}
						}
					} else {	// >0, annotated						
						String[] ss= tokens[useIDs[i].i].split(",\\s+");
						for (int j = 0; j < ss.length; j++) {
							if (prescan&& (!background)) {
								Integer id= null;
								try {
									id= Integer.parseInt(ss[j]);
								} catch (Exception e) {
									System.err.println(rowCtr);
									e.printStackTrace();
								}
								idLists[i].add(id);
							} else if (!background) {
								String id= useIDs[i].o+attrGrpIDseparator+ss[j];
								int p= Arrays.binarySearch(attrIDs, id);
								redTable[redRowCtr][p]= 1;
							} else { // background
								String id= useIDs[i].o+attrGrpIDseparator+ss[j];
								if (mapCountBg.get(id)== null) 
									mapCountBg.put(id, 1);
								else
									mapCountBg.put(id, mapCountBg.get(id)+1);
							}
							++cntPlus;
						}
						if (background) {	// count up for the group the annotated transcripts
							if (mapCountBg.get(useIDs[i].o)== null)
								mapCountBg.put((String) useIDs[i].o, 1);
							else
								mapCountBg.put((String) useIDs[i].o, 
										mapCountBg.get((String) useIDs[i].o)+ 1);
						}
					}
				}
				++redRowCtr;
			}
			System.err.println();
			
			if (prescan) {
				int sum= 0;
				for (int i = 0; i < idLists.length-1; i++) 
					sum+= idLists[i].size();
				System.err.println("[PUUH] found "+idLists[idLists.length-1].size()+" genes (" +
						+(noNameGeneIDCtr-1)+") genes without name, "
						+sum+" unique attributes in total, "
						+((System.currentTimeMillis()- t0)/ 1000)+ "sec.");
				return sum;
			} else if (!background) {
				long tabSize= ((long) rowIDs.length)* attrIDs.length;
				System.err.println("[YES!] redTable ("+rowIDs.length+" x "+attrIDs.length+", "
						+MyFile.humanReadableSize(tabSize)+") inited, "
						+cntPlus+" pos, "+(tabSize-cntPlus-cntMinus)+" neg, "+cntMinus+" not annotated, "
						+((System.currentTimeMillis()- t0)/1000)+" sec.");
			}
			
	
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return -1;
	}
}