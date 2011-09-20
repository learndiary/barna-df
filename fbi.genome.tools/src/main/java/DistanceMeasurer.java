import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Execute;
import fbi.genome.io.BufferedIteratorDisk;

/**
 * Takes an .aln file and computes distances between 
 * mapped reads according to custom criteria.
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class DistanceMeasurer {

	public static void main(String[] args) {
		Execute.initialize(2);
		File input= null; // new File("/Users/micha/projects/demassy/download/R209_K4Me3_IP2_2-2.aln");
		int pile= 1;
		boolean sameStrand= true;
		int uDist= 1000, dDist= 1000;
		if (sameStrand&& uDist!= 0) {
			System.err.println("same strand, set upstream distance to 0");
			uDist= 0;
		}
		
		PrintStream ps= null;
		try {
			ps= new PrintStream(new FileOutputStream(input.getParent()+
					File.separator+ input.getName().substring(0, input.getName().lastIndexOf("."))+
					(sameStrand?"_phase":"_dist")+ "_pile"+ pile+ "-"+ uDist+ "_"+dDist+ ".txt"));
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		

		
		getDistances(input, pile, sameStrand, uDist, dDist, ps);
		
		try {
			ps.flush();
			ps.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		Execute.shutdown();
	}
	
	// pile chr pos F/R
	static long getDistances(File f, int pile, boolean sameStrand, int uDistance, int dDistance, PrintStream ps) {
		long nrPiles= 0;
		BufferedIteratorDisk iter= new BufferedIteratorDisk(f);
		ByteArrayCharSequence cs= null;
		long t0= System.currentTimeMillis();
		int ctr= 0;
		boolean jumped= false;
		while(iter.hasNext()) {
			++ctr;
			if (ctr%100000== 0)
				System.err.println(ctr);
			// check pile
			cs= iter.next();
			if (jumped) {
//				if (!check(cs, ctr))
//					System.currentTimeMillis();
				jumped= false;
			}
			int p= cs.getTokenInt(0);
			if (p< pile)
				continue;
			++nrPiles;
			
			// get attributes
			ByteArrayCharSequence chr= cs.getToken(1).cloneCurrentSeq();
			int pos= -1;
			try {
				pos= cs.getTokenInt(2);
			} catch (Exception e) {
				System.currentTimeMillis();
			}
			boolean forward= false;
			if (cs.getToken(3).equals("F"))
				forward= true;
			
			if ((!sameStrand)&& (!forward)&& (uDistance== 0))
				continue;
			
			iter.mark();
			if (sameStrand) {
				getDistances(iter, pile, chr, pos, forward, dDistance, 1, 0, ps);
			} else {	// opposite strand
				getDistances(iter, pile, chr, pos, !forward, (forward?dDistance:uDistance), (forward?1:-1), (forward?26:-26), ps);
			}
			iter.reset();
			jumped= true;
			
		}
		
		System.err.println("found "+ nrPiles+ " piles, took "+ (System.currentTimeMillis()- t0)/ 1000+ " sec.");
		return nrPiles;
	}

	private static boolean check(ByteArrayCharSequence cs, int ctr) {
		
		
		try {
			BufferedReader buffy= new BufferedReader(new FileReader("/Users/micha/projects/demassy/download/IP5300109chrall_sorted.aln"));
			String s= null;
			for (int i = 0; i < ctr; ++i) 
				s= buffy.readLine();
			buffy.close();
			boolean val= cs.equals(s);
			
			if (!val)
				System.currentTimeMillis();
			return val;
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		
	}

	static private void getDistances(BufferedIteratorDisk iter, int pile, ByteArrayCharSequence chr, int pos,
			boolean forward, int dist, int factor, int summand, PrintStream ps) {
		
		while (iter.hasNext()) {
			ByteArrayCharSequence cs= iter.next();
			// get attributes
			int p= cs.getTokenInt(0);
			if (p< pile)
				continue;
			ByteArrayCharSequence chr2= cs.getToken(1).cloneCurrentSeq();
			if (!chr2.equals(chr))
				break;
			int pos2= cs.getTokenInt(2);
			int diff= (pos2- pos+ summand);
			if (diff> dist)
				break;
			diff*= factor;
			boolean forward2= false;
			if (cs.getToken(3).equals("F"))
				forward2= true;
			if (forward2!= forward)
				continue;
			
			// valid 
			ps.println(Integer.toString(diff));
		}
		
	}
	
	
}
