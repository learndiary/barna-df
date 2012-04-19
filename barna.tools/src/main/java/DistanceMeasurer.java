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

import barna.commons.ByteArrayCharSequence;
import barna.commons.Execute;
import barna.io.BufferedIteratorDisk;

import java.io.*;

/**
 * Takes an .aln file and computes distances between 
 * mapped reads according to custom criteria.
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class DistanceMeasurer {

	public static void main(String[] args) {
		
		Execute.initialize(2);
		
		// config
		// # /Users/micha/projects/demassy/download/R209_K4Me3_IP2_2-2_sorted.aln
		// # IP5300109chrall.aln
		File input= new File("/Users/micha/projects/demassy/download/IP5300109chrall_sorted.aln");
		int pile= 3;		
		boolean sameStrand= true;
		int distance= 1000; 
		
		// uDist, dDist= 1000;		
/*		if (sameStrand&& uDist!= 0) {
			System.err.println("same strand, set upstream distance to 0");
			uDist= 0;
		}
*/		
		
		PrintStream ps= null;
		try {
			ps= new PrintStream(new FileOutputStream(input.getParent()+
					File.separator+ input.getName().substring(0, input.getName().lastIndexOf("."))+
					(sameStrand?"_phase":"_dist")+ "_pile"+ pile+ "-"+ distance+ ".txt"));
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		

		
		getDistances(input, pile, sameStrand, distance, ps);
		
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
	static long getDistances_copy(File f, int pile, boolean sameStrand, int uDistance, int dDistance, PrintStream ps) {
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
				getDistances_copy(iter, pile, chr, pos, forward, dDistance, 1, 0, ps);
			} else {	// opposite strand
				getDistances_copy(iter, pile, chr, pos, !forward, (forward?dDistance:uDistance), (forward?1:-1), (forward?26:-26), ps);
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

	static private void getDistances_copy(BufferedIteratorDisk iter, int pile, ByteArrayCharSequence chr, int pos,
			boolean forward, int dist, int summand, int factor, PrintStream ps) {
		
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
			boolean forward2= false;
			if (cs.getToken(3).equals("F"))
				forward2= true;
			if (forward2!= forward)
				continue;
			
			// valid 
			ps.println(Integer.toString(diff));
		}
		
	}

	// pile chr pos F/R
		static long getDistances(File f, int pile, boolean sameStrand, int distance, PrintStream ps) {
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
				pos= cs.getTokenInt(2);
				
				boolean forward= false;
				if (cs.getToken(3).equals("F"))
					forward= true;
				
				// for distogram, consider only F
				// phasogram both
				if ((!sameStrand)&& (!forward))
					continue;
				
				iter.mark();
				if (sameStrand) {
					getDistances(iter, pile, chr, pos, 
							forward, distance, 0, ps);
				} else {	// opposite strand
					getDistances(iter, pile, chr, pos, 
							!forward, distance, (forward?26:-26), ps);	// can only be forward
				}
				iter.reset();
				jumped= true;
				
			}
			
			System.err.println("found "+ nrPiles+ " piles, took "+ (System.currentTimeMillis()- t0)/ 1000+ " sec.");
			return nrPiles;
		}

		static private void getDistances(BufferedIteratorDisk iter, int pile, ByteArrayCharSequence chr, int pos,
				boolean forward, int dist, int offset, PrintStream ps) {
			
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
				int diff= (pos2- pos+ offset);
				if (diff> dist)
					break;
				
				boolean forward2= false;
				if (cs.getToken(3).equals("F"))
					forward2= true;
				if (forward2!= forward)
					continue;
				
				// valid 
				ps.println(Integer.toString(Math.abs(diff)));
			}
			
		}
	
	
}
