package fbi.genome.reconstruction.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import fbi.genome.model.commons.MyHashMap;
import fbi.genome.reconstruction.closure.AliGraphClosure;
import fbi.genome.reconstruction.closure.Closure;


public class IsotigToGTF {

	static final String BLOCK_PLUS= ">>>>>";
	static final String BLOCK_MINUS= "<<<<<";
	
	static class Tig implements Comparable<Tig> {
		int start= -1;
		int end= -1;
		String id= null;
		boolean plus= true;
		
		public Tig(String id, int length) {
			this.id= id;
			this.start= 1;
			this.end= length;
		}
		public Tig(String id, int start, int end) {
			this.id= id;
			this.start= start;
			this.end= end;
		}
		public int length() {
			return (end- start+ 1);
		}
		@Override
		public int hashCode() {
			return id.hashCode();
		}
		@Override
		public boolean equals(Object paramObject) {
			return id.equals(((Tig) paramObject).id);
		}
		@Override
		public int compareTo(Tig paramT) {
			return id.compareTo(paramT.id);
		}
		@Override
		public String toString() {
			return id+"["+start+","+end+"]"+(plus?"+":"-");
		}
	}
	static class Group implements Iterable<Tig>, Iterator<Tig> {
		MyHashMap<Tig, Vector<Tig>> mapTix= null;
		MyHashMap<Tig, Vector<Tig>> mapCox= null;
		String gid= null;
		
		Vector<Tig> iterated= null;
		
		public Group(String gid) {
			this.gid= gid;
			mapTix= new MyHashMap<IsotigToGTF.Tig, Vector<Tig>>();
			mapCox= new MyHashMap<IsotigToGTF.Tig, Vector<Tig>>();
		}
		@Override
		public boolean hasNext() {
			return (iterated!= null&& iterated.size()< mapTix.size());
		}
		@Override
		public Tig next() {
			if (iterated== null|| iterated.size()== mapTix.size())
				return null;
			
			Tig current= null;
			if (iterated.size()== 0) {
				current= mapTix.keySet().iterator().next();	// random
			} else {
				for (int i = 0; current== null&& i < iterated.size(); i++) {
					Vector<Tig> cox= mapTix.get(iterated.elementAt(i));
					for (int j = 0; current== null&& j < cox.size(); j++) {
						Vector<Tig> tix= mapCox.get(cox.elementAt(j));
						for (int k = 0; current== null&& k < tix.size(); k++) {
							int p= iterated.indexOf(tix.elementAt(k));
							if (p< 0) 
								current= mapTix.getKey(tix.elementAt(k));	// baseTig
						}
					}
				}
				
				// disconnected components
				if (current== null) {
					Iterator<Tig> iter= mapTix.keySet().iterator();
					while(current== null&& iter.hasNext()) {
						Tig tig= iter.next();
						int p= iterated.indexOf(tig);
						if (p< 0) 
							current= mapTix.getKey(tig);	// baseTig
					}
				}
					
			}
			
			iterated.add(current);
			return current;
		}
		@Override
		public void remove() {
			if (iterated!= null&& iterated.size()> 0) {				
				mapTix.remove(iterated.lastElement());
				iterated.remove(iterated.size()- 1);
			}
		}
		@Override
		public Iterator<Tig> iterator() {
			iterated= new Vector<IsotigToGTF.Tig>(mapTix.size());
			return this;
		}
	}
	
	void parseIsotigLine(String line, int[] a) {
		int p= line.indexOf(' '), q= line.lastIndexOf(' ');
		for (int i= p+1; i<= q- 6; i+= 6) {
			int j= (i- p- 1)/ 6;
			String s= line.substring(i, i+ 6);
			if (s.equals(BLOCK_PLUS))
				a[j]= 1;
			else if (s.equals(BLOCK_PLUS))
				a[j]= -1;
			else
				a[j]= 0;
		}
	}
	
	int[][] parseIsogroupLine(String line) {
		String[] token= line.split("\\s");
		int nrIsos= Integer.parseInt(token[1].substring(token[1].indexOf('=')+ 1));
		int nrContig= Integer.parseInt(token[2].substring(token[2].indexOf('=')+ 1));
		int[][] a= new int[nrIsos][nrContig];
		return a;
	}
	
	/**
		 * Adds the isotig aligned blocks to the closure of the locus. If there are 
		 * any anti-sense alignment in blocks already used by other isotigs,the transcript
		 * is reversed. 
		 * @param vTix all contigs the isotig aligns to
		 * @param vCox the 
		 * @param g
		 * @param groupTx
		 * @param c
		 * @return whether the transcript was reversed.
		 */
		private static boolean addTx(Vector<Tig> vTix, Vector<Tig> vCox, Group g,
				Tig[] groupTx, Closure c) {
			
			assert(vTix.size()== vCox.size());
			
			// check tx
			boolean sense= true;
			for (int i = 0; sense&& i < vCox.size(); i++) {
				Tig tix= vTix.elementAt(i);
				Vector<Tig> tixOfCoxV= g.mapCox.get(vCox.elementAt(i));
				for (int j = 0; sense&& tixOfCoxV!= null&&
						j < tixOfCoxV.size(); j++) {
					if(tix.plus!= tixOfCoxV.elementAt(j).plus)
						sense= false;
				}
			}
			
			// swap
			if (!sense) {
				Vector<Tig> newTix= new Vector<IsotigToGTF.Tig>(),
							newCox= new Vector<IsotigToGTF.Tig>();
				for (int i = vTix.size()- 1; i >= 0; --i) {
					Tig h= vTix.elementAt(i);
					h.plus= !h.plus;
					int start= (newTix.size()== 0? 1: newTix.elementAt(newTix.size()- 1).end+ 1);
					int end= start+ h.length()- 1;
					h.start= start;
					h.end= end;
					newTix.add(h);
					
					h= vCox.elementAt(i);
					h.plus= !h.plus;
					newCox.add(h);
				}
				vTix= newTix;
				vCox= newCox;
			}
			
			// align
			for (int h = 0; h < vTix.size(); ++h) {
				
				Tig tix= vTix.elementAt(h);	// isotig00001[x,y]+
				Tig cox= vCox.elementAt(h); // contig00001[x,y]+
				
				Tig baseTix= g.mapTix.getKey(tix);
	//			if (!baseTix.plus)
	//				tix.plus= !tix.plus;
				Vector<Tig> tixAli= g.mapTix.get(baseTix);
				if (tixAli== null) 
					tixAli= new Vector<IsotigToGTF.Tig>();
				
				Tig baseCox= g.mapCox.getKey(cox);
				if (baseCox.start!= cox.start|| baseCox.end!= cox.end) {
					;//System.err.print("*");
					//continue;
					//System.err.println("Partial cox alignment ["+ cox.start+ ","+ cox.end+ "] <> ["+
					//		baseCox.start+ ","+ baseCox.end+ "]");
				}
				Vector<Tig> coxAli= g.mapCox.get(baseCox);
				if (coxAli== null) 
					coxAli= new Vector<IsotigToGTF.Tig>();
	
				// iterate tix already aligned to that cox
				int idx1= findIdx(tix, groupTx);
				for (int i = 0; i < coxAli.size(); i++) {
					Tig otherTx= coxAli.elementAt(i);
					int p= g.mapTix.get(otherTx).indexOf(baseCox);
					int idx2= findIdx(otherTx, groupTx);
					if (otherTx.plus== tix.plus) {
						// align the tix positions
						for (int j = tix.start, k= otherTx.start; j <= tix.end&& k<= otherTx.end; ++j, ++k) {
							if (!AliGraphClosure.alignablePositions(c, idx1, h+1, idx2, p+1)) {
							//if (!AliGraphClosure.alignablePositions(c, idx1, j+1, idx2, k+1)) {
								System.err.println("\nInconsistent alignment "+ tix);
								break;
							}
							AliGraphClosure.addAlignedPositions(c, idx1, h+1, idx2, p+1);
							AliGraphClosure.computeClosure(c);
						}
					} else {
						System.err.println("\nInconsistent reverse alignment "+ tix);
						for (int j = tix.start, k= otherTx.end; j <= tix.end&& k>= otherTx.start; ++j, --k) {
	//						if (!AliGraphClosure.alignablePositions(c, idx1, j, idx2, k)) {
	//							int lb= AliGraphClosure.getLowerBound(c, idx1, j, idx2);
	//							int ub= AliGraphClosure.getUpperBound(c, idx1, j, idx2);
	//							System.err.println("\nInconsistent reverse alignment "+ tix);
	//							break;
	//						}
	//						AliGraphClosure.addAlignedPositions(c, idx1, j, idx2, k);
						}
					}
				}
				
				
				coxAli.add(tix);
				tixAli.add(cox);
				g.mapCox.put(baseCox, coxAli);
				g.mapTix.put(baseTix, tixAli);
			}
			
			return (sense);
		}

	/**
		 * Adds the isotigs according to their connections to a newly
		 * created group.
		 * @param g
		 */
		static Group resolveGroup(Group g) {
			
			// alloc closure
			Tig[] groupTx= new Tig[g.mapTix.size()];
			g.mapTix.keySet().toArray(groupTx);
			Arrays.sort(groupTx);
			int[] len= new int[groupTx.length];
			for (int i = 0; i < len.length; i++) 
				len[i]= g.mapTix.get(groupTx[i]).size();
			Closure c= AliGraphClosure.newAligGraphClosure(len.length, len, 0, null);
	
			// prepare new group to add incrementally the isotigs to
			Group newGroup= new Group(g.gid);
			Iterator<Tig> iter= g.mapCox.keySet().iterator();
			while(iter.hasNext())
				newGroup.mapCox.put(iter.next(), null);
			iter= g.mapTix.keySet().iterator();
			while(iter.hasNext())
				newGroup.mapTix.put(iter.next(), null);
			
			// add transcripts in a good order
			int cntSense= 0, cntAsense= 0;
			Iterator<Tig> tigiter= g.iterator();
			while(tigiter.hasNext()) {
				Tig tig= tigiter.next();
				//System.err.println("adding "+ tig);
				
				Vector<Tig> vCox= g.mapTix.get(tig);
				Vector<Tig> vTix= new Vector<Tig>(vCox.size());
				for (int i = 0; i < vCox.size(); i++) {
					Vector<Tig> v= g.mapCox.get(vCox.elementAt(i));
					Tig t= v.get(v.indexOf(tig));
					vTix.add(t);
				}
				
				boolean sense= addTx(vTix, vCox, newGroup, groupTx, c);
				newGroup.mapTix.getKey(tig).plus= sense;
				if (sense)
					++cntSense;
				else
					++cntAsense;
			}
			
			// clean-up contigs, reverse the ones that break the linear order
	/*		iter= newGroup.mapTix.keySet().iterator();
			while(iter.hasNext()) {
				Tig baseTig= iter.next();	// knows transcript directionality
				Vector<Tig> coxV= newGroup.mapTix.get(baseTig);
				for (int i = 0; i < coxV.size(); i++) {
					if (coxV.elementAt(i).plus!= baseTig.plus)
						swapContig(newGroup, coxV.elementAt(i));
				}
			}
			// re-check
			iter= newGroup.mapTix.keySet().iterator();
			while(iter.hasNext()) {
				Tig baseTig= iter.next();	// knows transcript directionality
				Vector<Tig> coxV= newGroup.mapTix.get(baseTig);
				for (int i = 0; i < coxV.size(); i++) {
					if (coxV.elementAt(i).plus!= baseTig.plus)
						System.err.println("Inconsistent Contig "+ coxV.elementAt(i));
				}
			}
	*/		
			// turn around locus, s.t. there are more sense isotigs
			if (cntAsense> cntSense) {
				// TODO
			}		
			System.err.print("Group "+ newGroup.gid+ " "+ cntSense+ " sense, "
					+ cntAsense+ " anti-sense isotigs.");
			System.err.flush();
			
			// change names (+/-)
			for (int i = 0; i < groupTx.length; i++) 
				groupTx[i]= newGroup.mapTix.getKey(groupTx[i]);
			
			
			// output
			Vector<Integer> vLen= new Vector<Integer>();
			Vector<Boolean>[] vMap= new Vector[len.length];
			makeLayout(c, len, vLen, vMap);
//			printLayout(groupTx, vLen, vMap, System.err);
			
			Vector<Tig> v= getPlayground(newGroup, groupTx, len, vMap);
//			for (int i = 0; i < v.size(); i++) 
//				System.err.print("\t"+ v.elementAt(i));
//			System.err.println();
			
			File f= new File("/Users/micha/isotest.gtf");
			writeGTF(f, newGroup, groupTx, len, vMap, v);
			
			return newGroup;
		}

	private static void writeGTF(File f, Group newGroup, Tig[] groupTx, int[] len,
			Vector<Boolean>[] vMap, Vector<Tig> pground) {
		
		final String tab= "\t", exon= "exon", plus= "+", minus= "-", 
				dot= ".", space= " ", txID= "transcriptID", quot= "\"", 
				semic= ";", nl= "\n";
		
		try {
			BufferedWriter writer= new BufferedWriter(new FileWriter(f));
			for (int i = 0; i < vMap.length; i++) {
				Tig baseTix= newGroup.mapTix.getKey(groupTx[i]);
				int cntTxCox= 0;
				int cumuLen= 0;
				//Vector<Tig> vCox= newGroup.mapTix.get(baseTix);
				
				for (int j = 0; j < vMap[i].size(); j++) {
					if (!vMap[i].elementAt(j)) {
						cumuLen+= pground.elementAt(j).length();
						continue;
					}
					
					Tig cox= newGroup.mapTix.get(baseTix).elementAt(cntTxCox);
					Tig baseCox= newGroup.mapCox.getKey(cox);
					Vector<Tig> vTix= newGroup.mapCox.get(cox);
					Tig tix= vTix.elementAt(vTix.indexOf(baseTix));
					
					StringBuilder sb= new StringBuilder(newGroup.gid);
					sb.append(tab);
					sb.append(cox.id);
					// directionality of contig sequence when assembling the tx sequence
					sb.append(cox.plus^ baseTix.plus? minus: plus);
					
					sb.append(tab);
					sb.append(exon);
					sb.append(tab);
					
					sb.append(tab);
					sb.append(Integer.toString(cumuLen+ 1));
					sb.append(tab);					
					cumuLen+= baseCox.length();	// take base, alignments variable!
					sb.append(Integer.toString(cumuLen));
					
					sb.append(tab);
					// score: delta of alignment with respect to contig
					sb.append(Integer.toString(tix.length()- baseCox.length()));

					sb.append(tab);
					sb.append(baseTix.plus? plus: minus);
					
					sb.append(tab);
					sb.append(dot);

					sb.append(tab);
					sb.append(txID);
					sb.append(space);
					sb.append(quot);
					sb.append(baseTix.id);
					sb.append(quot);
					sb.append(semic);
					sb.append(nl);
					
					++cntTxCox;
					
					writer.write(sb.toString());
				}
			}
			
			writer.flush();
			writer.close();
		} catch (IOException e) {			
			e.printStackTrace();
		}
		

	}

	private static Vector<Tig> getPlayground(Group newGroup, Tig[] groupTx, int[] len,
			Vector<Boolean>[] vMap) {
		
		int[] p= new int[groupTx.length];
		Arrays.fill(p, 0);
		Vector<Tig> pground= new Vector<IsotigToGTF.Tig>(vMap[0].size());
		
		
		for (int i = 0; i < vMap[0].size(); i++) {
			Tig currentCox= null;
			for (int j = 0; j < vMap.length; j++) {
				if (!vMap[j].elementAt(i))
					continue;
				
				
				Vector<Tig> v= newGroup.mapTix.get(groupTx[j]);
				if (currentCox== null) 
					currentCox= v.elementAt(p[j]);
				else {
					if (!v.elementAt(p[j]).equals(currentCox))
						System.err.println("Inconsistent playground col "+ i+" in "+ groupTx[j]+": "+ 
								currentCox+ " <> "+ v.elementAt(p[j]));					
				}				
				++p[j];
			}
			if (currentCox== null)
				System.err.println("Null column "+ i);
			pground.add(newGroup.mapCox.getKey(currentCox));
		}
		
		for (int i = 0; i < p.length; i++) {
			if (p[i]!= len[i])
				System.err.println("Incomplete playground "+ i+ ": "+ p[i]+ "<>"+ len[i]);
		}
		
		return pground;
	}

	static void parseIsotig(File isotigFile, HashMap<String, Group> map) {
		
		// OBS: 5th always W, 7th always 1
		/* isotig00041	200	225	3	W	contig21293	1	26	+
		 * isotig00041	226	1308	4	W	contig07350	1	1083	-
		 */
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(isotigFile));
			Group lastG= null;
			Tig lastTx= null;
			int txInGroup= 0;
			long t= System.currentTimeMillis();
			
			for(String line; (line= buffy.readLine())!= null;) {
	
				if (line.startsWith("contig"))
					continue;	// TODO non-AS contigs.. (none in IsotigLayout)
				
				String[] s= line.split("\\s");
				Group g= map.get(s[0]);
				Tig cox= new Tig(s[5].substring(6), 
						Integer.parseInt(s[6]), 
						Integer.parseInt(s[7]));	// (-1)+ 1= 0
				Tig tix= new Tig(s[0], 
						Integer.parseInt(s[1]),
						Integer.parseInt(s[2]));
				tix.plus= cox.plus= s[8].equals("+");
				if (cox.length()!= tix.length())
					System.err.println("Inconsistent alignment length "+ cox.length()+"<>"+tix.length()+"\n"+line);

				if (g== lastG) {
					if (!s[0].equals(lastTx.id)) {	// new tx
						++txInGroup;
						Tig lastBaseTx= g.mapTix.getKey(lastTx);
						if(lastTx.end!= lastBaseTx.length())
							System.err.println("Check tx length "+ lastTx.id+ " "+ lastBaseTx.length()+ "<>"+ lastTx.end);
					}
				
				} else {	// new group
					if (lastG!= null) { 
						
						// check
						if (txInGroup!= lastG.mapTix.size())
							System.err.println("Check group size "+ lastG.gid+ " "+ lastG.mapTix.size()+"<>"+ txInGroup);
						txInGroup= 0;
						
						resolveGroup(lastG);
						System.err.println(" "+ ((System.currentTimeMillis()- t)/ 1000)+ " sec.");
						lastG= g;
						t= System.currentTimeMillis();
					}
					++txInGroup;
				}
				
				// add
				Vector<Tig> tigV= g.mapTix.get(tix); 
				if (tigV== null) {
					tigV= new Vector<IsotigToGTF.Tig>();
					g.mapTix.put(tix, tigV);
				}
				tigV.add(cox);
				Vector<Tig> cogV= g.mapCox.get(cox);
				if (cogV== null) {
					cogV= new Vector<IsotigToGTF.Tig>();
					g.mapCox.put(cox, cogV);
				}
				cogV.add(tix);

				lastG= g;
				lastTx= tix;
			}
			
			resolveGroup(lastG);
			System.err.println(" "+ (System.currentTimeMillis()- t)/ 1000);
							
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Adds the isotigs according to their connections to a newly
	 * created group.
	 * @param g
	 */
	static Group resolveGroup_nt(Group g) {
		
		// alloc closure
		Tig[] groupTx= new Tig[g.mapTix.size()];
		g.mapTix.keySet().toArray(groupTx);
		Arrays.sort(groupTx);
		int[] len= new int[groupTx.length];
		for (int i = 0; i < len.length; i++) 
			len[i]= groupTx[i].length();
		Closure c= AliGraphClosure.newAligGraphClosure(len.length, len, 0, null);

		// prepare new group to add incrementally the isotigs to
		Group newGroup= new Group(g.gid);
		Iterator<Tig> iter= g.mapCox.keySet().iterator();
		while(iter.hasNext())
			newGroup.mapCox.put(iter.next(), null);
		iter= g.mapTix.keySet().iterator();
		while(iter.hasNext())
			newGroup.mapTix.put(iter.next(), null);
		
		// add transcripts in a good order
		int cntSense= 0, cntAsense= 0;
		Iterator<Tig> tigiter= g.iterator();
		while(tigiter.hasNext()) {
			Tig tig= tigiter.next();
			System.err.println("adding "+ tig);
			
			Vector<Tig> vCox= g.mapTix.get(tig);
			Vector<Tig> vTix= new Vector<Tig>(vCox.size());
			for (int i = 0; i < vCox.size(); i++) {
				Vector<Tig> v= g.mapCox.get(vCox.elementAt(i));
				Tig t= v.get(v.indexOf(tig));
				vTix.add(t);
			}
			
			boolean sense= addTx(vTix, vCox, newGroup, groupTx, c);
			newGroup.mapTix.getKey(tig).plus= sense;
			if (sense)
				++cntSense;
			else
				++cntAsense;
		}
		
		// clean-up contigs, reverse the ones that break the linear order
/*		iter= newGroup.mapTix.keySet().iterator();
		while(iter.hasNext()) {
			Tig baseTig= iter.next();	// knows transcript directionality
			Vector<Tig> coxV= newGroup.mapTix.get(baseTig);
			for (int i = 0; i < coxV.size(); i++) {
				if (coxV.elementAt(i).plus!= baseTig.plus)
					swapContig(newGroup, coxV.elementAt(i));
			}
		}
		// re-check
		iter= newGroup.mapTix.keySet().iterator();
		while(iter.hasNext()) {
			Tig baseTig= iter.next();	// knows transcript directionality
			Vector<Tig> coxV= newGroup.mapTix.get(baseTig);
			for (int i = 0; i < coxV.size(); i++) {
				if (coxV.elementAt(i).plus!= baseTig.plus)
					System.err.println("Inconsistent Contig "+ coxV.elementAt(i));
			}
		}
*/		
		// turn around locus, s.t. there are more sense isotigs
		if (cntAsense> cntSense) {
			// TODO
		}		
		System.err.println("Group "+ newGroup.gid+ " "+ cntSense+ " sense, "
				+ cntAsense+ " anti-sense isotigs.");
		
		// output
		Vector<Integer> vLen= new Vector<Integer>();
		Vector<Boolean>[] vMap= new Vector[len.length];
		makeLayout(c, len, vLen, vMap);
		printLayout(groupTx, vLen, vMap, System.err);
		
		return newGroup;
	}

	private static void swapContig(Group g, Tig c) {

		Tig baseCox= g.mapCox.getKey(c); // enforce base
		baseCox.plus= !baseCox.plus;
		
		Vector<Tig> vTix= g.mapCox.get(c);
		for (int i = 0; i < vTix.size(); i++) {
			Tig tix= vTix.elementAt(i);
			tix.plus= !tix.plus; // turn around alignment isotig side
			
			Vector<Tig> vCox= g.mapTix.get(tix);
			Tig cox= vCox.get(vCox.indexOf(baseCox));
			cox.plus= !cox.plus;
		}
	}

	/**
	 * Adds the isotig aligned blocks to the closure of the locus. If there are 
	 * any anti-sense alignment in blocks already used by other isotigs,the transcript
	 * is reversed. 
	 * @param vTix all contigs the isotig aligns to
	 * @param vCox the 
	 * @param g
	 * @param groupTx
	 * @param c
	 * @return whether the transcript was reversed.
	 */
	private static boolean addTx_nt(Vector<Tig> vTix, Vector<Tig> vCox, Group g,
			Tig[] groupTx, Closure c) {
		
		assert(vTix.size()== vCox.size());
		
		// check tx
		boolean sense= true;
		int cnt= 0;
		for (int i = 0; /*sense&&*/ i < vCox.size(); i++) {
			Tig tix= vTix.elementAt(i);
			Vector<Tig> tixOfCoxV= g.mapCox.get(vCox.elementAt(i));
			for (int j = 0; /*sense&&*/ tixOfCoxV!= null&&
					j < tixOfCoxV.size(); j++) {
				if(tix.plus!= tixOfCoxV.elementAt(j).plus)
					sense= false;
				else
					++cnt;
			}
		}
		
		// swap
		if (!sense) {
			Vector<Tig> newTix= new Vector<IsotigToGTF.Tig>(),
						newCox= new Vector<IsotigToGTF.Tig>();
			for (int i = vTix.size()- 1; i >= 0; --i) {
				Tig h= vTix.elementAt(i);
				h.plus= !h.plus;
				int start= (newTix.size()== 0? 1: newTix.elementAt(newTix.size()- 1).end+ 1);
				int end= start+ h.length()- 1;
				h.start= start;
				h.end= end;
				newTix.add(h);
				
				h= vCox.elementAt(i);
				h.plus= !h.plus;
				newCox.add(h);
			}
			vTix= newTix;
			vCox= newCox;
		}
		
		// align
		for (int h = 0; h < vTix.size(); ++h) {
			
			Tig tix= vTix.elementAt(h);	// isotig00001[x,y]+
			Tig cox= vCox.elementAt(h); // contig00001[x,y]+
			
			Tig baseTix= g.mapTix.getKey(tix);
//			if (!baseTix.plus)
//				tix.plus= !tix.plus;
			Vector<Tig> tixAli= g.mapTix.get(baseTix);
			if (tixAli== null) 
				tixAli= new Vector<IsotigToGTF.Tig>();
			
			Tig baseCox= g.mapCox.getKey(cox);
			if (baseCox.start!= cox.start|| baseCox.end!= cox.end) {
				System.err.print("*");
				//continue;
				//System.err.println("Partial cox alignment ["+ cox.start+ ","+ cox.end+ "] <> ["+
				//		baseCox.start+ ","+ baseCox.end+ "]");
			}
			Vector<Tig> coxAli= g.mapCox.get(baseCox);
			if (coxAli== null) 
				coxAli= new Vector<IsotigToGTF.Tig>();

			// iterate tix already aligned to that cox
			int idx1= findIdx(tix, groupTx);
			for (int i = 0; i < coxAli.size(); i++) {
				Tig otherTx= coxAli.elementAt(i);
				int idx2= findIdx(otherTx, groupTx);
				if (otherTx.plus== tix.plus) {
					// align the tix positions
					for (int j = tix.start, k= otherTx.start; j <= tix.end&& k<= otherTx.end; ++j, ++k) {
						if (!AliGraphClosure.alignablePositions(c, idx1, j+1, idx2, k+1)) {
							System.err.println("\nInconsistent alignment "+ tix);
							break;
						}
						AliGraphClosure.addAlignedPositions(c, idx1, j+1, idx2, k+1);
						AliGraphClosure.computeClosure(c);
					}
				} else {
					System.err.println("\nInconsistent reverse alignment "+ tix);
					for (int j = tix.start, k= otherTx.end; j <= tix.end&& k>= otherTx.start; ++j, --k) {
//						if (!AliGraphClosure.alignablePositions(c, idx1, j, idx2, k)) {
//							int lb= AliGraphClosure.getLowerBound(c, idx1, j, idx2);
//							int ub= AliGraphClosure.getUpperBound(c, idx1, j, idx2);
//							System.err.println("\nInconsistent reverse alignment "+ tix);
//							break;
//						}
//						AliGraphClosure.addAlignedPositions(c, idx1, j, idx2, k);
					}
				}
			}
			
			
			coxAli.add(tix);
			tixAli.add(cox);
			g.mapCox.put(baseCox, coxAli);
			g.mapTix.put(baseTix, tixAli);
		}
		
		return (sense);
	}

	static HashMap<String, Group> parseLayout(File layoutFile) {
		/*
		 * >isogroup00013  numIsotigs=45  numContigs=63
		 *     Length : 123    45     423    ... (bp)
		 *     Contig : 123456 123456 123456 ... Total:
		 * isotig000001 <<<<<<        >>>>>> ... 1080
		 */
		try {			
			BufferedReader buffy= new BufferedReader(new FileReader(layoutFile));
			HashMap<String, Group> map= new HashMap<String, IsotigToGTF.Group>();
			System.err.print("Reading layouts ");
			System.err.flush();
			String line= buffy.readLine();
			int ctr= 0;
			long t0= System.currentTimeMillis();
			while (line!= null) {
				if (line.startsWith(">")) {
					
					++ctr;
					if (ctr%2500== 0) {
						System.err.print("*");
						System.err.flush();
					}
						
					
					// ID
					String[] s= line.trim().split("\\s+");
					String gid= s[0].substring(1);
					int nrTix= Integer.parseInt(s[1].substring(s[1].indexOf("=")+ 1));
					int nrCox= Integer.parseInt(s[2].substring(s[2].indexOf("=")+ 1));
					Group tixGroup= new Group(gid);
					
					// Length and Contig
					while(!line.substring(0,20).contains("Length"))
						line=buffy.readLine();						
					s= line.trim().split("\\s+");
					if (s.length!= nrCox+ 3)
						System.err.println("Wrong number of contigs "+ nrCox+ "<>"+ (s.length- 2));
					while(!line.substring(0,20).contains("Contig"))
						line=buffy.readLine();						
					String[] t= line.trim().split("\\s+");
					if (s.length!= nrCox+ 3) // Length, :, ..., (bp)
						System.err.println("Wrong number of contigs "+ nrCox+ "<>"+ (t.length- 2));
					for (int i = 2; i < s.length- 1; i++) {
						Tig tig= new Tig(t[i], Integer.parseInt(s[i]));
						tixGroup.mapCox.put(tig, null);	// vector for isotix
					}
	
					int cntTx= 0;
					while ((line= buffy.readLine())!= null&& !line.startsWith(">")) {
						if (line.trim().length()== 0) 
							continue;
						if (line.startsWith("NumContigsInIsogroup")) {
							line= null;
							break;
						}
						String tixID= line.substring(0,line.indexOf(' '));
						if (!tixID.startsWith("isotig"))
							System.err.println("Unexpected isotig line: "+ line);
						Tig tig= new Tig(tixID, Integer.parseInt(line.substring(line.lastIndexOf(' ')+ 1 )));
						tixGroup.mapTix.put(tig, null);
						map.put(tixID, tixGroup);
						++cntTx;
					}
					if (cntTx!= nrTix)
						System.err.println("Wront number of tix "+ nrTix+ "<>"+ cntTx);
						
	
				} else 
					line= buffy.readLine();
			}
			
			buffy.close();
			System.err.println("\nRead "+ map.size()+ " isotigs "+ ((System.currentTimeMillis()- t0)/1000)+ " sec.");
			return map;
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	static void parseIsotig_save(File isotigFile, HashMap<String, Group> map) {
		
		// OBS: 5th always W, 7th always 1
		/* isotig00041	200	225	3	W	contig21293	1	26	+
		 * isotig00041	226	1308	4	W	contig07350	1	1083	-
		 */
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(isotigFile));
			Group lastG= null;
			Tig lastTx= null;
			Closure c= null;
			int txInGroup= 0;
			Tig[] groupTx= null;
			
			for(String line; (line= buffy.readLine())!= null;) {

				String[] s= line.split("\\s");
				Group g= map.get(s[0]);
				Tig cox= new Tig(s[5].substring(6), 
						Integer.parseInt(s[6]), 
						Integer.parseInt(s[7]));	// (-1)+ 1= 0
				Tig tix= new Tig(s[0], 
						Integer.parseInt(s[1]),
						Integer.parseInt(s[2]));
				tix.plus= cox.plus= s[8].equals("+");
				if (cox.length()!= tix.length())
					System.err.println("Inconsistent alignment length "+ cox.length()+"<>"+tix.length()+"\n"+line);

				Tig baseTix= g.mapTix.getKey(tix);
				if (!baseTix.plus)
					tix.plus= !tix.plus;
				Vector<Tig> tixAli= g.mapTix.get(baseTix);
				if (tixAli== null) 
					tixAli= new Vector<IsotigToGTF.Tig>();
				
				Tig baseCox= g.mapCox.getKey(cox);
				if (baseCox.start!= cox.start|| baseCox.end!= cox.end) {
					System.err.print("*");
					//continue;
					//System.err.println("Partial cox alignment ["+ cox.start+ ","+ cox.end+ "] <> ["+
					//		baseCox.start+ ","+ baseCox.end+ "]");
				}
				Vector<Tig> coxAli= g.mapCox.get(baseCox);
				if (coxAli== null) 
					coxAli= new Vector<IsotigToGTF.Tig>();

				if (g== lastG) {
					++txInGroup;
					if (!s[0].equals(lastTx.id)) {
						Tig lastBaseTx= g.mapTix.getKey(lastTx);
						if(lastTx.end!= lastBaseTx.length())
							System.err.println("Check tx length "+ lastTx.id+ " "+ lastBaseTx.length()+ "<>"+ lastTx.end);
						lastTx= tix;
					}
				
				} else {	// new group
					if (lastG!= null&& txInGroup!= lastG.mapTix.size())
						System.err.println("Check group size "+ lastG.gid+ " "+ lastG.mapTix.size()+"<>"+ txInGroup);
					txInGroup= 0;
					
					// alloc closure
					groupTx= new Tig[g.mapTix.size()];
					g.mapTix.keySet().toArray(groupTx);
					Arrays.sort(groupTx);
					int[] len= new int[groupTx.length];
					for (int i = 0; i < len.length; i++) 
						len[i]= groupTx[i].length();
					c= AliGraphClosure.newAligGraphClosure(len.length, len, 0, null);
					
				}
				lastG= g;
				lastTx= tix;
				
				int idx1= findIdx(tix, groupTx);
				for (int i = 0; i < coxAli.size(); i++) {
					Tig otherTx= coxAli.elementAt(i);
					int idx2= findIdx(otherTx, groupTx);
					if (otherTx.plus== tix.plus) {
						for (int j = tix.start, k= otherTx.start; j <= tix.end&& k<= otherTx.end; ++j, ++k) {
							if (!AliGraphClosure.alignablePositions(c, idx1, j, idx2, k)) {
								System.err.println("\nInconsistent alignment "+ tix);
								break;
							}
							AliGraphClosure.addAlignedPositions(c, idx1, j, idx2, k);
						}
					} else {
						if (baseTix.plus) {
							if (tix.start== 1) {
								baseTix.plus= false; // turn around
							} else {
								System.err.println("Problem");
								continue;
							}
						}
						for (int j = tix.start, k= otherTx.end; j <= tix.end&& k>= otherTx.start; ++j, --k) {
							if (!AliGraphClosure.alignablePositions(c, idx1, j, idx2, k)) {
								int lb= AliGraphClosure.getLowerBound(c, idx1, j, idx2);
								int ub= AliGraphClosure.getUpperBound(c, idx1, j, idx2);
								System.err.println("\nInconsistent reverse alignment "+ tix);
								break;
							}
							AliGraphClosure.addAlignedPositions(c, idx1, j, idx2, k);
						}
					}
				}
				
				
				coxAli.add(tix);
				tixAli.add(cox);
				g.mapCox.put(baseCox, coxAli);
				g.mapTix.put(baseTix, tixAli);
				
			}
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	private static int findIdx(Tig tix, Tig[] groupTx) {
		int p= Arrays.binarySearch(groupTx, tix);
		if (p< 0)
			System.err.println("Not found "+ tix);
		return p;
	}

	static void makeLayout(Closure c, int[] seqLen, Vector<Integer> vLen, Vector<Boolean>[] vMap) {
		
		// trivial case, would deadlock
		if (seqLen.length== 1) {
			vMap= new Vector[1];
			vMap[0]= new Vector<Boolean>();
			for (int i = 0; i < seqLen[0]; i++) 
				vMap[0].add(true);
			return;
		}
			
		
		int[] p= new int[seqLen.length];
		Arrays.fill(p, 1);
		
		for (int i = 0; i < vMap.length; i++) 
			vMap[i]= new Vector<Boolean>();
		
		boolean left= true;
		while(left) {
			
			int[] lo= new int[seqLen.length];
			Arrays.fill(lo, Integer.MAX_VALUE);
			int[] hi= new int[seqLen.length];
			Arrays.fill(hi, Integer.MAX_VALUE);
			
			// find min blocks
			for (int i = 0; i < p.length; i++) {
				if (p[i]> seqLen[i])
					continue;
				for (int j = 0; j < hi.length; j++) {
					if (i== j|| p[j]> seqLen[j])
						continue;
					int lb= AliGraphClosure.predFrontier(c, i, p[i], j);
					int ub= AliGraphClosure.succFrontier(c, i, p[i], j);
					
					// TODO check this condition
					//if (lb<= lo[j]&& ub<= hi[j]) {
//					lo[j]= Math.min(lo[j], lb);
//					hi[j]= Math.min(hi[j], ub);
//					else
//					hi[j]= Math.min(lb, hi[j]);	// lb included in next

					if (lb== ub) { // matching
						if (lb== p[j]) {
							lo[j]= p[j];
							hi[j]= p[j];
						} else {
							hi[j]= Math.min(lb, hi[j]);	// exclude from potential range
						}
					} else { // range, possibly empty (lb< ub)
						if (!(lo[j]== hi[j]&& lo[j]== p[j])) {	// preserve match
							if (lb>= p[j])
								hi[j]= Math.min(lb+1, hi[j]);	// hope that lo[j] is ok
							else {
								lo[j]= p[j]- 1;
								hi[j]= Math.min(ub, hi[j]);
							}
						}
					}
					
				}
			}

			// move, top-down
			for (int i = 0; i < p.length; i++) {
				if (lo[i]== Integer.MAX_VALUE|| hi[i]== Integer.MAX_VALUE 
						|| lo[i]> hi[i]	// some sequence cannot align here
						|| lo[i]> p[i]|| hi[i]< p[i])
					continue; // lo= hi+ 1
				
				if (i>= 31)
					System.currentTimeMillis();
				
				if (lo[i]== hi[i]) {	// a match
					// find with it aligns
					// TODO extend: Vector<Integer> set= new Vector<Integer>();
					vLen.add(1);
					vMap[i].add(true);
					for (int j = 0; j < p.length; j++) {
						if (j== i)
							continue;
						if (lo[j]== Integer.MAX_VALUE) {
							vMap[j].add(false);
							continue;
						}
							
						int lb= AliGraphClosure.predFrontier(c, j, p[j], i);
						int ub= AliGraphClosure.succFrontier(c, j, p[j], i);

						// decide on aligned without/before converting coordinates
						vMap[j].add(lb== ub&& lb== lo[i]);
						if (vMap[j].lastElement()) {
							lo[j]= Integer.MAX_VALUE;	// skip
							++p[j];
						}
					}
					++p[i];
				} else {	// it's a block
					int len= hi[i]- lo[i]- 1;
					for (int k = 0; k < len; ++k) {
						vLen.add(1);
						for (int j = 0; j < vMap.length; j++)  
							vMap[j].add(i== j);
						p[i]+= 1;
					}
				}
				
			}
			
			// stop condition
			left= false;
			for (int i = 0; i < p.length; i++) 
				left|= (p[i]<= seqLen[i]);
		}
	}
	
	static boolean[][] joinPairwise(Closure c, int s1, int s2, int len1, int len2, boolean[][] a) {
		
		boolean[] a1= new boolean[len1+ len2];
		boolean[] a2= new boolean[len1+ len2];
		
		int totLen= 0, p1= 0, p2= 0;
		while(p1< len1|| p2< len2) {
			// check p1
			int[] b= AliGraphClosure.getTransitivityBounds(c, s1, p1, s2);
			--b[0]; // return val is 1-based
			--b[1];
			
			if (p1< len1) {
				if(b[0]== b[1]) {
					if (b[0]== p2) {
						a1[totLen]= true;
						a2[totLen++]= true;
						++p2;
						++p1;
						continue;
					}
				} else if (b[0]< b[1]) {
					if (p2> b[0]) {
						int blen= b[1]- b[0]- 1;
						for (int i = 0; i < blen; i++) {
							a1[totLen]= true;
							a2[totLen++]= false;
						}
						p1+= blen;
						continue;
					}
				} else
					System.out.println("check");
			}
			
			// check p2
			if (p2< len2) {
				b= AliGraphClosure.getTransitivityBounds(c, s2, p2, s1);
				--b[0];
				--b[1];
				if(b[0]== b[1]) {
					if (b[0]== p1) {
						a1[totLen]= true;
						a2[totLen++]= true;
						++p1;
						++p2;
						continue;
					}
				} else if (b[0]< b[1]) {
					if (p2> b[0]) {
						int blen= b[1]- b[0]- 1;
						for (int i = 0; i < blen; i++) {
							a1[totLen]= false;
							a2[totLen++]= true;
						}
						p2+= blen;
						continue;
					}
				} else
					System.out.println("check");
			}
		}
	
		return new boolean[][] {a1, a2};
		
	}

	static HashMap<String, Group> parseLayout_first(File layoutFile) {
		/*
		 * >isogroup00013  numIsotigs=45  numContigs=63
		 *     Length : 123    45     423    ... (bp)
		 *     Contig : 123456 123456 123456 ... Total:
		 * isotig000001 <<<<<<        >>>>>> ... 1080
		 */
		try {			
			BufferedReader buffy= new BufferedReader(new FileReader(layoutFile));
			HashMap<String, Group> map= new HashMap<String, IsotigToGTF.Group>();
			String line= buffy.readLine();
			while (line!= null) {
				if (line.startsWith(">")) {
					// ID
					String[] s= line.trim().split("\\s+");
					String gid= s[0].substring(1);
					int nrTix= Integer.parseInt(s[1].substring(s[1].indexOf("=")+ 1));
					int nrCox= Integer.parseInt(s[2].substring(s[2].indexOf("=")+ 1));
					Group tixGroup= new Group(gid);
					
					// Length and Contig
					while(!line.substring(0,20).contains("Length"))
						line=buffy.readLine();						
					s= line.trim().split("\\s+");
					if (s.length!= nrCox+ 3)
						System.err.println("Wrong number of contigs "+ nrCox+ "<>"+ (s.length- 2));
					while(!line.substring(0,20).contains("Contig"))
						line=buffy.readLine();						
					String[] t= line.trim().split("\\s+");
					if (s.length!= nrCox+ 3) // Length, :, ..., (bp)
						System.err.println("Wrong number of contigs "+ nrCox+ "<>"+ (t.length- 2));
					for (int i = 2; i < s.length- 1; i++) {
						Tig tig= new Tig(t[i], Integer.parseInt(s[i]));
						tixGroup.mapCox.put(tig, null);	// vector for isotix
					}

					int cntTx= 0;
					while ((line= buffy.readLine())!= null&& !line.startsWith(">")) {
						if (line.trim().length()== 0) 
							continue;
						String tixID= line.substring(0,line.indexOf(' '));
						if (!tixID.startsWith("isotig"))
							System.err.println("Unexpected isotig line: "+ line);
						Tig tig= new Tig(tixID, Integer.parseInt(line.substring(line.lastIndexOf(' ')+ 1 )));
						tixGroup.mapTix.put(tig, null);
						map.put(tixID, tixGroup);
						++cntTx;
					}
					if (cntTx!= nrTix)
						System.err.println("Wront number of tix "+ nrTix+ "<>"+ cntTx);
						

				} else 
					line= buffy.readLine();
			}
			
			buffy.close();
			System.err.println("Read "+ map.size()+ " isotigs.");
			return map;
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	static void parseIsotig_first(File isotigFile, HashMap<String, Group> map) {
		
		// OBS: 5th always W, 7th always 1
		/* isotig00041	200	225	3	W	contig21293	1	26	+
		 * isotig00041	226	1308	4	W	contig07350	1	1083	-
		 */
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(isotigFile));
			for(String line; (line= buffy.readLine())!= null;) {
				String[] s= line.split("\\s");
				Group g= map.get(s[0]);
				
				Tig cox= new Tig(s[5].substring(6), Integer.parseInt(s[7]));	// (-1)+ 1= 0
				Tig tix= new Tig(s[0], Integer.parseInt(s[2])- Integer.parseInt(s[1])+ 1);				
				if (cox.length()!= tix.length())
					System.err.println("Inconsistent alignment length "+ cox.length()+"<>"+tix.length()+"\n"+line);
				
				Tig baseCox= g.mapCox.getKey(cox);
				Vector<Tig> coxAli= g.mapCox.remove(baseCox);
				if (coxAli== null) 
					coxAli= new Vector<IsotigToGTF.Tig>();
				Tig baseTix= g.mapTix.getKey(tix);
				Vector<Tig> tixAli= g.mapTix.remove(baseTix);
				if (tixAli== null) 
					tixAli= new Vector<IsotigToGTF.Tig>();
				

				tix.plus= cox.plus= s[8].equals("+");
				if (tixAli.size()== 0|| tixAli.lastElement().plus== cox.plus) {	// tx happy
					baseTix.plus= cox.plus;
					tixAli.add(cox);
					coxAli.add(tix);
				} else {	// tx inconsistent
					if (coxAli.size()== 0) {	// cox not yet assigned to something
						baseCox.plus= false;	// ..turn it around!
						cox.plus= !cox.plus;
						tix.plus= !cox.plus;
						tixAli.add(cox);
						coxAli.add(tix);
					} else {	// cox is already occupied
						for (int i = 0; i < tixAli.size(); i++) {
							Vector<Tig> v= g.mapCox.get(tixAli.elementAt(i));
							if (v!= null&& v.size()> 0)
								System.err.println("Problem");
						}
					}
				}
				g.mapCox.put(baseCox, coxAli);
				g.mapTix.put(baseTix, tixAli);
				
			}
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	
	public static void main(String[] args) {
		
//		if (1== 1) {
//			closureTest();
//			System.exit(1);
//		}
		File layoutFile= new File("/Users/micha/Downloads/Files/454IsotigsLayout.txt");
		File contigFile= new File("/Users/micha/Downloads/Files/454Isotigs.txt");
		long t0= System.currentTimeMillis();
		HashMap<String, Group> map= parseLayout(layoutFile);
		parseIsotig(contigFile, map);
		System.err.println("Total "+ ((System.currentTimeMillis()- t0)/ 1000)+ " sec.");
	}
	
	static void closureTest() {
		
		int[] seqLen= new int[] {10,10,10};
		Closure c= AliGraphClosure.newAligGraphClosure(
				seqLen.length, 
				seqLen, 
				0, 
				null);
		AliGraphClosure.addAlignedPositions(c, 0, 3, 1, 3);
		//AliGraphClosure.computeClosure(c);
		AliGraphClosure.addAlignedPositions(c, 1, 5, 2, 5);
		AliGraphClosure.computeClosure(c);

//		Vector<Integer> vLen= new Vector<Integer>();
//		Vector<Boolean>[] vMap= new Vector[seqLen.length];
//		makeLayout(c, seqLen, vLen, vMap);
//		printLayout(vLen, vMap, System.out);
		
/*		boolean[][] a= joinPairwise(c, 0, 1, 10, 10, null);
		printLayout(a, System.out);
*/		
/*		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < 10; i++) {
				System.out.println("("+j+","+i+")");
				//AliGraphClosure.print_aligSets(c, j, i);
				//int[] a= AliGraphClosure.getTransitivityBounds(c, j, i, 1);
				//System.out.println(a[0]+","+a[1]);
				
			}
		}
*/		
	}
	
	static void printLayout(Tig[] groupTx, Vector<Integer> vLen, Vector<Boolean>[] vMap, PrintStream p) {
		StringBuilder[] sbs= new StringBuilder[vMap.length];
		for (int i = 0; i < sbs.length; i++) {
			sbs[i]= new StringBuilder(groupTx[i]+" ");
		}
		for (int i = 0; i < vLen.size(); i++) {
			String base= "["+ Integer.toString(i)+ "]";	// vLen.elementAt(i)
			String blank= "";
			for (int j = 0; j < base.length(); j++) 
				blank= blank+ " ";
			for (int j = 0; j < vMap.length; j++) {
				sbs[j].append(vMap[j].elementAt(i)?base:blank);
			}
		}
		for (int i = 0; i < sbs.length; i++) 
			p.println(sbs[i]);
		
	}
	
}
