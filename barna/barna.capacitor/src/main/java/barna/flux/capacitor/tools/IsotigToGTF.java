/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.flux.capacitor.tools;

import barna.flux.capacitor.closure.AliGraphClosure;
import barna.flux.capacitor.closure.Closure;
import barna.model.commons.IntVector;
import barna.model.commons.MyHashMap;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Vector;


public class IsotigToGTF {

	static final String BLOCK_PLUS= ">>>>>";
	static final String BLOCK_MINUS= "<<<<<";
	
	static class Tig implements Comparable<Tig> {
		int start= -1;
		int end= -1;
		String id= null;
		boolean plus= true;

		
		public Tig(Tig aTig) {
			start= aTig.start;
			end= aTig.end;
			id= aTig.id;
			plus= aTig.plus;
		}
		
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
		Closure c= null;
		int not= 0;
		
		Vector<Tig> iterated= null, deleted= null;
		
		public Group(String gid) {
			this.gid= gid;
			mapTix= new MyHashMap<IsotigToGTF.Tig, Vector<Tig>>();
			mapCox= new MyHashMap<IsotigToGTF.Tig, Vector<Tig>>();
		}
		
		public Group(Group aGroup) {
			this.gid= aGroup.gid;
			this.mapTix= clone(aGroup.mapTix);
			this.mapCox= clone(aGroup.mapCox);
			this.not= aGroup.not;
			if (aGroup.c!= null)
				this.c= new Closure(aGroup.c);
		}
		
		MyHashMap<IsotigToGTF.Tig, Vector<Tig>> clone(MyHashMap<IsotigToGTF.Tig, Vector<Tig>> map) {
			if (map== null)
				return null;
			MyHashMap<IsotigToGTF.Tig, Vector<Tig>> map2= new MyHashMap<IsotigToGTF.Tig, Vector<Tig>>(map.size());
			Iterator<Entry<Tig, Vector<Tig>>> iter= map.entrySet().iterator();
			while(iter.hasNext()) {
				Entry<Tig, Vector<Tig>> entry= iter.next();
				
				Tig key= null;
				if (entry.getKey()!= null)
					key= new Tig(entry.getKey());
				
				Vector<Tig> value= null;
				if (entry.getValue()!= null) {
					value= new Vector<Tig>(entry.getValue().size());
					for (int i = 0; i < entry.getValue().size(); i++) 
						value.add(new Tig(entry.getValue().elementAt(i)));
				}
				
				map2.put(key, value);
			}
			return map2;
		}
		@Override
		public boolean hasNext() {
			return (iterated!= null&& (deleted.size()+ iterated.size())< mapTix.size());
		}
		@Override
		public Tig next() {
			if (iterated== null|| (iterated.size()+ deleted.size())== mapTix.size())
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
							int p1= iterated.indexOf(tix.elementAt(k));
							int p2= deleted.indexOf(tix.elementAt(k));
							if (p1< 0&& p2< 0) 
								current= mapTix.getKey(tix.elementAt(k));	// baseTig
						}
					}
				}
				
				// disconnected components
				if (current== null) {
					Iterator<Tig> iter= mapTix.keySet().iterator();
					while(current== null&& iter.hasNext()) {
						Tig tig= iter.next();
						int p1= iterated.indexOf(tig);
						int p2= deleted.indexOf(tig);
						if (p1< 0&& p2< 0) 
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
				Tig lastTig= iterated.remove(iterated.size()- 1);
				deleted.add(lastTig);
			}
		}
		
		public void removeTig(Tig tig) {
			Vector<Tig> vTigs= mapTix.remove(tig);
			for (int i = 0; i < vTigs.size(); ++i) {
				Vector<Tig> vCog= mapCox.get(vTigs.elementAt(i));
				vCog.remove(tig);	// also remove from contig mapping
			}
		}
		
		@Override
		public Iterator<Tig> iterator() {
			iterated= new Vector<IsotigToGTF.Tig>(mapTix.size());
			deleted= new Vector<IsotigToGTF.Tig>();
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
	
	static void printLayout(Tig[] groupTx, Vector<Boolean>[] vMap, PrintStream p) {
		StringBuilder[] sbs= new StringBuilder[vMap.length];
		int maxLen= -1;
		for (int i = 0; i < groupTx.length; i++) { 
			if(groupTx[i]== null)
				continue;
			if (groupTx[i].toString().length()> maxLen)
				maxLen= groupTx[i].toString().length();
		}
		for (int i = 0; i < sbs.length; i++) {
			if(groupTx[i]== null)
				continue;
			sbs[i]= new StringBuilder(groupTx[i]+" ");
			for (int j = groupTx[i].toString().length(); j < maxLen; j++) {
				sbs[i].append(" ");
			}
		}
		for (int i = 0; i < vMap[0].size(); i++) {
			String base= "["+ Integer.toString(i)+ "]";	// vLen.elementAt(i)
			String blank= "";
			for (int j = 0; j < base.length(); j++) 
				blank= blank+ " ";
			for (int j = 0; j < vMap.length; j++) {
				if(groupTx[j]== null)
					continue;
				sbs[j].append(vMap[j].elementAt(i)?base:blank);
			}
		}
		for (int i = 0; i < sbs.length; i++) {
			if(groupTx[i]== null)
				continue;
			p.println(sbs[i]);
		}
		
	}

	static void makeLayout(Closure c, int[] seqLen, Vector<Boolean>[] vMap) {
			
			// trivial case, would deadlock
			if (seqLen.length== 1) {
				vMap[0]= new Vector<Boolean>();
				for (int i = 0; i < seqLen[0]; i++) 
					vMap[0].add(true);
				return;
			}
				
			// init
			int[] p= new int[seqLen.length];
			Arrays.fill(p, 1);
			for (int i = 0; i < p.length; i++) {	// mark unused
				if (seqLen[i]< 0)
					p[i]= -1;
			}
			for (int i = 0; i < vMap.length; i++) 
				vMap[i]= new Vector<Boolean>();

			// advance
			boolean left= true;
			while(left) {
				
				
				// find min blocks
				for (int i = 0; i < p.length; i++) {
					if (p[i]< 0|| p[i]> seqLen[i])
						continue;
					IntVector v= new IntVector();
					v.add(i);
					for (int j = 0; j < p.length; j++) {
						if (i== j|| p[j]< 0|| p[j]> seqLen[j])
							continue;
						int lb= AliGraphClosure.predFrontier(c, i, p[i], j);
						int ub= AliGraphClosure.succFrontier(c, i, p[i], j);
						
						if (lb== ub) { // matching
							if (lb== p[j]) {
								v.add(j);
							} else if (lb> p[j]){
								v= null;	// wait
								break;
								//hi[j]= Math.min(lb, hi[j]);	// exclude from potential range
							} else {	// match before, problem
								System.err.println("Problem, match before");
							}
						} else { // range, possibly empty (lb< ub)
							if (lb>= p[j]) {
								v= null;	// wait
								break;
							} else if (ub< p[j])	// align before, problem
								System.err.println("Problem, ali before");
						}
					}
					// advance
					if (v!= null) {
						int curLen= vMap[0].size();
						for (int j = 0; j < v.length; j++) {
							vMap[v.get(j)].add(true);
							++p[v.get(j)];
						}
						for (int j = 0; j < vMap.length; j++) {
							if (vMap[j].size()== curLen)	// p[j]> 0&& 
								vMap[j].add(false);
						}
						break;
					}
				}
				
				// stop condition
				left= false;
				for (int i = 0; i < p.length; i++) {
					if (p[i]> 0) {
//						if (vMap[i].lastElement())
//							++p[i];
						left|= (p[i]<= seqLen[i]);			// DEBUG !!! SHOULD READ <=
					}
				}
			}
		}

	/**
	 * Adds the isotigs according to their connections to a newly
	 * created group.
	 * @param g
	 */
	static Vector<Group> resolveGroup(Group g, File outputFile) {
		
		// init
		System.err.print(g.gid);
		Tig[] groupTx= new Tig[g.mapTix.size()];
		g.mapTix.keySet().toArray(groupTx);
		Arrays.sort(groupTx);
		int[] len= new int[groupTx.length];
		for (int i = 0; i < len.length; i++) 
			len[i]= g.mapTix.get(groupTx[i]).size();

		
		// build up groups
		Vector<Group> vGroups= new Vector<Group>();
		g.not= g.mapTix.size();
		for(int i= 0; g.not> 0; ++i) {
			//System.err.println(g.gid+"-"+(i+1)+", "+ g.mapTix.size()+ " transcripts left.");
			Group nuG= resolveGroup(g, groupTx, len);
			vGroups.add(nuG);
			if (g.not> 0) {
				for (int j = 0; j < g.iterated.size(); j++) {
					Tig tig= g.iterated.elementAt(j);
					g.removeTig(tig);
				}
			}
		}

		// output
		System.err.print(", split in "+vGroups.size()+ " subgroups: ");
		for (int i = 0; i < vGroups.size(); i++) {
			
			Group gx= vGroups.elementAt(i);
			System.err.print(gx.mapTix.size()+ " ");
			if (vGroups.size()> 1)
				gx.gid= gx.gid+ "-"+ Integer.toString(i+1);

			// truncate tx set NOO, closure based on complete
/*			Tig[] groupTx2= new Tig[gx.mapTix.size()];
			int[] len2= new int[groupTx2.length];
			for (int j = 0, p= 0; j < groupTx.length; j++) { 
				if (gx.mapTix.containsKey(groupTx[j])) {
					groupTx2[p]= groupTx[j];
					len2[p++]= len[j];
				}
			}
*/			
			// mask out tx not present in the group
			int[] len2= len.clone();
			Tig[] groupTx2= groupTx.clone();
			for (int j = 0; j < groupTx2.length; j++) {
				if (!gx.mapTix.containsKey(groupTx2[j])) {
					groupTx2[j]= null;
					len2[j]= -1;
				}
			}
			
			// layout
			Vector<Integer> vLen= new Vector<Integer>();
			Vector<Boolean>[] vMap= new Vector[len2.length];
			makeLayout(gx.c, len2, vMap);
			//System.err.println();
			//printLayout(groupTx2, vMap, System.err);
			
			Vector<Tig> v= getPlayground(gx, groupTx, len2, vMap);
			writeGTF(outputFile, gx, groupTx, len, vMap, v);
			
			PrintStream p;
			try {
				p = new PrintStream(new FileOutputStream(outputFile+".txt", true));
				p.println("=== Layout "+ gx.gid+ " === "+ v.size()+ " contigs, "+ gx.mapTix.size()+ " isotigs");
				for (int j = 0; j < v.size(); j++) 
					p.println("["+(j+1)+"]\t"+v.elementAt(j).id+"\t"+v.elementAt(j).length());
				p.println();
				printLayout(groupTx2, vMap, p);
				p.println();
				p.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			
		}
		
		return vGroups;
	}
				
			static Group resolveGroup(Group g, Tig[] groupTx, int[] len) {
				
				// prepare new group to add incrementally the isotigs to
				Group newGroup= new Group(g.gid);
				newGroup.c= AliGraphClosure.newAligGraphClosure(len.length, len, 0, null);
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
					
					// copy
					Group cGroup= new Group(newGroup);
					
					boolean added= addTx(tig, g, newGroup, groupTx);
					if (added) {
						--g.not;
						//g.removeTig(tig);	// substract from main group
					
					// revert
					} else {
						g.remove();	// remove only from iterator
						newGroup= cGroup;
						// purge tig from maps of this group
						Vector<Tig> cox= newGroup.mapTix.remove(tig);
						for (int i = 0; cox!= null&& i < cox.size(); i++) 
							newGroup.mapCox.get(cox.elementAt(i)).remove(tig);
					}
					
				}
				
				// TODO turn around locus, s.t. there are more sense isotigs
				if (cntAsense> cntSense) {
					;
				}		
				
				// change names (+/-)
				for (int i = 0; i < groupTx.length; i++) {
					if (newGroup.mapTix.containsKey(groupTx[i]))
						groupTx[i]= newGroup.mapTix.getKey(groupTx[i]);
				}
				
				return newGroup;
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
		 * @return whether the transcript was added
		 */
		private static boolean addTx(Tig tig, Group og, Group g, Tig[] groupTx) {
			
			//System.err.println("adding "+ tig);

			// get segments
			Vector<Tig> vCox= og.mapTix.get(tig);
			Vector<Tig> vTix= new Vector<Tig>(vCox.size());
			for (int i = 0; i < vCox.size(); i++) {
				Vector<Tig> v= og.mapCox.get(vCox.elementAt(i));
				Tig t= v.get(v.indexOf(tig));
				vTix.add(t);
			}
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
			
			// align transcript segments
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
						if (AliGraphClosure.alignablePositions(g.c, idx1, h+1, idx2, p+1)) {
							AliGraphClosure.addAlignedPositions(g.c, idx1, h+1, idx2, p+1);
							AliGraphClosure.computeClosure(g.c);
						} else {
							// System.err.println("\nInconsistent alignment "+ tix);
							return false;
						}
					} else {
						// should actually never incur now, by swapping before
						//System.err.println("\nInconsistent reverse alignment !!! "+ tix);
						return false;
    //					for (int j = tix.start, k= otherTx.end; j <= tix.end&& k>= otherTx.start; ++j, --k) {
	//						if (!AliGraphClosure.alignablePositions(c, idx1, j, idx2, k)) {
	//							int lb= AliGraphClosure.getLowerBound(c, idx1, j, idx2);
	//							int ub= AliGraphClosure.getUpperBound(c, idx1, j, idx2);
	//							System.err.println("\nInconsistent reverse alignment "+ tix);
	//							break;
	//						}
	//						AliGraphClosure.addAlignedPositions(c, idx1, j, idx2, k);
	//					}
					}
				}
				
				
				coxAli.add(tix);
				tixAli.add(cox);
				g.mapCox.put(baseCox, coxAli);
				g.mapTix.put(baseTix, tixAli);
			}
			

			g.mapTix.getKey(tig).plus= sense;
			// TODO update counter
			
			return true;			
		}

	private static void writeGTF(File f, Group newGroup, Tig[] groupTx, int[] len,
			Vector<Boolean>[] vMap, Vector<Tig> pground) {
		
		final String tab= "\t", exon= "exon", plus= "+", minus= "-", 
				dot= ".", space= " ", txID= "transcriptID", quot= "\"", 
				semic= ";", nl= "\n";
		
		try {
			BufferedWriter writer= new BufferedWriter(new FileWriter(f, true));
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
			if (len[i]> 0&& p[i]!= len[i])
				System.err.println("Incomplete playground "+ i+ ": "+ p[i]+ "<>"+ len[i]);
		}
		
		return pground;
	}

	static void parseIsotig(File isotigFile, File outputFile, HashMap<String, Group> map) {
		
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
						
						resolveGroup(lastG, outputFile);
						System.err.println("isotigs, "+ ((System.currentTimeMillis()- t)/ 1000)+ " sec.");
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
			
			resolveGroup(lastG, outputFile);
			System.err.println("isotigs, "+ (System.currentTimeMillis()- t)/ 1000+ " sec.");
							
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
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

	private static int findIdx(Tig tix, Tig[] groupTx) {
		int p= Arrays.binarySearch(groupTx, tix);
		if (p< 0)
			System.err.println("Not found "+ tix);
		return p;
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

	public static void main(String[] args) {
		
		// 454IsotigsLayout.txt
		// 454Isotigs.txt
		// Isogroup0038Layout.txt
		// Isogroup0038.txt
		File layoutFile= new File("/Users/micha/projects/pedro/Files/454IsotigsLayout.txt");
		File contigFile= new File("/Users/micha/projects/pedro/Files/454Isotigs.txt");
		File outputFile= new File("/Users/micha/projects/pedro/Files/isotest.gtf");
		if (outputFile.exists())
			outputFile.delete();
		File layOut= new File(outputFile.getAbsolutePath()+".txt");
		if (layOut.exists())
			layOut.delete();
		
		long t0= System.currentTimeMillis();
		HashMap<String, Group> map= parseLayout(layoutFile);
		parseIsotig(contigFile, outputFile, map);
		System.err.println("Total "+ ((System.currentTimeMillis()- t0)/ 1000)+ " sec.");
	}
	
}
