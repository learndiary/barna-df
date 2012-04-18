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

package barna.flux.capacitor.reconstruction;

import java.io.*;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

public class LParser {

	public static class ChainComparator implements Comparator<int[]> {
		public int compare(int[] a1, int[] a2) {
			int min= Math.min(a1.length, a2.length);
			for (int i = 0; i < min; i++) {
				if (a1[i]< a2[i])
					return -1;
				if (a1[i]> a2[i])
					return 1;
			}
			
			if (a1.length< a2.length)
				return -1;
			if (a1.length> a2.length)
				return 1;
			return 0;
		}
	}
	
	public static final ChainComparator DEFAULT_CCOMPARATOR= new ChainComparator();
	
	public static int[] chain2int(String chain) {
		String[] ss= chain.split("\\p{Punct}");
		int[] a= new int[ss.length];
		for (int i = 0; i < a.length; i++) {
			a[i]= Integer.parseInt(ss[i]);
		}
		return a;
	}
	
	public class LP {
		public int[] coords;
		public int[][] tstruct;
		public char[] types;
		public double[] vals; // [0]=OF, C1, ...
		public double[][] m;	// [0]=rhs, C1, ...
		
		public HashMap<Object, int[]> mapConstraints;
		String[] tids= null;
		int[][] chains= null;
		
		public String[] getTranscriptIDs() {
			
			if (tids == null) {
				Iterator iter= mapConstraints.keySet().iterator();
				Vector<String> v= new Vector<String>();
				while(iter.hasNext()) {
					Object o= iter.next();
					if (o instanceof String)
						v.add((String) o);
				}
				
				// sort according to constraint nb
//				for (int i = ss.length- 1; i < ss.length; i++) 
//					ss[i]= v.elementAt(i);
				int last= -1, idx= -1;
				String[] ss= new String[v.size()];
				for (int i = 0; i < ss.length; i++) {
					for (int j = 0; j < ss.length; j++) {
						if (v.elementAt(j)== null)
							continue;
						int x= mapConstraints.get(v.elementAt(j))[0];
						if (last== -1|| x< last) {
							last= x;
							idx= j;
						}
					}
					
					ss[i]= v.elementAt(idx);
					v.set(idx, null);
					last= -1; idx= -1; 
				}
				
				tids= ss;
			}

			return tids;
		}
		
		public int[][] getChains() {
			
			if (chains == null) {
				Iterator iter= mapConstraints.keySet().iterator();
				Vector<int[]> v= new Vector<int[]>();
				while(iter.hasNext()) {
					Object o= iter.next();
					if (o instanceof int[])
						v.add((int[]) o);
				}

				int[][] m= new int[v.size()][];
				for (int i = 0; i < m.length; i++) 
					m[i]= v.elementAt(i);
				
				Arrays.sort(m, DEFAULT_CCOMPARATOR);
				chains= m;
			}

			return chains;
		}
		
		public double getCorrection(int[] chain) {
			int[] cc= mapConstraints.get(chain);
			double res;
			if (vals[cc[0]]!= 0) {
				res= -vals[cc[0]];
				assert(vals[cc[1]]== 0);
			} else {
				res= vals[cc[1]];
				assert(vals[cc[0]]== 0);
			}
			return res;
		}
		
		public double[] getRestriction(int[] chain) {
			int[] cc= mapConstraints.get(chain);
			for (int i = 0; i < m.length; i++) {
				double[] a= m[i];
				if (m[i][cc[0]]!= 0&& m[i][cc[1]]!= 0)
					return m[i];
			}
			
			return null;
		}		
		
		public char getSiteSymbol(int coord) {
			int p= Arrays.binarySearch(coords, coord);
			return types[p];
		}
		
		public double getTranscriptExp(String id) {
			int x= getTranscriptNr(id);
			int c= vals.length- getTranscriptIDs().length+ x;
			return vals[c];
		}
		
		public int getTranscriptNr(String id) {
			for (int i = 0; i < tids.length; i++) {
				if (tids[i].equals(id))
					return i;
			}
			return -1;
		}
		
		
		public boolean transcriptPresent(String tid, int[] area) {
			int nr= getTranscriptNr(tid);
			for (int i = 0; i < getChains().length; i++) {
				double[] a= getRestriction(getChains()[i]);
				for (int j = 0; a!= null&& j < getChains()[i].length- 1; j++) {
					if (getChains()[i][j]>= area[1])
						break;
					if (getChains()[i][j]== area[0]&& getChains()[i][j+ 1]== area[1]&&
							a[a.length- getTranscriptIDs().length+ nr]!= 0) {
						return true;
					}
				}
			}
			return false;
		}
		
		public String toString() {
			StringBuffer sb= new StringBuffer();
			sb.append("Matrix:\n");
			for (int i = 0; i < m.length; i++) {
				for (int j = 1; j < m[i].length; j++) 
					sb.append(Double.toString(m[i][j])+ "\t");
				sb.append(m[i][0]+ "\n");
			}
			sb.append("\n");
			
			sb.append("Values:\n");
			sb.append("OF\t"+vals[0]+"\n");
			for (int i = 1; i < vals.length; i++) 
				sb.append("C"+ i+ "\t"+ vals[i]+ "\n");
			sb.append("\n");

			sb.append("Key:\n");
			Iterator iter= mapConstraints.keySet().iterator();
			while (iter.hasNext()) {
				Object o= iter.next();
				StringBuffer sb2= new StringBuffer();
				String id= null;
				if (o instanceof int[]) {
					int[] a= (int[]) o;
					for (int i = 0; i < a.length; i++) 
						sb2.append(a[i]+ " ");
					id= sb2.toString();
				} else
					id= o.toString();
				
				int[] cc= mapConstraints.get(o);
				sb.append(id+ "\t");
				for (int i = 0; i < cc.length; i++) 
					sb.append("C"+ cc[i]+ " ");
				sb.append("\n");
			}
			sb.append("\n");

			sb.append("Coordinates:\n");
			for (int i = 0; i < coords.length; i++) 
				sb.append(coords[i]+ "\t"+ types[i]+ "\n");
			sb.append("\n");
			
			return sb.toString();
		}
		
		
	}
	
	public static void main(String[] args) {
		
		File myZip= null;		
		try {
			myZip= new File("c:\\workspace\\Genome\\resources\\formats\\CME_W1_CI.8-PE_lp_mod.zip");
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		LParser myParser= new LParser(myZip);
		String name= myParser.getNames()[0];
		System.err.println("Reading "+ name);
		LP myLP= myParser.parseEntry(name);
		
		System.err.println(myLP);
	}
	
	File baseFile= null;
	public LParser(File base) {
		this.baseFile= base;
	}
	
	public String[] getNames() {
		
		if (baseFile.isDirectory())
			return baseFile.list();
		
		// else
		ZipFile zf= null;
		try {
			zf = new ZipFile(baseFile);
		} catch (ZipException e) {
			e.printStackTrace();
			return null;
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
		Enumeration entries = zf.entries();		
		String[] names= new String[zf.size()];
		for(int i= 0; i< names.length; ++i) {
			names[i]= ((ZipEntry) entries.nextElement()).getName();
		}
		return names;
	}
	
	public LP parseEntry(String name) {
		
		LP lp= new LP();
		try {
			InputStream ins= null;
			long size= 0;
			ZipFile zf= null;
			if (baseFile.isDirectory()) {
				File inFile= new File(baseFile+ File.separator+ name);
				size= inFile.length();
				ins= new FileInputStream(inFile);
			} else {
				try {
					zf = new ZipFile(baseFile);
				} catch (ZipException e) {
					e.printStackTrace();
					return null;
				} catch (IOException e) {
					e.printStackTrace();
					return null;
				}
				ZipEntry ze = zf.getEntry(name);
				size = ze.getSize();
				ins= zf.getInputStream(ze);
			}
			//System.out.println("Read "+ ze.getName());
			if (size > 0) {
				//System.out.println("Length is " + size);
				BufferedReader br = new BufferedReader(
						new InputStreamReader(ins));
				String line;
				while ((line = br.readLine()) != null&& (!line.startsWith("Model size:")));
				String[] ss= line.split("\\s+");
				int nbR= Integer.parseInt(ss[2]), nbC= Integer.parseInt(ss[4]);
				lp.m= new double[nbR][nbC+1];
				
				while ((line = br.readLine()) != null&& (!line.startsWith("Model name:")));
				br.readLine(); // C1       C2       C3       C4 ...
				br.readLine(); // Minimize         1        1        1
					// 0 0 1 -1 ... 0.0177016  =        0
				for (int i= 0; !(line= br.readLine()).startsWith("Type"); ++i) {
					ss= line.split("\\s+");
					lp.m[i][0]= Double.parseDouble(ss[ss.length- 1]);
					for (int j = 1; j < ss.length- 2; j++)  	// = 0
						lp.m[i][j]= Double.parseDouble(ss[j]);
				}
				
				// upbo.., lowbo..
				while ((line = br.readLine()) != null&& (!line.startsWith("Value of objective function:")));
				double vaLP= Double.parseDouble(line.split("\\s")[4]); // Value of objective function: 769.52828154 
				
				
				while (!(line = br.readLine()).startsWith("Actual values of"));	// Actual values of the variables:
				lp.vals= new double[nbC+ 1];
				lp.vals[0]= vaLP;
				for (int i = 1; i < lp.vals.length; i++) {
					line= br.readLine();
					ss= line.split("\\s+");
					lp.vals[i]= Double.parseDouble(ss[1]);
				}
				
				while (!(line = br.readLine()).startsWith("Settings:"));
				br.readLine();	// paired-end	true 
				br.readLine(); 	// lower bound 
				br.readLine(); 	// costfunc	linear 
				
				HashMap<Integer, Character> coordMap= new HashMap<Integer, Character>();
				String[] pp;
				lp.mapConstraints= new HashMap<Object, int[]>();
				while ((line = br.readLine()) != null) {
					ss= line.split("\\s+");
					if (ss.length< 2)
						break;
					ss[0]= (ss[0].endsWith("PE:")?ss[0].substring(0, ss[0].length()- 3):
						ss[0].substring(0, ss[0].length()- 1));
					int[] cc= new int[ss.length- 1];
					for (int i = 0; i < cc.length; i++) 
						cc[i]= Integer.parseInt(ss[i+1].substring(1));					
					int x= 0;
					for (; x < ss[0].length(); ++x) 
						if (Character.isLetter(ss[0].charAt(x)))
							break;
					if (x == ss[0].length()) {
						pp= ss[0].split("\\p{Punct}");
						for (int i = 0; i < pp.length; i++) {
//							if (!pp[i].contains("PE")) {
								char symbol= ss[0].charAt(ss[0].indexOf(pp[i])+ pp[i].length());
								coordMap.put(new Integer(pp[i]), symbol);
//							}
						}
						lp.mapConstraints.put(chain2int(ss[0]), cc);
					} else
						lp.mapConstraints.put(ss[0], cc);

				}

				while (!(line = br.readLine()).startsWith("Transcripts:"));
				lp.tstruct= new int[lp.getTranscriptIDs().length][];
				while ((line = br.readLine()) != null) {
					ss= line.split("\\s+");
					if (ss.length< 2)
						break;
					int idx= lp.getTranscriptNr(ss[0]);
					
					lp.tstruct[idx]= chain2int(ss[1]);
				}				
				
				br.close();
				if (zf!= null)
					zf.close();

				lp.coords= new int[coordMap.size()];
				lp.types= new char[lp.coords.length];
				Iterator<Integer> iter= coordMap.keySet().iterator();
				for (int i = 0; i < lp.coords.length; i++) 
					lp.coords[i]= iter.next();
				Arrays.sort(lp.coords);
				
				for (int i = 0; i < lp.types.length; i++) 
					lp.types[i]= coordMap.get(lp.coords[i]).charValue();
				
	      }
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return lp;
	}
	
}
