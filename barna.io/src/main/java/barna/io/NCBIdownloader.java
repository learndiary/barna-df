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

package barna.io;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Arrays;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class NCBIdownloader {

	class LeechThread2 extends Thread {
		Vector<String> v;
		public LeechThread2(Vector<String> vv) {
			this.v= vv;
		}
		
		@Override
		public void run() {
			
			++alive;
			
			for (int i = 0; i < v.size(); i++) {
				get(v.elementAt(i));
			}
			
			--alive;
		}
		
	}
	static class LeechThread extends Thread {
		int lo, hi;
		String pfx;
		char[] x;
		public LeechThread(String pfx, int lo, int hi, int width) {
			this.pfx= pfx;
			this.lo= lo;
			this.hi= hi;
			x= new char[width];
			Arrays.fill(x, '0');
		}
		
		@Override
		public void run() {
			
			++alive;
			
			for (int i = lo; i < hi; i++) {
				char[] a= Integer.toString(i).toCharArray();
				System.arraycopy(a, 0, x, x.length- a.length, a.length);
				getEST(pfx+ new String(x));
	
			}
			
			--alive;
		}
		
	}
	static int alive= 0;
	
	String database= null;
	File inFile= null;	
	int maxThread= 10;
	
	static boolean parse(NCBIdownloader dl, String[] args) {
		
		if (args== null|| args.length< 2) {
			System.err.println("usage: -jar leechNCBI.jar inputFile database [maxThread]");
			System.err.println("\tinputFile\tfile with each ID in one line");
			System.err.println("\tdatabase\t\'nucest\' (for EST), " +
					"\n\t\t\t\'nucleotide\' (for something else)");
			System.err.println("\t[maxThread]\tan (integer) number of simultaneous "
					+"\n\t\t\tconnections (default: 10)");
			System.err.println();
			System.err.println("by micha (2009), code adapted from the NCBI sissi tools");
			System.err.println("license: do whatever you want with my program, I am not\nresposible for anything.");
			System.err.println();
			System.err.println("Description: program uses the CGI interfaces\n");
			System.err.println("http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi");
			System.err.println("\nand\n");
			System.err.println("http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi");
			System.err.println("\nto get the fasta sequences with the IDs from \'inputFile\'.");
			System.err.println("IDs that generated an error during retrieval are written to");
			System.err.println("stderr, in case you try too aggressive connection settings");
			System.err.println("just catch them in a file and re-run the program. The");
			System.err.println("fetched fasta sequences come over stdout.");
			System.err.println();
			System.err.println("Example: java -jar leechNCBI.jar arabidopsis.txt nucest 30");
			System.err.println();
			System.err.println("cheers.");
			return false;
		}
		
		dl.inFile= new File(args[0]);
		if (!dl.inFile.exists()) { 
			System.err.println("This is not a file "+args[0]);
			return false;
		}
		
		dl.database= args[1];
		if (args.length> 2) {
			try{
				dl.maxThread= Integer.parseInt(args[2]);
			} catch (Exception e) {
				System.err.println("Ey man, an integer number of threads, not "+args[2]);
				return false;
			}
		}
		
		return true;
	}
	
	public static void main(String[] args) {
		//getArabidopsis();
		//getArabidopsis2();
		//getArabidopsis(new File("arabidopsis_EH_200K_missing.fasta"));	// arabidopsis_EH_200K_missing.fasta
		NCBIdownloader leecher= new NCBIdownloader();
		if (parse(leecher, args))
			leecher.run();
	}
	
	static String fName= "P:\\interspecies2\\weber\\arabidopsis_batch2_EL_1-341852.fasta";
	public static void getArabidopsis2() {
		int simu= 100;
		int lo= 795234, up= 995234;	// excl
		int delta= ((up- lo)/ simu);
		if ((up- lo)% simu!= 0)
			++delta;
		
		delta= 30;
		
		try {
			File f= new File(fName);
			if (f.exists())
				f.delete();
			
			for (int i = lo; i < up; i+= delta) {
				while(alive> delta) {
					try {
						Thread.currentThread().sleep(100);
					} catch (Exception e) {
						; // :)
					}
				}
				Thread t= new LeechThread("EH", i, Math.min(i+ delta, up), 6);
				t.start();
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	static BufferedWriter writer;
	public static void getArabidopsis() {
		try {
			writer= new BufferedWriter(new FileWriter("arabidopsis.fasta"));
			// EST sequence accession numbers in GenBank are 
			// EH795234 through EH995233 and 
			// EL000001 through EL341852.
			
			for (int i = 795234; i <= 995233; i++) {
				getEST("EH"+i/*, writer*/);
				writer.flush();
			}
			char[] x= new char[6];
			Arrays.fill(x, '0');
			for (int i = 1; i <= 341852; i++) {
				char[] a= Integer.toString(i).toCharArray();
				System.arraycopy(a, 0, x, x.length- a.length, a.length);
				getEST("EL"+ new String(x)/*, writer*/);
				writer.flush();
			}

			writer.flush();
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	static long wroteBytes= 0, lastTime= 0;
	public static synchronized void writeStream(InputStream stream) {
		++calls;
		if (calls% 10== 0) {
			System.out.println(Float.toString((float) (wroteBytes/ (double) (System.currentTimeMillis()- lastTime)))
					+ " Kbps, "+alive+" alives.");
			wroteBytes= 0;
			lastTime= System.currentTimeMillis();
		}
		int bWrite= 0;
		boolean tag= false;
		try {
			writer= new BufferedWriter(new FileWriter(fName, true));
			int c;
			while ((c= stream.read())!= -1) {				
				writer.write(c);
				if (c== '>') 
					tag= true;
				++bWrite;
			}
			writer.flush();
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		if (bWrite== 0|| (!tag))
			System.currentTimeMillis();
		wroteBytes+= bWrite;
	}
	
	static int calls= 0;
	public static void getEST(String ID) {
		// System.out.println(ID);
		// http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucest&retmax=1&usehistory=y&term=EH795234
		try {
			int count= -1, qkey= -1;
			String env= null;
			while (count< 0|| qkey< 0|| env== null) {
				URL query = new URL("http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" +
						"db=nucest&retmax=1&usehistory=y&term="+ID);
				HttpURLConnection con= (HttpURLConnection) query.openConnection();
				
				InputStream in= null;
				while(in== null)
					try {
						in= con.getInputStream();
					} catch (Exception e) {
						in= null;
						try {
							Thread.currentThread().sleep(500);
						} catch (Exception e2) {
							; //:)
						}
					}
				
				int c;
				StringBuilder sb= new StringBuilder();
				while ((c= in.read())!= -1) 
					sb.append((char) c);
				con.disconnect();
				//System.out.println(sb);
				try {
					Pattern p= Pattern.compile("<Count>(\\d+)</Count>.*<QueryKey>(\\d+)</QueryKey>.*<WebEnv>(\\S+)</WebEnv>");	//.*<QueryKey>(\\d+)</QueryKey>.*<WebEnv>(\\S+)</WebEnv>
					Matcher m= p.matcher(sb);
					m.find();

					count= Integer.parseInt(m.group(1));
					qkey= Integer.parseInt(m.group(2));
					env= m.group(3);
				} catch (Exception e) {
					; // :)
				}
			}
			
			if (count== 0)
				System.currentTimeMillis();
			
			int retmax= 3;
			for(int retstart = 0; retstart < count; retstart += retmax) {
				URL query= new URL("http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
						+ "rettype=fasta&retmode=text&retstart="+retstart+"&retmax="+retmax
						+ "&db=nucest&query_key="+qkey+"&WebEnv="+env);

				HttpURLConnection con= (HttpURLConnection) query.openConnection();
				InputStream in= null;
				while (in== null) {
					try {
						in= con.getInputStream();
					} catch (Exception e) {
						in= null;
						try {
							Thread.currentThread().sleep(500);
						} catch (Exception e2) {
							; //:)
						}
					}
				}
				
				writeStream(in);
				//sb= new StringBuilder();
//				while ((c= in.read())!= -1) {
//					//sb.append((char) c);
//					writer.write(c);
//				}
				con.disconnect();
			}
			//System.out.println(sb.toString());
			
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	public void run() {
		
		int ctr= 0;
		try {
			File f= new File(fName);
			if (f.exists())
				f.delete();
			
			BufferedReader buffy= new BufferedReader(new FileReader(inFile));
			int lines= 0;
			for (String s; (s= buffy.readLine())!= null; ++lines);
			buffy.close(); // =?
			int batchSize= lines/ maxThread;
			if (lines% maxThread!= 0)
				++batchSize;
			buffy= new BufferedReader(new FileReader(inFile));
			
			while (true) {				
				Vector<String> v= new Vector<String>();
				String s= null;
				for (int i = 0; i < batchSize&& (s= buffy.readLine())!= null; i++) { 
					v.add(s);
					++ctr;
				}
				
				while(alive> maxThread) {
					try {
						Thread.currentThread().sleep(100);
					} catch (Exception e) {
						; // :)
					}
				}
				
				Thread t= new LeechThread2(v);
				t.start();
				if (s== null)
					break;
			}
			
			buffy.close();
			System.out.flush();
			System.err.flush();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	
	void get(String ID) {
			// System.out.println(ID);
			// http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucest&retmax=1&usehistory=y&term=EH795234
			try {
				int count= -1, qkey= -1;
				String env= null;
				while (count< 0|| qkey< 0|| env== null) {
					URL query = new URL("http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" +
							"db="+database+"&retmax=1&usehistory=y&term="+ID);
					HttpURLConnection con= (HttpURLConnection) query.openConnection();
					
					InputStream in= null;
					while(in== null)
						try {
							in= con.getInputStream();
						} catch (Exception e) {
							in= null;
							try {
								Thread.currentThread().sleep(500);
							} catch (Exception e2) {
								; //:)
							}
						}
					
					int c;
					StringBuilder sb= new StringBuilder();
					while ((c= in.read())!= -1) 
						sb.append((char) c);
					con.disconnect();
					//System.out.println(sb);
					try {
						Pattern p= Pattern.compile("<Count>(\\d+)</Count>.*<QueryKey>(\\d+)</QueryKey>.*<WebEnv>(\\S+)</WebEnv>");	//.*<QueryKey>(\\d+)</QueryKey>.*<WebEnv>(\\S+)</WebEnv>
						Matcher m= p.matcher(sb);
						m.find();
	
						count= Integer.parseInt(m.group(1));
						qkey= Integer.parseInt(m.group(2));
						env= m.group(3);
					} catch (Exception e) {
						; // :)
					}
				}
				
				if (count== 0) {
					System.err.println(ID);
				}
				
				int retmax= 3;
				for(int retstart = 0; retstart < count; retstart += retmax) {
					URL query= new URL("http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
							+ "rettype=fasta&retmode=text&retstart="+retstart+"&retmax="+retmax
							+ "&db=nucest&query_key="+qkey+"&WebEnv="+env);
	
					HttpURLConnection con= (HttpURLConnection) query.openConnection();
					InputStream in= null;
					while (in== null) {
						try {
							in= con.getInputStream();
						} catch (Exception e) {
							in= null;
							try {
								Thread.currentThread().sleep(500);
							} catch (Exception e2) {
								; //:)
							}
						}
					}
					
					if (!write(in))
						System.err.println(ID);
					//sb= new StringBuilder();
	//				while ((c= in.read())!= -1) {
	//					//sb.append((char) c);
	//					writer.write(c);
	//				}
					con.disconnect();
				}
				//System.out.println(sb.toString());
				
			} catch (Exception e) {
				e.printStackTrace();
			}		
		}
	
	/**
	 * let the right one in:
	 * only one can enter to not mix output
	 * @param stream
	 * @return
	 */
	synchronized boolean write(InputStream stream) {
//		++calls;
//		if (calls% 10== 0) {
//			System.out.println(Float.toString((float) (wroteBytes/ (double) (System.currentTimeMillis()- lastTime)))
//					+ " Kbps, "+alive+" alives.");
//			wroteBytes= 0;
//			lastTime= System.currentTimeMillis();
//		}
		int bWrite= 0;
		boolean tag= false;
		try {
			if (fName!= null)
				writer= new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fName, true)));	//
			else
				writer= new BufferedWriter(new OutputStreamWriter(System.out));
			int c;
			while ((c= stream.read())!= -1) {				
				writer.write(c);
				if (c== '>') 
					tag= true;
				++bWrite;
			}
			writer.flush();
			if (fName!= null)
				writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	
		if (bWrite== 0|| (!tag))
			return false;
		wroteBytes+= bWrite;
		return true;
	}
	public static void getArabidopsis(File fx) {
		int simu= 100;
		int lo= 795234, up= 995234;	// excl
		int delta= ((up- lo)/ simu);
		if ((up- lo)% simu!= 0)
			++delta;
		
		delta= 30;
		
		int ctr= 0;
		try {
			File f= new File(fName);
			if (f.exists())
				f.delete();
			
			BufferedReader buffy= new BufferedReader(new FileReader(fx));
			while (true) {				
				Vector<String> v= new Vector<String>();
				String s= null;
				for (int i = 0; i < 30&& (s= buffy.readLine())!= null; i++) { 
					v.add(s);
					++ctr;
				}
				
				while(alive> delta) {
					try {
						Thread.currentThread().sleep(100);
					} catch (Exception e) {
						; // :)
					}
				}
				
//				Thread t= new LeechThread2(v);
//				t.start();
				if (s== null)
					break;
			}
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.out.println("demanded "+ctr);
	}
}
