package fbi.genome.sequencing.rnaseq.simulation;

import fbi.genome.model.constants.Constants;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Hashtable;
import java.util.Vector;

public class ThreadedFragFileReader extends Thread {
	
	private static int nr= 0;
	
	BufferedReader buffy;
	int lineLimit= 100000;
	Hashtable<String,int[][]> mapFrags;
	boolean finish= false;
	
	public ThreadedFragFileReader(Reader ffreader) {
		super();
		setName(this.getClass().getName()+"-"+(nr++));
		buffy= new BufferedReader(ffreader);
	}

	public void read() {
		Vector<int[]> v= new Vector<int[]>();
		String s, token[], currentID= null;
		finish= false;
		if (mapFrags== null)
			mapFrags= new Hashtable<String,int[][]>(lineLimit);
		try {
			for (int ctr= 1; (s= buffy.readLine())!= null; ++ctr) {
				
				token= s.split(Fragmenter.FRG_FILE_TAB);
				
				if (!token[2].equals(currentID)) {
					if (v.size()> 0) {
						int[][] m= new int[v.size()][];
						for (int i = 0; i < m.length; i++) 
							m[i]= v.elementAt(i);
						v.removeAllElements();
						mapFrags.put(currentID, m);
					}
					if (finish)
						break;
					currentID= token[2];
				}
				
				if ((!finish)&& ctr>= lineLimit) 
					finish= true;
				
				int[] p= new int[2];
				p[0]= Integer.parseInt(token[0]);
				p[1]= Integer.parseInt(token[1]);
				v.add(p);
			}
			if (s== null) 
				mapFrags= null;
			
		} catch (Exception e) {
			if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
				e.printStackTrace();			
		}
		
	}
	
	public void close() {
		try {
			buffy.close();
		} catch (IOException e) {
			if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
				e.printStackTrace();
		}
	}
	
	@Override
	public void run() {
		String[] token;
		Vector<int[]> v= new Vector<int[]>();
		while (true) {
			
			if (!finish) {
				String s= null;
				mapFrags= new Hashtable<String,int[][]>(lineLimit);
				
				try {
					int ctr= 1;
					String currentID= null;
					for (; (s= buffy.readLine())!= null; ++ctr) {
						
						token= s.split(Fragmenter.FRG_FILE_TAB);
						
						if (!token[2].equals(currentID)) {
							if (v.size()> 0) {
								int[][] m= new int[v.size()][];
								for (int i = 0; i < m.length; i++) 
									m[i]= v.elementAt(i);
								v.removeAllElements();
								mapFrags.put(currentID, m);
							}
							if (finish)
								break;
							currentID= token[2];
						}
						
						if ((!finish)&& ctr>= lineLimit) 
							finish= true;
						
						int[] p= new int[2];
						p[0]= Integer.parseInt(token[0]);
						p[1]= Integer.parseInt(token[1]);
						v.add(p);
					}
					
									
				} catch (IOException e) {
					if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
						e.printStackTrace();
				}
				
				if (s== null) {
					mapFrags= null;
					break;
				}
			}
			
			try {
				sleep(5000);
			} catch (InterruptedException e) {
				; // :)
			}
			
			
		}	// end while(true)
		
		try {
			buffy.close();
		} catch (IOException e) {
			if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
				e.printStackTrace();
		}
		
	}

	public synchronized void setFinish(boolean finish) {
		this.finish = finish;
	}

	public Hashtable<String, int[][]> getMapFrags() {
		return mapFrags;
	}
}
