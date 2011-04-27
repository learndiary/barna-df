package fbi.genome.io.bed;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.Progressable;
import fbi.commons.io.DevNullOutputStream;
import fbi.commons.thread.SyncIOHandler2;
import fbi.commons.tools.ArrayUtils;
import fbi.commons.tools.Interceptable;
import fbi.commons.tools.Sorter;
import fbi.genome.io.BufferedBACSReader;
import fbi.genome.io.DefaultIOWrapper;
import fbi.genome.io.ThreadedBufferedByteArrayStream;
import fbi.genome.io.rna.FMRD;
import fbi.genome.io.rna.ReadDescriptor;
import fbi.genome.io.rna.SolexaPairedEndDescriptor;
import fbi.genome.io.rna.UniversalReadDescriptor;
import fbi.genome.model.bed.BEDobject;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.commons.MyFile;
import fbi.genome.model.constants.Constants;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;
import java.util.concurrent.Future;

public class BEDwrapper extends DefaultIOWrapper {

	static void test() {
		System.out.println(((byte) -1)| (byte) 1);
		System.out.println(((byte) 2)& ((byte) -1));
		System.out.println(((byte) 2)& ((byte) 1));
		System.out.println(2&0);
		System.out.println(-1&Integer.MAX_VALUE);
	}
	
	BEDobject[] beds= null;
	
	File file;
	public BEDwrapper(String newFilePath) {
		super(newFilePath);
		file = new File(newFilePath);
		size = file.length();
		fileSep= guessFileSep();
	}
	
	HashSet<String> refIDset;
	
	private void addRefID(String ID) {
		if (refIDset== null)
			refIDset= new HashSet<String>();
		refIDset.add(ID);
	}
	
	private void clearRefIDs() {
		refIDset= null;
	}
	
	public boolean isApplicable() {
		try {
            Log.progressStart("checking");
			File f= new File(fPath+File.separator+fName);
			BufferedReader buffy= new BufferedReader(new FileReader(f), 10* 1024* 1024);
			int rowCtr= 0, perc= 0;
			long size= f.length(), bRead= 0;
			String lastC= null, nowC= null;
			int lastP= -1, nowP= -1;
			HashSet<String> setChr= new HashSet<String>();
			for (String s= buffy.readLine(); s!= null; s= buffy.readLine()) {
				++rowCtr;
				bRead+= s.length()+ 1;

                Log.progress(bRead, size);

				if (s.startsWith("track")|| s.startsWith("browser"))
					continue;

				int i= 0, len= s.length();
				char c= ' ';
				for (;i < len; i++) {
					c= s.charAt(i);
					if (c== ' '|| c== '\t')
						break;
				}
				nowC= s.substring(0, i);
				addRefID(nowC);
				while (i< len&& (c== ' '|| c== '\t'))
					c= s.charAt(++i);
				int p= i;
				for (;i < len; i++) {
					c= s.charAt(i);
					if (c== ' '|| c== '\t')
						break;
				}
				nowP= BEDobject.encodeInt(s, p, i);
				
				// not efficient
//				String[] sss= s.split("\\s");	// take whitspaces, for Ali M.
//				nowC= sss[0];
//				addRefID(nowC);
//				nowP= Integer.parseInt(sss[1]);
				
				boolean ascendingChr= false, chrYetRead= false;
				if (lastC!= null) {
					ascendingChr= lastC.compareTo(nowC)< 0;
					if (ascendingChr) {
						if (setChr.contains(nowC))
							chrYetRead= true;
						setChr.add(nowC);
					}
				}
				
				if ((!chrYetRead)&& (ascendingChr|| lastP<= nowP)) {					
					lastC= nowC;
					lastP= nowP;
					nowC= null;
				} else {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
						System.err.println("\n\tunsorted in line "+rowCtr+".");
					clearRefIDs();
					buffy.close();
					return false;
				}
			}
			buffy.close();
            Log.progressFinish(Constants.OK, true);
				
			nrUniqueLinesRead= rowCtr;
		} catch (Exception e) {
			e.printStackTrace();
            Log.progressFailed(" ERROR.");
		}
		
		
		return true;
	}

	public void read() throws Exception {
		read(0);
	}
	
	public void read(int bedLines) {

		if (bedLines== 0)
			bedLines= Integer.MAX_VALUE;
		Vector objV= new Vector();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(this.fPath+MyFile.separator+this.fName));
			for (String line= buffy.readLine();line!= null&& objV.size()< bedLines; line= buffy.readLine()) {
				if (line== null)
					break;
				line= line.trim();
				if (line.startsWith("browser")|| line.startsWith("track")|| line.length()< 1)
					continue;
				String[] tokens= line.split("\\s");	// not \\s+ for empty name
				if (tokens.length< 3) {
					System.out.println("WARNING: skipped incomplete line with "+tokens.length+" token");
					continue;
				}
				
				BEDobject bed= BEDobject.getRecycleObj(); // new BEDobject();
				bed.setChrom(tokens[0]);
				try {
					bed.setStart(Integer.parseInt(tokens[1]));
					bed.setEnd(Integer.parseInt(tokens[2]));
				} catch (NumberFormatException e) {
					System.currentTimeMillis();
					continue;
				}
				if (tokens.length> 3) 
					bed.setName(tokens[3]);
				if (tokens.length> 4) {
					try {
						bed.setScore(Integer.parseInt(tokens[4]));
					} catch (NumberFormatException e) {
						; //:) '.'
					}
				}
				
				if (tokens.length> 5) 
					bed.setStrand(tokens[5]);
				if (tokens.length> 6) 
					bed.setThickStart(Integer.parseInt(tokens[6]));
				if (tokens.length> 7) 
					bed.setThickEnd(Integer.parseInt(tokens[7]));
					
				if (tokens.length> 8) 
					bed.setCol(tokens[8]);

				if (tokens.length> 9) 
					bed.setBlockCount(Integer.parseInt(tokens[9]));
				if (tokens.length> 10) 
					bed.setBlockSizes(tokens[10]);
				if (tokens.length> 11)  
					bed.setBlockStarts(tokens[11]);

				objV.add(bed);
			}
			buffy.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		

		beds= (BEDobject[]) ArrayUtils.toField(objV);

	}

	HashMap<String,long[]> mapChr= new HashMap<String,long[]>(); // bytes and lines
	private ByteArrayCharSequence cs= new ByteArrayCharSequence(200);
	
	int nrUniqueLinesRead= -1;
	
	/**
	 * reads the rest of the lines from the reader and closes it.
	 */
	public void finish() {
		try {
			//ThreadedBufferedByteArrayStream buffy= getReader();
			BufferedBACSReader buffy= getReaderBACS();
			//for (cs= buffy.readLine(cs); cs.end!= 0; cs= buffy.readLine(cs)) {
			while (buffy.readLine(cs)> 0) {
				bytesRead+= cs.length()+ guessFileSep().length();
				++nrUniqueLinesRead;
			}
			close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		// DEBUG: output chr table
//		Object[] keys= mapChr.keySet().toArray();
//		for (int i = 0; i < keys.length; i++) {
//			if (mapChr.get(keys[i])!= null)
//				System.err.println(keys[i]+"\t"+mapChr.get(keys[i])[1]);
//		}
	}
	
	public int guessReadLen() {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fPath+ File.separator+ fName));
			String line= buffy.readLine();
			while(line.startsWith("track")|| line.startsWith("browser")) {
				line= buffy.readLine();
			}
			String[] ss= line.split("\\s");
			if (ss.length< 3) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
					System.err.println("\t[OOOH] Not a valid first BED line");
					System.err.println("\t"+line);
				}
				return -1;
			}
			if (ss.length< 11) 
				try {
					int x= Integer.parseInt(ss[1]);
					int y= Integer.parseInt(ss[2]);
					return (y- x);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						System.err.println("\t[OHNO] First BED line does not comply with specification");
						System.err.println("\t"+line);
					}
					return -1;
				}
			// else
			String[] sss= ss[10].split(",");
			int sum= 0;
			for (int i = 0; i < sss.length; i++) {
				try {
					sum+= Integer.parseInt(sss[i]);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						System.err.println("\t[OHNO] Invalid first BED line");
						System.err.println("\t"+line);
					}
					return -1;
				}
			}
			return sum;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return -1;
	}

	private int identTok= -1;
	private static final String TRACK= "track", BROWSER= "browser";
	private ByteArrayCharSequence lastLine= null;
	boolean reuse= true;
	
	private boolean addChr(String chrToki, long bytes, int lines) {
		if (mapChr.containsKey(chrToki))
			return false;
		else {

			//DEBUG
/*			String thisChr= chrToki.toString();
			if (true) {	// thisChr.equals("3RHet")
				String chk= null;
				try {
					// DEBUG
					File file = new File(this.fPath+MyFile.separator+this.fName);
					InputStream inputStream = new FileInputStream(file);
					inputStream.skip(bytes);	// must read next line 
					BufferedReader r2 = new BufferedReader(new InputStreamReader(inputStream));
					chk= r2.readLine();
					r2.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
				System.currentTimeMillis();
			}
*/			
			mapChr.put(chrToki,new long[] {bytes,lines});
			return true;
		}

	}
	
private BEDobject[] toObjectsOld(Vector<BEDobject> objV) {
		BEDobject[] beds= new BEDobject[objV.size()];
		for (int i = 0; i < beds.length; i++) 
			beds[i]= objV.elementAt(i);
		return beds;
	}

private BEDobject2[] toObjects(Vector<BEDobject2> objV) {
	BEDobject2[] beds= new BEDobject2[objV.size()];
	for (int i = 0; i < beds.length; i++) 
		beds[i]= objV.elementAt(i);
	return beds;
}

	private static final char TAB= '\t';
	
	public void write() throws Exception {
		write(false);
	}
	
	public int countLines(Progressable prog) {
		try {
			if (prog!= null)
				prog.start("progress ");
			int cnt= 0;
			File f= new File(this.fPath+MyFile.separator+this.fName);
			BufferedReader buffy= new BufferedReader(new FileReader(f));
			long bRead= 0, bTot= f.length();
			int perc= 0;
			for(String s; (s= buffy.readLine())!= null;++cnt,bRead+= s.length()+1) {
				if (bRead*10d/bTot> perc) {
					if (prog!= null)
						prog.progress();
					++perc;
				}
			}
			buffy.close();
			if (prog!= null)
				prog.finish();
			return cnt;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return -1;
	}
	
	int countAll;
	int countEntire;
	int countSplit;
	int countReads;
	public boolean checkReadDescriptor(UniversalReadDescriptor descriptor2) {
		
		ReadDescriptor descriptor= null;
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(getAbsFileName()));
			
			String s;
			while (((s= buffy.readLine())!= null)&&
					(s.trim().length()== 0
					|| s.startsWith(Constants.HASH)
					|| s.startsWith(BROWSER)
					|| s.startsWith(TRACK)));
					
			buffy.close();
			
			if (s== null)
				return false;
		
			String[] ss= s.split("\\s");
			if (ss.length< 4)
				return false;
			
			// check descriptor
			if (descriptor2.getAttributes(ss[3], null)== null) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("[OHNO] Read descriptor "+descriptor2+" not applicable for read ID\n\t"+ 
							ss[3]);
				return false;
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
			
		return true;
	}
	
	private int scanFileReadLines= 0;
	public boolean scanFile() {
        BufferedReader buffy = null;
        BufferedWriter tmpWriter = null;
        PipedInputStream in = null;
        PipedOutputStream out = null;
        Future sorterFuture = null;
		try {
			scanFileReadLines= 0;
			countAll= 0; countEntire= 0; countSplit= 0; countReads= 0;
			
			File f= new File(this.fPath+MyFile.separator+this.fName);
			buffy= new BufferedReader(new FileReader(f));
			int sepLen= guessFileSep().length();
			long bRead= 0, bTot= f.length();

			out = new PipedOutputStream();
			in = new PipedInputStream(out);
			tmpWriter= new BufferedWriter(new OutputStreamWriter(out));

            sorterFuture = Sorter.create(in, new DevNullOutputStream(), true)
                    .separator("\\s")
                    .field(0, false)
                    .addInterceptor(new Interceptable.Interceptor<String>() {
                        String lastLine= null;
                        public String intercept(String line) {
                            if (lastLine== null|| !line.equals(lastLine))
                                ++countReads;
                            lastLine= line;
                            return line;
                        }
                    })
                    .sortInBackground();

			final String COMA= ",";
			for(String s; (s= buffy.readLine())!= null;bRead+= s.length()+ sepLen) {
                if(!s.isEmpty())++countAll;
				++nrUniqueLinesRead;
				if (s.startsWith(BROWSER)|| s.startsWith(TRACK))
					continue;
				
				// last col tells you whether it is a 
				int p= s.length();
				while (p> 0&& Character.isWhitespace(s.charAt(--p))); // trailing ws
				while (p> 0&& !Character.isWhitespace(s.charAt(--p)));
				if (s.indexOf(COMA, p)>= 0)
					++countSplit;
				else if(p >= 0 && !s.isEmpty()){
					++countEntire;
                }

				// get ID
				p= 0;
				int cnt= 0;
				while(cnt< 3) {
					while (p< s.length()&& Character.isWhitespace(s.charAt(p++)));
					--p;
                    if (p >= 0){
                        while (p< s.length()&& !Character.isWhitespace(s.charAt(p++)));
                        if (p< s.length())
                            ++cnt;
                        else
                            break;
                    }else{
                        break;
                    }
				}
				int from= (cnt== 3)? p: -1, to= -1;
				while (p >= 0 && p< s.length()&& Character.isWhitespace(s.charAt(p++)));
				while (p >= 0 && p< s.length()&& !Character.isWhitespace(s.charAt(p++)));
				--p;
				if (p< s.length())
					to= p;
				
				if (from>= 0&& to>= 0) {
					String id= s.substring(from, to);
					tmpWriter.write(id);
					tmpWriter.write((int) '\n');
				}
			}
            tmpWriter.flush();
            tmpWriter.close();
            sorterFuture.get();
			return true;
			
		} catch (Exception e) {
            Log.error("Error while scanning BED file : " + e.getMessage(), e);
			return false;
		}finally {
            if(buffy != null)try {buffy.close();} catch (IOException e) {}
            if(tmpWriter != null)try {tmpWriter.flush();tmpWriter.close();} catch (IOException e) {}
            if(in != null)try {in.close();} catch (IOException e) {}
            if(out != null)try {out.close();} catch (IOException e) {}
            if(sorterFuture != null)sorterFuture.cancel(true);
        }
	}
	
	
	public void write(boolean append) {
		try {
			BufferedWriter buffy= new BufferedWriter(new FileWriter(this.fPath+MyFile.separator+this.fName, append));
			for (int i = 0; beds!= null&& i < beds.length&&beds[i]!= null; i++) {
				buffy.write(beds[i].toString()+"\n");
			}
			buffy.flush();
			buffy.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		test();
	}

	public BEDobject[] getBeds() {
		return beds;
	}

	public void setBeds(BEDobject[] beds) {
		this.beds = beds;
	}

	// not static for fileSep
	public File sortBED(File f) {
        InputStream in = null;
        OutputStream out = null;
        try {
            File outFile = File.createTempFile(f.getName() + "_", "_sorted");
            in = new BufferedInputStream(new FileInputStream(f));
            out = new BufferedOutputStream(new FileOutputStream(outFile));
            Sorter.create(in, out, true)
                    .separator("\\s")
                    .field(0, false)
                    .field(1, false)
                    .sort();
			return outFile;
			
        } catch (Exception e) {
            Log.progressFailed("ERROR");
            Log.error("Error while sorting file!", e);
        }finally {
            if(in != null ) try {in.close();} catch (IOException e) {}
            if(out != null ) try {out.close();} catch (IOException e) {}
        }
        return null;
	}

	public ByteArrayCharSequence sweepToChromosome(CharSequence chr) {
		try {
			long saveBytesRead= bytesRead;
			int saveUniqueLines= nrUniqueLinesRead;
			
			ByteArrayCharSequence cs= new ByteArrayCharSequence(this.cs.a.length);	// this.cs;
			BufferedBACSReader buffy= getReaderBACS();
			//for (cs= getReader().readLine(cs); cs.end!= 0; cs=getReader().readLine(cs)) {
			while (buffy.readLine(cs)> 0) {
				bytesRead+= cs.length()+guessFileSep().length();
				++nrUniqueLinesRead;
				if (cs.startsWith(chr))
					return cs;
			}
			
			reset(saveBytesRead, saveUniqueLines);
//			bytesRead= saveBytesRead;
//			nrUniqueLinesRead= saveUniqueLines;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return null;
	}
	
	public void reset() {
		reset(0,0);
	}
	
	public boolean close() {
		try {
			if (readerB!= null) {
				readerB.setCloseInputStream(true);
				readerB.close();
			}
			if (readerC!= null) {
				readerC.close();
			}
			return true;
		} catch (Exception e) {
			return false;
		}
	}
	
	public boolean reset(String chr) {
		// assert(mapChr.containsKey(chr)); 
		if (!mapChr.containsKey(chr))
			sweepToChromosome(chr);

		if (mapChr.get(chr)== null)
			return false;
		
		long[] bytesNlines= mapChr.get(chr);
		reset(bytesNlines[0], (int) bytesNlines[1]);
		return true;
	}
	
	public void reset(long bytes, int lines) {

		if (readerB!= null)
			try {
				readerB.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		if (readerC!= null)
			try {
				readerC.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		bytesRead= bytes;
		nrUniqueLinesRead= lines;
		if (readerB!= null&& readerB.isThreaded())
			readerB.setStop(true);
		if (reuse) {
			readerB= null;
			readerC= null;
			lastLine= null;
		}
	}
	
	File baseFile;
	File getFile() {
		 if (baseFile == null) {
			baseFile = new File(fPath+File.separator+fName);			
		}

		return baseFile;
	}

	public HashSet<String> getRefIDset() {
		return refIDset;
	}

	public int getNrLines() {
		return nrUniqueLinesRead;
	}

	public int getCountAll() {
		return countAll;
	}

	public int getCountEntire() {
		return countEntire;
	}

	public int getCountSplit() {
		return countSplit;
	}

	public int getCountReads() {
		return countReads;
	}

	BufferedBACSReader readerC;
	protected BufferedBACSReader getReaderBACS() {
		if (readerC == null) {
			try {
				InputStream inputStream = new FileInputStream(file);
				inputStream.skip(bytesRead);
				readerC= new BufferedBACSReader(inputStream);
				//new ThreadedBufferedByteArrayStream(10* 1024* 1024, inputStream, true, false);
				//readerB.setCloseInputStream(false);
				//reader = new BufferedReader(new InputStreamReader(inputStream));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	
		return readerC;
	}

	public int getScanFileReadLines() {
		return scanFileReadLines;
	}

	protected ThreadedBufferedByteArrayStream getReader() {
		if (readerB == null) {
			try {
				InputStream inputStream = new FileInputStream(file);
				inputStream.skip(bytesRead);
				readerB= new ThreadedBufferedByteArrayStream(10* 1024* 1024, inputStream, true, false);
				readerB.setCloseInputStream(false);
				//reader = new BufferedReader(new InputStreamReader(inputStream));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	
		return readerB;
	}

	public BEDobject[] read_old(String chr, int start, int end) {
			
			if (mapChr.containsKey(chr)) {
				if (mapChr.get(chr)== null)
					return null;
			}  
				
			
			--start;	// convert to bed coordinate
			
			Vector<BEDobject> objV= new Vector<BEDobject>();
			try {
				//BufferedReader buffy= getReader();
				ThreadedBufferedByteArrayStream buffy= getReader();
				ByteArrayCharSequence cs= this.cs; // new ByteArrayCharSequence(100);
				//String line;
				String lastChrRead= null;
				guessFileSep();
				boolean inited= false;
				if (reuse&& lastLine!= null) {
					cs= lastLine;
					lastLine= null;
					inited= true;
				}
	
				//for (cs= buffy.readLine(cs); cs.end!= 0; cs= buffy.readLine(cs)) {
				while (true) {
					
					long tmpBytes= bytesRead;			
	
					if (inited) {
						inited= false;
					} else {
						cs= buffy.readLine(cs);
						bytesRead+= cs.length()+ fileSep.length();
						++nrUniqueLinesRead;
					}
					if (cs.end== 0)
						break; // EOF
					
					// use this for debugging fpointer
	//				File file = new File(this.fPath+MyFile.separator+this.fName);
	//				InputStream inputStream = new FileInputStream(file);
	//				inputStream.skip(bytesRead);	// must read next line 
	//				BufferedReader r2 = new BufferedReader(new InputStreamReader(inputStream));
	//				String chk= r2.readLine();
	//				r2.close();
					
					if (cs.startsWith(TRACK)|| cs.startsWith(BROWSER))
						continue;
					
					// check if in range
					int bedStart= -1, bedEnd= -1;
					String chrToki= null;
					try {
						int toks= cs.countTokens(TAB);
						if (identTok< 0)
							identTok= toks;
						else
							if (identTok!= toks&& false) {	// can be now, we read reads and split reads
								if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
									System.err.println("\t\n[OHLALA] line "+nrUniqueLinesRead+" has not "+identTok+" elements as the lines before!");
									System.err.println("\t"+ cs.toString());
									System.err.println("\tcheck file "+ fName);
								}
							}
						if (toks< 3) {
							if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
								System.err.println("\t\n[OHNOO] line "+nrUniqueLinesRead+" has less than 3 token, I am skipping.");
								System.err.println("\t"+ cs.toString());
								System.err.println("\tcheck file "+ fName);
							}
							continue;
						}
						chrToki= cs.getToken(1, TAB).toString();
						try {
							bedStart= BEDobject.encodeInt(cs.getToken(2, TAB)); 
							bedEnd= BEDobject.encodeInt(cs.getToken(3, TAB));
						} catch (NumberFormatException e) {
							e.printStackTrace();
						}
					} catch (Exception e) {
						e.printStackTrace();
					}
					
					if (chrToki.compareTo(chr)> 0) {	// 090520 check whether reads too far
						if (reuse)
							lastLine= cs;
						else {
							bytesRead= tmpBytes; 
							--nrUniqueLinesRead;
						}
						if (objV== null)
							return null;
						BEDobject[] obj= toObjectsOld(objV);
						addChr(chrToki, bytesRead- cs.length()- guessFileSep().length(), nrUniqueLinesRead- 1);
						
						return obj;
					}
						
					if (lastChrRead== null) {	// first line read in this batch
						
						addChr(chrToki, bytesRead- cs.length()- guessFileSep().length(), nrUniqueLinesRead- 1);
						
						if (chr!= null&& !cs.subSequence(0, chr.length()).equals(chr)) {
							if (mapChr.containsKey(chr)) {
								if (mapChr.get(chr)== null) {
									if (reuse)
										lastLine= cs;
									else {
										bytesRead= tmpBytes;
										--nrUniqueLinesRead;
	//									buffy.setStop(true);	//close();
	//									readerB= null;
									}
									return null;
								} else {
									if (mapChr.get(chr)[0]> tmpBytes) { 	// only jump forward, never back
										long[] bytesNlines= mapChr.get(chr);
										reset(bytesNlines[0], (int) bytesNlines[1]);
										lastChrRead= chr.toString();
										buffy= getReader();
										cs= buffy.readLine(cs);
										//cs= new ByteArrayCharSequence(line);
										bytesRead+= cs.length()+ fileSep.length();
										++nrUniqueLinesRead;
										chrToki= cs.getToken(1, TAB).toString();
										try {
											bedStart= BEDobject.encodeInt(cs.getToken(2, TAB)); 
											bedEnd= BEDobject.encodeInt(cs.getToken(3, TAB));
										} catch (Exception e) {
											System.err.println("[ERROR] encodeint:\n"+ cs.toString());
											System.currentTimeMillis();
										}
									} else {
										reset(tmpBytes, nrUniqueLinesRead-1);
										return null;
									}
								}
							} else {
								ByteArrayCharSequence newCS= sweepToChromosome(chr);
								if (newCS== null) {	// not found
									mapChr.put(chr, null);	// BUG: not chrToki
									if (reuse) 
										lastLine= cs;
									else {
										bytesRead= tmpBytes;
										--nrUniqueLinesRead;
										if (buffy.isThreaded())
											buffy.setStop(true);	//close();
										buffy.close();
										readerB= null;
									}
									return null; 
								} else {
									
									this.cs= newCS;
									cs= newCS;
									chrToki= cs.getToken(1, TAB).toString();	
									addChr(chrToki,
											bytesRead- cs.length()- guessFileSep().length(), 
											nrUniqueLinesRead-1);
											
									buffy= getReader();
									lastChrRead= chr.toString();
									chrToki= cs.getToken(1, TAB).toString();
									bedStart= BEDobject.encodeInt(cs.getToken(2, TAB)); 
									bedEnd= BEDobject.encodeInt(cs.getToken(3, TAB));
								}
							}
						}
					}
	
					//line= line.trim();
					if (cs.startsWith("browser")|| cs.startsWith("track")|| cs.length()< 1)
						continue;
					
					if (lastChrRead== null)
						lastChrRead= chrToki;
					else {
						if (!lastChrRead.equals(chrToki)) {	// changes chr
							if (reuse)
								lastLine= cs;
							else {
								bytesRead= tmpBytes;
								--nrUniqueLinesRead;
							}
							addChr(chrToki, bytesRead- cs.length()- guessFileSep().length(), nrUniqueLinesRead);
							break;
						} 
					}
					
	
					boolean stop= false, continues= false;
					
					if (start>= 0&& bedEnd< start)
						continues= true;
					else if (end>= 0&& bedStart> end)
						stop= true;
					
					if (continues)
						continue;
					if (stop) {	// not found on this chr
						if (reuse)
							lastLine= cs;
						else {
							bytesRead= tmpBytes;
							--nrUniqueLinesRead;
						}
						if (objV.size()== 0)
							return null;
						return toObjectsOld(objV);
					}
					
					
					// create object
					//String[] tokens= line.split("\\s");	// not \\s+ for empty name
					BEDobject bed= BEDobject.createBEDobject(cs, 
							chrToki, bedStart, bedEnd); //.getRecycleObj();
					objV.add(bed);
	
				}
				if (buffy.isThreaded())
					buffy.setStop(true);
				
				//readerB= null;	// always reset, for jumping around
	
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			return toObjectsOld(objV);
		}

	public int get(String chr, int start, int end, SyncIOHandler2 handler, OutputStream ostream) {
			
			int count= 0;
			
			if (mapChr.containsKey(chr)) {
				if (mapChr.get(chr)== null)
					return count;
			}  
				
			
			--start;	// convert to bed coordinate
			try {
				BufferedBACSReader buffy= getReaderBACS();
				ByteArrayCharSequence cs= this.cs; // new ByteArrayCharSequence(100);
				String lastChrRead= null;
				guessFileSep();
				boolean inited= false;
				if (reuse&& lastLine!= null) {
					cs= lastLine;
					lastLine= null;
					inited= true;
				}
	
				//for (cs= buffy.readLine(cs); cs.end!= 0; cs= buffy.readLine(cs)) {
				while (true) {
					
					long tmpBytes= bytesRead;			
	
					if (inited) {
						inited= false;
					} else {
						if (buffy.readLine(cs)<= 0)
							break;	// EOF
						bytesRead+= cs.length()+ fileSep.length();
						++nrUniqueLinesRead;
					}
					
					// use this for debugging fpointer
	//				File file = new File(this.fPath+MyFile.separator+this.fName);
	//				InputStream inputStream = new FileInputStream(file);
	//				inputStream.skip(bytesRead);	// must read next line 
	//				BufferedReader r2 = new BufferedReader(new InputStreamReader(inputStream));
	//				String chk= r2.readLine();
	//				r2.close();
					
					if (cs.startsWith(TRACK)|| cs.startsWith(BROWSER))
						continue;
					
					// check if in range
					int bedStart= -1, bedEnd= -1;
					String chrToki= null;
					try {
						int toks= cs.countTokens(TAB);
						if (identTok< 0)
							identTok= toks;
						else
							if (identTok!= toks&& false) {	// can be now, we read reads and split reads
								if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
									System.err.println("\t\n[OHLALA] line "+nrUniqueLinesRead+" has not "+identTok+" elements as the lines before!");
									System.err.println("\t"+ cs.toString());
									System.err.println("\tcheck file "+ fName);
								}
							}
						if (toks< 3) {
							if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
								System.err.println("\t\n[OHNOO] line "+nrUniqueLinesRead+" has less than 3 token, I am skipping.");
								System.err.println("\t"+ cs.toString());
								System.err.println("\tcheck file "+ fName);
							}
							continue;
						}
						chrToki= cs.getToken(1, TAB).toString();
						try {
							bedStart= BEDobject.encodeInt(cs.getToken(2, TAB)); 
							bedEnd= BEDobject.encodeInt(cs.getToken(3, TAB));
						} catch (NumberFormatException e) {
							e.printStackTrace();
						}
					} catch (Exception e) {
						e.printStackTrace();
					}
					
					if (chrToki.compareTo(chr)> 0) {	// 090520 check whether reads too far
						if (reuse)
							lastLine= cs;
						else {
							bytesRead= tmpBytes; 
							--nrUniqueLinesRead;
						}
						addChr(chrToki, bytesRead- cs.length()- guessFileSep().length(), nrUniqueLinesRead- 1);
						return count;
					}
						
					if (lastChrRead== null) {	// first line read in this batch
						
						addChr(chrToki, bytesRead- cs.length()- guessFileSep().length(), nrUniqueLinesRead- 1);
						
						if (chr!= null&& !cs.subSequence(0, chr.length()).equals(chr)) {
							if (mapChr.containsKey(chr)) {
								if (mapChr.get(chr)== null) {
									if (reuse)
										lastLine= cs;
									else {
										bytesRead= tmpBytes;
										--nrUniqueLinesRead;
									}
									return 0;
								} else {
									if (mapChr.get(chr)[0]> tmpBytes) { 	// only jump forward, never back
										long[] bytesNlines= mapChr.get(chr);
										reset(bytesNlines[0], (int) bytesNlines[1]);
										lastChrRead= chr.toString();
										buffy= getReaderBACS();
										if (buffy.readLine(cs)<= 0)
											break;
										bytesRead+= cs.length()+ fileSep.length();
										++nrUniqueLinesRead;
										chrToki= cs.getToken(1, TAB).toString();
										try {
											bedStart= BEDobject.encodeInt(cs.getToken(2, TAB)); 
											bedEnd= BEDobject.encodeInt(cs.getToken(3, TAB));
										} catch (Exception e) {
											System.err.println("[ERROR] encodeint:\n"+ cs.toString());
											System.currentTimeMillis();
										}
									} else {
										reset(tmpBytes, nrUniqueLinesRead-1);
										return 0;
									}
								}
							} else {
								ByteArrayCharSequence newCS= sweepToChromosome(chr);
								if (newCS== null) {	// not found
									mapChr.put(chr, null);	// BUG: not chrToki
									if (reuse) 
										lastLine= cs;
									else {
										bytesRead= tmpBytes;
										--nrUniqueLinesRead;
										readerB= null;
									}
									return 0; 
								} else {
									
									this.cs= newCS;
									cs= newCS;
									chrToki= cs.getToken(1, TAB).toString();	
									addChr(chrToki,
											bytesRead- cs.length()- guessFileSep().length(), 
											nrUniqueLinesRead-1);
											
									buffy= getReaderBACS();
									lastChrRead= chr.toString();
									chrToki= cs.getToken(1, TAB).toString();
									bedStart= BEDobject.encodeInt(cs.getToken(2, TAB)); 
									bedEnd= BEDobject.encodeInt(cs.getToken(3, TAB));
								}
							}
						}
					}
	
					//line= line.trim();
					if (cs.startsWith("browser")|| cs.startsWith("track")|| cs.length()< 1)
						continue;
					
					if (lastChrRead== null)
						lastChrRead= chrToki;
					else {
						if (!lastChrRead.equals(chrToki)) {	// changes chr
							if (reuse)
								lastLine= cs;
							else {
								bytesRead= tmpBytes;
								--nrUniqueLinesRead;
							}
							addChr(chrToki, bytesRead- cs.length()- guessFileSep().length(), nrUniqueLinesRead);
							break;
						} 
					}
					
	
					boolean stop= false, continues= false;
					
					if (start>= 0&& bedEnd< start)
						continues= true;
					else if (end>= 0&& bedStart> end)
						stop= true;
					
					if (continues)
						continue;
					if (stop) {	// not found on this chr
						if (reuse)
							lastLine= cs;
						else {
							bytesRead= tmpBytes;
							--nrUniqueLinesRead;
						}
						return count;
					}
					
					
					// write line
					handler.writeLine(cs, ostream);
					++count;
//					BEDobject2 bed= new BEDobject2(cs); 
//					objV.add(bed);
					
	
				}
				
				//readerB= null;	// always reset, for jumping around
	
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			return count;
		}

	public BEDobject2[] read(String chr, int start, int end) {
				
				if (mapChr.containsKey(chr)) {
					if (mapChr.get(chr)== null)
						return null;
				}  
					
				
				--start;	// convert to bed coordinate
				
				Vector<BEDobject2> objV= new Vector<BEDobject2>();
				try {
					//BufferedReader buffy= getReader();
					//ThreadedBufferedByteArrayStream buffy= getReader();
					BufferedBACSReader buffy= getReaderBACS();
					ByteArrayCharSequence cs= this.cs; // new ByteArrayCharSequence(100);
					//String line;
					String lastChrRead= null;
					guessFileSep();
					boolean inited= false;
					if (reuse&& lastLine!= null) {
						cs= lastLine;
						lastLine= null;
						inited= true;
					}
		
					//for (cs= buffy.readLine(cs); cs.end!= 0; cs= buffy.readLine(cs)) {
					//ByteArrayCharSequence lastLine= null;
					while (true) {
						
						long tmpBytes= bytesRead;			
		
						if (inited) {
							inited= false;
						} else {
							lastLine= cs.cloneCurrentSeq();
							if (buffy.readLine(cs)<= 0)
								break;	// EOF
							cs.resetFind();
							bytesRead+= cs.length()+ fileSep.length();
							++nrUniqueLinesRead;
						}
						
						// use this for debugging fpointer
		//				File file = new File(this.fPath+MyFile.separator+this.fName);
		//				InputStream inputStream = new FileInputStream(file);
		//				inputStream.skip(bytesRead);	// must read next line 
		//				BufferedReader r2 = new BufferedReader(new InputStreamReader(inputStream));
		//				String chk= r2.readLine();
		//				r2.close();
						
						if (cs.startsWith(TRACK)|| cs.startsWith(BROWSER))
							continue;
						
						// check if in range
						int bedStart= -1, bedEnd= -1;
						String chrToki= null;
						try {
							int toks= cs.countTokens(TAB);
							if (identTok< 0)
								identTok= toks;
							else
								if (identTok!= toks&& false) {	// can be now, we read reads and split reads
									if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
										System.err.println("\t\n[OHLALA] line "+nrUniqueLinesRead+" has not "+identTok+" elements as the lines before!");
										System.err.println("\t"+ cs.toString());
										System.err.println("\tcheck file "+ fName);
									}
								}
							if (toks< 3) {
								if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
									System.err.println("\t\n[OHNOO] line "+nrUniqueLinesRead+" has less than 3 token, I am skipping.");
									System.err.println("\t"+ cs.toString());
									System.err.println("\tcheck file "+ fName);
								}
								continue;
							}
							
							cs.resetFind();	// DEBUG 100827: somehow getToken(0) results in NullPointer otherwise 
							chrToki= cs.getToken(0).toString();	// getToken(1, TAB)
//							try {
								bedStart= cs.getTokenInt(1); 	//BEDobject.encodeInt(cs.getToken(2, TAB)); 
								bedEnd= cs.getTokenInt(2);		//BEDobject.encodeInt(cs.getToken(3, TAB));
//							} catch (NumberFormatException e) {
//								e.printStackTrace();
//							}
						} catch (Exception e) {
							e.printStackTrace();
						}
						
						if (chrToki.compareTo(chr)> 0) {	// 090520 check whether reads too far
							if (reuse)
								lastLine= cs;
							else {
								bytesRead= tmpBytes; 
								--nrUniqueLinesRead;
							}
							if (objV== null)
								return null;
							BEDobject2[] obj= toObjects(objV);
							addChr(chrToki, bytesRead- cs.length()- guessFileSep().length(), nrUniqueLinesRead- 1);
							
							return obj;
						}
							
						if (lastChrRead== null) {	// first line read in this batch
							
							addChr(chrToki, bytesRead- cs.length()- guessFileSep().length(), nrUniqueLinesRead- 1);
							
							if (chr!= null&& !cs.subSequence(0, chr.length()).equals(chr)) {
								if (mapChr.containsKey(chr)) {
									if (mapChr.get(chr)== null) {
										if (reuse)
											lastLine= cs;
										else {
											bytesRead= tmpBytes;
											--nrUniqueLinesRead;
		//									buffy.setStop(true);	//close();
		//									readerB= null;
										}
										return null;
									} else {
										if (mapChr.get(chr)[0]> tmpBytes) { 	// only jump forward, never back
											long[] bytesNlines= mapChr.get(chr);
											reset(bytesNlines[0], (int) bytesNlines[1]);
											lastChrRead= chr.toString();
											buffy= getReaderBACS();
											if (buffy.readLine(cs)<= 0)
												break;
											//cs= new ByteArrayCharSequence(line);
											bytesRead+= cs.length()+ fileSep.length();
											++nrUniqueLinesRead;
											chrToki= cs.getToken(1, TAB).toString();
											try {
												bedStart= BEDobject.encodeInt(cs.getToken(2, TAB)); 
												bedEnd= BEDobject.encodeInt(cs.getToken(3, TAB));
											} catch (Exception e) {
												System.err.println("[ERROR] encodeint:\n"+ cs.toString());
											}
										} else {
											reset(tmpBytes, nrUniqueLinesRead-1);
											return null;
										}
									}
								} else {
									ByteArrayCharSequence newCS= sweepToChromosome(chr);
									if (newCS== null) {	// not found
										mapChr.put(chr, null);	// BUG: not chrToki
										if (reuse) 
											lastLine= cs;
										else {
											bytesRead= tmpBytes;
											--nrUniqueLinesRead;
											readerB= null;
										}
										return null; 
									} else {
										
										this.cs= newCS;
										cs= newCS;
										chrToki= cs.getToken(1, TAB).toString();	
										addChr(chrToki,
												bytesRead- cs.length()- guessFileSep().length(), 
												nrUniqueLinesRead-1);
												
										buffy= getReaderBACS();
										lastChrRead= chr.toString();
										chrToki= cs.getToken(1, TAB).toString();
										bedStart= BEDobject.encodeInt(cs.getToken(2, TAB)); 
										bedEnd= BEDobject.encodeInt(cs.getToken(3, TAB));
									}
								}
							}
						}
		
						//line= line.trim();
						if (cs.startsWith("browser")|| cs.startsWith("track")|| cs.length()< 1)
							continue;
						
						if (lastChrRead== null)
							lastChrRead= chrToki;
						else {
							if (!lastChrRead.equals(chrToki)) {	// changes chr
								if (reuse)
									lastLine= cs;
								else {
									bytesRead= tmpBytes;
									--nrUniqueLinesRead;
								}
								addChr(chrToki, bytesRead- cs.length()- guessFileSep().length(), nrUniqueLinesRead);
								break;
							} 
						}
						
		
						boolean stop= false, continues= false;
						
						if (start>= 0&& bedEnd< start)
							continues= true;
						else if (end>= 0&& bedStart> end)
							stop= true;
						
						if (continues)
							continue;
						if (stop) {	// not found on this chr
							if (reuse)
								lastLine= cs;
							else {
								bytesRead= tmpBytes;
								--nrUniqueLinesRead;
							}
							if (objV.size()== 0)
								return null;
							return toObjects(objV);
						}
						
						
						// create object
						//String[] tokens= line.split("\\s");	// not \\s+ for empty name
						BEDobject2 bed= new BEDobject2(cs); 
	//						BEDobject.createBEDobject(cs, 
	//							chrToki, bedStart, bedEnd); //.getRecycleObj();
						objV.add(bed);
		
					}
					
					//readerB= null;	// always reset, for jumping around
		
				} catch (Exception e) {					
					e.printStackTrace();
				}
				
				return toObjects(objV);
			}

	public ReadDescriptor checkReadDescriptor(boolean pairedEnd) {
		
		ReadDescriptor descriptor= null;
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(getAbsFileName()));
			
			String s;
			while (((s= buffy.readLine())!= null)&&
					(s.trim().length()== 0
					|| s.startsWith(Constants.HASH)
					|| s.startsWith(BROWSER)
					|| s.startsWith(TRACK)));
					
			buffy.close();
			
			if (s== null)
				return null;
		
			String[] ss= s.split("\\s");
			if (ss.length< 4)
				return null;
			
			// check descriptor
			descriptor= new SolexaPairedEndDescriptor();
			if (pairedEnd) {
				if (!descriptor.isPairedEnd(ss[3])) {
					// descriptor.getPairedEndInformation(bedWrapper.getBeds()[0].getName()))== 0)
					descriptor= new FMRD();
					if (!descriptor.isPairedEnd(ss[3])) {
						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
							System.err.println("[OHNO] Could not detect the format of read descriptor:\n\t"+ 
									s);
						return null;
					}
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
			
		return descriptor;
	}

}
