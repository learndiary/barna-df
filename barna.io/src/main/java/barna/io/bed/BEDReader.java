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

package barna.io.bed;

import barna.commons.ByteArrayCharSequence;
import barna.commons.Progressable;
import barna.commons.io.DevNullOutputStream;
import barna.commons.log.Log;
import barna.commons.system.OSChecker;
import barna.commons.utils.ArrayUtils;
import barna.commons.utils.Interceptable;
import barna.commons.utils.LineComparator;
import barna.io.*;
import barna.model.Mapping;
import barna.model.bed.BEDMapping;
import barna.model.bed.BEDobject;
import barna.model.constants.Constants;
import barna.model.rna.UniversalReadDescriptor;

import java.io.*;
import java.util.*;
import java.util.concurrent.Future;

public class BEDReader extends AbstractFileIOWrapper implements MappingReader {

    private boolean paired;

    static void test() {
		System.out.println(((byte) -1)| (byte) 1);
		System.out.println(((byte) 2)& ((byte) -1));
		System.out.println(((byte) 2)& ((byte) 1));
		System.out.println(2&0);
		System.out.println(-1&Integer.MAX_VALUE);
	}
	
	/**
	 * @deprecated
	 */
	ThreadedBufferedByteArrayStream readerB= null;
	
	BEDMapping[] beds= null;
    File tmpDir = null;
    boolean sortInRam = false;
    UniversalReadDescriptor descriptor;
    MSIterator<BEDMapping> currentIter;
    private int scoreFilter = -1;
	/**
	 * Creates an instance using a given file and
	 * line comparator.
	 * @param inputFile file to use
	 * @param comparator comparator describing required file
	 * sorting 
	 */
	public BEDReader(File inputFile, LineComparator<CharSequence> comparator, boolean sortInRam, UniversalReadDescriptor descriptor, File tmpDir) {
        this(inputFile, comparator, sortInRam, descriptor, tmpDir, -1);
	}

	public BEDReader(File inputFile, LineComparator<CharSequence> comparator, boolean sortInRam, UniversalReadDescriptor descriptor, File tmpDir, int scoreFilter) {
		super(inputFile);
		this.comparator= (comparator== null? COMPARATOR_DEFAULT: comparator);
        this.sortInRam = sortInRam;
        this.descriptor = descriptor;
        this.scoreFilter = scoreFilter;
	}

	/**
	 * Creates an instance using a specific file 
	 * and the default comparator.
	 * @param inputFile
	 */
	public BEDReader(File inputFile, Boolean sortInRam, UniversalReadDescriptor descriptor, File tmpDir) {
		this(inputFile, sortInRam, descriptor, tmpDir, -1);
	}

	public BEDReader(File inputFile, Boolean sortInRam, UniversalReadDescriptor descriptor, File tmpDir, int scoreFilter) {
		this(inputFile, COMPARATOR_DEFAULT, sortInRam, descriptor, tmpDir, scoreFilter);
	}

	/**
	 * Creates an instance using a specific path to a file 
	 * and the default line comparator.
	 * @param absolutePath path to the file the wrapper is based on
	 */
	public BEDReader(String absolutePath, Boolean sortInRam, UniversalReadDescriptor descriptor, File tmpDir) {
		this(new File(absolutePath), sortInRam,  descriptor, tmpDir);
	}

	    public BEDReader(String absolutePath) {
        this(new File(absolutePath), false,  null, null);   // TODO Not the proper way!!!
    }

    public BEDReader(File inputFile) {
        this(inputFile, false,  null, null);   // TODO Not the proper way!!!
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
	
	@Override
	public boolean isApplicable() {
		long lines= isApplicable(inputFile);
		if (lines< 0)
			return false;
		return true;
	}
	
	/**
	 * Checks for correct sorting, returns number of lines read (&lt;0 if not applicable).
	 * @param inputFile file from which is read
	 * @return number of lines read, or -(number of lines read) up to the unsorted
	 * entry
	 */
	public long isApplicable(File inputFile) {
		
		FileInputStream fis= null;
		try {
			fis= new FileInputStream(inputFile);
			long linesOK= -1;
			if (comparator== COMPARATOR_DEFAULT)
				linesOK= isApplicableDefault(fis, FileHelper.getSize(inputFile));
			else
				linesOK= isApplicable(fis, FileHelper.getSize(inputFile));
			return linesOK;
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			if (fis!= null)
				try {
					fis.close();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
		}
	}
	
	/**
	 * Checks for correct sorting, returns number of lines read (&lt;0 if not applicable).
	 * <b>Note:</b> does not close the given stream.
	 * @param inputStream stream from which is read
	 * @param size total size of data in the stream, if known, otherwise &le; 0
	 * @return number of lines read, or -(number of lines read) up to the unsorted
	 * entry
	 */
	public long isApplicable(InputStream inputStream, long size) {

        LineComparator clone = new LineComparator(comparator);
        try {
            Log.progressStart("checking");
			BufferedReader buffy= new BufferedReader(new InputStreamReader(inputStream), 10* 1024* 1024);
			long rowCtr= 0;
			long bRead= 0;
			String lastRow= null;
			
			for (String s= buffy.readLine(); s!= null; s= buffy.readLine()) {
				++rowCtr;
				bRead+= s.length()+ 1;

				if (size> 0)
					Log.progress(bRead, size);

				if (s.startsWith("track")|| s.startsWith("browser"))
					continue;

				if (lastRow== null) {
					lastRow= s;
					continue;
				}
				
				if (clone.compare(lastRow, s)> 0) {
					Log.info("\n\tunsorted in line "+rowCtr+".");
					buffy.close();
					return (-rowCtr);
				}
				
				lastRow= s;
				
			}
			buffy.close();
            Log.progressFinish(Constants.OK, true);
				
			return rowCtr;
		} catch (Exception e) {
            Log.progressFailed(" ERROR.");
            throw new RuntimeException(e);
		}
	}

	/**
	 * @deprecated based on BedObject
	 */
	public void read() {
		read(0);
	}
	
	/**
	 * @deprecated based on BedObject
	 * @param bedLines
	 */
	public void read(int bedLines) {

		if (bedLines== 0)
			bedLines= Integer.MAX_VALUE;
		Vector objV= new Vector();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(getInputFile()));
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
			throw new RuntimeException(e);
		}
		

		beds= (BEDMapping[]) ArrayUtils.toField(objV);

	}

	HashMap<String,long[]> mapChr= new HashMap<String,long[]>(); // bytes and lines
	private ByteArrayCharSequence cs= new ByteArrayCharSequence(200);
	
	int nrUniqueLinesRead= 0;
	
	/**
	 * reads the rest of the lines from the reader and closes it.
	 */
	public void finish() {
		try {
			//ThreadedBufferedByteArrayStream buffy= getReader();
			BufferedBACSReader buffy= getReaderBACS();
			//for (cs= buffy.readLine(cs); cs.end!= 0; cs= buffy.readLine(cs)) {
			while (buffy.readLine(cs)!= null) {
				bytesRead+= cs.length()+ getLineSeparator().length();
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
			BufferedReader buffy= new BufferedReader(new FileReader(getInputFile()));
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
	/**
	 * Default comparator, sort (1) chromosome, (2) position.
	 */
	public static final LineComparator<CharSequence> COMPARATOR_DEFAULT=
		new LineComparator<CharSequence>(false, "\t", 0)
                .addComparator(new LineComparator<CharSequence>(true, "\t", 1))
				.addComparator(new LineComparator<CharSequence>(true, "\t", 2))
				.addComparator(new LineComparator<CharSequence>(new Comparator<CharSequence>() {
                    @Override
                    public int compare(CharSequence o1, CharSequence o2) {
                        int n1 = o1.length(), n2 = o2.length();
                        for (int i1 = 0, i2 = 0; i1 < n1 && i2 < n2; i1++, i2++) {
                            char c1 = o1.charAt(i1);
                            char c2 = o2.charAt(i2);
                            if (c1 != c2) {
                                return c1 - c2;
                            }
                        }
                        return n1 - n2;
                    }
                }));


	/**
	 * Default comparator for read pairing, sort (1) chromosome, (2) name,
	 * (3) position.
	 */
	public static final LineComparator<CharSequence> COMPARATOR_PAIRED_END=
		new LineComparator<CharSequence>(false, "\t", 0)
			.addComparator(new LineComparator<CharSequence>(false, "\t", 3))
					.addComparator(new LineComparator<CharSequence>(true, "\t", 1));
	
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

private BEDMapping[] toObjects(Vector<BEDMapping> objV) {
	if (objV.size()== 0)
		return null;
	BEDMapping[] beds= new BEDMapping[objV.size()];
	for (int i = 0; i < beds.length; i++) 
		beds[i]= objV.elementAt(i);
	return beds;
}

	private static final char TAB= '\t';
	
	public void write() {
		write(false);
	}
	
	public int countLines(Progressable prog) {
		try {
			if (prog!= null)
				prog.start("progress ");
			int cnt= 0;
			BufferedReader buffy= new BufferedReader(new FileReader(getInputFile()));
			long bRead= 0, bTot= getInputSize();
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
	
	/**
	 * Sorting rules for that file.
	 */
	LineComparator<CharSequence> comparator;
	
	int countAll;
	int countEntire;
	int countSplit;
	int countReads;
	public boolean checkReadDescriptor(UniversalReadDescriptor descriptor) {

		BufferedReader buffy= null;
		try {
			buffy= new BufferedReader(new FileReader(getInputFile()));
			
			String s;
			while (((s= buffy.readLine())!= null)&&
					(s.trim().length()== 0
					|| s.startsWith(Constants.HASH)
					|| s.startsWith(BROWSER)
					|| s.startsWith(TRACK)));
					
			if (s== null)
				return false;
		
			String[] ss= s.split("\\t");
			if (ss.length< 4)
				return false;
			
			// check descriptor
			if (descriptor.getAttributes(ss[3], null)== null) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("[OHNO] Read descriptor "+descriptor+" not applicable for read ID\n\t"+ 
							ss[3]);
				return false;
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (buffy!= null) 
				try {
					buffy.close();
				} catch (Exception e) {
					throw new RuntimeException(e);
				}
		}
			
		return true;
	}
	
	private int scanFileReadLines= 0;
	public void scanFile() {
		
		warnFirstSkip= true;
		skippedLines= 0;
		
        BufferedReader buffy = null;
        BufferedWriter tmpWriter = null;
        PipedInputStream in = null;
        PipedOutputStream out = null;
        Future sorterFuture = null;
		try {
			scanFileReadLines= 0;
			countAll= 0; countEntire= 0; countSplit= 0; countReads= 0;
			
			buffy= new BufferedReader(new FileReader(getInputFile()));
			int sepLen= getLineSeparator().length();
			long bRead= 0, bTot= getInputSize();

			out = new PipedOutputStream();
			in = new PipedInputStream(out);
			tmpWriter= new BufferedWriter(new OutputStreamWriter(out));

            sorterFuture = Sorter.create(in, new DevNullOutputStream(), true, "\t")
                    .field(0, false)
                    .addInterceptor(new Interceptable.Interceptor<String>() {
                        String lastLine = null;

                        public String intercept(String line) {
                            if (lastLine == null || !line.equals(lastLine)) {
                                ++countReads;
                            }
                            lastLine = line;
                            return line;
                        }
                    })
                    .sortInBackground();

			final String COMA= ",";
			for(String s; (s= buffy.readLine())!= null;bRead+= s.length()+ sepLen) {
                if(!s.isEmpty())
                    ++countAll;
				++nrUniqueLinesRead;
				if (s.startsWith(BROWSER)|| s.startsWith(TRACK)) {
					++skippedLines;
					if (warnFirstSkip) {
						Log.warn("Skipped control line: "+ cs.toString());
						warnFirstSkip= false;
					}
					continue;
				}
				
				// last col tells you whether it is a split-mappings
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
				while (p >= 0 && p< s.length()&& s.charAt(p++)!='\t');
				--p;
				if (p< s.length())
					to= p;
				
				if (from>= 0&& to>= 0) {
					String id= s.substring(from, to);
					tmpWriter.write(id);
					tmpWriter.write(OSChecker.NEW_LINE);
				} else {
					++skippedLines;
					if (warnFirstSkip) {
						Log.warn("Skipped too shortline: "+ s);
						warnFirstSkip= false;
					}
				}
			}
            tmpWriter.flush();
            tmpWriter.close();
            sorterFuture.get();
			
		} catch (Exception e) {
			throw new RuntimeException(e);
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
			BufferedWriter buffy= new BufferedWriter(new FileWriter(getInputFile(), append));
			for (int i = 0; beds!= null&& i < beds.length&&beds[i]!= null; i++) {
				buffy.write(beds[i].toString()+barna.commons.system.OSChecker.NEW_LINE);
			}
			buffy.flush();
			buffy.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	public File getSortedFile(File tmpFile, LineComparator<CharSequence> comparator) {
	
		//BEDReader wrapper= new BEDReader(inputFile, comparator);
		if (!this.isApplicable()) {
			if (tmpFile== null)
				try {
					tmpFile= FileHelper.createTempFile(
							FileHelper.stripExtension(inputFile.getName()), 
							FileHelper.getExtension(inputFile));
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			this.sort(tmpFile);
			
			inputFile= tmpFile;
		}

		return inputFile;
	}


	public BEDobject[] getBeds() {
		return new BEDobject[beds.length];     //TODO not working
	}

	public void setBeds(BEDMapping[] beds) {
		this.beds = beds;
	}

	public Sorter getSorter(InputStream in, OutputStream out) {
		Sorter sorter= Sorter.create(in, out, true, "\t")
			.field(comparator);
			
		return sorter;
	}
	
	@Override
	public void sort(OutputStream outputStream) {
        InputStream in = null;
        OutputStream out = null;
        try {
            FileInputStream iStream = new FileInputStream(getInputFile());
            getSorter(iStream, outputStream).sort();
			return;
			
        } catch (Exception e) {
            Log.progressFailed("ERROR");
            Log.error("Error while sorting file!", e);
        }finally {
            if(in != null ) try {in.close();} catch (IOException e) {}
            if(out != null ) try {out.close();} catch (IOException e) {}
        }
	}

	public ByteArrayCharSequence sweepToChromosome(CharSequence chr) {
		try {
			long saveBytesRead= bytesRead;
			int saveUniqueLines= nrUniqueLinesRead;
			
			ByteArrayCharSequence cs= new ByteArrayCharSequence(this.cs.chars.length);	// this.cs;
			BufferedBACSReader buffy= getReaderBACS();
			//for (cs= getReader().readLine(cs); cs.end!= 0; cs=getReader().readLine(cs)) {
			while (buffy.readLine(cs)!= null) {
                bytesRead+= cs.length()+getLineSeparator().length();
                ++nrUniqueLinesRead;
                if (cs.getToken(0).equals(chr))
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

    @Override
    public boolean isPaired() {
        return descriptor.isPaired();
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
            if (nrUniqueLinesRead!=1) //only if nothing read or not at the first line
                lastLine= null;
		}
		
	}
	
	public HashSet<String> getRefIDset() {
		return refIDset;
	}

	public int getNrLines() {
		return nrUniqueLinesRead;
	}

	public int getCountMappings() {
		return countAll;
	}

	public int getCountContinuousMappings() {
		return countEntire;
	}

	public int getCountSplitMappings() {
		return countSplit;
	}

	public int getCountReads() {
		return countReads;
	}

	BufferedBACSReader readerC;
	protected BufferedBACSReader getReaderBACS() {
		if (readerC == null) {
			try {
				InputStream inputStream = new FileInputStream(inputFile);
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



	/**
     * Comparator for comparing read identifiers according to the provided descriptor.
     */
    BEDMappingComparator comp= null;

    /**
     * Returns an instance for comparing read identifiers according to the provided descriptor.
     * @return instance for comparing read identifiers according to the provided descriptor
     */
    private Comparator<? super Mapping> getDescriptorComparator() {
        if (comp == null) {
            comp = new BEDMappingComparator(
                    this.descriptor);
        }

        return comp;
    }

    /**
     * Creates a temporary file in the location provided, iff write access is
     * available there. Otherwise the file is created in the custom or system
     * temporary directory.
     *
     * @param location a file in the target directory or the directory itself,
     * may be <code>null</code>
     * @param name prefix of the file to be created, class name is appended
     * at the beginning
     * @param extension (optional) suffix of the temporary file that is created
     * @param deleteOnExit flag for calling the <code>deleteOnExit()</code>
     * method for the file
     * @return a temporary file according to the specifications
     */
    public File createTempFile(File location, String name, String extension, boolean deleteOnExit) {

        // get location
        if (location== null)
            location= this.tmpDir;
        else {
            if (!location.isDirectory())
                location= location.getParentFile();
            if (!location.canWrite())
                location= this.tmpDir;
        }

        // get name
        if (name== null)
            name= getClass().getSimpleName();
        else
            name= getClass().getSimpleName()+ "_"+ name;

        File f= null;
        try {
            f= FileHelper.createTempFile(name, extension, location);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return createFile(f, deleteOnExit);
    }

    /**
     * Control gateway for file creation from the main class,
     * adds a hook for delete on exit in case.
     *
     * @param f the file that has been created
     * @param deleteOnExit flag to mark for deletion on exit
     * @return
     */
    protected File createFile(File f, boolean deleteOnExit) {
        if (deleteOnExit)
            f.deleteOnExit();

        return f;
    }

    /**
     * Retrieves all mappings in a certain region from a BED input file.
     *
     * @param chromosome the locus for which reads are to be read
     * @param from start coordinate on chromosome
     * @param to end coordinate on chromosome
     * @return an iterator instance that enumerates all mappings in the specified region
     */
    public MSIterator read(String chromosome, int from, int to) {
        return read(chromosome, from, to, 0, 1);
    }

    /**
     * Out-of-memory-proof method that retrieves all mappings in a certain region from a BED input file.
     * The strategy is try-and-see, first it is attempted to try to load all requested reads into memory (RAM);
     * if latter attempt fails, the iterator is initialized on disk. Method retries if disk/filesystem blocks.
     * in the latter case.
     *
     * @param chromosome the locus for which reads are to be read
     * @param from start coordinate on chromosome
     * @param to end coordinate on chromosome
     * @param retryCount number of retries that are attempted in the case of disk/filesystem temporarily unreachable
     * @param timeInSeconds time between retries
     * @return an iterator instance that enumerates all mappings in the specified region
     */
    private MSIterator read(String chromosome, int from, int to, int retryCount, long timeInSeconds) {
        if (this.sortInRam) {
            try{
                return readRAM(chromosome, from, to);
            }catch (OutOfMemoryError memoryError){
                System.gc();
                Thread.yield();
                Log.warn("Not enough memory to sort BED entries in RAM. Switching to disk sorting. This run is NOT failed!\n " +
                        "You can increase the amount of memory used " +
                        "by the capacitor using the FLUX_MEM environment variable. For example: export FLUX_MEM=\"6G\"; flux-capacitor ... to use" +
                        "6 GB of memory.");
                return readDisk(chromosome, from, to, retryCount, timeInSeconds);
            }
        }else{
            return readDisk(chromosome, from, to, retryCount, timeInSeconds);
        }
    }


    /**
     * Loads all mappings in the respective region into RAM.
     * @param chromosome the locus for which reads are to be read
     * @param from start coordinate on chromosome
     * @param to end coordinate on chromosome
     * @return an iterator instance that enumerates elements of an array stored in RAM
     */
    private MSIterator<BEDMapping> readRAM(String chromosome, int from, int to) {

        if (from> to|| from< 0|| to< 0)
            throw new RuntimeException("BED reading range error: "+from+" -> "+to);
        // init iterator
        currentIter= null;
        // memory
        MappingReaderState state= readState(chromosome, from, to);
        if (state.result== null)
            return null;
        ArrayList<BEDMapping> beds= /*(BEDMapping[])*/ state.result;//TODO move to Mapping
        Collections.sort(beds, getDescriptorComparator());
        currentIter= new BEDMappingIterator(beds, this.descriptor);

        return currentIter;

    }

    /**
     * Writes all mappings in the respective region to disk, retries if disk/filesystem blocks.
     * @param chromosome the locus for which reads are to be read
     * @param from start coordinate on chromosome
     * @param to end coordinate on chromosome
     * @return an iterator instance that enumerates elements of an array stored in RAM
     */
    private MSIterator<BEDMapping> readDisk(String chromosome, int from, int to, int retryCount, long timeInSeconds) {

        if (from> to|| from< 0|| to< 0)
            throw new RuntimeException("BED reading range error: "+from+" -> "+to);

        // init iterator
        currentIter= null;

        try {
            // read, maintain main thread
            PipedInputStream  pin= new PipedInputStream();
            PipedOutputStream pout= new PipedOutputStream(pin);
            Comparator<CharSequence> c= new BEDDescriptorComparator(this.descriptor);
            File tmpFile= createTempFile(null, chromosome+ ":"+ from+ "-"+ to+ ".", "bed", true);
            BEDMappingIteratorDisk biter= new BEDMappingIteratorDisk(pin, tmpFile, c, this.descriptor);
            biter.init();
            currentIter= biter;
            MappingReaderState state= readState(pout, chromosome, from, to);
            pout.flush();
            pout.close();
            if (state.count== 0)
                return null;
        } catch (IOException e) {
            /*
             * "Resource temporarily unavailable"
             * Catch this exception and try again after sleeping for a while
             */
            if(e.getMessage().contains("Resource temporarily unavailable")){
                if(retryCount < 6){
                    Log.warn("Filesystem reports : 'Resource temporarily unavailable', I am retrying ("+(retryCount+1)+")");
                    try {Thread.sleep(1000 * (timeInSeconds));} catch (InterruptedException e1) {}
                    return readDisk(chromosome, from, to, retryCount + 1, timeInSeconds * 6);
                }
            }
            throw new RuntimeException(
                    "Could not get reads for locus "+ chromosome+ ":"+ from+ "-"+ to +", retried " + retryCount + " times", e);
        }

        return currentIter;

    }


	/**
	 * Reads BED objects from the underlying <code>InputStream</code> 
	 * overlapping the area specified by <code>chr</code>, 
	 * <code>start</code> and <code>end</code> and returns them as 
	 * an array.
	 * @param chr chromosome name of the specified area
	 * @param start start position of the specified area
	 * @param end end position of the specified area
	 */
	public MappingReaderState readState(String chr, int start, int end) {
		MappingReaderState state= new MappingReaderState();
		state.result= new ArrayList<BEDMapping>();
		state.count= 0;
		readState(chr, start, end, null, state);
		if (state.count== 0)
			state.result= null;
		else
			state.result= new ArrayList<BEDMapping>(state.result);
		return state;
	}
	
	/**
	 * Reads the BED lines from the underlying <code>InputStream</code> 
	 * overlapping the area specified by <code>chr</code>, 
	 * <code>start</code> and <code>end</code> into the 
	 * provided <code>OutputStream</code>.
	 * @param os stream to which the output is written
	 * @param chr chromosome name of the specified area
	 * @param start start position of the specified area
	 * @param end end position of the specified area
	 */
	public MappingReaderState readState(OutputStream os, String chr, int start, int end) {
		return readState(chr, start, end, os, null);
	}
	
	int skippedLines= 0;
	boolean warnFirstSkip= true;
	protected MappingReaderState readState(String chr, int start, int end, OutputStream os, MappingReaderState state) {
	
				if (state== null)
					state= new MappingReaderState();
				// else 
				state.count= 0l;
				state.state= MappingReaderState.STATE_OK;
				state.result= null;
				state.nextChr= null;
		
				if (bytesRead== 0) {
					warnFirstSkip= true;
					skippedLines= 0;
				}
		
				state.count= 0l;
				if (mapChr.containsKey(chr)) {
					if (mapChr.get(chr)== null) {
						state.state= MappingReaderState.STATE_CHROMOSOME_NOT_FOUND;
						return state;
					}
				}  
					
				--start;	// convert to bed coordinate
				try {
					BufferedBACSReader buffy= getReaderBACS();
					ByteArrayCharSequence cs= this.cs; // new ByteArrayCharSequence(100);
					String lastChrRead= null;
					getLineSeparator();
					boolean inited= false;
					if (reuse&& lastLine!= null) {
						cs= lastLine;
						lastLine= null;
						inited= true;
					}
		
					while (true) {
						
						long tmpBytes= bytesRead;			
		
						if (inited) {
							inited= false;
						} else {
							lastLine= cs.cloneCurrentSeq();
							if (buffy.readLine(cs)== null) {
								state.state= MappingReaderState.STATE_END_OF_FILE;
								return state;	// EOF
							}
							cs.resetFind();
							bytesRead+= cs.length()+ lineSeparator.length();
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
							int toks= cs.countTokens();
							if (identTok< 0)
								identTok= toks;
//							else
//								if (identTok!= toks&& false) {	// can be now, we read reads and split reads
//									if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
//										System.err.println("\t\n[OHLALA] line "+nrUniqueLinesRead+" has not "+identTok+" elements as the lines before!");
//										System.err.println("\t"+ cs.toString());
//										System.err.println("\tcheck file "+ getInputFile().getName());
//									}
//								}
							if (toks< 3) {
								++skippedLines;
								if (warnFirstSkip) {
									Log.warn("\t\n[OHNOO] line "+nrUniqueLinesRead+" has less than 3 token, I am skipping."
											+ "\t"+ cs.toString()+ "\tcheck file "+ getInputFile().getName());
									warnFirstSkip= false;
								}
								continue;
							}
							
							cs.resetFind();	// DEBUG 100827: somehow getToken(0) results in NullPointer otherwise 
							chrToki= cs.getToken(0).toString();	// getToken(1, TAB)
							bedStart= cs.getTokenInt(1); 	//BEDobject.encodeInt(cs.getToken(2, TAB)); 
							bedEnd= cs.getTokenInt(2);		//BEDobject.encodeInt(cs.getToken(3, TAB));
						} catch (Exception e) {
							e.printStackTrace();
						}

						// limit on objects that are read
						if (state.count>= maxBEDobjects) {
							state.nextChr= chrToki;
							if (reuse)
								lastLine= cs;
							else {
								bytesRead= tmpBytes; 
								--nrUniqueLinesRead;
							}
							return state;
						}
						

						if (chrToki.compareTo(chr)> 0) {	// end of chromosome reached
							if (reuse)
								lastLine= cs;
							else {
								bytesRead= tmpBytes; 
								--nrUniqueLinesRead;
							}
							
							addChr(chrToki, bytesRead- cs.length()- getLineSeparator().length(), nrUniqueLinesRead- 1);
							
							state.state= MappingReaderState.STATE_END_OF_CHROMOSOME;
							state.nextChr= chrToki;
							return state;
						}
							
						if (lastChrRead== null) {	// first line read in this batch
							
							addChr(chrToki, bytesRead- cs.length()- getLineSeparator().length(), nrUniqueLinesRead- 1);
							
							if (chr!= null&& !cs.subSequence(0, chr.length()).equals(chr)) {
								if (mapChr.containsKey(chr)) {
									if (mapChr.get(chr)== null) {
										if (reuse)
											lastLine= cs;
										else {
											bytesRead= tmpBytes;
											--nrUniqueLinesRead;
										}
										state.state= MappingReaderState.STATE_CHROMOSOME_NOT_FOUND;
										return state;
									} else {
										if (mapChr.get(chr)[0]> tmpBytes) { 	// only jump forward, never back
											long[] bytesNlines= mapChr.get(chr);
											reset(bytesNlines[0], (int) bytesNlines[1]);
											lastChrRead= chr.toString();
											buffy= getReaderBACS();
											if (buffy.readLine(cs)== null)
												break;
											//cs= new ByteArrayCharSequence(line);
											bytesRead+= cs.length()+ lineSeparator.length();
											++nrUniqueLinesRead;
											chrToki= cs.getToken(0).toString();
											bedStart= BEDobject.encodeInt(cs.getToken(1));
											bedEnd= BEDobject.encodeInt(cs.getToken(2));
										} else {
											reset(tmpBytes, nrUniqueLinesRead-1);
											state.state= MappingReaderState.STATE_CHROMOSOME_NOT_FOUND;	// ?
											state.nextChr= chrToki;
											return state;
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
										state.state= MappingReaderState.STATE_CHROMOSOME_NOT_FOUND;
										state.nextChr= chrToki;
										return state; 
									} else {
										
										this.cs= newCS;
										cs= newCS;
										chrToki= cs.getToken(0).toString();
										addChr(chrToki,
												bytesRead- cs.length()- getLineSeparator().length(), 
												nrUniqueLinesRead-1);
												
										buffy= getReaderBACS();
										lastChrRead= chr.toString();
										chrToki= cs.getToken(0).toString();
										bedStart= BEDobject.encodeInt(cs.getToken(1));
										bedEnd= BEDobject.encodeInt(cs.getToken(2));
									}
								}
							}
						}
		
						//line= line.trim();
						if (cs.startsWith("browser")|| cs.startsWith("track")|| cs.length()< 1) {
							++skippedLines;
							if (warnFirstSkip) {
								Log.warn("Skipped control line: "+ cs.toString());
								warnFirstSkip= false;
							}
							continue;
						}
						
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
								addChr(chrToki, bytesRead- cs.length()- getLineSeparator().length(), nrUniqueLinesRead);
								break;
							} 
						}
						
		
						boolean stop= false, continues= false;

                        //put the "=" to exclude reads on the boundaries
						if (start>= 0&& bedEnd<= start)
							continues= true;
						else if (end>= 0&& bedStart>= end)
							stop= true;

                        //take only reads entirely contained into the locus
                        if (!continues && !stop) {
                            if (bedStart<start||bedEnd>end)
                                continue;
                        }

						if (continues)
							continue;
						if (stop) {	// not found on this chr
							if (reuse)
								lastLine= cs;
							else {
								bytesRead= tmpBytes;
								--nrUniqueLinesRead;
							}
							state.state= MappingReaderState.STATE_OK; // ? more info ? out of range ?
							state.nextChr= chrToki;
							return state;
						}
						
						// create object
						if (os== null) {
							BEDMapping bed= new BEDMapping(cs, descriptor);
//							/*((Vector<Mapping>)*/ state.result.add(bed);      //using generic
//							++state.count;
                            state.addResult(bed);
						} else {
							os.write(cs.chars, cs.start, cs.length());
							os.write(Constants.NL);
							//os.flush();
							state.count+= cs.chars.length+ 1; 
						}
		
					}
					
				} catch (Exception e) {					
					e.printStackTrace();
				}

				return state;
			}

    //TODO check this method
	/*public ReadDescriptor checkReadDescriptor(boolean pairedEnd) {
		
		ReadDescriptor descriptor= null;
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(getInputFile()));
			
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
					descriptor= new UniversalReadDescriptor();
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
	}*/

	@Override
	public int getNrInvalidLines() {
		return skippedLines;
	}

	public HashMap<String, long[]> getMapChr() {
		return mapChr;
	}

	/**
	 * @deprecated maybe a bit more efficient, assumes tradtional BED
	 * sorting hierachy of fields 1,2,3
	 */
	private long isApplicableDefault(InputStream inputStream, long size) {
		try {
	        Log.progressStart("checking");
			BufferedReader buffy= new BufferedReader(new InputStreamReader(inputStream), 10* 1024* 1024);
			long rowCtr= 0;
			long bRead= 0;
			String lastC= null, nowC= null;
			int lastP= -1, nowP= -1;
			HashSet<String> setChr= new HashSet<String>();
			for (String s= buffy.readLine(); s!= null; s= buffy.readLine()) {
				++rowCtr;
				bRead+= s.length()+ 1;
	
				if (size> 0)
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
				//addRefID(nowC);
				while (i< len&& (c== ' '|| c== '\t'))
					c= s.charAt(++i);
				int p= i;
				for (;i < len; i++) {
					c= s.charAt(i);
					if (c== ' '|| c== '\t')
						break;
				}
				nowP= BEDobject.encodeInt(s, p, i);
				
				// avoid String.split(), not efficient
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
					Log.info("\n\tunsorted in line "+rowCtr+".");
					//clearRefIDs();
					buffy.close();
					return (-rowCtr);
				}
			}
			buffy.close();
	        Log.progressFinish(Constants.OK, true);
				
			return rowCtr;
		} catch (Exception e) {
	        Log.progressFailed(" ERROR.");
	        throw new RuntimeException(e);
		}
	}

	int maxBEDobjects= Integer.MAX_VALUE;
	
	public void setMaxBEDObjects(int i) {
		maxBEDobjects= i;
	}

	@Override
	public boolean isApplicable(UniversalReadDescriptor descriptor) {
		
		// TODO re-implement line-wise reading
		BufferedBACSReader buffy= null;
		try {
			buffy= new BufferedBACSReader(new FileInputStream(inputFile));
			ByteArrayCharSequence b= buffy.readLine(null);
			
			return descriptor.isApplicable(b);
			
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		} finally {
			if (buffy!= null)
				try {
					buffy.close();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
		}
	}

    @Override
    public Iterator<Mapping> iterator() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public class MappingReaderState {

        public static final int STATE_OK= 0;

        public static final int STATE_END_OF_CHROMOSOME= 1;

        public static final int STATE_END_OF_FILE= 2;

        public static final int STATE_CHROMOSOME_NOT_FOUND= 3;

        private long count= 0l;
        private ArrayList<BEDMapping> result= null;
        private byte state= STATE_OK;
        private String nextChr= null;

        public MappingReaderState() {
            reset();
        }

        public void reset() {
            count= 0l;
            result= null;
            state= STATE_OK;
            nextChr= null;
        }

        public void addResult(BEDMapping bed) {
            if (result == null)
                result = new ArrayList<BEDMapping>();
            if(scoreFilter < 0 || bed.getScore() >= scoreFilter){
                result.add(bed);
                count++;
            }
        }

        public ArrayList<BEDMapping> getResults() {
            return result;
        }

        public byte getState() {
            return state;
        }

        public String getNextChromosome() {
            return nextChr;
        }

        public void resetResults() {
            result = null;
            count = 0l;
        }
    }
}
