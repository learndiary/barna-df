package fbi.genome.io;

//import io.gff.GTFSorter;
//import io.gff.ThreadedBufferedReader;
//import io.gff.GTFSorter.ByteArrayCharSequence;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.file.FileHelper;
import fbi.genome.model.commons.MyArrays;
import fbi.genome.model.commons.MyFile;
import fbi.genome.model.constants.Constants;

import java.io.*;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Vector;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.regex.Pattern;

public class UnixStreamSort2 extends Thread {
	public static final long MAX_FILE_SIZE_BYTES= 50000000;
	public static final String SFX_SORTED= "_sorted";
	File baseFile= null;
	boolean stable= true;
	DesignatedHierarchicalFieldComparator comp= null;
	long memSizeBytes= 70000000;		// Runtime.maxMemory()
	boolean silent= true;
	int tidField= -1;
	int[] fNrs;
	boolean[] fNum;

	private boolean debug= false;
	int multiProcessors= 4;
	Vector<File> tmpFileV= new Vector<File>();
	int avgLineLen= 40;
	long size= -1;
	private long lineCnt= 0;
	
	BufferedInputStream inStream= null;
	PipedOutputStream outStream= null;
	PipedInputStream outInStream= null;
	
	
	String lineSep= System.getProperty("line.separator");
	
	public static File createFileSfxSorted(File original) {
		return new File(MyFile.append(original.getAbsolutePath(), SFX_SORTED));
	}
	
	public static File createFileSfxSorted(String path) {
		return new File(MyFile.append(path, SFX_SORTED));
	}
	
	static class DesignatedFieldsLine {
		Comparable[] fields;
		int[] fIndex;
		ByteArrayCharSequence cs;
		static Pattern p= Pattern.compile("\\s");
		
		public DesignatedFieldsLine(ByteArrayCharSequence charSeq, int[] fieldNrs, boolean[] numeric) {
			this.cs= charSeq;
			// sort?
			for (int i = 1; i < fieldNrs.length; i++) {
				if (fieldNrs[i-1]< fieldNrs[i]) {
					int[] fieldNrSorted= fieldNrs.clone();
					Arrays.sort(fieldNrSorted);
					boolean[] numericSorted= new boolean[numeric.length];
					for (int j = 0; j < numericSorted.length; j++) {
						int p= Arrays.binarySearch(fieldNrSorted, fieldNrs[j]);
						numericSorted[p]= numeric[i];
					}
					numeric= numericSorted;
					fieldNrs= fieldNrSorted;
					break;
				}
			}
			
			init(charSeq, fieldNrs, numeric);
		}
		
		public Comparable getValue(int fieldNr) {
			int p= Arrays.binarySearch(fIndex, fieldNr);	
			return fields[p];
		}
		public String toString() {
			return cs.toString();
		}
		
		void init(ByteArrayCharSequence input, int[] fieldNrs, boolean[] numeric) {
			 
			fIndex= fieldNrs;
			fields= new Comparable[fieldNrs.length];
			
			// get values
			cs.resetFind();
			for (int i = 0; i < fieldNrs.length; i++) 
				fields[i]= numeric[i]? cs.getTokenInt(fieldNrs[i]- 1): cs.getToken(fieldNrs[i]- 1);
		}
				
	}
	
	public static class HierarchicalFieldComparator implements Comparator<String> {
		boolean numerical;
		int fieldNr= -1;
		HierarchicalFieldComparator subComparator= null;		

		public HierarchicalFieldComparator() {
			super();
		}
		
		public HierarchicalFieldComparator(int field, boolean numerical) {
			this();
			setFieldNr(field);
			setNumerical(numerical);
		}
		
		//@Override
		public int compare(String o1, String o2) {
			
			int val= 0;
			if (fieldNr>= 0) {
				String[] line1= o1.split("\\s");
				String[] line2= o2.split("\\s");
				if (numerical) {
					try {
						val= new Integer(line1[fieldNr]).compareTo(new Integer(line2[fieldNr]));
					} catch (Exception e) {
						try {
							val= new Double(line1[fieldNr]).compareTo(new Double(line2[fieldNr]));
						} catch (Exception ex) {
							System.err.println("field "+fieldNr+" not numeric as declared "+line1[fieldNr]+","+line2[fieldNr]); 
							return 0;
						}
					}
				} else
					val= line1[fieldNr].compareTo(line2[fieldNr]);					
			
			} else
				val= o1.compareTo(o2);
			
			if (val== 0 && subComparator!= null)
				return subComparator.compare(o1, o2);
			return val;
		}

		public boolean isNumerical() {
			return numerical;
		}

		public void setNumerical(boolean numerical) {
			this.numerical = numerical;
		}

		public int getFieldNr() {
			return (fieldNr+ 1);
		}

		public void setFieldNr(int fieldNr) {
			this.fieldNr = fieldNr-1;
		}

		public HierarchicalFieldComparator getSubComparator() {
			return subComparator;
		}

		public HierarchicalFieldComparator setSubComparator(HierarchicalFieldComparator subComparator) {
			this.subComparator = subComparator;
			return this;
		}
	}
	
	int ctrTmpFiles= 0;
	public static String PFX_SORT= "sort";
	public static final String MERGER_TID= "Merger-";
	static int mergerCtr= 0;
	
	/**
	 * constant mem, therefore only 2 files
	 * @author msammeth
	 *
	 */
	class MergerThread extends Thread {
		
		File inFileA= null, inFileB= null;
		File outFile= null;		
		
		public MergerThread(File oneFile) {
			outFile= oneFile;
		}
		public MergerThread(File inFileA, File inFileB) {
			super(MERGER_TID+ (++mergerCtr));
			this.inFileA= inFileA;
			this.inFileB= inFileB;
		}
		
		DesignatedFieldsLine getLine(InputStream istream, BufferedByteArrayReader reader, ByteArrayCharSequence cs) {
			cs= reader.readLine(istream, cs);
			if (cs== null)
				return null;
			return new DesignatedFieldsLine(cs, fNrs, fNum);
		}
		
		
		@Override
		public void run() {
			
			if (inFileA== null|| inFileB== null)
				return;
			
			BufferedByteArrayReader buffyA= new BufferedByteArrayReader(),
									buffyB= new BufferedByteArrayReader();
			try {
				outFile= File.createTempFile("sort", "merged");
				//System.err.println("created "+ outFile.getName());
				BufferedOutputStream ostream= new BufferedOutputStream(new FileOutputStream(outFile));
				InputStream istreamA= new FileInputStream(inFileA), 
					istreamB= new FileInputStream(inFileB);
				ByteArrayCharSequence csA= new ByteArrayCharSequence(1024),
										csB= new ByteArrayCharSequence(1024);
				DesignatedFieldsLine lineA= getLine(istreamA, buffyA, csA), 
									lineB= getLine(istreamB, buffyB, csB);
				int ctrA= 0, ctrB= 0;
				while (lineA!= null&& lineB!= null) {
					int val= comp.compare(lineA, lineB);
					if (val< 0) {
						ostream.write(lineA.cs.a, lineA.cs.start, lineA.cs.length());
						ostream.write('\n');
						lineA= getLine(istreamA, buffyA, csA);
						++ctrA;
					} else if (val> 0) {
						ostream.write(lineB.cs.a, lineB.cs.start, lineB.cs.length());
						ostream.write('\n');
						lineB= getLine(istreamB, buffyB, csB);
						++ctrB;
					} else {
						ostream.write(lineA.cs.a, lineA.cs.start, lineA.cs.length());
						ostream.write('\n');
						ostream.write(lineB.cs.a, lineB.cs.start, lineB.cs.length());
						ostream.write('\n');
						lineA= getLine(istreamA, buffyA, csA);
						lineB= getLine(istreamB, buffyB, csB);
						++ctrA;
						++ctrB;
					}
				}
				
//				if (ctrA== 0|| ctrB== 0)
//					System.currentTimeMillis();
				
				if (lineA!= null) {
					ostream.write(lineA.cs.a, lineA.cs.start, lineA.cs.length());
					ostream.write('\n');
					while ((csA= buffyA.readLine(istreamA, csA))!= null) {
						ostream.write(csA.a, csA.start, csA.length());
						ostream.write('\n');
					}
				}
				istreamA.close();
				if (!inFileA.delete())
					inFileA.delete();
//				else
//					System.err.println("deleted "+ inFileA.getName());
				
				if (lineB!= null) {
					ostream.write(lineB.cs.a, lineB.cs.start, lineB.cs.length());
					ostream.write('\n');
					while ((csB= buffyB.readLine(istreamB, csB))!= null) {
						ostream.write(csB.a, csB.start, csB.length());
						ostream.write('\n');
					}
				}
				istreamB.close();
				if (!inFileB.delete())
					inFileB.delete();
//				else
//					System.err.println("deleted "+ inFileB.getName());
				
				ostream.flush();
				ostream.close();
				
//				if (outFile.length()== 0)
//					System.currentTimeMillis();
				
			} catch (Exception e) {
				e.printStackTrace();
			}
			
		}
		
		public File getOutFile() {
			while (isAlive())
				try {
					join();
				} catch (InterruptedException e) {
					; // :)
				}
			return outFile;
		}
	}
	
	class SortThread extends Thread {
		DesignatedFieldsLine[] lines;
		File file;
		
		@Override
		public void run() {
			try {
				if (stable)
					Arrays.sort(lines, comp);	
				else {
					System.out.println("TODO: Implement Quicksort");
					Arrays.sort(lines, comp);	
				}
				
				// File.createTempFile("sort_", "_"+lines.length);
				file= UnixStreamSort2.this.createTempFile(Integer.toString(ctrTmpFiles++), null); 
				tmpFileV.add(file);
				
//				try {
//					file.deleteOnExit();
//				} catch (NullPointerException e) {
//					; // http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6526376
//				}
				BufferedOutputStream bossy= new BufferedOutputStream(new FileOutputStream(file));
				char[] eol= System.getProperty("line.separator").toCharArray();
				for (int i = 0; i < lines.length; i++) {
					bossy.write(lines[i].cs.a, lines[i].cs.start, lines[i].cs.length());
					for (int j = 0; j < eol.length; j++) 
						bossy.write((byte) eol[j]);
				}
				bossy.flush();
				bossy.close();
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	public static class DesignatedHierarchicalFieldComparator implements Comparator<DesignatedFieldsLine> {
		int fieldNr= -1;
		boolean numeric= false;
		DesignatedHierarchicalFieldComparator subComparator= null;		
	
		public DesignatedHierarchicalFieldComparator() {
			super();
		}
		
		public DesignatedHierarchicalFieldComparator(int field) {
			this();
			setFieldNr(field);
		}
		
		public DesignatedHierarchicalFieldComparator(int field, boolean numeric) {
			this(field);
			this.numeric= numeric;
		}
		
		//@Override
		public int compare(DesignatedFieldsLine o1, DesignatedFieldsLine o2) {
			
			int val= 0;
			if (fieldNr>= 0) 
				val= o1.getValue(fieldNr).compareTo(o2.getValue(fieldNr));
			else
				val= o1.toString().compareTo(o2.toString());
			
			if (val== 0 && subComparator!= null)
				return subComparator.compare(o1, o2);
//			if (o1.fields[0].toString().startsWith("chr6_")^o2.fields[0].toString().startsWith("chr6_"))
//				System.out.println(o1.fields[0]+"\n"+o2.fields[0]+"\n"+val);
			return val;
		}
	
		public int getFieldNr() {
			return (fieldNr);
		}
	
		public void setFieldNr(int fieldNr) {
			this.fieldNr = fieldNr;
		}
	
		public DesignatedHierarchicalFieldComparator getSubComparator() {
			return subComparator;
		}
	
		public DesignatedHierarchicalFieldComparator setSubComparator(DesignatedHierarchicalFieldComparator subComparator) {
			this.subComparator = subComparator;
			return this;
		}
	}

	public UnixStreamSort2(InputStream inStream) {
		this();
		this.inStream= new BufferedInputStream(inStream);
		this.comp= new DesignatedHierarchicalFieldComparator();
	}
	
	public PipedOutputStream getOutputStream() {
		if (outStream == null) {
			connectPipedOutput();
		}
		return outStream;
	}
	
	public PipedInputStream getOutInStream() {
		if (outInStream == null) {
			connectPipedOutput();
		}
		return outInStream;
	}
	
	private void connectPipedOutput() {
		outStream = new PipedOutputStream();
		try {
			outInStream= new PipedInputStream(outStream);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public UnixStreamSort2(InputStream inStream, int fieldNr, boolean numerical) {
		this(inStream);
		comp.setFieldNr(fieldNr);
		//comp.setNumerical(numerical);
	}
	
	public UnixStreamSort2(InputStream inStream, DesignatedHierarchicalFieldComparator comp) {
		this(inStream);
		this.comp= comp;
		
		initFNum();
	}
	
	public UnixStreamSort2(File file, DesignatedHierarchicalFieldComparator comp) {
		this(file);
		this.comp= comp;
		
		initFNum();
	}
	
	void initFNum() {
		int cnt= 1;
		DesignatedHierarchicalFieldComparator cc= comp.getSubComparator();
		for (;cc!= null;cc= cc.getSubComparator()) {
			++cnt;
		}
	
//		fNrs= new int[] {1,2,5,8,tidField};
//		fNum= new boolean[] {true, false, true, false, false};
		
		fNrs= new int[cnt];	// TODO 081217 we have to sort it here
		for (int i = 0; i < fNrs.length; i++) 
			fNrs[i]= 0;
		fNum= new boolean[cnt];
		cnt= 0;
		cc= comp;
		for (;cc!= null;++cnt,cc=cc.getSubComparator()) {
//			int p= Arrays.binarySearch(fNrs, 0, cnt, cc.fieldNr);	// TODO jdk16
//			assert(p<0);
//			p= -(p+1);
			int p= 0;
			for (; p < cnt; p++) 
				if (fNrs[p]> cc.fieldNr)
					break;			
			if (p< cnt) {
				System.arraycopy(fNrs, p, fNrs, p+1, cnt-p);
				System.arraycopy(fNum, p, fNum, p+1, cnt-p);
			}
			fNrs[p]= cc.fieldNr;
			fNum[p]= cc.numeric;
		}
	}

	
	public static void main(String[] args) {
		long t0= System.currentTimeMillis();
		File file1= new File("/home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_UNSORT.gtf");
//		UnixSort mySort1= new UnixSort(baseFile, 4, true);
//		mySort1.setSilent(false);
//		File file1= mySort1.sort();
		UnixStreamSort2 mySort2= new UnixStreamSort2(file1, 1, false);
		mySort2.setSilent(false);
		mySort2.debug= true;
		mySort2.start();
		
		fbi.genome.model.commons.DevNullReaderThread dev0= new fbi.genome.model.commons.DevNullReaderThread();
		dev0.setIn(mySort2.getOutInStream());
		dev0.start();
	}
	
	BufferedByteArrayReader buffy;
	class ThreadQueue extends Thread {
		ArrayBlockingQueue<MergerThread> inq= new ArrayBlockingQueue<MergerThread>(10);
		ThreadQueue outQT;
		boolean end= false;
		int ctr= 0, nr= 0;
		
		public ThreadQueue(int nr) {
			super("ThreadQueue-"+nr);
			this.nr= nr;
		}
		
		@Override
		public void run() {
			while (!end) {
				synchronized (inq) {
					while (inq.size()< 2&& !end)
						try {
							inq.wait();
						} catch (InterruptedException e) {
							; // :)
						}
					if (inq.size()>= 2) {
						
						if (outQT== null) {
							outQT= new ThreadQueue(nr+1);
							outQT.start();
						}
						MergerThread t1= inq.poll(), t2= inq.poll();
						while (t1.isAlive()|| t2.isAlive())
							try {
								t1.join();
								t2.join();
							} catch (InterruptedException e1) {
								;	// :)
							}
						MergerThread t= new MergerThread(t1.getOutFile(), t2.getOutFile());
						++ctr;
						t.start();
						boolean succeed= false;
						while (!succeed)
							try {
								synchronized (outQT.inq) {
									outQT.inq.put(t);
									outQT.inq.notifyAll();
								}
								succeed= true;
							} catch (InterruptedException e) {
								; // :)
							}
					}
				}
			}
			if (inq.size()> 0) {
				assert(inq.size()== 1);	// Assertion error, s. getFinalFile()
				if (ctr> 0) {
					if (outQT== null) {
						outQT= new ThreadQueue(nr+1);
						outQT.start();
					}
					MergerThread t= new MergerThread(inq.poll().getOutFile());
					++ctr;
					t.start();
					boolean succeed= false;
					while (!succeed)
						try {
							synchronized (outQT.inq) {
								outQT.inq.put(t);
								outQT.inq.notifyAll();
							}
							succeed= true;
						} catch (InterruptedException e) {
							; // :)
						}
				}
			}
			if (outQT!= null) {
				outQT.end= true;
				outQT.interrupt();
			}
			
		}
		
		public File getFinalFile() {
			while(isAlive())
				try {
					join();
				} catch (InterruptedException e) {
					; // :)
				}
			if (ctr== 1) {
				MergerThread t= outQT.inq.peek();	// otherwise Assertionerror in ThreadQueue.run()
				while (t.isAlive())
					try {
						t.join();
					} catch (InterruptedException e) {
						; // :)
					}
				return t.outFile;
			} else {
				if (outQT== null) {
					//System.err.println(ctr+ " "+inq.size());
					return inq.peek().getOutFile();
				}
				return outQT.getFinalFile();
			}
		}
	}
	
	File outFile= null;
	public File getOutFile() {
		return outFile;
	}

	public void setOutFile(File outFile) {
		this.outFile = outFile;
	}

	public void run() {
		long t0= System.currentTimeMillis();
		buffy= new BufferedByteArrayReader();
		ThreadQueue tq= new ThreadQueue(1);
		tq.start();
		try {
			InputStream iStream= baseFile== null? inStream: new FileInputStream(baseFile);
			File x;
			while ((x= nextSortedFile(iStream))!= null) {
				MergerThread t= new MergerThread(x);
				t.start();
				boolean succeed= false;
				while (!succeed)
					try {
						synchronized (tq.inq) {
							tq.inq.put(t);
							tq.inq.notifyAll();
						}
						succeed= true;
					} catch (InterruptedException e) {
						; // :)
					}
			}
			tq.end= true;
			tq.interrupt();

			// stream
			outFile= tq.getFinalFile();
			if (baseFile!= null) {
				File newOutFile= createFileSfxSorted(baseFile); // File.createTempFile(baseFile.getName(), "sorted");
				if (FileHelper.move(outFile, newOutFile, null))
					outFile= newOutFile;
			} else if (inStream!= null) {
				BufferedOutputStream bossy= new BufferedOutputStream(getOutputStream());
				BufferedInputStream fin= new BufferedInputStream(new FileInputStream(outFile));
				int n= 0;
				byte[] b= new byte[1024];
				try {
					while ((n= fin.read(b))> 0) 
						bossy.write(b, 0, n);
					fin.close();
					if (!outFile.delete())
						outFile.delete();
					bossy.flush();
					bossy.close();
				} catch (Exception e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
						e.printStackTrace();
				}
			}
			
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				e.printStackTrace();
			return;
		}
		
//		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP&& !silent) {
//			System.err.println("\tcreated "+ctrTmpFiles+" temporary files, " +
//					"removed "+cnt+" files, failed to remove "+(tmpFileV.size()- cnt)+" files.");
//		}
	}
	
	public void setSilent(boolean silent) {
		this.silent = silent;
	}

	public boolean isStable() {
		return stable;
	}

	public void setStable(boolean stable) {
		this.stable = stable;
		if (!stable)
			System.out.println("TODO: Implement Quicksort");
	}

	
	long maxSortedFileSize= 100000;
	File nextSortedFile(InputStream iStream) {
		try {
			long t0= System.currentTimeMillis();
			if (!silent) {
				if (Constants.progress!= null) {
					Constants.progress.setString("dividing");
					Constants.progress.setValue(0);
				} else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
					System.err.print("\tdividing ");
					System.err.flush();
				}
			}

			long bytes= 0;
			ByteArrayCharSequence buf= new ByteArrayCharSequence(1024);
			DesignatedFieldsLine[] v= new DesignatedFieldsLine[(int) (maxSortedFileSize/ buf.a.length)];			
			int lctr= 0, incr= Math.max(10, (v.length/ 10));
			while(bytes<  maxSortedFileSize&& (buf=buffy.readLine(iStream, buf))!= null) {
				bytes+= buf.length();
				if (lctr== v.length) {
					DesignatedFieldsLine[] w= new DesignatedFieldsLine[v.length+ incr];
					System.arraycopy(v, 0, w, 0, v.length);
					v= w;
				}
				v[lctr++]= new DesignatedFieldsLine(buf.cloneCurrentSeq(), fNrs, fNum); 
			}
			if (buf== null&& lctr== 0)
				return null;
			if (stable)
				Arrays.sort(v, 0, lctr, comp);	
			else {
				System.out.println("TODO: Implement Quicksort");
				Arrays.sort(v, 0, lctr, comp);	
			}
			File out= File.createTempFile("sort", "sorted");
			BufferedOutputStream ostream= new BufferedOutputStream(new FileOutputStream(out));
			for (int i = 0; i < lctr; i++) {
				ostream.write(v[i].cs.a, v[i].cs.start, v[i].cs.length());
				ostream.write('\n');
			}
			ostream.flush();
			ostream.close();
			return out;
			
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				e.printStackTrace();
			return null;
		}
		//return null;
	}

	private String id= null;
	private UnixStreamSort2() {
		setName("StreamSort "+hashCode());
		id= Constants.getTimestampID();
	}
	public UnixStreamSort2(File file) {
		this();
		this.baseFile= file;
		this.size= file.length();
		this.comp= new DesignatedHierarchicalFieldComparator();
	}

	public UnixStreamSort2(File file, int fieldNr, boolean numerical) {
		this(file);
		comp.setFieldNr(fieldNr);
		//comp.setNumerical(numerical);
	}

	private File createTempFile(String name, String ext) {
		File f= new File(System.getProperty(Constants.PROPERTY_TMPDIR)
				+ File.separator
				+ (Constants.globalPfx== null? PFX_SORT: Constants.globalPfx+ "_"+ PFX_SORT)
				+ "."
				+ this.id
				+ "."+ name
				+ (ext== null?"": "."+ ext));
		f.deleteOnExit();
		return f;
	}
	public long getSize() {
		return size;
	}

	public void setSize(long size) {
		this.size = size;
	}

	public long getMemSizeBytes() {
		return memSizeBytes;
	}

	public void setMemSizeBytes(long memSizeBytes) {
		this.memSizeBytes = memSizeBytes;
	}

	public int getMultiProcessors() {
		return multiProcessors;
	}

	public void setMultiProcessors(int multiProcessors) {
		this.multiProcessors = multiProcessors;
	}

	public int getTidField() {
		return tidField;
	}

	public void setTidField(int tidField) {
		this.tidField = tidField;
	}

	public String getLineSep() {
		return lineSep;
	}

	public void setLineSep(String lineSep) {
		this.lineSep = lineSep;
	}

	void divideAndSort() {
			try {
				long t0= System.currentTimeMillis();
				if (!silent) {
					if (Constants.progress!= null) {
						Constants.progress.setString("dividing");
						Constants.progress.setValue(0);
					} else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						System.err.print("\tdividing ");
						System.err.flush();
					}
				}
	
				ThreadedBufferedByteArrayStream buffy= null;
				if (baseFile!= null) {
					buffy= new ThreadedBufferedByteArrayStream(10000, baseFile, true);
					id= baseFile.getName();
				} else if (inStream!= null) {
					buffy= new ThreadedBufferedByteArrayStream(10000, inStream, true);
					id= "stream";
				}
				
				long bytes= 0;
				Vector<SortThread> vt= new Vector<SortThread>(multiProcessors);
				Vector<DesignatedFieldsLine> lineV= null;
				int lineCtr= 0, batchCtr= 0;
				long sharedMemSize= memSizeBytes/ multiProcessors;				
				long bytesRead= 0;
				int lastPerc= 0;
				ByteArrayCharSequence cs= new ByteArrayCharSequence(10000);
				
				long maxChunkSizeBytes= Math.min(sharedMemSize, MAX_FILE_SIZE_BYTES);
				if(debug)
					System.err.println("parallel "+multiProcessors+" x "+maxChunkSizeBytes);
				for (ByteArrayCharSequence line = buffy.readLine(cs); line.end!= 0; line = buffy
					.readLine(cs)) {
					
					bytesRead+= line.length()+ lineSep.length();
					++lineCnt;
					if (size> 0) {
						int perc = (int) ((bytesRead * 10d) / size);
						if (lastPerc!= 9&& perc > lastPerc) {	// ;)
							++lastPerc;
							if (!silent) {
								if (Constants.progress!= null)
									Constants.progress.progress();	// setValue(perc)
								else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
									System.err.print("*");
									System.err.flush();
								}
							}
						}
					}
					
					long currMem= Runtime.getRuntime().totalMemory(), currFreeMem= Runtime.getRuntime().freeMemory(); 
					if (lineV!= null&& (bytes+ line.length()+lineSep.length()> maxChunkSizeBytes||
							(currFreeMem< (currMem* 0.25)))) {	// 090522: deliberative mem bound, maxMem not well suited as it updates and size is not exactly related to Xmx
	
						if (debug)
							System.err.println("sorter started "+currFreeMem+" / "+currMem);
						// substart sort
						batchCtr++;
						while (vt.size()== multiProcessors) {
							vt.elementAt(0).join();
							//tmpFileV.add(vt.elementAt(0).file);
							vt.remove(0);
						}
						SortThread t= new SortThread();
						t.setName("Sorter "+batchCtr);
						t.lines= (DesignatedFieldsLine[]) MyArrays.toField(lineV.toArray());
						//t.start();
						t.run();	// // 090522: bound memory !!!!
						vt.add(t);
						
						lineCtr= 0;
						bytes= 0;
						lineV.removeAllElements();
						System.gc();
					}
					bytes+= line.length()+lineSep.length();
					if (lineV== null) {
						avgLineLen= line.length();
						//System.err.println("avg line len "+avgLineLen);
						lineV= new Vector<DesignatedFieldsLine>((int) ((float) maxChunkSizeBytes/ (avgLineLen+ lineSep.length())),
								10);
					}				
					lineV.add(new DesignatedFieldsLine(line.cloneCurrentSeq(), fNrs, fNum)); 
					++lineCtr;
				}
	//			buffy.close();
				
				if (lineV.size()> 0) {
					// substart sort
					batchCtr++;
					while (vt.size()== multiProcessors) {
						vt.elementAt(0).join();
						//tmpFileV.add(vt.elementAt(0).file);
						vt.remove(0);
					}
					SortThread t= new SortThread();
					t.setName("Sorter "+batchCtr);
					t.lines= (DesignatedFieldsLine[]) MyArrays.toField(lineV.toArray());
					t.start();		
					vt.add(t);			
				}
				
				while (!vt.isEmpty()) {
					vt.elementAt(0).join();
					//tmpFileV.add(vt.elementAt(0).file);
					vt.remove(0);
				}
				
				if (!silent) {
					if (Constants.progress!= null)
						Constants.progress.finish(Constants.OK, System.currentTimeMillis()- t0);
					else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						if (lastPerc< 9)
							System.err.print("*");
						System.err.println(" "+(System.currentTimeMillis()- t0)/ 1000+" sec.");
					}
				}
		
			} catch (Exception e) {
				e.printStackTrace();
			}
				
		}
	
}
