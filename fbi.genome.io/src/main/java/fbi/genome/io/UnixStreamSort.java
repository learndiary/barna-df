package fbi.genome.io;

//import io.gff.GTFSorter;
//import io.gff.ThreadedBufferedReader;
//import io.gff.GTFSorter.ByteArrayCharSequence;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Vector;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import commons.ByteArrayCharSequence;
import fbi.genome.model.commons.MyArrays;
import fbi.genome.model.commons.MyFile;
import fbi.genome.model.constants.Constants;

public class UnixStreamSort extends Thread {
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
			init(charSeq, fieldNrs, numeric);
		}
		
		public Comparable getValue(int fieldNr) {
			int p= Arrays.binarySearch(fIndex, fieldNr);	
//			if (p< 0)
//				System.currentTimeMillis();
			return fields[p];
		}
		public String toString() {
			return cs.toString();
		}
		
		void init(ByteArrayCharSequence input, int[] fieldNrs, boolean[] numeric) {
			 
			 int index = 0;
			 fields= new Comparable[fieldNrs.length];
			 fIndex= fieldNrs;
			 int mCtr= 0, fCtr= 0;;
			Matcher m= p.matcher(input);
	        while(m.find()) {
	        	++mCtr;
	        	if (fCtr< fieldNrs.length&& mCtr== fieldNrs[fCtr]) {
	        		ByteArrayCharSequence match = (ByteArrayCharSequence) input.subSequence(index, m.start());
	        		if (numeric[fCtr])
	        			fields[fCtr]= new Integer(match.parseInt());
	        		else
	        			fields[fCtr]= match;
	        		
	        		++fCtr;
	        	}
	            index = m.end();
	        }
	        // last
	        if (fCtr< fieldNrs.length&& fieldNrs[fCtr]== ++mCtr) {
        		ByteArrayCharSequence match = (ByteArrayCharSequence) input.subSequence(index, input.length());
        		if (numeric[fCtr])
        			fields[fCtr]= new Integer(match.parseInt());
        		else
        			fields[fCtr]= match;
	        }
		}
		
		
	}
	
	class ShutdownHooker extends Thread {
		@Override
		public void run() {
			for (int i = 0; i < tmpFileV.size(); i++) {
				tmpFileV.elementAt(i).delete();
			}
			System.err.println("Exited clean.");
		}
	}
	
	/**
	 * Gets a stream and sends it to nirvana..
	 * @author msammeth
	 *
	 */
	// this one is for you, Tobi
	static class DevNullReaderThread extends Thread {
		InputStream in= null;
		
		@Override
		public void run() {
			while (true) {
				try {
					int c= 0;
					for (int i = 0; i < in.available(); i++) {
						c= in.read();
						if (c== -1)
							break;
					}
					if (c== -1)
						break;
				} catch (Exception e) {
					e.printStackTrace();
				}
				try {
					sleep(500);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	
	class PriorityReaderThreadComparator implements Comparator<PriorityReaderThread2> {
		
		Comparator<DesignatedFieldsLine> baseComp= null;
		
		public PriorityReaderThreadComparator(Comparator<DesignatedFieldsLine> comp) {
			this.baseComp= comp;
		}
		
		//@Override
		public int compare(PriorityReaderThread2 o1, PriorityReaderThread2 o2) {
			
			int c= baseComp.compare(o1.peek(), o2.peek());
//			// chr1 + AA313056 @ 103827836 after AW583633 @ 103968366
//			String line1= o1.peek().toString();
//			String line2= o2.peek().toString();
//			if (line1.toString().contains("AA313056")||
//					line1.toString().contains("AW583633")||
//					line2.toString().contains("AA313056")||
//					line2.toString().contains("AW583633"))
//				System.currentTimeMillis();
			if (c== 0)
				return (o1.getPriority()- o2.getPriority()); 
			return c;
		}
	}
	

	
	class PriorityReaderThread2 extends Thread {
				
		
		public static final boolean debug= true;
		public static final int READ_LINES_AHEAD_LIMIT= 1000;
		int priority= -1;
		ThreadedBufferedByteArrayStream buffy= null;
		ArrayBlockingQueue<DesignatedFieldsLine> lineBuffer= null;
		//CountingConcurrentLinkedQueue lineBuffer= null;
		int readLinesAheadLimit= -1, minReadLinesAheadLimit= -1;
		long readBytesAheadLimit= -1, minReadBytesAheadLimit= -1;
		float loadFactor= 0.5f;
		long byteCtr= 0;
		File sortedFile= null;

		ByteArrayCharSequence cs;
		public PriorityReaderThread2(int priority, File file) {
			this(priority, file, 10000);
		}
		public PriorityReaderThread2(int priority, File file, int size) {
			super();
			sortedFile= file;
			this.priority= priority;			
			setName("Reader Thread, priority "+ priority+": "+file.getName());
			this.buffy= new ThreadedBufferedByteArrayStream(10000,sortedFile, true);
			//System.out.println("cap "+lineBuffer.remainingCapacity());
			cs= new ByteArrayCharSequence(10000);
			

			lineBuffer= new ArrayBlockingQueue<DesignatedFieldsLine>(10);
		}
		
		public boolean close() {
			return this.buffy.close();
		}
		
		public DesignatedFieldsLine remove() {
			
			if (lineBuffer.isEmpty()) 
				fill();
			
			if (lineBuffer.isEmpty())
				return null;

			return  lineBuffer.remove();
		}

		public DesignatedFieldsLine peek() {

			if (lineBuffer.isEmpty()) 
				fill();
			
			if (lineBuffer.isEmpty())
				return null;

			return  lineBuffer.peek();
		}
		
		private void fill() {
			try {
				ByteArrayCharSequence line = buffy.readLine(cs); 
				if (line.end!= 0)
					lineBuffer.add(new DesignatedFieldsLine(line.cloneCurrentSeq(), fNrs, fNum));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		

		public boolean isEmpty() {
			if (lineBuffer.isEmpty()) 
				fill();
			
			return lineBuffer.isEmpty();
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
	
	public class VectorQueue<E> extends ArrayBlockingQueue<E> {
		public VectorQueue(int capacity, int increment, boolean fair) {
			super(capacity, fair);
		}
	}
	
	public class CountingConcurrentLinkedQueue extends ConcurrentLinkedQueue<String> {
		int elementCount= 0;
		
		@Override
		public synchronized boolean add(String e) {
			boolean truhu= super.add(e);
			++elementCount;
			return truhu;
		}

		@Override
		public synchronized String remove() {
			String s= super.remove();
			--elementCount;
			return s;
		}
		
		@Override
		public synchronized int size() {
			return elementCount;
		}
	}
	
	int ctrTmpFiles= 0;
	public static String PFX_SORT= "sort";
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
				file= UnixStreamSort.this.createTempFile(Integer.toString(ctrTmpFiles++), null); 
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
			
			if (o1== null || o2== null) {
				return 0;
			}
			
			int val= 0;
			if (fieldNr>= 0) { 
				Comparable c1= o1.getValue(fieldNr), c2= o2.getValue(fieldNr);
				if (c1== null|| c2== null)
					return 0;
				val= c1.compareTo(c2);
			} else
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

	public UnixStreamSort(InputStream inStream) {
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
	
	public UnixStreamSort(InputStream inStream, int fieldNr, boolean numerical) {
		this(inStream);
		comp.setFieldNr(fieldNr);
		//comp.setNumerical(numerical);
	}
	
	public UnixStreamSort(InputStream inStream, DesignatedHierarchicalFieldComparator comp) {
		this(inStream);
		this.comp= comp;
		
		initFNum();
	}
	
	public UnixStreamSort(File file, DesignatedHierarchicalFieldComparator comp) {
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
		UnixStreamSort mySort2= new UnixStreamSort(file1, 1, false);
		mySort2.setSilent(false);
		mySort2.debug= true;
		mySort2.start();
		
		fbi.genome.model.commons.DevNullReaderThread dev0= new fbi.genome.model.commons.DevNullReaderThread();
		dev0.setIn(mySort2.getOutInStream());
		dev0.start();
	}
	
	File mergefps() {
		try {
			long t0= System.currentTimeMillis();
			if (!silent) {
				if (Constants.progress!= null) {
					Constants.progress.setString("merging");
					Constants.progress.setValue(0);
				} else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
					System.err.print("\tmerging ");
					System.err.flush();
				}

			}						
			// The head of this queue is the least element with respect to the specified ordering..
			PriorityReaderThreadComparator comp2= new PriorityReaderThreadComparator(comp);
			PriorityQueue<PriorityReaderThread2> queue= new PriorityQueue<PriorityReaderThread2>(tmpFileV.size(), comp2);
			for (int i = 0; i < tmpFileV.size(); i++) { 
				PriorityReaderThread2 prThread= new PriorityReaderThread2(i, tmpFileV.elementAt(i));
				//System.err.println(memSizeBytes+" "+tempFilesV.size()+" "+tempFilesV.elementAt(i).length());
//				prThread.setReadBytesAheadLimit(
//						Math.min(distributableMemSize/ tempFilesV.size(), tempFilesV.elementAt(i).length()));				
				//prThread.start();
				//prThread.waitForInit();
				queue.add(prThread);
			}
			
			BufferedOutputStream bossy= null;
			File outFile= null;
			if (baseFile!= null) {
				outFile= createFileSfxSorted(baseFile); // File.createTempFile(baseFile.getName(), "sorted");
				//tmpFiles.add(outFile);
				bossy= new BufferedOutputStream(new FileOutputStream(outFile));
			} else if (inStream!= null)
				bossy= new BufferedOutputStream(getOutputStream());
			int lastPerc= 0, linesRead= 0;
			char[] eol= System.getProperty("line.separator").toCharArray();
			while(!queue.isEmpty()) {
				if (!silent) {	// && size> 0
					int perc = (int) ((linesRead * 10d) / lineCnt);
					if (perc > lastPerc) { // ;) lastPerc!= 9&& 
						++lastPerc;
						if (!silent) {
							if (Constants.progress!= null)
								Constants.progress.progress();	// setValue(perc)						lastPerc= perc;
							else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
								System.err.print("*");
								System.err.flush();
							}
						}
					}
				}
				++linesRead;

				PriorityReaderThread2 prThread= queue.poll();
				DesignatedFieldsLine line= prThread.remove();
//				if (line.cs.end> line.cs.a.length)
//					System.currentTimeMillis();
				bossy.write(line.cs.a, line.cs.start, line.cs.length());
				for (int i = 0; i < eol.length; i++) 
					bossy.write((byte) eol[i]);
				if (prThread.isEmpty())
					prThread.close();
				else
					queue.add(prThread);
			}
			bossy.flush();
			bossy.close();
			Iterator<PriorityReaderThread2> iter= queue.iterator();
			while (iter.hasNext())
				iter.next().close();
			
			if (!silent) {
				if (Constants.progress!= null)
					Constants.progress.finish(Constants.OK, System.currentTimeMillis()- t0);
				else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
					if (lastPerc!= 9)	// ;)
						System.err.print("*");
					System.err.println(" "+(System.currentTimeMillis()- t0)/ 1000+" sec.");
				}
			}
			
			return outFile;
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return null;
	}

	public void run() {
		long t0= System.currentTimeMillis();
//		System.out.println(memSizeBytes+" bytes");
		
		// side FX
//		Thread t= new ShutdownHooker();		
//		t.setName("Shutdown Hooker");
//		Runtime.getRuntime().addShutdownHook(t);
		
		
		divideAndSort();
//		if (!silent)
//			System.out.println("dividing/1sorting "+ (System.currentTimeMillis()- t0)/1000 +" sec");
//		if (1== 1) {
//			for (int i = 0; i < tempFilesV.size(); i++) 
//				System.out.println(tempFilesV.elementAt(i).getName());
//			try {
//				System.in.read();
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		}

		long t1= System.currentTimeMillis();
		File f= mergefps();
//		if (!silent) {
//			System.out.println("merging "+ (System.currentTimeMillis()- t1)/1000 +" sec");
//			System.out.println("total "+ (System.currentTimeMillis()- t0)/1000 +" sec");
//			if (f!= null)
//				System.out.println(f.getAbsolutePath());
//		}

		int cnt= 0;
		for (int i = 0; i < tmpFileV.size(); i++) {
			if (tmpFileV.elementAt(i).delete())
				++cnt;
		}

		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP&& !silent) {
			System.err.println("\tcreated "+ctrTmpFiles+" temporary files, " +
					"removed "+cnt+" files, failed to remove "+(tmpFileV.size()- cnt)+" files.");
		}
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

	private String id= null;
	private UnixStreamSort() {
		setName("StreamSort "+hashCode());
		id= Constants.getTimestampID();
	}
	public UnixStreamSort(File file) {
		this();
		this.baseFile= file;
		this.size= file.length();
		this.comp= new DesignatedHierarchicalFieldComparator();
	}

	public UnixStreamSort(File file, int fieldNr, boolean numerical) {
		this(file);
		comp.setFieldNr(fieldNr);
		//comp.setNumerical(numerical);
	}

	private File createTempFile(String name, String ext) {
		File f= new File(System.getProperty(Constants.PROPERTY_TMPDIR)
				+ File.separator
				//+ (Constants.globalPfx== null? PFX_SORT: Constants.globalPfx+ "_"+ PFX_SORT)
				+ Constants.getGlobalPfx()+ "_"+ PFX_SORT
				+ "."
				//+ this.id
				+ Constants.getTimestampID()
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

	File mergefps_old() {
			try {
				long t0= System.currentTimeMillis();
				if (!silent) {
					if (Constants.progress!= null) {
						Constants.progress.setString("merging");
						Constants.progress.setValue(0);
					} else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						System.err.print("\tmerging ");
						System.err.flush();
					}
	
				}						
				// The head of this queue is the least element with respect to the specified ordering..
				PriorityReaderThreadComparator comp2= new PriorityReaderThreadComparator(comp);
				PriorityQueue<PriorityReaderThread2> queue= new PriorityQueue<PriorityReaderThread2>(tmpFileV.size(), comp2);
				for (int i = 0; i < tmpFileV.size(); i++) { 
					PriorityReaderThread2 prThread= new PriorityReaderThread2(i, tmpFileV.elementAt(i));
					//System.err.println(memSizeBytes+" "+tempFilesV.size()+" "+tempFilesV.elementAt(i).length());
	//				prThread.setReadBytesAheadLimit(
	//						Math.min(distributableMemSize/ tempFilesV.size(), tempFilesV.elementAt(i).length()));				
					//prThread.start();
					//prThread.waitForInit();
					queue.add(prThread);
				}
				
				BufferedOutputStream bossy= null;
				File outFile= null;
				if (baseFile!= null) {
					outFile= createFileSfxSorted(baseFile); // File.createTempFile(baseFile.getName(), "sorted");
					//tmpFiles.add(outFile);
					bossy= new BufferedOutputStream(new FileOutputStream(outFile));
				} else if (inStream!= null)
					bossy= new BufferedOutputStream(getOutputStream());
				int lastPerc= 0, linesRead= 0;
				char[] eol= System.getProperty("line.separator").toCharArray();
				while(!queue.isEmpty()) {
					if (!silent) {	// && size> 0
						int perc = (int) ((linesRead * 10d) / lineCnt);
						if (perc > lastPerc) { // ;) lastPerc!= 9&& 
							++lastPerc;
							if (!silent) {
								if (Constants.progress!= null)
									Constants.progress.progress();	// setValue(perc)						lastPerc= perc;
								else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
									System.err.print("*");
									System.err.flush();
								}
							}
						}
					}
					++linesRead;
	
					PriorityReaderThread2 prThread= queue.poll();
					DesignatedFieldsLine line= prThread.remove();
	//				if (line.cs.end> line.cs.a.length)
	//					System.currentTimeMillis();
					bossy.write(line.cs.a, line.cs.start, line.cs.length());
					for (int i = 0; i < eol.length; i++) 
						bossy.write((byte) eol[i]);
					if (prThread.isEmpty())
						prThread.close();
					else
						queue.add(prThread);
				}
				bossy.flush();
				bossy.close();
				Iterator<PriorityReaderThread2> iter= queue.iterator();
				while (iter.hasNext())
					iter.next().close();
				
				if (!silent) {
					if (Constants.progress!= null)
						Constants.progress.finish(Constants.OK, System.currentTimeMillis()- t0);
					else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						if (lastPerc!= 9)	// ;)
							System.err.print("*");
						System.err.println(" "+(System.currentTimeMillis()- t0)/ 1000+" sec.");
					}
				}
				
				return outFile;
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			return null;
		}
	
}
