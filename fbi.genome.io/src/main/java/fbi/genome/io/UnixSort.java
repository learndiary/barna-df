package fbi.genome.io;

import fbi.genome.io.UnixSort.PriorityReaderThread.PriorityReaderThreadComparator;
import fbi.genome.model.commons.MyArrays;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Vector;
import java.util.concurrent.ArrayBlockingQueue;

public class UnixSort {
	File baseFile= null;
	boolean stable= true;
	HierarchicalFieldComparator comp= null;
	long memSizeBytes= 1000000;
	boolean silent= true;
	private boolean debug= false;
	
	
	String lineSep= System.getProperty("line.separator");
	
	static class PriorityReaderThread extends Thread {
		
		static class PriorityReaderThreadComparator implements Comparator<PriorityReaderThread> {
			
			Comparator<String> baseComp= null;
			
			public PriorityReaderThreadComparator(Comparator<String> comp) {
				this.baseComp= comp;
			}
			
			//@Override
			public int compare(PriorityReaderThread o1, PriorityReaderThread o2) {
				
				int c= baseComp.compare(o1.peek(), o2.peek());
				if (c== 0)
					return (o1.getPriority()- o2.getPriority()); 
				return c;
			}
		}
		
		
		
		public static final boolean debug= true;
		public static final int READ_LINES_AHEAD_LIMIT= 1000;
		int priority= -1;
		BufferedReader buffy= null;
		ArrayBlockingQueue<String> lineBuffer= null;
		int readLinesAheadLimit= 1000;
		int minReadLinesAheadLimit= 0;
		float loadFactor= 0.5f;
		
		public PriorityReaderThread(int priority, File file) {
			super();
			this.priority= priority;			
			setName("Reader Thread, priority "+ priority+": "+file.getName());
			setReadLinesAheadLimit(READ_LINES_AHEAD_LIMIT);
			try {
				this.buffy= new BufferedReader(new FileReader(file));
				while (buffy.ready()&& lineBuffer.size()< readLinesAheadLimit) 
					lineBuffer.add(buffy.readLine());
					
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		public String remove() {
			
			if (lineBuffer.isEmpty()&& (emergency())) {
				if (debug)
					System.out.println("SOS, no data!");
				return null;
			}
			
			String line= lineBuffer.poll();
			
			if (lineBuffer.size()<= minReadLinesAheadLimit&& isReady())
				interrupt();
			
			return line;
		}

		private boolean emergency() {
			while (lineBuffer.isEmpty()) {
				if (isReady()) {
					interrupt();
					try {
						sleep(10);
					} catch (InterruptedException e) {
						if (debug)
							System.currentTimeMillis(); // :)
					}
//					if (debug&& lineBuffer.isEmpty())
//						System.out.println("w8ing for data "+ getName());
				} else {
					return true; // SOS
				}
			}
			return false;
		}
		
		public String peek() {

			if (lineBuffer.isEmpty()&& (emergency())) {
				if (debug)
					System.out.println("SOS, no data!");
				return null;
			}

			return  lineBuffer.peek();
		}
		

		public boolean isReady() {
			try {
				return 	buffy.ready();
			} catch (Exception e) {
				; // :)  java.io.IOException: Stream closed
			}
			return false;
		}
		
		public boolean isEmpty() {
			if (lineBuffer.isEmpty()) {
				return emergency();
			}
			return false;
		}
		
//		@Override
//		public void interrupt() {
//			if (!lineBuffer.isEmpty())
//				System.currentTimeMillis();
//			super.interrupt();
//		}
		
		@Override
		public void run() {
			
			try {
				int ctr= 0;
				while(buffy.ready()) {
					++ctr;
//					if (debug&& lineBuffer.size()== readLinesAheadLimit)
//						System.out.println("queue full");
					while (buffy.ready()&& lineBuffer.size()< readLinesAheadLimit) 
						lineBuffer.add(buffy.readLine());
					
					long t0= System.currentTimeMillis();
					while(buffy.ready()&& lineBuffer.size()> minReadLinesAheadLimit)	// sometimes called when buffer full..
						try {
							sleep(5000);
						} catch (InterruptedException ex) {
//							if (debug&& (!lineBuffer.isEmpty())) {
//								System.out.println("interrupted "+lineBuffer.size());
//								System.currentTimeMillis(); // :)
//							}
						}
//					System.out.println("interrupted "+lineBuffer.size()+" "+(System.currentTimeMillis()- t0)+ " x "+ctr);
				}
				buffy.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		public int getReadLinesAheadLimit() {
			return readLinesAheadLimit;
		}

		public void setReadLinesAheadLimit(int readLinesAheadLimit) {
			this.readLinesAheadLimit = readLinesAheadLimit;
			this.minReadLinesAheadLimit = (int) (readLinesAheadLimit* loadFactor);
			ArrayBlockingQueue<String> lineBufferNew= new ArrayBlockingQueue<String>(this.readLinesAheadLimit, true);
			while (lineBuffer!= null&& (!lineBuffer.isEmpty()))
				lineBufferNew.add(lineBufferNew.poll());
			this.lineBuffer= lineBufferNew;
		}
		
		
	}
	
	public class HierarchicalFieldComparator implements Comparator<String> {
		boolean numerical;
		int fieldNr= -1;
		HierarchicalFieldComparator subComparator= null;		
		
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
				}
				
				val= line1[fieldNr].compareTo(line2[fieldNr]);					
			}
			
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

		public void setSubComparator(HierarchicalFieldComparator subComparator) {
			this.subComparator = subComparator;
		}
	}
	
	public UnixSort(File file) {
		this.baseFile= file;
		this.comp= new HierarchicalFieldComparator();
	}
	
	public UnixSort(File file, int fieldNr, boolean numerical) {
		this(file);
		comp.setFieldNr(fieldNr);
		comp.setNumerical(numerical);
	}
	
	public static void main(String[] args) {
		long t0= System.currentTimeMillis();
		File file1= new File("/home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_UNSORT.gtf");
//		UnixSort mySort1= new UnixSort(baseFile, 4, true);
//		mySort1.setSilent(false);
//		File file1= mySort1.sort();
		UnixSort mySort2= new UnixSort(file1, 1, false);
		mySort2.setSilent(false);
		mySort2.debug= true;
		File file2= mySort2.sort();
		System.out.println((System.currentTimeMillis()- t0)/ 1000+ " sec");
		System.out.println(file2.getName());
	}
	
	void sort(File file) {
		String[] name= file.getName().split("@");
		int nbLines= new Integer(name[1]).intValue();
		String[] lines= new String[nbLines];
		int lineCtr= 0;
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(file));
			while (buffy.ready()) 
				lines[lineCtr++]= buffy.readLine();
			if (lineCtr!=  lines.length)
				System.err.println("error1");
			buffy.close();
			
			if (stable)
				Arrays.sort(lines, comp);
			else {
				System.out.println("TODO: Implement Quicksort");
				Arrays.sort(lines, comp);
			}
			
			BufferedWriter writer= new BufferedWriter(new FileWriter(file));
			for (int i = 0; i < lines.length; i++) {
				writer.write(lines[i]);
				writer.write(System.getProperty("line.separator"));
			}
			writer.flush();
			writer.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	File mergefps(File[] files) {
		try {
			// The head of this queue is the least element with respect to the specified ordering..
			PriorityReaderThreadComparator comp2= new PriorityReaderThreadComparator(comp);
			PriorityQueue<PriorityReaderThread> queue= new PriorityQueue<PriorityReaderThread>(files.length, comp2);
			for (int i = 0; i < files.length; i++) { 
				PriorityReaderThread prThread= new PriorityReaderThread(i, files[i]);
				prThread.start();
				queue.add(prThread);
			}
			
			File outFile= File.createTempFile(baseFile.getName(), "sorted");
			BufferedWriter writer= new BufferedWriter(new FileWriter(outFile));
			while(!queue.isEmpty()) {
				PriorityReaderThread prThread= queue.poll();
				writer.write(prThread.remove());
				writer.write(System.getProperty("line.separator"));
				if (!prThread.isEmpty())
					queue.add(prThread);
			}
			writer.flush();
			writer.close();
			
			
			return outFile;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return null;
	}

	public File sort() {
		long t0= System.currentTimeMillis();
		File[] files= divide(baseFile);
		if (!silent)
			System.out.println("dividing "+ (System.currentTimeMillis()- t0)/1000 +" sec");
//		if (debug) {
//			for (int i = 0; i < files.length; i++) 
//				System.out.println(files[i].getName());
//		}
		
//		t0= System.currentTimeMillis();
//		for (int i = 0; i < files.length; i++) 
//			sort(files[i]);
//		if (!silent)
//			System.out.println("sorting "+ (System.currentTimeMillis()- t0)/1000 +" sec");

		t0= System.currentTimeMillis();
		File f= mergefps(files);
		if (!silent)
			System.out.println("merging "+ (System.currentTimeMillis()- t0)/1000 +" sec");

		return f;
	}
	
	File[] divide(File file) {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(file));
			long bytes= 0;
			Vector<File> v= new Vector<File>();
			File outFile= File.createTempFile(baseFile.getName(), "sort");
			BufferedWriter writer= new BufferedWriter(new FileWriter(outFile));
			int lineCtr= 0;
			while (buffy.ready()) {
				String line= buffy.readLine();
				if (bytes+ line.length()+lineSep.length()> memSizeBytes) {
					writer.flush();
					writer.close();					
					File rFile= new File(outFile.getAbsolutePath()+"@"+lineCtr);
					outFile.renameTo(rFile);
					outFile= File.createTempFile(baseFile.getName(), "sort");
					v.add(rFile);
					writer= new BufferedWriter(new FileWriter(outFile));
					lineCtr= 0;
					bytes= 0;
				}
				bytes+= line.length()+lineSep.length();
				writer.write(line);
				writer.write(lineSep);
				++lineCtr;
			}
			buffy.close();
			writer.flush();
			writer.close();
			File rFile= new File(outFile.getAbsolutePath()+"@"+lineCtr);
			outFile.renameTo(rFile);
			v.add(rFile);
			
			return (File[]) MyArrays.toField(v);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return null;
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
	
}
