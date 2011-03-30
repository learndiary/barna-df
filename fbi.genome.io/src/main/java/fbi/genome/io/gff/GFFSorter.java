package fbi.genome.io.gff;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.file.FileHelper;
import fbi.genome.io.UnixStreamSort;
import fbi.genome.io.UnixStreamSort.DesignatedHierarchicalFieldComparator;
import fbi.genome.model.commons.MyArrayHashMap;
import fbi.genome.model.commons.MyFile;
import fbi.genome.model.constants.Constants;
import fbi.genome.model.gff.GFFObject;

import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class GFFSorter {

	public static  class ThreadedBufferedByteArrayStream {
		static int btCtr= 0;
		BufferThread bt;
		boolean reader= true;

		byte[] cb;
		protected int nChars;
		protected int nextChar;
		/**
		 * If the next character is a line feed, skip it 
		 */
		protected boolean skipLF = false;
		
		class BufferThread extends Thread {
			byte[] buf;
			int len;
			File f;
			boolean stop= false;
			InputStream inStream;

			public BufferThread(int bufSize, File f) {
				this(bufSize);
				this.f= f;
			}
			
			private BufferThread(int bufSize) {
				setName("byte[] buffer - "+btCtr++);
				buf= new byte[bufSize];
				len= bufSize+1;	// ready to overwrite
			}
			
			public BufferThread(int bufSize, InputStream inStream) {
				this(bufSize);
				this.inStream= inStream;
			}
			

			
			public synchronized int read(byte[] tgt, int req) {
				try {
					while (len== 0|| len> buf.length) {					
						//System.out.println("not read: shouldnt happen, synchronized.");
							wait();
					}
				} catch (InterruptedException e) {
					return 0;
				}
				int nb= Math.min(req, len);
				if (nb== -1) 
					return nb;
				System.arraycopy(buf, 0, tgt, 0, nb);
				interrupt();
				len= buf.length+1;	// mark read
				return nb;
			}

			public synchronized int write(byte[] tgt, int req) {
				while (len<= buf.length) {					
					//System.out.println("not read: shouldnt happen, synchronized.");
					try {
						wait();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
				int nb= Math.min(req, len);
				if (nb== -1) 
					return nb;
				System.arraycopy(tgt, 0, buf, 0, nb);
				interrupt();
				len= buf.length+1;	// mark read
				return nb;
			}
			
			@Override
			public void run() {
				try {
					if (reader) {
						if (f!= null)
							inStream= new FileInputStream(f);
						int stat= 0;
						while (stat!= -1) {
							synchronized (this) {
								if (len> buf.length) {
									int inter= 0;
									while (inter>= 0)
										try {
											stat= inStream.read(buf, 0+ inter, buf.length- inter);
											inter= -1;
										} catch (InterruptedIOException ioe) {
											inter= ioe.bytesTransferred;
										}
									len= stat;
									notify();
								}
							}
							try {
								sleep(5000);
							} catch (InterruptedException e) {
								;	// :) 
							}
						}
						if (f!= null)
							inStream.close();
						
						
					} else {
						FileOutputStream stream=  new FileOutputStream(f);
						while (!stop) {
							synchronized (this) {
								if (len> 0) {
									stream.write(buf, 0, len);
									len= buf.length+ 1;
									notify();
								}
							}
							try {
								wait();
							} catch (InterruptedException e) {
								;	// :) 
							}
						}
						stream.flush();
						stream.close();

					}
					
				} catch (Exception e) {
					e.printStackTrace();
				}
			}

		}
		
		public void setStop(boolean stop) {
			bt.stop = stop;
			bt.interrupt();
			while(bt.isAlive())
			try {
				bt.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		public ThreadedBufferedByteArrayStream(int size, File f, boolean reader) {
			this.reader= reader;
			bt= new BufferThread(size, f);
			bt.start();
			this.cb= new byte[size];
			nChars= bt.read(cb, cb.length); // get more into buffer
			nextChar= 0;
		}
		
		public ThreadedBufferedByteArrayStream(int size, InputStream inStream, boolean reader) {
			this.reader= reader;
			bt= new BufferThread(size, inStream);
			bt.start();
			this.cb= new byte[size];
			nChars= bt.read(cb, cb.length); // get more into buffer
			nextChar= 0;
		}

		public ByteArrayCharSequence readLine(ByteArrayCharSequence cs) throws IOException {
			//StringBuffer s = null;
			int startChar;
			cs.end= cs.start;
			
//			synchronized (lock) {
//				ensureOpen();
//				boolean omitLF = skipLF; // ignoreLF || skipLF;
	
				bufferLoop: for (;;) {
	
					if (nextChar >= nChars) {
						nChars= bt.read(cb, cb.length); // get more into buffer
						nextChar= 0;
					}
					if (nextChar >= nChars) { /* EOF */
						return cs;
					}
	
					boolean eol = false;
					byte c = 0;
					int i;
	
					// Skip a leftover '\n', if necessary, 081114 has to be while for windows separator!!
					while (skipLF && (cb[nextChar] == '\n'|| cb[nextChar] == '\r')) {
						nextChar++;
						if (nextChar >= nChars) {
							nChars= bt.read(cb, cb.length); // get more into buffer
							nextChar= 0;
						}
						if (nextChar >= nChars) { /* EOF */
							return cs;
						}

					}
					skipLF = false;
//					omitLF = false;
	
					charLoop: for (i = nextChar; i < nChars; i++) {
						c= cb[i];
//						System.out.print(c+" ");
//						System.out.flush();
//						if (cs.end>= cs.a.length)
//							System.currentTimeMillis();
						if ((c == '\n') || (c == '\r')) {
							eol = true;
							break charLoop;
						}
						if (cs.end== cs.a.length) {
							cs.extend();
						}
						cs.a[cs.end++] = (byte) cb[i];
					}
	
					startChar = nextChar;
					nextChar = i;
	
					if (eol) {
//						nextChar++;
//						if (c == '\r') {
//							skipLF = true;
//						}
						skipLF= true;
						return cs;
					}
	
//					for (int j = 0; j < (i-startChar); j++) {
//						cs.a[cs.end++]= (byte) cb[startChar+j];
//					}
				}
//			}
		}

		public void writeLine(ByteArrayCharSequence cs) {
			;//
		}
		
	}

	public static boolean debug = false;

	public static int MAX_ID_LEN = 1000, MAX_POS_LEN = 24, MAX_CHR_LEN = 100;

	static boolean silent = false;

	/**
	 * Ein Thread, der eine InputStream ausliest und die erhaltenen Daten ins
	 * Nirvana schickt (bzw. nach /dev/null f??r den Linux-user ;-)). Eigentlich
	 * werden sie nur dem GarbageCollector ?berlassen. <br>
	 * <br>
	 * Wird ein Name f??r den Stream ??bergeben (z.B. "standard_out"), so wird
	 * der Inhalt des Streams unter Angabe des Namens nach System.out
	 * geschrieben.<br>
	 * Der Thread beendet sich selbst, sobald der stream geschlossen wird oder
	 * ein .interrupt() aufgerufen wird.<br>
	 * Erstellungsdatum: (03.04.01 21:55:05)
	 * 
	 * @author Tobias Vogele
	 */
	// Clip Out ColumnS
	public static class Cocs extends Thread {

		// needed for run()
		// inherited from interface Constants
		protected boolean DEBUG = false;

		/**
		 * The Inputstream, to be read.
		 */
		protected InputStream in;

		/**
		 * The Outputstream, to be written.
		 */
		private OutputStream out;

		String sepChar = "\t";

		int[] skipFields = new int[0];

		/**
		 * DevNullThread - constructs a new DevPipeReader with the given
		 * streams.
		 */
		public Cocs(InputStream in, OutputStream out) {
			this();
			this.in = in;
			this.out = out;
		}

		public Cocs() {
			setName("Cocs "+hashCode());
		}

		public InputStream getIn() {
			return in;
		}

		public void setIn(InputStream in) {
			this.in = in;
		}

		public OutputStream getOut() {
			return out;
		}

		/**
		 * The main method. Reads continously from the input stream(s) as long
		 * as there is something to read and interrupt() is not called
		 */
		public void run() {

			try {
				BufferedReader buffy = new BufferedReader(
						new InputStreamReader(in));
				BufferedWriter writer = new BufferedWriter(
						new OutputStreamWriter(out));
				String eol = System.getProperty("line.separator");
				char sepOut = ' ';
				if (sepChar.length() == 1)
					sepOut = sepChar.charAt(0);
				for (String line = buffy.readLine(); line != null; line = buffy
						.readLine()) {
					String[] tokens = line.split(sepChar);
					int ctr = 0;
					for (int i = 0; i < tokens.length; i++) {
						if (ctr < skipFields.length && skipFields[ctr] == i) {
							++ctr;
							continue;
						}
						writer.write(tokens[i]);
						if (i+ 1< tokens.length)
							writer.write(sepOut);
					}
					writer.write(eol);
				}
				buffy.close();
				writer.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		public String getSepChar() {
			return sepChar;
		}

		public void setSepChar(String sepChar) {
			this.sepChar = sepChar;
		}

		public int[] getSkipFields() {
			return skipFields;
		}

		public void setSkipFields(int[] skipFields) {
			for (int i = 0; i < skipFields.length; i++)
				--skipFields[i];
			Arrays.sort(skipFields);
			this.skipFields = skipFields;
		}

	}

	public class ArrayCharSequence implements CharSequence {
		public static final int defaultSize = 0;
	
		char[] a;
	
		int start, end;
	
		public ArrayCharSequence(char[] value) {
			a = value;
			this.start = 0;
			this.end = a.length;
		}
	
		public ArrayCharSequence(char[] value, int start, int end) {
			a = value;
			this.start = start;
			this.end = end;
		}
	
		//@Override
		public char charAt(int index) {
			return (char) a[start + index];
		}
	
		//@Override
		public int length() {
			return (end - start);
		}
	
		/**
		 * The resulting CharSequence operates on the same byte[] !
		 */
		//@Override
		public CharSequence subSequence(int start, int end) {
			// byte[] b= new byte[end-start];
			// System.arraycopy(a, start, b, 0, b.length);
			// return new ByteArrayCharSequence(b);
			return new ArrayCharSequence(a, this.start + start, this.start
					+ end);
	
		}
	
		public char[] getArray() {
			return a;
		}
	}

	public static void main(String[] args) {
		silent = false;
		File f= sort(new File(args[0]));
		System.out.println(f);
		try {
			GFFReader reader= new GFFReader(f.getAbsolutePath());
			reader.isApplicable();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	static byte[] encode(String tid, String chr, char mapKeySepChar) {
		byte[] key = new byte[tid.length() + 1 + chr.length()];
		for (int i = 0; i < tid.length(); i++)
			key[i] = (byte) tid.charAt(i);
		key[tid.length()] = (byte) mapKeySepChar;
		int offset = tid.length() + 1;
		for (int i = 0; i < chr.length(); i++)
			key[offset + i] = (byte) chr.charAt(i);

		return key;
	}

	static byte[] encode(ByteArrayCharSequence tid, ByteArrayCharSequence chr, char mapKeySepChar) {
		byte[] key = new byte[tid.length() + 1 + chr.length()];
		System.arraycopy(tid.a, tid.start, key, 0, tid.length());
		key[tid.length()] = (byte) mapKeySepChar;
		System.arraycopy(chr.a, chr.start, key, tid.length(), chr.length());
		return key;
	}

	static int find(Pattern p, CharSequence input, Pattern q) {
		 int index = 0;
		Matcher m= p.matcher(input);
	    for(int mCtr= 0;m.find();index= m.end()) {
	    	++mCtr;
	    	if (mCtr< 9) 
	    		continue;	    	
	    	CharSequence match = input.subSequence(index, m.start());
	    	Matcher m2= q.matcher(match);
	    	if (m2.find())
	    		return mCtr;
	    }
	    // cannot occur
	    // result[resCtr++]= input.subSequence(index, input.length());
	    
	    return -1;
	}

	static ByteArrayCharSequence[] find(Pattern p, CharSequence input, int[] fieldNrs) {
		 int index = 0;
		 ByteArrayCharSequence[] result= new ByteArrayCharSequence[fieldNrs.length];
		 int mCtr= 0, fCtr= 0;;
		Matcher m= p.matcher(input);
        while(m.find()) {
        	++mCtr;
//        	ByteArrayCharSequence s=(ByteArrayCharSequence) input.subSequence(index, m.start());
        	if (fCtr< fieldNrs.length&& mCtr== fieldNrs[fCtr]) {
        		ByteArrayCharSequence match = (ByteArrayCharSequence) input.subSequence(index, m.start());
        		result[fCtr++]= match;
        	}
            index = m.end();
        }
        // last
        if (fCtr< fieldNrs.length&& fieldNrs[fCtr]== ++mCtr) 
        	result[fCtr++]= (ByteArrayCharSequence) input.subSequence(index, input.length());
        
        return result;
	}
	
	
	public static File sort(File f) {
		try {
			String eol = FileHelper.guessFileSep(f);	//System.getProperty("line.separator");
			long t0 = System.currentTimeMillis();
			long size = f.length();
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[SORT] I am sorting the file now.");
			long estLineCount = (size / 100);
			int estIDCount = 0;
			if (estLineCount <= 50000) // reference annotation
				estIDCount = (int) (estLineCount / 40); // 10 exons per
														// transcript, 2 for CDS
														// and exon line
			else if (estLineCount <= 500000) // mRNA collection
				estIDCount = (int) (estLineCount / 6); // 7 exons per trpt, no
														// CDS
			else
				// ESTs, 25 mio lines -> 4 mio IDs
				estIDCount = (int) (estLineCount / 6);
			int hashMapSize = (int) (estIDCount);
			int incrSize = (int) (estIDCount * 0.3);
			MyArrayHashMap<byte[], Integer> map = new MyArrayHashMap<byte[], Integer>(
					hashMapSize); // for reference annot
			map.setIncrementSize(incrSize);
//			MyArrayHashMap<byte[], Integer> map2 = new MyArrayHashMap<byte[], Integer>(hashMapSize/2);
//			map2.setIncrementSize(incrSize);
//			if (!silent)
//				System.err.println("basehash created " + hashMapSize + ":"
//						+ incrSize);
//			if (debug)
//				System.in.read();
	
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.print("\tscanning IDs ");
				System.err.flush();
			}
			ThreadedBufferedByteArrayStream buffy = new ThreadedBufferedByteArrayStream(10000, f, true);
			int tidPos = -1;
			char mapKeySepChar = '_';
			long bytesRead = 0;
			int lastPerc = 0;
			long sumIDsize = 0, sumPossize = 0, sumCnt = 0;
			int lineCnt = 0;
			ByteArrayCharSequence cs= new ByteArrayCharSequence(1000);
			Pattern patty= Pattern.compile("\\s"), putty= Pattern.compile(GFFObject.TRANSCRIPT_ID_TAG), putty2= Pattern.compile(GFFObject.GENE_ID_TAG);			
//			String s= "chr1\thg18_intronEst\texon\t231869009\t231869096\t0.000000\t+\t.\ttranscript_id\t\"AA001030\";\tgene_id\t\"AA001030\"";
//			CharSequence[] r= find(patty, s, new int[] {1,4,7});
			// attentionAttention, the nrs have to be sorted - look in find()
			int[] fieldNrs= new int[]{1,4,7,-1};	// ,-1 for gene
			for (ByteArrayCharSequence line = buffy.readLine(cs); line.end!= 0; line = buffy
					.readLine(cs)) {
	
				++lineCnt;
//				if (lineCnt % 1000 == 0) {
//					System.gc(); // *** DEBUG ***
//					Thread.currentThread().yield();
//				}
	
				bytesRead += line.length() + eol.length();
				if (!silent) {
					int perc = (int) ((bytesRead * 10d) / size);
					if (perc > lastPerc && Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						++lastPerc;
						System.err.print("*");
						System.err.flush();
					}
				}
				
				if (cs.charAt(0)== '#')
					continue;
				
				if (fieldNrs[3]< 0) {
					fieldNrs[3]= find(patty, cs, putty)+1;
					if (fieldNrs[3]<1) {
						System.err.println("No transcript ID found.\n"+line);
						System.exit(1);
					}
				}
				// also gene, not needed
//				if (fieldNrs[3]< 0) {
//					fieldNrs[3]= find(patty, cs, putty2)+1;
//					if (fieldNrs[3]<1) {
//						System.err.println("No gene ID found.\n"+line);
//						System.exit(1);
//					}
//				}

				ByteArrayCharSequence[] fields= find(patty, cs, fieldNrs);	// 1,4,tid,gid
				boolean found= true;
				for (int i = 0; i < fields.length; i++) {
					if (fields[i]== null) {
						found= false;
						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
							System.err.println("\n\t[WARNING] I could not find field number "+fieldNrs[i]+" in line");
							System.err.println(cs.toString());
							System.err.println("\tline skipped check format of GTF file");
							System.err.println("\t(first 8 fields and transcript_id in the same column!)");
						}
						break;
					}
				}
				if (!found)
					System.exit(-1);
				
					// transcript start
				int start = fields[1].parseInt();
				byte[] key= encode(fields[3], fields[0], mapKeySepChar);
				byte strand= GFFObject.parseStrand(fields[2]);
				strand=(strand== 0)?1:strand;	// convert unknown to 1, for multiplication afterwards map.put() 
				
				Integer oldVal = map.get(key);
//				String DEBUG= key.toString();
				// 090901: consistency check of aligned transcripts, 
				// must be to the same strand on the same chromosome
				if (oldVal!= null&& !oldVal.equals(Integer.MIN_VALUE)) {
					if (strand* oldVal< 0) {
						oldVal= Integer.MIN_VALUE;
						map.put(key, oldVal);
						continue;
					}
				}
				if (oldVal == null || Math.abs(oldVal) > start) {
					map.put(key, (strand* start));
					if (oldVal == null) {
						sumPossize += fields[1].length();
						sumIDsize += fields[2].length();
						++sumCnt;
					}
				}
				// debug strand 0
//				int xxx= map.get(key);
//				System.currentTimeMillis();
//				// chr1 + AA313056 @ 103827836 after AW583633 @ 103968366
//				if (fields[3].toString().contains("AW583633"))
//					System.currentTimeMillis();
				
				// gene start, no we dont need that
/*				key= encode(fields[3], fields[0], mapKeySepChar);
				oldVal = map2.get(key);
				if (oldVal == null || oldVal > start) {
					map.put(key, start);
					if (oldVal == null) {
						sumPossize += fields[1].length();
						sumIDsize += fields[2].length();
						++sumCnt;
					}
				}*/
			}
			buffy.setStop(true);
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println(" " + (System.currentTimeMillis() - t0)
						/ 1000 + " sec.");
			
			System.gc(); 
			Thread.currentThread().yield();
			
			if (debug) {
				System.err.println("map attributes: " + map.size()
						+ " elements in " + map.getCapacity() + " table.");
				//System.in.read();
			}

			
			// tpos tid $1 $2 ...
			// sort chr, strand, tpos, tid, pos, (feat)
			// sort 3, 9, 1num, 2, 6num, (5)
//			DesignatedHierarchicalFieldComparator comp = new UnixStreamSort.DesignatedHierarchicalFieldComparator(
//					2)
//					.setSubComparator(new UnixStreamSort.DesignatedHierarchicalFieldComparator(
//							4)
//							.setSubComparator(new UnixStreamSort.DesignatedHierarchicalFieldComparator(
//									1)
//									.setSubComparator(new UnixStreamSort.DesignatedHierarchicalFieldComparator(
//											5)
//											.setSubComparator(new UnixStreamSort.DesignatedHierarchicalFieldComparator(
//													3)))));
//			fNrs= new int[] {1,2,5,8,tidField};
//			fNum= new boolean[] {true, false, true, false, false};
			DesignatedHierarchicalFieldComparator comp = new UnixStreamSort.DesignatedHierarchicalFieldComparator(
					2, false)	// chr
					.setSubComparator(new UnixStreamSort.DesignatedHierarchicalFieldComparator(
							8, false)	// strand 	
							.setSubComparator(new UnixStreamSort.DesignatedHierarchicalFieldComparator(
									1, true)	// tpos
									.setSubComparator(new UnixStreamSort.DesignatedHierarchicalFieldComparator(
											fieldNrs[3]+1, false)	// tid
											.setSubComparator(new UnixStreamSort.DesignatedHierarchicalFieldComparator(
													5, true)))));	// pos
			
	
			buffy = new ThreadedBufferedByteArrayStream(10000, f, true);			
			PipedOutputStream out = new PipedOutputStream();
			PipedInputStream in = new PipedInputStream(out);
			BufferedWriter writer= new BufferedWriter(new OutputStreamWriter(
					out));
			UnixStreamSort sorter = new UnixStreamSort(in, comp);
			sorter.setLineSep(eol);
			sorter.setTidField(fieldNrs[3]+1);	// for add col
			double avgIDsize = ((double) sumIDsize) / sumCnt, avgPosSize = ((double) sumPossize / sumCnt);
			if (debug)
				System.err.println("avg ID size " + avgIDsize
						+ ", avg Pos size " + avgPosSize);
			sorter.setSize(size
					+ (long) (map.size() * (avgIDsize + avgPosSize + 2)));
			sorter.setSilent(false);
			long sorterMemSize = (long) (Runtime.getRuntime().freeMemory() * 0.75);
			sorter.setMemSizeBytes(sorterMemSize);
			if (debug)
				System.err.println("sorter memory: " + sorterMemSize);
			sorter.setMultiProcessors(1);
			OutputStream outStr = System.err;
			File outFile = null;
			if (true) {
				outFile = new File(MyFile.append(f.getAbsolutePath(), "_sorted")); 
					//File.createTempFile(f.getName() + "_", "_sorted");
				outStr = new FileOutputStream(outFile);
			}
			Cocs pipe = new Cocs(sorter.getOutInStream(), outStr);
			pipe.setSkipFields(new int[] { 1});
			pipe.setSepChar("\t");
			pipe.start();
			sorter.start();
			
			
			fieldNrs= new int[]{fieldNrs[0], fieldNrs[3]};	// 090901: start pos from map
			HashSet<String> setInvalidTx= new HashSet<String>();
			int nrInvalidLines= 0;
			for (ByteArrayCharSequence line = buffy.readLine(cs); line.end!= 0; line = buffy
					.readLine(cs)) {
	
				bytesRead += line.length() + eol.length();
				if (!silent) {
					int perc = (int) ((bytesRead * 10d) / size);
					if (false&& perc > lastPerc && Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						++lastPerc;
						System.err.print("*");
						System.err.flush();
					}
				}
								
				ByteArrayCharSequence[] fields= find(patty, cs, fieldNrs);	// 1,4,-1
				
				byte[] key= encode(fields[1], fields[0], mapKeySepChar);
				int val= map.get(key);
				if (val== Integer.MIN_VALUE) {
					++nrInvalidLines;
					setInvalidTx.add(fields[1].toString());
					continue;
				} else
					System.currentTimeMillis();
					
//				if (fields[1].toString().contains("mod(mdg4)-RZ"))
//					System.currentTimeMillis();
				writer.write(Integer.toString(Math.abs(val))); // tpos
				writer.write("\t");
				writer.write(line.toString());
				writer.write(eol);
			}
			map= null;
			buffy.setStop(true);
			writer.flush();
			writer.close();
			System.gc();
			try {
				//pipe.join(); didnt work out
				sorter.join();	// wait 4 sorting
			} catch (InterruptedException e) {
				return null; 
			}
	
			
			if (nrInvalidLines> 0)
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
					System.err.println("\tskipped "+nrInvalidLines+" lines in "
							+ setInvalidTx.size()+ " transcripts");
					System.err.print("\t");
					Iterator<String> iter= setInvalidTx.iterator();
					while (iter.hasNext())
						System.err.print(iter.next()+ " ");
					System.err.println();
				}
			
			return outFile;
	
		} catch (Exception e) {
			e.printStackTrace();
		}
	
		return null;
	}

}
