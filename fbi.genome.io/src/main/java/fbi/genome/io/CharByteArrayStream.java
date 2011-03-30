package fbi.genome.io;

import fbi.commons.ByteArrayCharSequence;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;


public class CharByteArrayStream {
		static int btCtr= 0;
		BufferThread bt;
		boolean reader= true, closeInputStream= false;

		byte[] cb;
		protected int nChars;
		protected int nextChar;
		/**
		 * If the next character is a line feed, skip it 
		 */
		protected boolean skipLF = false;
		public ByteArrayCharSequence lastLine= null;
		
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
				while (len== 0|| len> buf.length) {					
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
						while (stat!= -1&& !stop) {
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
								sleep(100);
							} catch (InterruptedException e) {
								;	// :) 
							}
						}
						if (f!= null|| closeInputStream)
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
			System.gc();
			while(bt.isAlive())
				try {
					bt.interrupt();
					bt.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
		}
		
		public CharByteArrayStream(int size, File f, boolean reader) {
			this.reader= reader;
			bt= new BufferThread(size, f);
			bt.start();
			this.cb= new byte[size];
			nChars= bt.read(cb, cb.length); // get more into buffer
			nextChar= 0;
		}
		
		public CharByteArrayStream(int size, InputStream inStream, boolean reader) {
			this.reader= reader;
			bt= new BufferThread(size, inStream);
			bt.start();
			this.cb= new byte[size];
			nChars= bt.read(cb, cb.length); // get more into buffer
			nextChar= 0;
		}

		
		
		public ByteArrayCharSequence readLine(ByteArrayCharSequence cs) throws IOException {
			
			final ReadableByteChannel inputChannel = Channels.newChannel(System.in);
			ByteBuffer buf= ByteBuffer.allocateDirect(16 * 1024);			
			
			if (lastLine!= null) {
				ByteArrayCharSequence l= lastLine;
				lastLine= null;
				return l;
			}
			
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
						if (cs.end>= cs.a.length)
							cs.extend();
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
		
		public boolean close() {
			setStop(true);
			return true;
		}

		public boolean isCloseInputStream() {
			return closeInputStream;
		}

		public void setCloseInputStream(boolean closeInputStream) {
			this.closeInputStream = closeInputStream;
		}
		
	}

