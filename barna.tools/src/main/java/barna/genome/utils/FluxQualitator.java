package barna.genome.utils;


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


//import io.ThreadedBufferedByteArrayStream;

import java.io.*;


public class FluxQualitator extends Thread {
	
	public static int OFFSET_SOLEXA= 64, OFFSET_PHRED= 33,
		MIN_PHRED= 0, MAX_SOLEXA= 40, MAX_PHRED= 93;
	
	static String CLI_SHORT_IS_SOLEXA= "issolexa",
		CLI_SHORT_IS_PHRED= "isphred",
		CLI_SHORT_IS_NUMERIC= "isnumeric",
		CLI_SHORT_TO_SOLEXA= "tosolexa",
		CLI_SHORT_TO_PHRED= "tophred";
	
	static File fileIn= null;	
	static PipedInputStream pipeIn;
	static PipedOutputStream pipeOut;
	
	boolean check= false, convert= false;
	static int isNumeric= 0, isQuality= 0, toNumeric= -1, toQuality= -1;	// solexaOld= 1, solexaNew= 2, phred= -1
	static long ctrLines= 0, ctrQual= 0;
	
	static void parse(String[] args) {
		for (int i = 0; i < args.length; i++) {
			boolean cliLong= false;
			String s= null;
			if (args[i].startsWith("--")) {
				s= args[i].substring(2);
				cliLong=true;
			} else if (args[i].startsWith("-")) 
				s= args[i].substring(1);
			if (s== null) {
				fileIn= new File(args[i]);
				continue;
			}
			
			s= s.toLowerCase();
			if (s.equals(CLI_SHORT_IS_SOLEXA))
				isQuality= 1;
			else if (s.equals(CLI_SHORT_IS_PHRED))
				isQuality= -1;
			else if (s.equals(CLI_SHORT_IS_NUMERIC))
				isNumeric= 1;
			else if (s.equals(CLI_SHORT_TO_SOLEXA))
				toQuality= 1;
			else if (s.equals(CLI_SHORT_TO_PHRED))
				toQuality= -1;
		}
		
	}
	
	static boolean doCheck= false, doConvert= false;
	static PipedInputStream dummyIn;
	public static void main(String[] args) {
		
		if (1== 2) {	// debug stdin
			dummyIn= new PipedInputStream();
			System.setIn(dummyIn);
			Thread t= new Thread(){
				@Override
				public void run() {
					try {
						PipedOutputStream pipo= new PipedOutputStream(dummyIn);
						BufferedReader buffy= new BufferedReader(new FileReader("resources\\formats\\080815_S2-DRSC-Untreated-1_s_1_sequence.txt"));
						byte[] bb= new byte[1024];
						char[] cc= new char[1024];
						int len= -1;
						while ((len= buffy.read(cc))!= -1) {
							for (int i = 0; i < len; i++) 
								bb[i]= (byte) cc[i];
							pipo.write(bb, 0, len);
						}
						
					} catch (Exception e) {
						e.printStackTrace();
					}
					
				}
			};
			t.start();
			try {
				Thread.currentThread().sleep(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		long t0= System.currentTimeMillis();
		
		parse(args);
		doCheck= (isNumeric== 0|| isQuality== 0);
		doConvert= (toNumeric!= 0&& toQuality!= 0&& (isNumeric!= toNumeric|| isQuality!= toQuality));
		
		System.err.println("\nHowdy, I am the FluxQualitator!\n");
		int ready= 0;
		if (fileIn== null)
			try {
				ready= System.in.available();
			} catch (IOException ex) {
			}
		if (fileIn== null&& ready== 0) {
			System.err.println("This is how I work:");
			System.err.println("Input comes via stdin or from a file.");
			System.err.println("Output comes on stdout");
			System.err.println("Optional parameters");
			System.err.println("-isPhred\tinput are Phred qualities");
			System.err.println("-isSolexa\tinput are Solexa qualities");
			System.err.println("-isNumeric\tinput qualities are given as numbers");
			System.err.println("\t\t\trather than ASCII symbols (default)");
			//System.err.println("-toPhred\toutput are Phred qualities");
			System.err.println("-toSolexa\toutput are Solexa qualities (default is Phred)");
			System.err.println("\nIf the input parameters are incomplete, I will try to guess them.");
			System.err.println("\nHave a continued nice day.");
			System.exit(1);
		}
		
		if (!(doCheck|| doConvert)) {
			if (toQuality== 0|| toNumeric== 0)
				System.err.println("no output format specified");
			if (isNumeric== toNumeric&& isQuality== toQuality)
				System.err.println("input format is the same as desired output");
			System.err.println("There is nothing I can do for you.");
			System.exit(1);
		}
		
		System.err.println("Let's see what's up:");
		System.err.println("You told me to ");
		if (doCheck)
			System.err.println("\tfind out what format the input is");
		else {
			System.err.print("\ttake a file with ");
			if (isNumeric> 0)
				System.err.print("numeric ");
			else if (isNumeric< 0)
				System.err.print("ASCII ");
			if (isQuality> 0)
				System.err.print("Solexa ");
			else if (isQuality< 0)
				System.err.print("Phred ");
			System.err.println("qualities");
		}
		if (doConvert) {
			System.err.print("\tand convert it to ");
			if (toNumeric> 0)
				System.err.print("numeric ");
			else if (toNumeric< 0)
				System.err.print("ASCII ");
			if (toQuality> 0)
				System.err.print("Solexa ");
			else
				System.err.print("Phred ");
			System.err.println("qualities");
		}
		
		FluxQualitator checker= null, converter= null;
		if (doCheck) {
			if (doConvert) {
				pipeOut= new PipedOutputStream();
				try {
					pipeIn= new PipedInputStream(pipeOut);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			checker= new FluxQualitator();
			checker.check= true;
			checker.start();
		} 
		if (doConvert) {
			converter= new FluxQualitator();
			converter.convert= true;			
			converter.start();
			try {
				converter.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			System.err.println("I converted "+ctrLines+" lines, "+ctrQual+" qualities\nand finished flux in "+(System.currentTimeMillis()- t0)/ 1000+" sec.");
			System.err.println("\tCheers!");
		}
		
	}
	
	boolean solexa2phred= false, numeric2phred= false, numeric2solexa= false;	
	
	@Override
	public void run() {
		if (check)
			check();
		else
			convert();
	}
	
	public static int solexa2phredQuality(int qualSol) {
		// 10*log(10**(solexa_quality/10.0) + 1, 10)

		int q= (int) (10 * 
				Math.log(1 + Math.pow(10, (qualSol / 10d))) 
						/ Math.log(10));
		assert(q>= 0);
		return q;
	}
	
	public static int phred2solexaQuality(int qualPhred) {
		//10*log(10**(phred_quality/10.0) - 1, 10)
		int q= (int) (10 * 
				Math.log(Math.pow(10, (qualPhred / 10d))- 1)
						/ Math.log(10));
		return q;
	}
	
	public static byte solexa2phredASCII(byte asciiSol) {
		int q= solexa2phredQuality(asciiSol- OFFSET_SOLEXA);
		return (byte) (q+ OFFSET_PHRED);
	}
	
	public static byte phred2solexaASCII(byte asciiSol) {
		int q= phred2solexaQuality(asciiSol- OFFSET_PHRED);
		return (byte) (q+ OFFSET_SOLEXA);
	}
	
	
	static final byte sixtyfour= (byte) 64, thirtythree= (byte) 33;
	
	void convert() {
		try {			
			int bufSize= 1024, maxReadLen= 1000, p= 0, pOrd= 0, qStart= -1;
			boolean preQual= false, lastCR= false;
			byte[] buf= new byte[bufSize], ord= null, ords= null;

			while(isNumeric== 0|| isQuality== 0)
				try {
					int avail= System.in.available();					
					Thread.currentThread().sleep(100);
					if (isNumeric== 0&& isQuality== 0&& fileIn== null&& avail== 0) {
						System.err.println("[FATAL] I failed to determine the input format");
						System.err.println("\tcannot perform conversion");
						System.exit(1);
					}
				} catch (InterruptedException e) {
					; // :)
				}

			System.err.println("[VEGETARIAN] Converting.. ");
			System.err.println("\tto format:\t"+ (toNumeric>0?"numeric":"ascii"));
			System.err.println("\tto base:\t"+ (toQuality>0?"solexa":"phred"));
			if (isNumeric== 1&& toNumeric== -1) {
				ord= new byte[3];
				ords= new byte[maxReadLen];				
			}

			BufferedInputStream in= new BufferedInputStream(getInputStream(), bufSize);
			BufferedOutputStream buffy= new BufferedOutputStream(System.out, bufSize); 
			for (int len= 0; (len= in.read(buf))!= -1; qStart=Math.min(qStart, 0)) {
				for (int i = 0; i < len; ++i) {
					boolean cr= buf[i]== '\r'|| buf[i]== '\n';
					if (cr) {
						if (lastCR) {
							if (qStart>= 0)
								++qStart;
							continue;
						}
						++ctrLines;
					}
					if (qStart>= 0) {
						if (isNumeric== -1) {
							assert(isQuality!= toQuality);
							if (cr)
								qStart= -1;
							else { 
								buf[i]= (toQuality== -1?solexa2phredASCII(buf[i])
										:phred2solexaASCII(buf[i]));
								++ctrQual;
							}
						} else { // isNumeric== 1
							if (cr|| buf[i]== ' ') {
								if (p> 0) {
									byte asc= (byte) (parseByte(ord, p)+ (isQuality> 0?sixtyfour:thirtythree));									
									ords[pOrd++]= (toQuality> 0? phred2solexaASCII(asc):solexa2phredASCII(asc));
									++ctrQual;
								}
								p= 0;								
								if (cr) {
									len-= i;
									System.arraycopy(buf, i, buf, 0, len);
									buffy.write(ords, 0, pOrd);
									//buffy.flush();
									pOrd= 0;
									qStart= -1;
									i= 0;
								}
							} else 
								ord[p++]= (byte) ((buf[i]>= 48&& buf[i]<= 57)?buf[i]- 48:buf[i]);
							
						
						} 
						
					} else if (preQual&& cr) {
						if (isNumeric== 1) {
							int wLen= i+ 1;
							buffy.write(buf, 0, wLen);
							//buffy.flush();
							len-= wLen;
							System.arraycopy(buf, wLen, buf, 0, len);
							i= -1;
							qStart= 0;
						} else
							qStart= i;
						preQual= false;
					} else if (lastCR&& buf[i]== '+')
						preQual= true;

					lastCR= cr;
				}
				if (!(isNumeric== 1&& qStart>= 0)) {
					buffy.write(buf, 0, len);
					//buffy.flush();
				}
			}
			in.close();
			buffy.flush();
			buffy.close();
				
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public InputStream getInputStream() {
		if (fileIn== null) {
			if (check|| (doConvert&& !doCheck)) {
				return System.in;
			} else {
				return pipeIn;
			}
		}
		try {
			return new FileInputStream(fileIn);
		} catch (Exception e) {
			return null;
		}
	}
	
	int min= Integer.MAX_VALUE, max= Integer.MIN_VALUE;
	private int setMinMaxQual(int x) {
		if (x> max) {
			max= x;
			if (max> MAX_SOLEXA)
				return -1;
		}
		if (x< min) {
			min= x;
			if (min< MIN_PHRED)
				return 1;
		}
		return 0;
	}
	
	// check quality checking:
	// file with IIIIIII estimated solexa1.0, not phred!
	void check() {
		try {
			System.err.println("[WHATSUP] Checking format.. ");
			BufferedReader buffy= new BufferedReader(new InputStreamReader(getInputStream()));
			BufferedWriter bos= null;
			if (pipeOut!= null)
				bos= new BufferedWriter(new OutputStreamWriter(pipeOut));
			String s;
			boolean qual= false, finished= false;
			int quals= 0, maxQuals= 1000;
			while ((!finished)&& (s= buffy.readLine())!= null) {
				if (s.charAt(0)== '+')
					qual= true;
				else if (qual) {
					if (isNumeric>= 0) {
						String[] ss= s.split("\\s");
						if (ss.length== 1)
							isNumeric= -1;
						else 
							isNumeric= 1;
						for (int i = 0; isQuality== 0&& isNumeric>= 0&& i < ss.length; i++) {
							if (ss[i].equals(""))
								continue;
							try {
								int x= Integer.parseInt(ss[i]);
								isQuality= setMinMaxQual(x);
								if (x< MIN_PHRED)
									isQuality= 1;
								if (x> MAX_SOLEXA)
									isQuality= -1;
								++quals;
							} catch (Exception e) {
								isNumeric= -1;
							}
						}
					}
					if (isNumeric<= 0) {
						for (int i = 0; i < s.length(); i++) {
							int x= s.charAt(i);
							setMinMaxQual(x);
							if (x- OFFSET_PHRED< MIN_PHRED)
								isQuality= 1;
							if (x- OFFSET_SOLEXA> MAX_SOLEXA) 
								isQuality= -1;
							++quals;
						}
					}
					qual= false;
				}
				if (bos!= null) {
					bos.write(s+ System.getProperty("line.separator"));
					if (quals>= maxQuals)
						finished= true;
				}
				if (isNumeric!= 0&& isQuality!= 0)
					finished= true;
			}
			//buffy.close(); // System.in
			
			assert(isNumeric!= 0);
			if (isQuality== 0) {
				if (isNumeric> 0) {
					if (max< MAX_SOLEXA)
						isQuality= 1;
					else
						isQuality= -1;
				} else {
					if ((max- min)> MAX_SOLEXA)
						isQuality= -1;
					else
						isQuality= 1;
				}
				if (max- OFFSET_SOLEXA== MAX_SOLEXA) 
					isQuality= 1;
				if (max- OFFSET_PHRED== MAX_PHRED)
					isQuality= -1;
			}
			if (isQuality== 1&& (min- OFFSET_SOLEXA>= MIN_PHRED))
				isQuality= 2;

			System.err.println("\tis format:\t"+(isNumeric==0?"undecided":(isNumeric>0?"numeric":"ascii")));
			System.err.println("\tis base:\t"+(isQuality==0?"undecided":(isQuality>0?(isQuality==1?"solexa1.0":"solexa1.3"):"phred")));
			
			if (finished&& bos!= null) 
				while ((s= buffy.readLine())!= null) {
					bos.write(s+ System.getProperty("line.separator"));
				}
			if (bos!= null) {
				bos.flush();
				bos.close();
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	byte parseByte(byte[] a, int len) {
		byte r= 0;
		for (int i = 0; i < len; i++) {
			if (a[len- 1- i]== '-')
				r*= -1;
			else
				r+= a[len- 1- i]* Math.pow(10, i);
		}
		return r;
	}
	


	
	
}
