package fbi.commons.file;


import fbi.commons.Log;
import fbi.commons.StringUtils;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.nio.channels.WritableByteChannel;
import java.util.*;
import java.util.zip.*;

/**
 * File Utilities
 *
 */
public class FileHelper {
    /**
     * Indicates file compression
     */
	public static final byte COMPRESSION_NONE= 0, COMPRESSION_ZIP= 1, COMPRESSION_GZIP= 2;

    /**
     * Reads the file until the first new line character appears. This looks for
     * Unix (\n) and Windows (\r) separators. If none is found, an empty string is returned.
     *
     * @param f the file
     * @return separator the file separator or empty string
     */
    public static String guessFileSep(File f) {
        String fileSep= "";
        BufferedReader buffy = null;
        try {
            buffy= new BufferedReader(new FileReader(f));
            char[] b= new char[1];
            while(buffy.read(b)!= -1){
                if(b[0] =='\n' || b[0] == '\r'){
                    // found a separator
                    char first = b[0];
                    fileSep += b[0];
                    // check if this is followed by another new line \r or \n
                    if(buffy.read(b) != -1){
                        if(b[0] =='\n' || b[0] == '\r'){
                            if (first == '\r' || b[0] == '\r')
                                fileSep += b[0];
                        }
                    }
                    //return fileSep;
                    break;
                }
            }
        }catch (IOException e){
            Log.error("Unable to identify newline character in " + f.getAbsolutePath(), e);
        }finally {
            if(buffy != null){
                try {
                    buffy.close();
                } catch (IOException ignore) {
                    // ignore
                }
            }
        }
        return fileSep;
    }

	public static byte getCompression(File f) {
		try {
			GZIPInputStream s= new GZIPInputStream(new FileInputStream(f));
			s.close();
			return COMPRESSION_GZIP;
		} catch (IOException ex) {
			try {
				new ZipFile(f);
				return COMPRESSION_ZIP;
			} catch (Exception e) {	// actually only ZipException, but what shalls
				return COMPRESSION_NONE;
			}
		}
	}

	public static final String[] COMPRESSION_KEYWORDS= new String[] {"none", "zip", "gzip"};
	public static byte getCompression(String compression) {
		compression= compression.toLowerCase();
		for (int i = 0; i < COMPRESSION_KEYWORDS.length; i++) 
			if (compression.equals(COMPRESSION_KEYWORDS[i]))
				return (byte) i;
		return COMPRESSION_NONE;
	}
	
	public static String getCompressionString(byte compression) {
		if (compression== COMPRESSION_NONE)
			return "none";
		if (compression== COMPRESSION_ZIP)
			return SFX_ZIP;
		if (compression== COMPRESSION_GZIP)
			return SFX_GZIP;
		return "";
	}
			
	public static void inflate(File src, File dest, byte compression) throws Exception {
		InflaterInputStream in= null;
		if (compression== COMPRESSION_ZIP) 
			in= new ZipInputStream(new FileInputStream(src));
		else if (compression== COMPRESSION_GZIP)
			in= new MultiMemberGZIPInputStream(new FileInputStream(src)); // GZIPInputStream
		if (in== null)
			return;

		int bufSize= 65536;
		BufferedOutputStream buffy= new BufferedOutputStream(new FileOutputStream(dest), bufSize);
		byte[] buf= new byte[bufSize];
		int rec= 0, perc= 0;
		long read= 0, max= src.length()* 10;

        Log.progressStart("\tinflating");

		while((rec= in.read(buf, 0, buf.length))!= -1) {	// in.available()!= 0 , fails for concat
			// not the problem
//			rec= 0;
//			while (rec== 0) {
//				rec= in.read(buf, 0, buf.length);
//				if (rec== 0)
//					try {in.wait(500);} catch (InterruptedException ex) {} // :)
//			}
			read+= rec;
            Log.progress(read, max);
			if (rec> 0)
				buffy.write(buf, 0, rec);
		}
		in.close();
		buffy.flush();
		buffy.close();
        Log.progressFinish(StringUtils.OK, true);
	}
	
	public static void deflate(File src, File dest, byte compression) throws Exception {
		DeflaterOutputStream out= null;
		if (compression== COMPRESSION_ZIP) {
			ZipOutputStream zout= new ZipOutputStream(new FileOutputStream(dest));
			zout.putNextEntry(new ZipEntry(src.getName()));
			out= zout;
		} else if (compression== COMPRESSION_GZIP)
			out= new GZIPOutputStream(new FileOutputStream(dest));
		if (out== null)
			return;

		int bufSize= 65536;
		BufferedInputStream buffy= new BufferedInputStream(new FileInputStream(src), bufSize);
		byte[] buf= new byte[bufSize];
		int rec= 0, perc= 0;
		long read= 0, max= src.length();

        Log.progressStart("deflating");
		while((rec= buffy.read(buf, 0, buf.length))> 0) {
			read+= rec;
            Log.progress(read, max);
			out.write(buf, 0, rec);
		}
		buffy.close();
		if (out instanceof ZipOutputStream)
			((ZipOutputStream) out).closeEntry();
		out.flush();
		out.close();
        Log.progressFinish();
	}
	
	static boolean silent= true;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) 
	{
		String dir = "C:\\";
		String fileName = args[0];
		int lineCount = Integer.parseInt(args[1]);
		FileHelper.tail(dir + fileName, lineCount);		
	}
 
	public static Vector tail(String fileName, int lineCount)
	{
		return tail(fileName, lineCount, 2000);
	}
	
	public static String getExtension(File f) {
		return getExtension(f.getName());
	}
	public static String getExtension(String s) {
		
		int p= -1;
		for (int i = s.length()- 1; i >= 0; --i) {
			if (s.charAt(i)== '.') {
				p= i;
				break;
			}
			if (s.charAt(i)== File.separatorChar) 
				break;
		}
		if (p< 0)
			return "";
		return s.substring(p+1);
	}
 
	/**
	 * Given a byte array this method:
	 * a. creates a String out of it
	 * b. reverses the string
	 * c. extracts the lines
	 * d. characters in extracted line will be in reverse order, 
	 *    so it reverses the line just before storing in Vector.
	 *     
	 *  On extracting required numer of lines, this method returns TRUE, 
	 *  Else it returns FALSE.
	 *   
	 * @param bytearray
	 * @param lineCount
	 * @param lastNlines
	 * @return
	 */
	private static boolean parseLinesFromLast(byte[] bytearray, int lineCount, Vector lastNlines)
	{
		String lastNChars = new String (bytearray);
		StringBuffer sb = new StringBuffer(lastNChars);
		lastNChars = sb.reverse().toString();
		StringTokenizer tokens= new StringTokenizer(lastNChars,"\n");
		while(tokens.hasMoreTokens())
		{
			StringBuffer sbLine = new StringBuffer((String)tokens.nextToken());			
			lastNlines.add(sbLine.reverse().toString());
			if(lastNlines.size() == lineCount)
			{
				return true;//indicates we got 'lineCount' lines
			}
		}
		return false; //indicates didn't read 'lineCount' lines
	}
	
	
	
	public static boolean canWrite(File f) {
		if (f== null)
			return false;
		if (!f.exists()) {			
			f= f.getParentFile();
			if (f== null|| !f.exists()) {
				return false;
			}
			return canWriteToDir(f);
		} else {
			if (f.isDirectory())
				return canWriteToDir(f);
			else {
				try {
					System.getSecurityManager().checkWrite(f.getAbsolutePath());
					return true;
				} catch (SecurityException e) {
					return false;
				}
				//return f.canWrite();
			}
		}
	
		
	}
	
	public static boolean canWriteToDir(File f) {
		if (f.canWrite())
			return true;
		else {
			/*	Bug ID:   	 4939819
			 *  When I investigated this more closely, it turns out
			 *  those two directories had "read-only" attribute bit set,
			 *  so in a way JFileChooser (and File) was working
			 *  correctly (though read-only attribute has no effect on
			 *  Windows directories).
			 */
			// do it always
//			if (SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_WINNT
//					|| SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_VISTA) {
				try {
					File delme= File.createTempFile("SystemInspector", "canWrite", f);
					delme.delete();
					return true;
				} catch (Exception e) {
					return false;
				}
//			} else
//				return false;
		}
	}

	/**
	 * Reads last N lines from the given file. File reading is done in chunks.
	 * 
	 * Constraints:
	 * 1 Minimize the number of file reads -- Avoid reading the complete file
	 * to get last few lines.
	 * 2 Minimize the JVM in-memory usage -- Avoid storing the complete file 
	 * info in in-memory.
	 *
	 * Approach: Read a chunk of characters from end of file. One chunk should
	 * contain multiple lines. Reverse this chunk and extract the lines. 
	 * Repeat this until you get required number of last N lines. In this way 
	 * we read and store only the required part of the file.
	 * 
	 * 1 Create a RandomAccessFile.
	 * 2 Get the position of last character using (i.e length-1). Let this be curPos.
	 * 3 Move the cursor to fromPos = (curPos - chunkSize). Use seek().
	 * 4 If fromPos is less than or equal to ZERO then go to step-5. Else go to step-6
	 * 5 Read characters from beginning of file to curPos. Go to step-9.
	 * 6 Read 'chunksize' characters from fromPos.
	 * 7 Extract the lines. On reading required N lines go to step-9.
	 * 8 Repeat step 3 to 7 until 
	 *			a. N lines are read.
	 *		OR
	 *			b. All lines are read when num of lines in file is less than N. 
	 * Last line may be a incomplete, so discard it. Modify curPos appropriately.
	 * 9 Exit. Got N lines or less than that.
	 *
	 * @param fileName
	 * @param lineCount
	 * @param chunkSize
	 * @return
	 */
	public static Vector tail(String fileName, int lineCount, int chunkSize)
	{
		try
		{
			RandomAccessFile raf = new RandomAccessFile(fileName,"r");
			Vector lastNlines = new Vector();			
			int delta=0;
			long curPos = raf.length() - 1;
			long fromPos;
			byte[] bytearray;
			while(true)
			{				
				fromPos = curPos - chunkSize;
				System.out.println(curPos);
				System.out.println(fromPos);				
				if(fromPos <= 0)
				{
					raf.seek(0);
					bytearray = new byte[(int)curPos];
					raf.readFully(bytearray);
					parseLinesFromLast(bytearray, lineCount, lastNlines);
					break;
				}
				else
				{					
					raf.seek(fromPos);
					bytearray = new byte[chunkSize];
					raf.readFully(bytearray);
					if(parseLinesFromLast(bytearray, lineCount, lastNlines))
					{
						break;
					}
					delta = ((String)lastNlines.get(lastNlines.size()-1)).length();
					lastNlines.remove(lastNlines.size()-1);
					curPos = fromPos + delta;	
				}
			}
			Enumeration e = lastNlines.elements();
			while(e.hasMoreElements())
			{
				System.out.println(e.nextElement());
			}			
			return lastNlines;
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * @deprecated
	 * @param fileName
	 * @return
	 */
	public static int countLines_bug(String fileName) {
		final String delim="\n\r";
		try {
			int cntLines= 0;
			BufferedReader buffy= new BufferedReader(new FileReader(fileName));
			char[] cbuf= new char[80];
			boolean endNL= false;
			for (int x ; (x= buffy.read(cbuf)) != -1;) {
				for (int i = 0; i < x; i++) {
					if (delim.indexOf(cbuf[i])>= 0) {
						if (i== 0&& endNL)
							endNL= false;
						else {
							++cntLines;
							if (i== x-1)
								endNL= true;
							else {
								endNL= false;
								if (delim.indexOf(cbuf[i+1])>= 0)
									++i;
							}
						}
					}
				}
			}
			buffy.close();
			return cntLines;	// NOT +1
		} catch (Exception e) {
			; // :)
		}
		return -1;
	}

    /**
     * Reads the complete file and counts the line numbers. Error are catched and this returns -1 if
     * the number of lines could not be counted successfully.
     *
     * @param file the filename
     * @return lines the number of lines or -1 in case of any errors
     */
	public static int countLines(File file) {
        BufferedReader buffy = null;
		try {
			int cntLines= 0;
			buffy= new BufferedReader(new FileReader(file));
			for(String s;(s= buffy.readLine())!= null;++cntLines);
			return cntLines;
		} catch (Exception e) {
			; // :)
		}finally {
            if(buffy != null){
                try {
                    buffy.close();
                } catch (IOException ignore) {
                    // ignore
                }
            }
        }
		return -1;
    }
    /**
     * Reads the complete file and counts the line numbers. Error are catched and this returns -1 if
     * the number of lines could not be counted successfully.
     *
     * @param fileName the filename
     * @return lines the number of lines or -1 in case of any errors
     */
	public static int countLines(String fileName) {
        return countLines(new File(fileName));
	}
	
	public static boolean move(File from, File to) {
		return move(from, to, null);
	}
	
	/**
	 * <code>File.renameTo()</code> fails on UNIX to move files between different FSs.
	 * Here is an workaround.
	 * @param from
	 * @param to
	 * @return
	 */
	public static boolean move(File from, File to, String msg) {

        Log.progressStart(msg != null ? msg : "moving");

		if (to.exists())
			to.delete();
		
		boolean ok= from.renameTo(to);
		if (ok)
			return true;
		try { 
			fastChannelCopy(from, to, false);
		} catch (Exception e) {
            Log.error("[FATAL] couldn't copy: "+ e.getMessage(), e);
			return false;
		}
		Log.progressFinish();
		return from.delete();
	}
	
	public static boolean copy(File from, File to) {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(from));
			BufferedWriter wright= new BufferedWriter(new FileWriter(to));
			char[] cbuf= new char[1024];
			long bytesRead= 0, bytesTotal= from.length();
			int perc= 0;
			for (int c= 0, nb= 0; (nb= buffy.read(cbuf))>= 0; ++c) {
				bytesRead+= nb;
				if (!silent)
					StringUtils.printPercentage(perc, bytesRead, bytesTotal, System.err);
				wright.write(cbuf, 0, nb);
				if (c% 10== 0)
					wright.flush();
			}
			buffy.close();
			wright.flush();
			wright.close();
			return true;
			
		} catch (Exception e) {
			return false;
		}
	}

	public static void setSilent(boolean silent) {
		FileHelper.silent = silent;
	}
	
	public static File replaceSfx(File org, String newSfx) {
		String s= org.getAbsolutePath();
		int p= s.lastIndexOf('.');
		if (p< 0)
			return new File(s+ newSfx);
		return new File(s.substring(0, p)+ newSfx);
	}

	public static boolean checkForOverwrite(PrintStream p, File f) {
		if (!f.exists())
			return true;
		p.println("Confirm overwriting file "+f+" (y/n)");
		int b= 'n';
		try {
			b= System.in.read();
		} catch (Exception e) {
			p.println("Could not confirm for overwriting output file "+f);
			return false;	// eg, no stdin available
		}
		if (b== 'y'|| b== 'Y') {
			if (f.isDirectory())
				rmDir(f);
			else
				f.delete();
			return true;
		}
		return false;
	}

	public static boolean rmDir(File dir) {
		boolean success= true;
		if (dir.isDirectory()) {
			String[] files= dir.list();
			for (int i = 0; i < files.length; i++) 
				success&= rmDir(new File(dir.getAbsolutePath()+File.separator+files[i]));
		}
		success&= dir.delete();
		return success;
	}
	
	/**
	 * break a path down into individual elements and add to a list.
	 * example : if a path is /a/b/c/d.txt, the breakdown will be [d.txt,c,b,a]
	 * @param f input file
	 * @return a List collection with the individual elements of the path in
reverse order
	 */
	private static List getPathList(File f) {
		List l = new ArrayList();
		File r=null;
		try {
			r = f.getCanonicalFile();
			while(r != null) {
				l.add(r.getName());
				r = r.getParentFile();
			}
		}
		catch (IOException e) {
			e.printStackTrace();
			l = null;
		}
		return l;
	}

	/**
	 * figure out a string representing the relative path of
	 * 'f' with respect to 'r'
	 * @param r home path
	 * @param f path of file
	 */
	private static String matchPathLists(List r,List f) {
		int i;
		int j;
		String s;
		// start at the beginning of the lists
		// iterate while both lists are equal
		s = "";
		i = r.size()-1;
		j = f.size()-1;

		// first eliminate common root
		while((i >= 0)&&(j >= 0)&&(r.get(i).equals(f.get(j)))) {
			i--;
			j--;
		}

		// for each remaining level in the home path, add a ..
		for(;i>0;i--) {
			s += ".." + File.separator;
		}

		// for each level in the file path, add the path
		for(;j>0;j--) {
			s += f.get(j) + File.separator;
		}

		// file name
		s += f.get(j);
		return s;
	}

	/**
	 * get relative path of File 'f' with respect to 'home' directory
	 * example : home = /a/b/c
	 *           f    = /a/d/e/x.txt
	 *           s = getRelativePath(home,f) = ../../d/e/x.txt
	 * @param home base path, should be a directory, not a file, or it doesn't
make sense
	 * @param f file to generate path for
	 * @return path from home to f as a string
	 */
	public static String getRelativePath(File home,File f){
		File r;
		List homelist;
		List filelist;
		String s;

		homelist = getPathList(home);
		filelist = getPathList(f);
		
		s = matchPathLists(homelist,filelist);

		return s;
	}


    /**
     * Converts a given File as a Unix path relative to the given directory. If directory is
     * null, the current working directory is used
     *
     * @param file the source file
     * @param directory the directory (null permitted)
     * @return path to the file relative to the directory
     */
    public static String toRelative(File file, File directory) {
        String absFile = file.getAbsolutePath();
        String absDir = directory.getAbsolutePath();
        absFile = absFile.substring(absDir.length() + 1);
        //absFile = absFile.replace('\\', '/');
        return absFile;
    }

    /**
     * Resolves a file from a relative path. If the given base directory is null, the current working
     * directory is used
     *
     * @param rel the relative path
     * @param dir the directory to which the path is relative to
     * @return file referenced by the relative path
     */
    public static File fromRelative(String rel, File dir) {
        if(isAbsolute(rel)) return new File(rel);

        // find the right splitter
        String sep = System.getProperty("file.separator");
        String matcher = sep.equals("/") ? "\\/" : "\\\\";
        String[] split = rel.split(matcher);
        File f;
        if(dir != null){
            f = new File(dir.getAbsolutePath());
        }else{
            f = new File(".");
        }
        for (String p : split) {
            f = new File(f, p);
        }
        return f;
    }

    /**
     * Returns true if the given path is absolute
     *
     * @param path the path
     * @return absolute true if path is absolute
     */
    public static boolean isAbsolute(String path){
        if(path.startsWith("/")) return true;
        if (path.matches("^.:\\\\.*")) return true;
        return false;
    }


	public static void fastChannelCopy(File inputFile, File outputFile, boolean append) throws IOException {
		final InputStream input = new FileInputStream(inputFile);  
		final OutputStream output = new FileOutputStream(outputFile, append);  
			// get an channel from the stream  
		final ReadableByteChannel inputChannel = Channels.newChannel(input);  
		final WritableByteChannel outputChannel = Channels.newChannel(output);  
		fastChannelCopy(inputChannel, outputChannel, inputFile.length());
			// closing the channels  
		inputChannel.close();  
		outputChannel.close();
	}
	
	public static void fastChannelCopy(final ReadableByteChannel src, final WritableByteChannel dest, double max) throws IOException {
			final ByteBuffer buffer = ByteBuffer.allocateDirect(16 * 1024);
			int perc= 0, x= -1;
			long tot= 0;
			while ((x= src.read(buffer)) != -1) {
				tot+= x;
                Log.progress(tot, (long) max);
				// prepare the buffer to be drained
				buffer.flip();  
				// write to the channel, may block  
				dest.write(buffer);  
				// If partial transfer, shift remainder down  
				// If buffer is empty, same as doing clear()  
				buffer.compact();  
			}  
			// EOF will leave buffer in fill state  
			buffer.flip();  
			// make sure the buffer is fully drained.  
			while (buffer.hasRemaining()) {  
				dest.write(buffer);  
			}  
		}

	public static int cleanup(String nameDir, String pfx,
			String sfx, PrintStream stream) {
		
		File dir= new File(nameDir);
		if ((!dir.exists())|| (!dir.isDirectory()))
			return -1;
		
		if (stream!= null)
			System.err.println("[MRPROPER] cleaning all files " +
					(pfx== null?"":"with prefix "+ pfx)+ 
					(sfx== null?"":"with suffix "+ sfx)+
					" in folder:\n\t"+ nameDir);
		String[] fNames= dir.list();
		int cnt= 0, cntFail= 0;
		for (int i = 0; i < fNames.length; i++) {
			if ((pfx!= null&& !fNames[i].startsWith(pfx))
				|| (sfx!= null&& !fNames[i].endsWith(sfx)))
				continue;
			File f= new File(dir+ File.separator+ fNames[i]);
			//System.err.println(f.getName());
			boolean failed= false;
			if (f.isDirectory()) {
				if (!rmDir(f)) 
					failed= true;
			} else if (!f.delete()) {
				failed= true;
			}
			if (failed) {
				++cntFail;
//				if (stream!= null)
//					System.err.println("\tfailed to remove "+fNames[i]);
			} else {
				++cnt;
//				if (stream!= null)
//					System.err.println("\tremoved "+fNames[i]);
			}
		}
		if (stream!= null)
			System.err.println("\tremoved "+cnt+" files, failed to remove "+
					cntFail+ " files.");

		return cntFail;
	}

	private static boolean zipRecursive(ZipOutputStream zos, String pfx, File in, int depth) {
		if (depth== 0&& !silent) {
            Log.progressStart("zipping");
		}
		
		if (in.isDirectory()) {
			
			pfx= pfx== null? in.getName(): pfx+ "/"+ in.getName();
			File[] files= in.listFiles();
			int perc= 0;
			for (int i = 0; i < files.length; i++) {
				if (depth== 0&& !silent) {
                    Log.progress(i, files.length);
				}
				if (!zipRecursive(zos, pfx, files[i], depth+ 1))
					return false;
			}
			return true;
		} else {
			try {
				String name= pfx== null?in.getName():pfx+"/"+in.getName();
				zos.putNextEntry(new ZipEntry(name));
				BufferedReader buffy= new BufferedReader(new FileReader(in));
				// TODO read in buffer and write at once?
				int perc= 0;
				long size= in.length(), bread= 0l;
				for (String s= null; (s = buffy.readLine())!= null; ) {
					bread+= s.length()+ 1;
					if (depth== 0&& !silent) {
                        Log.progress(bread, size);
					}
					zos.write(s.getBytes());
					zos.write(System.getProperty("line.separator").getBytes());
				}
				zos.closeEntry();
				buffy.close();
				if (depth== 0&& !silent) {
                    Log.progressFinish(StringUtils.OK, true);
				}
				return true;
			} catch (Exception e) {
				return false;
			}
		}
	}
	
	public static boolean zip(File in, File out) {
		
		if (in== null|| !in.exists())
			return false;

		boolean returnVal= true;
		try {
		    ZipOutputStream zos = new ZipOutputStream(
		    		new BufferedOutputStream(new FileOutputStream(out)));
		    returnVal= zipRecursive(zos, null, in, 0);
		    zos.flush();
		    zos.close();
			System.currentTimeMillis();
		} catch (Exception e) {
			return false;
		}
		
		return returnVal;
	}

	public static final String SFX_ZIP= "zip", SFX_GZIP= "gz";
	public static String getCompressionExtension(byte compression) {
		if (compression== COMPRESSION_NONE)
			return null;
		if (compression== COMPRESSION_ZIP)
			return SFX_ZIP;
		if (compression== COMPRESSION_GZIP)
			return SFX_GZIP;
		return null;
	}

	public static boolean zipFilesInDir(File dir, File dst) {
		
		try {
            if (!silent){
                Log.progressStart("zipping");
            }
			File[] ff= dir.listFiles();
			ZipOutputStream zos = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(dst)));


			for (int i = 0; i < ff.length; i++) {
				if (!silent) {
                    Log.progress(i, ff.length);
				}
				if (ff[i].isDirectory()) 
				    if (!zipRecursive(zos, null, ff[i], 0))
				    	return false;
				String name= ff[i].getName();
				zos.putNextEntry(new ZipEntry(name));
				BufferedReader buffy= new BufferedReader(new FileReader(ff[i]));
				// TODO read in buffer and write at once?
	//			int perc= 0;
	//			long size= ff[i].length(), bread= 0l;
				for (String s= null; (s = buffy.readLine())!= null; ) {
	/*				bread+= s.length()+ 1;
					if (!silent) {
						if (bread* 10f/ size> perc) {
							++perc;
							if (progress!= null)
								progress.progress();
							else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
								System.err.print(Constants.STAR);
								System.err.flush();
							}
						}
					}
	*/				
					zos.write(s.getBytes());
					zos.write(System.getProperty("line.separator").getBytes());
				}
				zos.closeEntry();
				buffy.close();
			}
		    zos.flush();
		    zos.close();
            if(!silent){
                Log.progressFinish();
            }
		} catch (Exception e) {
            if(!silent){
                Log.progressFailed("ERROR");
            }
            Log.error("Error while zipping directory: " + e.getMessage(),e);
			return false;
		}
		return true;
	}

    public static String append(String s, String sfx) {
        return append(s, sfx, false, null);
    }

    public static String append(String s, String sfx, boolean stripSfx, String newSfx) {
        if (newSfx== null)
            newSfx= getExtension(s);
        newSfx= '.'+ newSfx;
        int p= s.lastIndexOf('.');
        String nuFname= (p>= 0)?
            s.substring(0, p)+ sfx+ (stripSfx?"":(newSfx== null? "": newSfx)):
            s+ sfx+ newSfx;

        return nuFname;
    }

    /**
     * Creates a temp file using the name as prefix and appending an optional extension
     *
     * @param name the name
     * @param ext the extension (null permitted)
     * @return file temp file
     * @throws IOException in case of any errors
     */
    public static File createTempFile(String name, String ext) throws IOException {
        if(ext != null && ext.length() > 0 && !ext.startsWith(".")) ext = "."+ext;
        return File.createTempFile(name, ext != null && ext.length() > 0 ? ext : "");
    }
}
