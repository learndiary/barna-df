package barna.genome.io;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.nio.channels.WritableByteChannel;
import java.util.concurrent.Callable;

/**
 * Class for capsulating file copy processes: stream-to-file, 
 * file-to-stream, and stream-to-stream. The actual process of 
 * copying data from the source to the target can be executed
 * in parallel by implementing the <code>Callable</code>
 * interface.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 */
public class Copier implements Callable<Void> {
	
	/**
	 * Default size of the buffer used for copying data from 
	 * the source to the target.
	 */
	public static int DEFAULT_BUFFER_SIZE= 100;
	
	/**
	 * A source stream.
	 */
	InputStream inputStream;
	/**
	 * A target stream.
	 */
	OutputStream outputStream;
	/**
	 * A source file.
	 */
	File inputFile;
	/**
	 * A target file.
	 */
	File outputFile;
	/**
	 * The actual size of the buffer that is used to copy from
	 * the source to the target.
	 */
	int bufferSize;
	
	/**
	 * Constructor to copy from a source stream to a target stream,
	 * employing a default buffer size.
	 * @param inputStream the source stream
	 * @param outputStream the target stream
	 */
	public Copier(InputStream inputStream, OutputStream outputStream) {
		this(inputStream, outputStream, DEFAULT_BUFFER_SIZE);
	}
	/**
	 * Constructor to copy from a source file to a target stream,
	 * employing a default buffer size.
	 * @param inputFile the source file
	 * @param outputStream the target stream
	 */
	public Copier(File inputFile, OutputStream outputStream) {
		this(inputFile, outputStream, DEFAULT_BUFFER_SIZE);
	}
	/**
	 * Constructor to copy from a source stream to a target file,
	 * employing a default buffer size.
	 * @param inputStream the source stream
	 * @param outputFile the target file
	 */
	public Copier(InputStream inputStream, File outputFile) {
		this(inputStream, outputFile, DEFAULT_BUFFER_SIZE);
	}
	/**
	 * Constructor to copy from a source to a target file,
	 * employing a default buffer size.
	 * @param inputFile the source file
	 * @param outputFile the target file
	 */
	public Copier(File inputFile, File outputFile) {
		this(inputFile, outputFile, DEFAULT_BUFFER_SIZE);
	}

	/**
	 * Constructor to copy from a source stream to a target stream,
	 * employing a custom buffer size.
	 * @param inputStream the source stream
	 * @param outputStream the target stream
	 * @param bufferSize the size of the buffer used for copying
	 */
	public Copier(InputStream inputStream, OutputStream outputStream, int bufferSize) {
		this.inputStream= inputStream;
		this.outputStream= outputStream;
		this.bufferSize= bufferSize;
	}
	
	/**
	 * Constructor to copy from a source file to a target stream,
	 * employing a custom buffer size.
	 * @param inputFile the source file
	 * @param outputStream the target stream
	 * @param bufferSize the size of the buffer used for copying
	 */
	public Copier(File inputFile, OutputStream outputStream, int bufferSize) {
		this.inputFile= inputFile;
		this.outputStream= outputStream;
		this.bufferSize= bufferSize;
	}
	
	/**
	 * Constructor to copy from a source stream to a target file,
	 * employing a custom buffer size.
	 * @param inputStream the source stream
	 * @param outputFile the target file
	 * @param bufferSize the size of the buffer used for copying
	 */
	public Copier(InputStream inputStream, File outputFile, int bufferSize) {
		this.inputStream= inputStream;
		this.outputFile= outputFile;
		this.bufferSize= bufferSize;
	}
	
	/**
	 * Constructor to copy from a sourceto a target file,
	 * employing a custom buffer size.
	 * @param inputFile the source file
	 * @param outputFile the target file
	 * @param bufferSize the size of the buffer used for copying
	 */
	public Copier(File inputFile, File outputFile, int bufferSize) {
		this.inputFile= inputFile;
		this.outputFile= outputFile;
		this.bufferSize= bufferSize;
	}
	
	/**
	 * Determines the nature of source and target (i.e., whether 
	 * streams or files have been provided), and implements the
	 * copy process as reading/writing blocks of 
	 * <code>bufferSize</code> from the source to the target, 
	 * respectively.
	 * @throws Exception
	 */
	@Override
	public Void call() throws Exception {
		
		if (inputFile!= null) 
			inputStream= new FileInputStream(inputFile);
		if (outputFile!= null) 
			outputStream= new FileOutputStream(outputFile);
		
//		byte[] b= new byte[bufferSize];
//		while(inputStream.read(b)>= 0) {
//			int k= inputStream.read(b);
//			outputStream.write(b, 0, k);
//		}
		// more efficient
		final ReadableByteChannel src= Channels.newChannel(inputStream);
		final WritableByteChannel dest = Channels.newChannel(outputStream);
		int x= -1;
		long tot= 0;
		final ByteBuffer buffer = ByteBuffer.allocateDirect(bufferSize);
		while ((x = src.read(buffer)) != -1) {
            tot += x;
            //Log.progress(tot, (long) max);
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
        src.close();
        dest.close();
        
        // close if necessary
		outputStream.flush();
		if (inputFile!= null) 
			inputStream.close();
		if (outputFile!= null) 
			outputStream.close();
		
		return null;
	}
}