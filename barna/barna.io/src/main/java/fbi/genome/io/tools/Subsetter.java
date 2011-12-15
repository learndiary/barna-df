package fbi.genome.io.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.Random;

import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import fbi.commons.Log;
import fbi.commons.flux.FluxTool;
import fbi.commons.options.HelpPrinter;
import fbi.genome.io.FileHelper;

/**
 * Class to create a subset of line from a file. 
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
@Cli(name = "subsetter", description = "Extract a random subset of lines from a file")
public class Subsetter implements FluxTool<Void> {

	/**
	 * Input file from which is read.
	 */
	File input;

	/**
	 * Output file to which the subset is written.
	 */
	File output;
	
	/**
	 * Number of lines in the input file.
	 */
	int inputLines= -1;

	/**
	 * Number of lines to be subset.
	 */
	int numberLines= -1;

	/**
	 * Set the input file from which is read.
	 * @param input input file from which is read
	 */
	@Option(name = "i", longName = "input", description = "Input File", required = true)
	public void setInput(File input) {
		this.input= input;
	}
	
	/**
	 * Set the (approximate) number of lines to be extracted from the input file 
	 * @param number of lines in the subset to be retrieved 
	 */
	@Option(name = "n", longName = "number", description = "Number of lines to be subset", required = true)
	public void setNumber(int number) {
		this.numberLines= number;
	}
	
	/**
	 * Set the output file to which is written.
	 * @param output output file to which is written
	 */
	@Option(name = "o", longName = "output", description = "Output File", required = false)
	public void setOutput(File output) {
		this.output= output;
	}
	
	/**
	 * Set the number of lines in the input file
	 * @param inputLines number of lines in the input file
	 */
	@Option(name = "l", longName = "lines", description = "Number of lines in the input file", required = false)
	public void setLines(int inputLines) {
		this.inputLines= inputLines;
	}

	/**
	 * Subset the file.
	 */
	@Override
	public Void call() throws Exception {
		
		// init output, if file
		if (output!= null) 
			Log.outputStream= new PrintStream(output);

		// doit
		subset(input, inputLines, Log.outputStream, numberLines);
		
		// close output, if file
		if (output!= null)
			Log.outputStream.close();
		
		return null;
	}
	
	/**
	 * Subsets the number of lines in a file to extract (about) a desired number.
	 * @param inputFile file from which line superset are read
	 * @param inputLines number of liens in the input file, set to <=0 if unknown 
	 * @param ostream stream to which the subset of lines retrieved from the input
	 * is written
	 * @param numberLines number of lines to be extracted from the input, has to be 
	 * strictly <= inputLines
	 */
	public static void subset(File inputFile, long inputLines, OutputStream ostream, int numberLines) {
		
		// get total line count
		if (inputLines<= 0) 
			inputLines= FileHelper.countLines(inputFile);
		if (numberLines>= inputLines) 
			throw new RuntimeException("Number of lines in subset has to be strictly less " +
					"than number of lines in the input ("+ numberLines+ " >= "+ inputLines+ ")");

		// subset
		BufferedReader buffy= null;
		BufferedWriter writer= null;
		try {
			
			buffy= new BufferedReader(new FileReader(inputFile));
			writer= new BufferedWriter(new OutputStreamWriter(ostream));
			double p= numberLines/ (double) inputLines;
			Random random= new Random();
			for (String s= null; (s= buffy.readLine())!= null;) {
				double t= random.nextDouble();
				if (t> p) 
					continue;
				
				writer.write(s+ "\n");
			}
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			if (buffy!= null)
				try {
					buffy.close();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			if (writer!= null)
				try {
					writer.flush();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
		}
	}

	/**
	 * Check parameters.
	 */
	@Override
	public boolean validateParameters(HelpPrinter printer,
			ArgumentProcessor toolArguments) {
		if (!input.exists()) {
			Log.error("Cannot find input file "+ input.getAbsolutePath());
			return false;
		}
		
		if (!input.canRead()) {
			Log.error("Cannot read from input file "+ input.getAbsolutePath());
			return false;
		}
		
		if (output!= null&& output.exists()) {
			Log.error("Output file exists "+ output.getAbsolutePath());
			return false;
		}
		
		if (output!= null&& !output.getParentFile().exists()) {
			Log.error("Parent folder for output does not exist "+ output.getParentFile().getAbsolutePath());
			return false;
		}
		
		if (output!= null&& output.canWrite()) {
			Log.error("Cannot write to "+ output.getAbsolutePath());
			return false;
		}
		
		return true;
	}

	/**
	 * Subsets the number of lines in a file to extract (about) a desired number.
	 * NOTE: Ensures read pairs are retrieved, requires sorting according to some
	 * read identifier scheme.
	 * @deprecated works, but make it general before production
	 * @param inputFile file from which line superset are read
	 * @param inputLines number of liens in the input file, set to <=0 if unknown 
	 * @param ostream stream to which the subset of lines retrieved from the input
	 * is written
	 * @param numberLines number of lines to be extracted from the input, has to be 
	 * strictly <= inputLines
	 */
	public static void subsetPairs(File inputFile, long inputLines, OutputStream ostream, int numberLines) {
		
		// get total line count
		if (inputLines<= 0) 
			inputLines= FileHelper.countLines(inputFile);
		if (numberLines>= inputLines) 
			throw new RuntimeException("Number of lines in subset has to be strictly less " +
					"than number of lines in the input ("+ numberLines+ " >= "+ inputLines+ ")");
	
		// subset
		BufferedReader buffy= null;
		BufferedWriter writer= null;
		try {
			
			buffy= new BufferedReader(new FileReader(inputFile));
			writer= new BufferedWriter(new OutputStreamWriter(ostream));
			double p= numberLines/ (2* (double) inputLines);
			Random random= new Random();
			for (String s= null, last= null; (s= buffy.readLine())!= null;last= s) {
				double t= random.nextDouble();
				if (t> p) 
					continue;
				
				String[] ss= s.split("\\s");
				if (ss[3].charAt(ss[3].length()- 3)== 'S') {
					writer.write(s+ "\n");
					s= buffy.readLine();
					writer.write(s+ "\n");
				} else if (last!= null) {
					writer.write(last+ "\n");
					writer.write(s+ "\n");
				}
			}
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			if (buffy!= null)
				try {
					buffy.close();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			if (writer!= null)
				try {
					writer.flush();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
		}
	}
	
}
