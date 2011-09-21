package fbi.genome.io.tools;

import java.io.File;
import java.io.PrintStream;

import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import fbi.commons.Execute;
import fbi.commons.Log;
import fbi.commons.flux.FluxTool;
import fbi.commons.options.HelpPrinter;

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
	int number= -1;

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
		this.number= number;
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
		Execute.initialize(1);
		
		// init output 
		if (output!= null) 
			Log.outputStream= new PrintStream(output);

		// get total line count
		if (inputLines<= 0) 
			;
			
		Execute.shutdown();
		return null;
	}

	@Override
	public boolean validateParameters(HelpPrinter printer,
			ArgumentProcessor toolArguments) {
		// Log.error
		return true;
	}
	
}
