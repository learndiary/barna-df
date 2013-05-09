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

package barna.io.tools;

import barna.commons.RandomFactory;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.Tool;
import barna.commons.log.Log;
import barna.io.FileHelper;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Class to create a subset of line from a file. 
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */

public class Subsetter implements Tool<Void> {

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
	public void setInput(File input) {
		this.input= input;
	}
	
	/**
	 * Set the (approximate) number of lines to be extracted from the input file 
	 * @param number of lines in the subset to be retrieved 
	 */
	public void setNumber(int number) {
		this.numberLines= number;
	}
	
	/**
	 * Set the output file to which is written.
	 * @param output output file to which is written
	 */
	public void setOutput(File output) {
		this.output= output;
	}
	
	/**
	 * Set the number of lines in the input file
	 * @param inputLines number of lines in the input file
	 */
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
			Random random= RandomFactory.get();
			for (String s= null; (s= buffy.readLine())!= null;) {
				double t= random.nextDouble();
				if (t> p) 
					continue;
				
				writer.write(s+ barna.commons.system.OSChecker.NEW_LINE);
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


    @Override
    public String getName() {
        return "subsetter";
    }

    @Override
    public String getDescription() {
        return "Extract a random subset of lines from a file";
    }
    @Override
    public String getLongDescription() {
        return null;
    }


    @Override
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        parameters.add(JSAPParameters.flaggedParameter("input", 'i').type(File.class).help("Input File").valueName("input").required().get());
        parameters.add(JSAPParameters.flaggedParameter("output", 'o').type(File.class).help("Output File").valueName("output").required().get());
        parameters.add(JSAPParameters.flaggedParameter("number", 'n').help("Number of lines to be subset").valueName("subset").get());
        parameters.add(JSAPParameters.flaggedParameter("lines", 'l').help("Number of lines in the input file").valueName("lines").get());
        return parameters;
    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        setInput(args.getFile("input"));
        setOutput(args.getFile("output"));
        if(args.userSpecified("number")) setNumber(Integer.parseInt(args.getString("number")));
        if(args.userSpecified("lines")) setLines(Integer.parseInt(args.getString("lines")));

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
			Log.error("Cannot write to " + output.getAbsolutePath());
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
			Random random= RandomFactory.get();
			for (String s= null, last= null; (s= buffy.readLine())!= null;last= s) {
				double t= random.nextDouble();
				if (t> p) 
					continue;
				
				String[] ss= s.split("\\s");
				if (ss[3].charAt(ss[3].length()- 3)== 'S') {
					writer.write(s+ barna.commons.system.OSChecker.NEW_LINE);
					s= buffy.readLine();
					writer.write(s+ barna.commons.system.OSChecker.NEW_LINE);
				} else if (last!= null) {
					writer.write(last+ barna.commons.system.OSChecker.NEW_LINE);
					writer.write(s+ barna.commons.system.OSChecker.NEW_LINE);
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
