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

package barna.astalavista;

import org.apache.commons.cli.*;

public class AttributeExtractor {
    public static final byte VERBOSE_NORMAL = 1;
    public static final byte VERBOSE_SHUTUP = 0;

    public static byte verboseLevel= VERBOSE_NORMAL;

	private static Options options;
	private static Options getOptions() {
		if (options == null) {
			options = new Options();
			options.addOption("c", "coding", false, "extract coding transcripts");
			options.addOption("n", "noncoding", false, "extract non-coding transcripts");
		}

		return options;
	}
	
	public static void main(String[] args) {
		Parser parser= new PosixParser();
		CommandLine line= null;
		try {
			line= parser.parse(getOptions(), args);
		} catch (ParseException e) {
			if (verboseLevel> VERBOSE_SHUTUP)
				System.err.println("[UPS] "+e.getMessage());
			System.exit(-1);
		}
		
		AttributeExtractor myExtractor= new AttributeExtractor();
		myExtractor.init(line);
		myExtractor.run();
	}
	
	public void init(CommandLine line) {
		; //
	}
	
	public void run() {
		; //
	}
	
	
}
