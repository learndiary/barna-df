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

import barna.io.SpliceGraphIO;
import barna.model.IntronModel;

import java.io.File;

public class SJextractor {

	
	public static void main(String[] args) {
		
		int usFlank= 0, dsFlank= 0;
		IntronModel iModel= null;
		File inFile= null, outFile= null;
		try {
			usFlank= Integer.parseInt(args[0]);
			dsFlank= Integer.parseInt(args[1]);
			for (int i = 2; i < args.length; i++) {
				if (args[i].equalsIgnoreCase("--iModel")) {
					iModel= new IntronModel();
					iModel.read(new File(args[++i]));
				} else if (args[i].equalsIgnoreCase("--ann")) {
					inFile= new File(args[++i]);
				} else if (args[i].equalsIgnoreCase("--genome")) {
					barna.model.Graph.overrideSequenceDirPath= args[++i];
				} else if (args[i].equalsIgnoreCase("--out")) {
					outFile= new File(args[++i]);
				}
			}
		} catch (Exception e) {
			System.err.println("[NONO] Use this: asta -c extractSJ usFlankLen dsFlankLen --ann <file> --genome <path> --iModel <file> [--out <file>]");
			System.err.println("\twhere");
			System.err.println("usLength\tis the length of the upstream flank of the splice junction");
			System.err.println("dsLength\tis the length of the downstream flank of the splice junction");
			System.err.println("--ann\tthe annotation file (GTF format)");
			System.err.println("--genome\tis the path to the genomic sequence (i.e., the chromFa directory in UCSC convention)");
			System.err.println("--iModel\ta file with the intron model");
			System.err.println("--out\tspecifies an optional output file (default: stdout)");
			System.exit(-1);
		}
		
		SpliceGraphIO.extractSpliceJunctions(usFlank, dsFlank, iModel, inFile, outFile);
	}
}
