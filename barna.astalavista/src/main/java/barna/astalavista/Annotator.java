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

public class Annotator extends Thread {

	public static final String[] GENE_LOCALIZATION= new String[] {"cdsRef", "cdsMax", "cdsAll", "cdsEvent"};  
	public static final int  GENE_LOC_NO= -1, GENE_LOC_REFCDS= 0, GENE_LOC_MAXCDS= 1, GENE_LOC_ALLCDS= 2, GENE_LOC_EVENTCDS=3;
	
	int geneLoc= GENE_LOC_REFCDS;
	
	
	static void printUsage() {
		System.err.println("Howdy, I am the AStalavista Event Annotator.\n");
		System.err.println("I understand:\n\tAnnotator event.gtf annotation.gtf [options]\n");
		System.err.println("where [options] can be one or more of the following things:\n");
		
		// cds
		System.err.print("\t[");
		for (int i = 0; i < GENE_LOCALIZATION.length; i++) {
			System.err.print("-"+ GENE_LOCALIZATION[i]);
			if (i< GENE_LOCALIZATION.length- 1)
				System.err.print("|");
			else
				System.err.print("]\t");
		}
		System.err.println("annotate the cds");
	}
	
	
	
}
