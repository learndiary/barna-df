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

package barna.model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Vector;

public class IntronModel {

	public static final String TAG_MODEL= "#MODEL";
	
	Vector<String[]> seqCombi= null; // 0 sin 
	int minLen= -1, maxLen= -1, maxDonLen= -1, maxAccLen= -1;
	String name= null;
	
	/**
	 * creates default GT/AG intron model
	 */
	public IntronModel() {
		seqCombi= new Vector<String[]>(1);
		seqCombi.add(new String[]{"GT","AG"});
		minLen= 69;
		maxLen= 696969;
		name= "default";
	}
	
	public int getMaxAcceptorLength() {
		if (maxAccLen < 0) {
			for (int i = 0; i < seqCombi.size(); i++) 
				maxAccLen = Math.max(maxAccLen, seqCombi.elementAt(i)[1].length()); 
		}

		return maxAccLen;
	}
	
	public int getMaxDonorLength() {
		if (maxDonLen < 0) {
			for (int j = 0; j < seqCombi.size(); j++) 
				maxDonLen = Math.max(maxDonLen, seqCombi.elementAt(j)[0].length());
		}

		return maxDonLen;
	}

	public boolean isValid(int len) {
		if ((minLen>= 0&& len<minLen)|| (maxLen>= 0&& len>maxLen))
			return false;
		return true;
	}
	
	public boolean isValid(String don, String acc) {
		
		don= don.toUpperCase();
		acc= acc.toUpperCase();
		for (int i = 0; i < seqCombi.size(); i++) {
			if (don.startsWith(seqCombi.elementAt(i)[0])&& acc.startsWith(seqCombi.elementAt(i)[1]))
				return true;
		}
		return false;
	}
	
	public boolean read(File f) {
		seqCombi= null;
		minLen= -1;
		maxLen= -1;
		name= f.getName();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(f));
			seqCombi= new Vector<String[]>();
			for (String s= null; (s= buffy.readLine())!= null;) {
				if (s.startsWith(TAG_MODEL)) {
					parseHeaderLine(s);
					continue;
				}
				String[] ss= s.split("\\s");
				ss[0]= ss[0].toUpperCase();
				ss[1]= ss[1].toUpperCase();
				assert(ss.length== 2);
				seqCombi.add(ss);
			}
			return (seqCombi.size()> 0);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return false;
	}
	
	private void parseHeaderLine(String line) {
		String[] ss= line.split("\\s");
		minLen= Integer.parseInt(ss[1]);
		maxLen= Integer.parseInt(ss[2]);
	}

	public String getName() {
		return name;
	}
}
