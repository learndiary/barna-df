/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.genome.model;

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
