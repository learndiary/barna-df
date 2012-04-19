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

/*
 * Created on May 6, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package barna.model.commons;

import barna.model.Exon;
import barna.model.Gene;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * 
 * 
 * @author micha
 */
public class Toolbox {

	public static int seqIdentity(String seq1, String seq2) {
		int shorter= java.lang.Math.min(seq1.length(), seq2.length());
		int id= 0;
		for (int i = 0; i < shorter; i++) {
			if (Character.toUpperCase(seq1.charAt(i))== 
				Character.toUpperCase(seq2.charAt(i)))
				++id;
		}
		return id;
	}
	
	public static String getAbsFileName(String fName) {
		return new File(fName).getAbsolutePath();
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
	 * 
	 * @param fName
	 * @return absolute file name
	 */
	public static String checkFileExists(String fName) {
		File f= new File(fName);
		if (f.exists()) {
			System.out.println("File "+fName+" exists. Specify new name or press <CR> to overwrite.");
			BufferedReader r= new BufferedReader(new InputStreamReader(System.in));
			String s= "";
			try {
				s= r.readLine();
			} catch (IOException e) {
				e.printStackTrace();
			}
			s= s.trim();
			if (s.length()< 1) {
				System.out.println("Overwritten.");
				if (f.isDirectory()) 
					rmDir(f);
				f.delete();
				return getAbsFileName(fName);
			}
			String p= "";
			int pos= fName.lastIndexOf(File.separator);
			if (pos>= 0)
				p= fName.substring(0, pos);
			if (p.length()> 0)
				s= p+ File.separator+ s;
			System.out.println("Redirected to file "+ s);
			return getAbsFileName(s);
		}
		System.out.println("Output in file "+ fName);
		return getAbsFileName(fName);
	}

	/**
	 * limit are 100 files the same day
	 * @param fName
	 * @return absolute file name
	 */
	public static String getDateStampIncr(String fName) {
		String baseName= barna.model.commons.MyFile.getFileNameOnly(fName);
		String path= fName.substring(0, fName.lastIndexOf(File.separator));
		String timestamp= MyTime.getHexDate();
		int ctr= 1;
		String s= (ctr<10)?"0"+java.lang.Integer.toString(ctr):java.lang.Integer.toString(ctr);		
		fName= timestamp+"-"+s+"_"+baseName;		
		File f= new File(path+File.separator+fName);
		while (f.exists()) {
			++ctr;
			s= (ctr<10)?"0"+java.lang.Integer.toString(ctr):java.lang.Integer.toString(ctr);
			fName= timestamp+"-"+s+"_"+baseName;
			f= new File(path+File.separator+fName);
		}
		System.out.println("Output in file "+ fName);
		return getAbsFileName(path+File.separator+fName);
	}
	
	static public String[] constraintFormater(String[] sequences, Point[][] constraints) {
		
			// clone (not loosing original)
		StringBuffer[] seqs= new StringBuffer[sequences.length];
		for (int i = 0; i < sequences.length; i++) 
			seqs[i]= new StringBuffer(sequences[i].toLowerCase());
		
			// constraints to uppercase
		for (int i = 0; i < constraints.length; i++) {
			for (int j = 0; j < constraints[i].length; j++) {
				seqs[i].replace(constraints[i][j].x, (constraints[i][j].y+ 1), 
						seqs[i].substring(constraints[i][j].x, (constraints[i][j].y+ 1)).toUpperCase());
			}
		}
		
			// convert result
		String[] result= new String[seqs.length];
		for (int i = 0; i < seqs.length; i++) 
			result[i]= seqs[i].toString();
		
		return result;
	}
	
	static public Point[] extractConstraints(Gene g, int offset) {
				
		Exon[] ex= g.getExons();
		Point[] pa= new Point[ex.length];
		for (int i = 0; i < ex.length; i++) 
			pa[i]= new Point(ex[i].getStart()- offset, ex[i].getEnd()- offset);		
		
		return pa;
	}
}
