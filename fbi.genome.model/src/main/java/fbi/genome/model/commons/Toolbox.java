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

/*
 * Created on May 6, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package fbi.genome.model.commons;

import java.awt.Point;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

import fbi.genome.model.Exon;
import fbi.genome.model.Gene;

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
		String baseName= fbi.genome.model.commons.MyFile.getFileNameOnly(fName);
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
