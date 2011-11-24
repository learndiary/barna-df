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

package fbi.genome.model.commons;


import java.io.PrintStream;

import fbi.commons.StringUtils;

public class MyFile extends java.io.File {
	
	public static String humanReadableSize(long size) {

		if (size< 1000)
			return Long.toString(size);
		
		char[] dimC= new char[] {'K','M','G','T'}; 
		
		float frac= 0;
		int i = dimC.length-1;
		for (; i >= 0; --i) {
			long bas= (long) Math.pow(10, 3* (i+1));
			if (size> bas) {
				frac= ((float) size)/ bas;
				break;
			}
		}
		if (frac< 10)
			return StringUtils.fprint(frac, 1)+ dimC[i];
		else
			return Integer.toString(Math.round(frac))+ dimC[i];
	}
	
	
	public static String toUnixFileSeparators(String inFName) {
		StringBuffer sb= new StringBuffer(inFName);
		for (int i = 0; i < sb.length(); i++) {
			if (sb.charAt(i)== '\\')
				sb.replace(i,i+1,"/");
		}
		return sb.toString();
	}
	
	public static boolean checkForOverwrite(PrintStream p, MyFile f) {
		if (!f.exists())
			return true;
		p.println("Confirm overwriting file "+f+" (y/n)");
		int b= 'n';
		try {
			b= System.in.read();
		} catch (Exception e) {
			System.err.println("Could not confirm for overwriting output file "+f);
			return false;	// eg, no stdin available
		}
		if (b== 'y'|| b== 'Y') {
			if (f.isDirectory())
				Toolbox.rmDir(f);
			else
				f.delete();
			return true;
		}
		return false;
	}

	public static final String SFX_GZ= "gz";
	public static final String SFX_ZIP= "zip";
	public static final String[] SFX_COMPRESSED= new String[] {SFX_GZ, SFX_ZIP};
	
	public MyFile(String name) {
		super(name); 
	}
	
	public static MyFile createNewFile(String fName, int posDate) {
		
		String timestamp= MyTime.getHexDate();
		
		int ctr= 1;
		String s= (ctr<10)?"0"+java.lang.Integer.toString(ctr):java.lang.Integer.toString(ctr);
		String ctrStr= timestamp+"-"+s;
		String pfx= fName.substring(0,posDate);
		String sfx= fName.substring(posDate, fName.length());
		
		fName= pfx+ctrStr+sfx;		
		MyFile f= new MyFile(fName);
		while (f.exists()) {
			++ctr;
			s= (ctr<10)?"0"+java.lang.Integer.toString(ctr):java.lang.Integer.toString(ctr);
			ctrStr= timestamp+"-"+s;
			fName= pfx+ctrStr+sfx;
			f= new MyFile(fName);
		}

		return new MyFile(fName);
	}
	
	public String getFileNameOnly() {
		return getFileNameOnly(getAbsolutePath());
	}
	
	public static String getFileNameOnly(String absolutePath) {
		int pos= absolutePath.lastIndexOf(MyFile.separator);
		if (pos< 0)
			return absolutePath;
		return absolutePath.substring(pos+1);
	}

	public String getExtension() {
		int pos= getFileNameOnly().lastIndexOf('.');
		if (pos< 0)
			return null;
		return getFileNameOnly().substring(pos+1);
	}
	
	public static String getExtension(String absPath) {
		String fnOnly= getFileNameOnly(absPath);
		int pos= fnOnly.lastIndexOf('.');
		if (pos< 0)
			return null;
		return fnOnly.substring(pos+1);
	}
	
	public String getFileNameWithoutExtension() {
		int pos= getFileNameOnly().lastIndexOf('.');
		if (pos< 0)
			return getFileNameOnly();
		return getFileNameOnly().substring(0, pos);
	}
		
	public static boolean hasCompressedSuffix(String fName) {
		String sfx= getExtension(fName);
		for (int i = 0; i < MyFile.SFX_COMPRESSED.length; i++) {
			if (sfx.equalsIgnoreCase(MyFile.SFX_COMPRESSED[i]))
				return true;
		}
		return false;
	}


	public static String append(String s, String sfx) {
		return append(s, sfx, false, null);
	}	

	public static String append(String s, String sfx, boolean stripSfx, String newSfx) {
		if (newSfx== null)
			newSfx= MyFile.getExtension(s);
		newSfx= '.'+ newSfx;
		int p= s.lastIndexOf('.');
		String nuFname= (p>= 0)?
			s.substring(0, p)+ sfx+ (stripSfx?"":(newSfx== null? "": newSfx)):
			s+ sfx+ newSfx;
		
		return nuFname;
	}


	public static String stripExtension(String fileName) {
		for (int i = fileName.length()- 1; i >= 0; --i) {
			if (fileName.charAt(i)== '.')
				return fileName.substring(0,i);
		}
		return fileName;
	}

}
