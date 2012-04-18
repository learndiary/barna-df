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

package barna.model.commons;

public class MyInteger implements Comparable<MyInteger> {
	int value;
	
	public MyInteger(int val) {
		this.value= val;	
	}
	
	//@Override
	public int compareTo(MyInteger o) {
		return value- ((MyInteger) o).value;
	}

	public int getValue() {
		return value;
	}

	public void setValue(int value) {
		this.value = value;
	}
	
	@Override
	public String toString() {
		return Integer.toString(value);
	}
	
	public int intValue() {
		return value;
	}
}

