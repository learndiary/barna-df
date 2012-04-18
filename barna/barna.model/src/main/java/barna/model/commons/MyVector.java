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

import java.util.Collection;

public class MyVector<E> extends java.util.Vector<E> {
	public MyVector() {
		super();
		// TODO Auto-generated constructor stub
	}

	public MyVector(Collection<? extends E> c) {
		super(c);
		// TODO Auto-generated constructor stub
	}

	public MyVector(int initialCapacity, int capacityIncrement) {
		super(initialCapacity, capacityIncrement);
		// TODO Auto-generated constructor stub
	}

	public MyVector(int initialCapacity) {
		super(initialCapacity);
		// TODO Auto-generated constructor stub
	}

	public Object[] getElementData() {
		return elementData;
	}
	
	public void setCapacityIncrement(int capacityIncrement) {
		this.capacityIncrement = capacityIncrement;
	}
}
