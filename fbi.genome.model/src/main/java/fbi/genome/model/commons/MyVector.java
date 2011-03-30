package fbi.genome.model.commons;

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
