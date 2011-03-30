package fbi.genome.model.commons;

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

