package fbi.genome.model.commons;

public class DoubleVector {

	static final int DEFAULT_LOAD_FACTOR= 13;
	double[] vector;
	public int length= -1;
	int incrementSize;
	
	public DoubleVector() {
		this(DEFAULT_LOAD_FACTOR);
	}
	
	public DoubleVector(int initialCapacity) {
		this(initialCapacity, DEFAULT_LOAD_FACTOR);
	}
	
	public DoubleVector(int initialCapacity, int loadFactor) {
		vector= new double[initialCapacity];
		incrementSize= loadFactor;
		length= 0;
	}
	
	public double getValue(int pos) {
		if (pos< length)
			return vector[pos];
		return 0;
	}
	
	public void add(double x) {
		if (length== vector.length)
			extendVector();
		vector[length++]= x;
	}
	
	void extendVector() {
		double[] newVector= new double[vector.length+ incrementSize];
		for (int i = 0; i < vector.length; i++) 
			newVector[i]= vector[i];
		vector= newVector;
	}
	
	public int size() {
		return length;
	}
	
	public double remove(int pos) {
		double result= vector[pos];
		--length;
		for (int i = pos; i < length; i++) 
			vector[pos]= vector[pos+1];
		return result;
	}
	
	public double get(int pos) {
		return vector[pos];
	}
	
	public void set(int pos, double val) {
		if (pos< vector.length)
			vector[pos]= val;
	}
	
	public double[] toDoubleArray() {
		double[] result= new double[length];
		for (int i = 0; i < length; i++) 
			result[i]= vector[i];
		return result;
	}
	
	public String toString() {
		StringBuffer sb= new StringBuffer("[");
		for (int i = 0; i < length; i++) {
			sb.append(java.lang.Double.toString(vector[i]));
			if (i< length-1)
				sb.append(",");
		}
		sb.append("]");
		return sb.toString();
	}

	public void insert(double val, int p) {
		
		if (p< 0)
			p= (p+1)* (-1);
		
		double[] newA= new double[vector.length+ 1];
		for (int i = 0; i < p; i++) 
			newA[i]= vector[i];
		newA[p]= val;
		for (int i = p+1; i < newA.length; i++) 
			newA[i]= vector[i-1];
		
		vector= newA;
		length= newA.length;
	}

	public void removeAll() {
		this.length= 0;
	}
}
