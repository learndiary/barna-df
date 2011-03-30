package fbi.genome.sequencing.rnaseq.reconstruction;

public interface Matrix {

	public void add(int p, int len, byte dir);
	public void add(int p, byte dir);
	public int get(int p1, int p2, int readLen, int[] insertMinMax, byte dir);
	public void add(int p1, int p2, int len);
	public int get(int p1, int p2, int p3, int p4, int readLen);
	public int getSum();
	public int getLength();
	public void merge(Matrix n, int readLen);
	public void fill(int[] insertSize, int readLen);
	public byte[] toByteArray();
	public int project(int[][] b);
}
