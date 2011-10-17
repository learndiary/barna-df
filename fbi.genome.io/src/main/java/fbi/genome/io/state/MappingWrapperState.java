package fbi.genome.io.state;

public class MappingWrapperState {

	public static final int STATE_OK= 0;
	
	public static final int STATE_END_OF_CHROMOSOME= 1;
	
	public static final int STATE_END_OF_FILE= 2;
	
	public static final int STATE_CHROMOSOME_NOT_FOUND= 3;

	public long count= 0l;
	public Object result= null;
	public byte state= STATE_OK;
	public String nextChr= null;
	
	public MappingWrapperState() {
		reset();
	}
	
	public void reset() {
		count= 0l;
		result= null;
		state= STATE_OK;
		nextChr= null;
	}
}
