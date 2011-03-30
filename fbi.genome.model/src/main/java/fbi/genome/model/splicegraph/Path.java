package fbi.genome.model.splicegraph;

import java.util.Iterator;
import java.util.Stack;

/**
 * feed forward, only add possible.
 * @author msammeth
 *
 */

public class Path {
	long[] transcripts= null;
	Stack<Edge> edgeStack;
	int length= 0;
	
	protected Path clonePath() {
		Path p= new Path((Stack<Edge>) this.edgeStack.clone());
		p.length= this.length;
		//p.sourceEdge= this.sourceEdge;
		p.transcripts= this.transcripts;
		return p;
	}
	
	public boolean isEmpty() {
		return edgeStack.isEmpty();
	}
	
	private Path(Stack<Edge> edgeM) {
		edgeStack= edgeM;
	}

	public Path() {
		edgeStack= new Stack<Edge>();
	}

	public long[] getTranscripts() {
		return transcripts;
	}

	public void setTranscripts(long[] transcripts) {
		this.transcripts = transcripts;
	}
	
	public String toString() {
		Iterator<Edge> iter= edgeStack.iterator();
		StringBuffer sb= new StringBuffer();
		while (iter.hasNext()) {
			Edge e= iter.next();
			if (sb.length()== 0) {
				sb.append(e.getTail());
				sb.append(",");
			}
			sb.append(e.getHead());
			sb.append(",");
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append(":");
		for (int i = 0; i < transcripts.length; i++) 
			sb.append(transcripts[i]);		
		
		return sb.toString();
	}
	
	public Edge getSinkEdge() {
		return (edgeStack== null)?null:edgeStack.peek();
	}
	
	public Edge removeSinkEdge() {
		return (edgeStack== null)?null:edgeStack.pop();
	}
	
	public Node getSinkNode() {
		return (edgeStack== null)?null:edgeStack.peek().head;
	}	

	public void addEdge(Edge newEdge, int newLen) {
		edgeStack.push(newEdge);
		length+= newLen;
		if (transcripts== null)
			transcripts= newEdge.transcripts;
		else
			transcripts= Graph.intersect(transcripts, newEdge.transcripts);
	}
	
	public int length() {
		return length;
	}
	
}
