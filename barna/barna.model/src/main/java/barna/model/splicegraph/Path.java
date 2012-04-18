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

package barna.model.splicegraph;

import java.util.Iterator;
import java.util.Stack;

/**
 * feed forward, only add possible.
 * @author msammeth
 *
 */

public class Path {
	long[] transcripts= null;
	Stack<AbstractEdge> edgeStack;
	int length= 0;
	
	protected Path clonePath() {
		Path p= new Path((Stack<AbstractEdge>) this.edgeStack.clone());
		p.length= this.length;
		//p.sourceEdge= this.sourceEdge;
		p.transcripts= this.transcripts;
		return p;
	}
	
	public boolean isEmpty() {
		return edgeStack.isEmpty();
	}
	
	private Path(Stack<AbstractEdge> edgeM) {
		edgeStack= edgeM;
	}

	public Path() {
		edgeStack= new Stack<AbstractEdge>();
	}

	public long[] getTranscripts() {
		return transcripts;
	}

	public void setTranscripts(long[] transcripts) {
		this.transcripts = transcripts;
	}
	
	public String toString() {
		Iterator<AbstractEdge> iter= edgeStack.iterator();
		StringBuffer sb= new StringBuffer();
		while (iter.hasNext()) {
			AbstractEdge e= iter.next();
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
	
	public AbstractEdge getSinkEdge() {
		return (edgeStack== null)?null:edgeStack.peek();
	}
	
	public AbstractEdge removeSinkEdge() {
		return (edgeStack== null)?null:edgeStack.pop();
	}
	
	public Node getSinkNode() {
		return (edgeStack== null)?null:edgeStack.peek().head;
	}	

	public void addEdge(AbstractEdge newEdge, int newLen) {
		edgeStack.push(newEdge);
		length+= newLen;
		if (transcripts== null)
			transcripts= newEdge.transcripts;
		else
			transcripts= SplicingGraph.intersect(transcripts, newEdge.transcripts);
	}
	
	public int length() {
		return length;
	}
	
}
