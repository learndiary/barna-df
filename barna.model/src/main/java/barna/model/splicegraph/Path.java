/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
