package fbi.genome.sequencing.rnaseq.simulation;

import java.io.File;
import java.io.FileInputStream;
import java.util.Random;

import org.junit.Assert;
import org.junit.Test;
import static junit.framework.Assert.*;

import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.Gene;
import fbi.genome.model.Transcript;
import fbi.genome.model.bed.BEDobject2;

public class SequencerTest {

	@Test
	public void testReadGeneration() throws Exception {
		
		int nrTests= 100; // nr. of tests per transcript
		int maxBoundary= 100; // nt added at beginning/end of tx
		Random rand= new Random();
		
		GFFReader reader= new GFFReader(new File(
				getClass().getResource("/sacCer2_sgdGene_fromUCSC110329.gtf").getFile()).getAbsolutePath());
		Gene[] genes;
		for (reader.read(); (genes=reader.getGenes())!= null; reader.read()) {
			for (int i = 0; i < genes.length; i++) {
				Gene g= genes[i];
				int n= g.getTranscriptCount();
				for (int j = 0; j < n; j++) {
					Transcript t= g.getTranscripts()[j];
					int tlen= t.getExonicLength();						
					for (int k = 0; k < nrTests; k++) {
						int maxLeft= maxBoundary;
						if (t.getStrand()>= 0)
							maxLeft= Math.min(maxLeft, t.get5PrimeEdge());
						int maxRight= maxBoundary;
						if (t.getStrand()< 0)
							maxRight= Math.min(maxRight, Math.abs(t.get3PrimeEdge()));
						int boundLeft= (maxLeft== 0? 0: rand.nextInt(maxLeft));
						int boundRight= (maxRight== 0? 0: rand.nextInt(maxRight));
						
						byte absDir= t.getStrand();
						int fragStart= rand.nextInt(tlen+ boundLeft+ boundRight)- boundLeft;
						int fragEnd= fragStart+ rand.nextInt(tlen+ boundRight- fragStart);
						int fragLen= fragEnd- fragStart+ 1;
						if (fragLen<= 1)
							continue;
						int readLen= 1+ rand.nextInt(fragLen- 1);
						boolean left= rand.nextBoolean();
						int start= (left? fragStart: fragEnd- readLen+ 1);
						int end= (left? fragStart+ readLen- 1: fragEnd);
						assert(start<= end);
						BEDobject2 obj= new BEDobject2();
						int polyA= Sequencer.createRead(obj, start, end, t, k, 
								absDir, fragStart, fragEnd, left, 1);
						
						assertEquals(readLen, obj.getLength());
						int bedLeft= t.getGenomicPosition(start);
						int bedRight= t.getGenomicPosition(end);
						if (absDir< 0) {
							int h= bedLeft;
							bedLeft= -bedRight;
							bedRight= -h;
						}
						--bedLeft;
						
						if (polyA> 0) {
							if (polyA== readLen) {
								assertEquals(0, obj.getStart());
								assertEquals(readLen, obj.getEnd());
							}
						} else {
							assertEquals(bedLeft, obj.getStart());
							assertEquals(bedRight, obj.getEnd());
						}						
						obj.getLength();
					}
				}
			}
		}
	}
	
}
