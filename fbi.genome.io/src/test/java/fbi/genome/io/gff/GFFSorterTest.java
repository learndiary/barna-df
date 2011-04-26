package fbi.genome.io.gff;

import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import static junit.framework.Assert.assertEquals;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class GFFSorterTest {

    private static File unsorted;
    private static File sorted;

    @BeforeClass
    public static void setUp(){
        unsorted = new File(GFFSorterTest.class.getResource("/testGtfSort.gtf").getFile());
        sorted = new File(GFFSorterTest.class.getResource("/testGtfSort.gtf.sorted").getFile());
    }


    @Test
    public void testGFFSort() throws IOException {
        File sortedFile = GFFSorter.sort(unsorted);

        BufferedReader o = new BufferedReader(new FileReader(sorted));
        BufferedReader s = new BufferedReader(new FileReader(sortedFile));
        String l1, l2;
        try{
            int c = 0;
            while( (l1 = o.readLine()) != null && (l2 = s.readLine()) != null){
                assertEquals("Missmatch in line " + c,l1, l2);
                c++;
            }
        }finally {
            o.close();
            s.close();
        }
    }
}
