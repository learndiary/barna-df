
package test

import org.junit.Test
import static junit.framework.Assert.assertEquals

/**
 * 
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */

class SimpleTest{

    @Test
    public void testMe(){
        assert 1==1
        assertEquals("eins", "eins")

        String parameter = "-la"
        // start process
        Process process= "ls ${parameter}".execute()
        process.waitFor()

        println "EXIT: ${process.exitValue()}"
        println "STDERR: ${process.err.text}"
        println "STDOUT: ${process.in.text}"

        // alternative to combine wait and stdout
        println "ls ${parameter}".execute().text
    }
}