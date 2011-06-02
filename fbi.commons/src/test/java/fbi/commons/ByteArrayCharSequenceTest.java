package fbi.commons;

import org.junit.Test;

import static junit.framework.Assert.*;
import static org.junit.Assert.assertFalse;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ByteArrayCharSequenceTest {


    @Test
    public void testConstructors() throws Exception {
        CharSequence s = new ByteArrayCharSequence(100);
        assertEquals(0, s.length());

        s = new ByteArrayCharSequence("ABC");
        assertEquals(3, s.length());
        assertEquals("ABC", s.toString());

        byte[] bytes = "ABC".getBytes();
        s = new ByteArrayCharSequence(bytes);
        assertEquals(3, s.length());
        assertEquals("ABC", s.toString());

        s = new ByteArrayCharSequence(bytes, 1,2);
        assertEquals(1, s.length());
        assertEquals("B", s.toString());
    }

    @Test
    public void testSet() throws Exception {
        ByteArrayCharSequence s = new ByteArrayCharSequence("ABC");
        assertEquals(3, s.length());
        assertEquals("ABC", s.toString());

        s.set("A");
        assertEquals(1, s.length());
        assertEquals("A", s.toString());

        s.set("ABCDEF");
        assertEquals(6, s.length());
        assertEquals("ABCDEF", s.toString());

        s = new ByteArrayCharSequence("ABC");
        s.set("ABCDEF");
        assertEquals(6, s.length());
        assertEquals("ABCDEF", s.toString());

    }

    @Test
    public void testCharByteAt() throws Exception {
        ByteArrayCharSequence s = new ByteArrayCharSequence("ABC");
        assertEquals('A', s.charAt(0));
        assertEquals('B', s.charAt(1));
        assertEquals('C', s.charAt(2));

        assertEquals('A', s.byteAt(0));
        assertEquals('B', s.byteAt(1));
        assertEquals('C', s.byteAt(2));

        try{
            assertEquals('C', s.byteAt(3));
            fail();
        }catch (Exception ignore){
            // expected
        }
    }

    @Test
    public void testSubSequence() throws Exception {
        ByteArrayCharSequence s = new ByteArrayCharSequence("ABC");
        ByteArrayCharSequence sub = s.subSequence(0, 2);
        assertEquals(2, sub.length());
        assertEquals("AB", sub.toString());

        s.chars[0] = 'X';
        assertEquals("XBC", s.toString());
        assertEquals("XB", sub.toString());
    }


    @Test
    public void testTrim(){
        ByteArrayCharSequence ss = new ByteArrayCharSequence("   abc  ");
        ss.trim();
        assertEquals("abc", ss.toString());

        ss = new ByteArrayCharSequence("   abc  \t ");
        ss.trim();
        assertEquals("abc", ss.toString());
    }

    @Test
    public void testFind() throws Exception {
        ByteArrayCharSequence s = new ByteArrayCharSequence("ABC\tDEF\t123");

        assertTrue(s.find(0));
        assertTrue(s.find(1));
        assertTrue(s.find(2));
        assertFalse(s.find(3));
        assertFalse(s.find(-1));
        assertEquals("ABC", s.getToken(0).toString());
        assertEquals("DEF", s.getToken(1).toString());
        assertEquals("123", s.getToken(2).toString());
        assertNull(s.getToken(3));
        assertNull(s.getToken(-1));


    }

    @Test
    public void testFindWithEmptyField() throws Exception {
        ByteArrayCharSequence s = new ByteArrayCharSequence("ABC\t\t123");

        assertTrue(s.find(0));
        assertTrue(s.find(1));
        assertTrue(s.find(2));
        assertFalse(s.find(3));
        assertFalse(s.find(-1));
        assertEquals("ABC", s.getToken(0).toString());
        assertEquals("", s.getToken(1).toString());
        assertEquals("123", s.getToken(2).toString());
        assertNull(s.getToken(3));
        assertNull(s.getToken(-1));
    }

    @Test
    public void testGetTokenInt() throws Exception {
        ByteArrayCharSequence s = new ByteArrayCharSequence("ABC\t\t123");

        assertEquals(123, s.getTokenInt(2));
        try{
            s.getTokenInt(3);
            fail();
        }catch (Exception ignore){
            // expected
        }

        try{
            s.getTokenInt(0);
            fail();
        }catch (Exception ignore){
            // expected
        }

        try{
            s.getTokenInt(1);
            fail();
        }catch (Exception ignore){
            // expected
        }
    }

    @Test
    public void testReplaceInteger() throws Exception {
        ByteArrayCharSequence s = new ByteArrayCharSequence("ABC\t\t123");
        s.replace(0, 2);
        assertEquals("2\t\t123", s.toString());

        s.replace(1, 2);
        assertEquals("2\t2\t123", s.toString());

        s.replace(2, 2);
        assertEquals("2\t2\t2", s.toString());

        s.replace(2, 100);
        assertEquals("2\t2\t100", s.toString());
    }
    @Test
    public void testCountTokens() throws Exception {
        assertEquals(3, new ByteArrayCharSequence("ABC\t\t123").countTokens());
        assertEquals(2, new ByteArrayCharSequence("ABC\t123").countTokens());
        assertEquals(1, new ByteArrayCharSequence("ABC").countTokens());
        assertEquals(0, new ByteArrayCharSequence("").countTokens());
    }

    @Test
    public void testAppend() throws Exception {
        ByteArrayCharSequence s = new ByteArrayCharSequence("ABC");
        s.append("DEF");
        assertEquals("ABCDEF", s.toString());

        s.append(100);
        assertEquals("ABCDEF100", s.toString());

        s.append((byte)'c');
        assertEquals("ABCDEF100c", s.toString());

    }

    @Test
    public void testStartsWith() throws Exception {
        ByteArrayCharSequence s = new ByteArrayCharSequence("ABC");
        assertTrue(s.startsWith("A"));
        assertTrue(s.startsWith("AB"));
        assertTrue(s.startsWith("ABC"));
        assertFalse(s.startsWith("X"));
        assertFalse(s.startsWith("a"));
    }
    @Test
    public void testEndsWith() throws Exception {
        ByteArrayCharSequence s = new ByteArrayCharSequence("ABC");
        assertTrue(s.endsWith("C"));
        assertTrue(s.endsWith("BC"));
        assertTrue(s.endsWith("ABC"));
        assertFalse(s.endsWith("AB"));
        assertFalse(s.endsWith("a"));
    }

    @Test
    public void testCountCharacters() throws Exception {
        assertEquals(1, ByteArrayCharSequence.countCharacters(0));
        assertEquals(1, ByteArrayCharSequence.countCharacters(1));
        assertEquals(2, ByteArrayCharSequence.countCharacters(-1));
        assertEquals(2, ByteArrayCharSequence.countCharacters(10));
        assertEquals(3, ByteArrayCharSequence.countCharacters(100));
        assertEquals(4, ByteArrayCharSequence.countCharacters(1000));
        assertEquals(5, ByteArrayCharSequence.countCharacters(10000));
        assertEquals(6, ByteArrayCharSequence.countCharacters(100000));
        assertEquals(7, ByteArrayCharSequence.countCharacters(1000000));
        assertEquals(8, ByteArrayCharSequence.countCharacters(10000000));
        assertEquals(9, ByteArrayCharSequence.countCharacters(100000000));
        assertEquals(10, ByteArrayCharSequence.countCharacters(1000000000));
    }

    @Test
    public void testGetChars() throws Exception {
        byte[] b = new byte[10];

        ByteArrayCharSequence.insertAsCharacters(100, 3, b);
        assertEquals("100", new String(b, 0, 3));

        ByteArrayCharSequence.insertAsCharacters(100, 8, b);
        assertEquals("100", new String(b, 5, 3));
    }
}
